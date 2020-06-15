(use-modules
  (ice-9 rdelim)
  (opencog)
  (opencog bioscience)
  (opencog exec)
  (opencog pln)
  (opencog ure))

;; --- Utils
(load "bio-as-utils.scm")

(define (gather-cardinality)
  (define alist (list))
  (for-each
    (lambda (go) (set! alist (assoc-set! alist (cog-name go) (get-cardinality go))))
    (get-go-categories))
  alist)

;; --- Load knowledge bases
(primitive-load "kbs/GO_2020-04-01.scm")
(primitive-load "kbs/GO_annotation_gene-level_2020-04-01.scm")
(primitive-load "kbs/Go-Plus-GO_2020-05-04.scm")

;; --- Step 1: Turn the "positively/negatively regulates" relations represented in EvaluationLinks into MemberLinks
(format #t "--- Turning EvaluationLinks into MemberLinks...\n")

(pln-load 'empty)
; TODO: To be replaced by an actual PLN rule
(pln-load-from-path "pln/rules/wip/evaluation-to-member.scm")

(define (evaluation-to-member pred)
  (Bind
    (VariableList
      (Variable "$A")
      (Variable "$B"))
    (Evaluation
      pred
      (List
        (Variable "$A")
        (Variable "$B")))
    (ExecutionOutput
      (GroundedSchema "scm: evaluation-to-member-2-formula")
      (List
        (Member
          (Variable "$A")
          (SatisfyingSetScope
            (Variable "$X")
            (Evaluation
              pred
              (List
                (Variable "$X")
                (Variable "$B")))))
        (Member
          (Variable "$B")
          (SatisfyingSetScope
            (Variable "$Y")
            (Evaluation
              pred
              (List
                (Variable "$A")
                (Variable "$Y")))))
        (Evaluation
          pred
          (List
            (Variable "$A")
            (Variable "$B")))))))

(define members-sat
  (append-map
    cog-outgoing-set
    (append-map
      (lambda (p)
        (cog-outgoing-set (cog-execute! (evaluation-to-member (Predicate p)))))
      (list "GO_positively_regulates" "GO_negatively_regulates"))))

; Then turn the SatisfyingSetScopeLinks created above into ConceptNodes
(format #t "--- Turning SatisfyingSetScopeLinks into ConceptNodes...\n")

; TODO: To be replaced by an actual PLN rule
(define (sat-set-to-concept s)
  (Concept
    (string-join
      (list (cog-name (gadr s)) (cog-name (gaddr s)) (cog-name (gdddr s)))
      "-")))
(define step-1-results
  (cog-outgoing-set
    (cog-execute!
      (Bind
        (VariableSet
          (TypedVariable (Variable "$X") (Type "ConceptNode"))
          (TypedVariable (Variable "$Y") (Type "SatisfyingSetScopeLink")))
        (Present (Member (Variable "$X") (Variable "$Y")))
        (Member
          (Variable "$X")
          (ExecutionOutput
            (GroundedSchema "scm: sat-set-to-concept")
            (List (Variable "$Y"))))))))

(for-each cog-extract-recursive (cog-get-atoms 'SatisfyingSetScopeLink))

(write-atoms-to-file "results/pln-step-1.scm" step-1-results)

;; --- Step 2: Turn the MemberLinks created in step 1 into InheritanceLinks
(format #t "--- Turning MemberLinks into InheritanceLinks...\n")
(define step-2-results (map (lambda (x) (Inheritance (gar x) (gdr x))) step-1-results))

(write-atoms-to-file "results/pln-step-2.scm" step-2-results)

;; --- Step 3: Infer new members
(format #t "--- Inferring new members...\n")

(pln-load 'empty)
(pln-load-from-path "rules/translation.scm")
(pln-load-from-path "rules/transitivity.scm")
; (Inheritance C1 C2) |- (Subset C1 C2)
(pln-add-rule-by-name "present-inheritance-to-subset-translation-rule")
; (Subset C1 C2) (Subset C2 C3) |- (Subset C1 C3)
(pln-add-rule-by-name "present-subset-transitivity-rule")
; (Member G C1) (Subset C1 C2) |- (Member G C2)
(pln-add-rule-by-name "present-mixed-member-subset-transitivity-rule")

(define step-3-results
  (cog-outgoing-set
    (pln-fc
      (Inheritance (Variable "$X") (Variable "$Y"))
      #:vardecl
        (VariableSet
          (TypedVariable (Variable "$X") (Type "ConceptNode"))
          (TypedVariable (Variable "$Y") (Type "ConceptNode")))
      #:maximum-iterations 12
      #:fc-full-rule-application #t)))

(write-atoms-to-file "results/pln-step-3.scm" step-3-results)

(define inferred-go-cardinality-alist (gather-cardinality))

;; --- Step 4: Calculate and/or assign TVs
(format #t "--- Assigning TVs to all the members/subsets...\n")
(for-each
  (lambda (x) (cog-set-tv! x (stv 1 1)))
  (append (cog-get-atoms 'MemberLink) (cog-get-atoms 'SubsetLink)))
(define genes (get-genes))
(define go-terms (get-go-categories))
(define go-related-terms
  (filter
    (lambda (c)
      (or (string-prefix? "GO_positively_regulates" (cog-name c))
          (string-prefix? "GO_negatively_regulates" (cog-name c))))
    (cog-get-atoms 'ConceptNode)))
(define usize (length genes))
(define (concept-mean x) (exact->inexact (/ (get-cardinality x) usize)))
(define (concept-tv x) (stv (concept-mean x) (count->confidence usize)))
(define go-terms-with-tvs
  (map (lambda (x) (cog-set-tv! x (concept-tv x))) (append go-terms go-related-terms)))

(write-atoms-to-file "results/pln-step-4.scm" go-terms-with-tvs)

;; --- Step 5: Infer inverse SubsetLinks
(format #t "--- Inferring inverse SubsetLinks...\n")
; (Subset A B) -> (Subset B A)
(define inversed-go-subsets (map true-subset-inverse (cog-get-atoms 'SubsetLink)))

(write-atoms-to-file "results/pln-step-5.scm" inversed-go-subsets)

;; --- Step 6: Infer all AttractionLinks
(format #t "--- Inferring all AttractionLinks...\n")

(pln-load 'empty)
; (Subset A B) |- (Subset (Not A) B)
(pln-add-rule-by-name "subset-condition-negation-rule")
; (Subset A B) (Subset (Not A) B) |- (Attraction A B)
(pln-add-rule-by-name "subset-attraction-introduction-rule")

(define step-6-results
  (cog-outgoing-set
    (pln-bc
      (Attraction (Variable "$X") (Variable "$Y"))
      #:vardecl
        (VariableSet
          (TypedVariable (Variable "$X") (Type "ConceptNode"))
          (TypedVariable (Variable "$Y") (Type "ConceptNode")))
      #:maximum-iterations 12
      #:complexity-penalty 10)))

(write-atoms-to-file "results/pln-step-6.scm" step-6-results)

#!
; ----- Step 7: Estimate the intensional similarity between GOs
(format #t "--- Inferring all IntensionalSimilarityLinks...\n")
(pln-load 'empty)
(pln-add-rule-by-name "intensional-similarity-direct-introduction-rule")
(define step-7-results
  (cog-outgoing-set
    (pln-bc (IntensionalSimilarity (Variable "$X") (Variable "$Y")))))

(write-atoms-to-file "results/pln-step-7.scm" step-7-results)

; ----- Step 8: Estimate the intensional difference between GOs
(format #t "--- Inferring all IntensionalDifferenceLinks...\n")
(pln-load 'empty)
(pln-add-rule-by-name "intensional-difference-direct-introduction-rule")
(define step-8-results
  (cog-outgoing-set
    (pln-bc (IntensionalDifference (Variable "$X") (Variable "$Y")))))

(write-atoms-to-file "results/step-8-results.scm" step-8-results))
!#

;; --- Export everything that would be useful
(write-atoms-to-file "results/pln-morphism.scm"
  (append
    genes
    go-terms
    go-related-terms
    (cog-get-atoms 'EvaluationLink)
    (cog-get-atoms 'MemberLink)
    (cog-get-atoms 'InheritanceLink)
    (cog-get-atoms 'SubsetLink)
    (cog-get-atoms 'AttractionLink)
    (cog-get-atoms 'IntensionalSimilarityLink)
    (cog-get-atoms 'IntensionalDifferenceLink)))

#!
;; ----- Misc ----- ;;
(define go-aging (list
  "GO:0016853" "GO:0000423" "GO:0016236" "GO:0003756"
  "GO:0016805" "GO:0006090" "GO:0015908" "GO:0006120"
  "GO:0045190" "GO:0006412" "GO:0009060" "GO:0006642"
))

(define go-biological-process
  (map
    (lambda (g) (cog-name (gadr g)))
    (filter
      (lambda (e) (string=? "biological_process" (cog-name (gddr e))))
      (cog-incoming-by-type (Predicate "GO_namespace") 'EvaluationLink))))

(define go-biological-process-with-members
  (filter
    (lambda (g) (> (assoc-ref inferred-go-cardinality-alist g) 0))
    go-biological-process))

(define (write-list-to-file filename lst)
  (define fp (open-output-file filename))
  (for-each
    (lambda (x) (display (string-append x "\n") fp))
    lst)
  (close-port fp))

(write-list-to-file "go-aging.txt" go-aging)
(write-list-to-file "go-biological-process.txt" go-biological-process)
(write-list-to-file "go-biological-process-with-members.txt" go-biological-process-with-members)

(define target-go-list go-biological-process-with-members)
(define (shuffle original-lst shuffled-lst)
  (define n (random (length original-lst) (random-state-from-platform)))
  (define g (list-ref original-lst n))
  (if (= (length original-lst) 1)
    (cons (car original-lst) shuffled-lst)
    (shuffle (delete g original-lst) (cons g shuffled-lst))))
(define shuffled-target-go-list (shuffle target-go-list (list)))
(set! target-go-list shuffled-target-go-list)

(define gp-fp (open-output-file "go-pairs.txt"))
(do ((i 0 (+ i 2)))
    ((>= (+ i 3) (length target-go-list)))
  (display
    (string-append (list-ref target-go-list i) "," (list-ref target-go-list (+ i 1)) "\n")
    gp-fp))
(close-port gp-fp)

(pln-load 'empty)
(pln-add-rule-by-name "intensional-difference-direct-introduction-rule")
(define intensional-difference-results (list))
(define intensional-difference-csv-fp (open-output-file "results/pln-intensional-differences.csv"))
(call-with-input-file "go-pairs.txt"
  (lambda (fp)
    (let ((line (read-line fp)))
      (while (not (eof-object? line))
        (let* ((pairs (string-split line #\,))
               (go1-str (list-ref pairs 0))
               (go2-str (list-ref pairs 1))
               (go1 (Concept go1-str))
               (go2 (Concept go2-str))
               (result (cog-outgoing-set (pln-bc (IntensionalDifference go1 go2))))
               (tv-mean (cog-mean (car result)))
               (tv-conf (cog-confidence (car result))))
          (format #t "--- Done calculating intensional difference between ~a and ~a\n"
            go1-str go2-str)
          (display
            (string-append line "," (number->string tv-mean) "," (number->string tv-conf) "\n")
            intensional-difference-csv-fp)
          (set! intensional-difference-results (append intensional-difference-results result)))
        (set! line (read-line fp))))))
(write-atoms-to-file "results/pln-intensional-differences.scm" intensional-difference-results)
(close-port intensional-difference-csv-fp)

(define go-cardinality-fp (open-output-file "results/go-cardinality.csv"))
(for-each
  (lambda (x)
    (display (string-append (car x) "," (number->string (cdr x)) "\n") go-cardinality-fp))
  inferred-go-cardinality-alist)
(close-port go-cardinality-fp)
!#
