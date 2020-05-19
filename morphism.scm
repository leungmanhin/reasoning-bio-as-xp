(use-modules
  (ice-9 rdelim)
  (opencog)
  (opencog bioscience)
  (opencog exec)
  (opencog pln)
  (opencog ure))

(load "bio-as-utils.scm")

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

(define step-1-results
  (append-map
    cog-outgoing-set
    (append-map
      (lambda (p)
        (cog-outgoing-set (cog-execute! (evaluation-to-member (Predicate p)))))
      (list "GO_positively_regulates" "GO_negatively_regulates"))))

(write-atoms-to-file "results/pln-step-1.scm" step-1-results)

;; --- Step 2: Turn the SatisfyingSetScopeLinks created in step 1 into ConceptNodes
(format #t "--- Turning SatisfyingSetScopeLinks into ConceptNodes...\n")

; TODO: To be replaced by an actual PLN rule
(define (sat-set-to-concept s)
  (Concept
    (string-join
      (list (cog-name (gadr s)) (cog-name (gaddr s)) (cog-name (gdddr s)))
      "-")))
(define step-2-results
  (cog-outgoing-set
    (cog-execute!
      (Bind
        (VariableSet
          (TypedVariable (Variable "$X") (TypeChoice (Type "ConceptNode") (Type "GeneNode")))
          (TypedVariable (Variable "$Y") (Type "SatisfyingSetScopeLink")))
        (Present (Member (Variable "$X") (Variable "$Y")))
        (Member
          (Variable "$X")
          (ExecutionOutput
            (GroundedSchema "scm: sat-set-to-concept")
            (List (Variable "$Y"))))))))

(for-each cog-extract-recursive (cog-get-atoms 'SatisfyingSetScopeLink))

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

;; --- Step 4: Calculate and/or assign TVs
(format #t "--- Assigning TVs to all the members/subsets...\n")
(for-each
  (lambda (x) (cog-set-tv! x (stv 1 1)))
  (append (cog-get-atoms 'MemberLink) (cog-get-atoms 'SubsetLink)))
(define genes (get-genes))
(define go-categories (get-go-categories))
(define usize (length genes))
(define (concept-mean x) (exact->inexact (/ (get-cardinality x) usize)))
(define (concept-tv x) (stv (concept-mean x) (count->confidence usize)))
(define go-categories-with-tvs
  (map (lambda (x) (cog-set-tv! x (concept-tv x))) go-categories))

(write-atoms-to-file "results/pln-step-4.scm" go-categories-with-tvs)

;; --- Step 5: Infer inverse SubsetLinks
(format #t "--- Inferring inverse SubsetLinks...\n")
(define go-subsets (get-go-subsets))
; (Subset A B) -> (Subset B A)
(define inversed-go-subsets (map true-subset-inverse go-subsets))

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
    go-categories
    (cog-get-atoms 'MemberLink)
    (cog-get-atoms 'SubsetLink)
    (cog-get-atoms 'AttractionLink)
    (cog-get-atoms 'IntensionalSimilarityLink)
    (cog-get-atoms 'IntensionalDifferenceLink)))
