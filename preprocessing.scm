(use-modules (opencog) (opencog bioscience) (opencog exec) (opencog ure) (opencog pln))

(load "bio-as-utils.scm")

(primitive-load "kbs/GO_2020-04-01.scm")
(primitive-load "kbs/GO_annotation_gene-level_2020-04-01.scm")
(primitive-load "kbs/Go-Plus.scm")
(primitive-load "kbs/string_ggi_2020-03-31.scm")
(primitive-load "kbs/NCBI2Reactome_PE_Pathway.txt_2020-04-01.scm")
(primitive-load "kbs/reactome_2020-04-01.scm")
(primitive-load "kbs/smpdb_gene_2020-04-03.scm")

(pln-load 'empty)

;; --- Step 1: Turn the relations represented in EvaluationLinks into MemberLinks
(format #t "--- Turning EvaluationLinks into MemberLinks...\n")
(define go-preds (list "GO_positively_regulates" "GO_negatively_regulates"))
(define string-preds (list "reaction" "catalysis" "inhibition" "ptmod" "expression" "binding" "activation"))

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
      (append go-preds string-preds))))

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

(for-each cog-extract-recursive (cog-get-atoms 'SatisfyingSetScopeLinks))

(write-atoms-to-file "results/pln-step-2.scm" step-2-results)

;; --- Step 3: Generate SubsetLinks for the above members
(format #t "--- Directly introducing SubsetLinks for the concepts...\n")

; TODO: To be replaced by an actual PLN rule
(pln-load-from-path "pln/rules/extensional/subset-direct-introduction.scm")
(define step-3-results
  (append-map
    cog-outgoing-set
    (cog-outgoing-set
      (cog-execute!
        (Bind
          (VariableSet
            (TypedVariable (Variable "$X") (TypeChoice (Type "ConceptNode") (Type "GeneNode")))
            (TypedVariable (Variable "$Y") (Type "ConceptNode")))
          (Present
            (Variable "$X")
            (Variable "$Y"))
          (ExecutionOutput
            (GroundedSchema "scm: subset-direct-introduction")
            (List
              ;; Conclusion
              (Subset (Variable "$X") (Variable "$Y"))
              ;; Premises
              (Variable "$X")
              (Variable "$Y"))))))))

(write-atoms-to-file "results/pln-step-3.scm" step-3-results)

;; --- Step 4: Infer new members
(format #t "--- Inferring new members...\n")
(pln-load-from-path "rules/translation.scm")
(pln-load-from-path "rules/transitivity.scm")
(pln-add-rule-by-name "present-inheritance-to-subset-translation-rule")
(pln-add-rule-by-name "present-subset-transitivity-rule")
(pln-add-rule-by-name "present-mixed-member-subset-transitivity-rule")
(define step-4-results
  (cog-outgoing-set
    (pln-fc
      (Inheritance (Variable "$X") (Variable "$Y"))
      #:vardecl (VariableSet (TypedVariable (Variable "$X") (Type "ConceptNode")) (TypedVariable (Variable "$Y") (Type "ConceptNode")))
      #:maximum-iterations 12
      #:fc-full-rule-application #t)))

(write-atoms-to-file "results/pln-step-4.scm" step-4-results)

;; --- Step 5: Calculate and/or assign TVs
(format #t "--- Assigning TVs to all the members...\n")
(define genes (get-genes))
(define go-categories (get-go-categories))
(define usize (length genes))
(define (concept-mean x) (exact->inexact (/ (get-cardinality x) usize)))
(define (concept-tv x) (stv (concept-mean x) (count->confidence usize)))
(define go-categories-with-tvs
  (map (lambda (x) (cog-set-tv! x (concept-tv x))) go-categories))

;; --- Step 6: Infer inverse SubsetLinks
(format #t "--- Inferring inverse SubsetLinks...\n")
(define go-subsets (get-go-subsets))
(define inversed-go-subsets (map true-subset-inverse go-subsets))
(define inversed-go-subsets-with-pos-tvs
  (filter gt-zero-mean-and-confidence? inversed-go-subsets))

;; --- Step 7: Infer all AttractionLinks
(format #t "--- Inferring all AttractionLinks...\n")
(pln-add-rule-by-name "subset-condition-negation-rule")
(pln-add-rule-by-name "subset-attraction-introduction-rule")
(define step-7-results
  (cog-outgoing-set
    (pln-bc
      (Attraction (Variable "$X") (Variable "$Y"))
      #:vardecl
        (VariableSet
          (TypedVariable (Variable "$X") (Type "ConceptNode"))
          (TypedVariable (Variable "$Y") (Type "ConceptNode")))
      #:maximum-iterations 12
      #:complexity-penalty 10)))

(write-atoms-to-file "results/pln-step-7.scm" step-7-results)

;; --- Step 8: Output AttractionLinks with null mean
(format #t "--- Filtering out AttractionLinks...\n")
(define non-null-attractions (filter all-nodes-non-null-mean? (cog-outgoing-set step-7-results)))

(write-atoms-to-file "results/pln-attractions.scm" non-null-attractions)
