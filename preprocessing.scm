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
(pln-load-from-path "pln/rules/wip/evaluation-to-member.scm")
; TODO: To be replaced by a PLN rule
(define (evaluation-to-member pred)
  (BindLink
    (VariableList
      (VariableNode "$A")
      (VariableNode "$B"))
    (EvaluationLink
      pred
      (ListLink
        (VariableNode "$A")
        (VariableNode "$B")))
    (ExecutionOutputLink
      (GroundedSchemaNode "scm: evaluation-to-member-2-formula")
        (ListLink
          (MemberLink
            (VariableNode "$A")
            (SatisfyingSetScopeLink
              (VariableNode "$X")
              (EvaluationLink
                pred
                (ListLink
                  (VariableNode "$X")
                  (VariableNode "$B")))))
          (MemberLink
            (VariableNode "$B")
            (SatisfyingSetScopeLink
              (VariableNode "$Y")
              (EvaluationLink
                pred
                (ListLink
                  (VariableNode "$A")
                  (VariableNode "$Y")))))
          (EvaluationLink
            pred
            (ListLink
              (VariableNode "$A")
              (VariableNode "$B")))))))

(define step-1-results
  (append-map
    cog-outgoing-set
    (append-map
      (lambda (p)
        (cog-outgoing-set (cog-execute! (evaluation-to-member (Predicate p)))))
      (append go-preds string-preds))))

(write-atoms-to-file "results/pln-step-1.scm" step-1-results)

;; --- Step 2: Turn the SatisfyingSetScopeLinks created in step 1 into ConceptNodes
; TODO: To be replaced by a PLN rule
(format #t "--- Turning SatisfyingSetScopeLinks into ConceptNodes...\n")
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

(write-atoms-to-file "results/pln-step-2.scm" step-2-results)

;; --- Step 3: Generate SubsetLinks for the above members
(format #t "--- Directly introducing SubsetLinks for the concepts...\n")
(define (inh-from-pair p) (Inheritance (gar p) (gdr p)))
(define (rev-inh-from-pair p) (Inheritance (gdr p) (gar p)))
(define all-members
  (append (get-member-links 'GeneNode 'ConceptNode) (get-member-links 'ConceptNode 'ConceptNode)))
(define inhs (map inh-from-pair all-members))
(define rev-inhs (map rev-inh-from-pair all-members))
(define all-inhs (append inhs rev-inhs))
(pln-add-rule-by-name "subset-direct-introduction-rule")
(define step-3-results (append (cog-get-atoms 'SubsetLink) (append-map cog-outgoing-set (map pln-bc all-inhs))))

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
(define results-lst-with-tvs
  (map
    (lambda (x) (cog-set-tv! x (stv 1 1)))
    (append
      step-4-results
      ; These are the existing members in the kbs
      (get-member-links 'GeneNode 'ConceptNode))))
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

;; --- Step 8: Filter out relationships involving GO concepts with null mean
(format #t "--- Filtering out relationships...\n")
(define non-null-go-categories-with-tvs
  (filter non-null-mean? go-categories-with-tvs))
(define non-null-results-lst-with-tvs
  (filter all-nodes-non-null-mean? results-lst-with-tvs))
(define non-null-inversed-go-subsets-with-pos-tvs
  (filter all-nodes-non-null-mean? inversed-go-subsets-with-pos-tvs))
(define non-null-attractions
  (filter all-nodes-non-null-mean? (cog-outgoing-set step-7-results)))

;; --- Step 9: Write results in file
(format #t "--- Writing final results...\n")
(define all-results
  (append
    non-null-go-categories-with-tvs
    non-null-results-lst-with-tvs
    non-null-inversed-go-subsets-with-pos-tvs
    non-null-attractions))

(write-atoms-to-file "results/pln-stage-1.scm" all-results)
