(use-modules (opencog) (opencog bioscience) (opencog exec))

(define (get-attraction source)
  (define pm-result
    (cog-execute!
      (Get
        (TypedVariable
          (Variable "$x")
          (Signature
            (Attraction
              source
              (Type "ConceptNode")
            )
          )
        )
        (Variable "$x")
      )
    )
  )

  (define results (cog-outgoing-set pm-result))

  (cog-extract pm-result)

  results
)

(define bd (ConceptNode "GO:0030183"))
(define bp (ConceptNode "GO:0042098"))
(define td (ConceptNode "GO:0030183"))
(define tp (ConceptNode "GO:0030217"))

#!
(primitive-load "kbs/biogrid_gene_gene_3.5.177.scm")
(primitive-load "kbs/GO.scm")
(primitive-load "kbs/Go-Plus.scm")
(primitive-load "kbs/GO_annotation.scm")
(primitive-load "kbs/NCBI2Reactome_PE_Pathway.txt.scm")
(primitive-load "kbs/reactome.scm")
(primitive-load "kbs/smpdb_gene.scm")

(define a1 (get-attraction bd))
(define a2 (get-attraction bp))
(define a3 (get-attraction td))
(define a4 (get-attraction tp))
(format #t "========== Looking at: ~a\n~a\n--- Length: ~a\n\n" "GO:0042098" a1 (length a1))
(format #t "========== Looking at: ~a\n~a\n--- Length: ~a\n\n" "GO:0042100" a2 (length a2))
(format #t "========== Looking at: ~a\n~a\n--- Length: ~a\n\n" "GO:0030217" a3 (length a3))
(format #t "========== Looking at: ~a\n~a\n--- Length: ~a\n\n" "GO:0030183" a4 (length a4))
!#
