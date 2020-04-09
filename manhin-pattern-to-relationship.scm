;; Given files with mined patterns (conjunctions), reason to turn them
;; into Similarity or Inheritance relationships.

;; Parameters
(define rs 0)                           ; Random seed
(define ss 1)                         ; Subsampled portion of the KBs

;; Filename containing the mined patterns
(define mp-filename "results/mine-tb.scm")

;; Load modules & utils
(use-modules (srfi srfi-1))
(use-modules (opencog randgen))
(use-modules (opencog logger))
(use-modules (opencog ure))
(use-modules (opencog miner))
(use-modules (opencog bioscience))
(use-modules (opencog pln))
(load "bio-as-utils.scm")

(define log-filename
  (string-append "log/inheritance-links-tb"
                 "-rs=" (number->string rs)
                 "-ss=" (number->string ss)
                 "-mp-filename=" (basename mp-filename ".scm")
                 ".log"))

;; (cog-logger-set-timestamp! #f)
;; (cog-logger-set-sync! #t)
(cog-logger-set-level! "debug")
(cog-logger-set-filename! log-filename)
;; (ure-logger-set-timestamp! #f)
;; (ure-logger-set-sync! #t)
(ure-logger-set-level! "debug")
(ure-logger-set-filename! log-filename)

;; Load KBs to reason on
(define db-lst (load-kbs (list
; Debby
"kbs/biogrid_gene_gene_3.5.177.scm"
"kbs/GO.scm"
"kbs/Go-Plus.scm"
"kbs/GO_annotation.scm"
"kbs/NCBI2Reactome_PE_Pathway.txt.scm"
"kbs/reactome.scm"
"kbs/smpdb_gene.scm"
;                             "kbs/GO.scm"
;                             "kbs/Go-Plus.scm"
;                             "kbs/GO_annotation.scm"
;                             "kbs/reactome.scm"
;                             "kbs/NCBI2Reactome_PE_Pathway.txt.scm"
;                             "kbs/ChEBI2Reactome_PE_Pathway.txt.scm"
;                             "kbs/smpdb_chebi.scm"
;                             "kbs/smpdb_gene.scm"
;                             "kbs/biogrid_gene_gene_3.5.177.scm"
;                             "kbs/hagr.scm"
;                             "kbs/GTEx_median_tissue_expression.scm"
;                             "kbs/moses-results-992.scm"
                         )
                         #:subsmp ss))

(cog-logger-debug "db-lst:\n~a" db-lst)

;; Load patterns
(define pattern-evaluations
  (load-pattern-evaluations mp-filename))

;; Create inheritance links for each pair and their reverse
(define (inh-from-pair p) (Inheritance (car p) (cadr p)))
(define (rev-inh-from-pair p) (Inheritance (cadr p) (car p)))
(define inhs (map inh-from-pair pattern-evaluations))
(define rev-inhs (map rev-inh-from-pair pattern-evaluations))
(define all-inhs (append inhs rev-inhs))

;; Load pln rules of interest
(pln-load #:rule-base 'empty)
(pln-add-rule-by-name "inheritance-direct-introduction-rule")

;; For each pair infer their inheritances
(cog-logger-debug "main (cog-atomspace) = ~a" (cog-atomspace))
(define all-inhs-tv (map gar (map pln-bc all-inhs)))

(write-atoms-to-file "inheritance-links-tb.scm" all-inhs-tv)

;; Get strengths (useful for looking at their distribution, see
;; histogram-inheritance-strengths.gp)
(define strengths (map cog-mean all-inhs-tv))
