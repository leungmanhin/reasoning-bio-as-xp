;; Parameters
(define jb 8)                           ; Number of jobs for the task
(define rs 0)                           ; Random seed
; JJJ
(define ss 0)                         ; Subsampled portion of the KBs
(define ms 100)                           ; Minimum support (ignored if mf is positive)
(define mf 0.001)                          ; Minimum frequency (ignored if negative)
(define mi 1000)                          ; Maximum number of iterations
(define mc 4)                           ; Maximum number of conjuncts
(define mv 3)                           ; Maximum number of variables
(define su 'nisurp)                       ; Surprisingness measure

;; Load modules & utils
(use-modules (srfi srfi-1))
(use-modules (opencog randgen))
(use-modules (opencog logger))
(use-modules (opencog ure))
(use-modules (opencog miner))
(use-modules (opencog bioscience))
(load "bio-as-utils.scm")

;; Set random seed
(cog-randgen-set-seed! rs)

;; Set loggers
(define log-filename (string-append "log/manhin-mining.log"))

(cog-logger-set-level! "debug")
(cog-logger-set-filename! log-filename)

(ure-logger-set-level! "debug")
(ure-logger-set-filename! log-filename)

;; Load preprocessed KBs, get the list of trees to mine
(define db-lst (load-kbs (list
                           ;; Debby
                           ; "kbs/biogrid_gene_gene_3.5.177.scm"
                           ; "kbs/GO.scm"
                           ; "kbs/Go-Plus.scm"
                           ; "kbs/GO_annotation.scm"
                           ; "kbs/NCBI2Reactome_PE_Pathway.txt.scm"
                           ; "kbs/reactome.scm"
                           ; "kbs/smpdb_gene.scm"
                           "kbs/GO.scm"
                           "kbs/Go-Plus.scm"
                           "kbs/GO_annotation.scm"
                           "kbs/reactome.scm"
                           "kbs/NCBI2Reactome_PE_Pathway.txt.scm"
                           "kbs/ChEBI2Reactome_PE_Pathway.txt.scm"
                           "kbs/smpdb_chebi.scm"
                           "kbs/smpdb_gene.scm"
                           "kbs/biogrid_gene_gene_3.5.177.scm"
                           "kbs/hagr.scm"
                           "kbs/GTEx_median_tissue_expression.scm"
                           "kbs/moses-results-992.scm"
                           "kbs/string_ggi.scm"
                           "kbs/string_ppi.scm"
                         )
                         #:subsmp ss))

;; Post-process by adding extra knowledge
(define db-lst (append db-lst (add-extra-kb)))

;; Call pattern miner
(define results (cog-mine db-lst
                          #:jobs jb
                          #:minimum-support ms
                          #:minimum-frequency mf
                          #:maximum-iterations mi
                          #:conjunction-expansion #f
                          #:maximum-conjuncts mc
                          #:maximum-variables mv
                          #:maximum-spcial-conjuncts 4
                          #:surprisingness su))

(cog-logger-debug "Final results:\n~a" results)

;; Write results in a file
(define miner-results-filename "results/manhin-mining.scm")
(write-atoms-to-file miner-results-filename results)
