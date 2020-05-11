(use-modules (opencog) (opencog bioscience) (opencog ure) (opencog pln) (ice-9 rdelim))

(load "bio-as-utils.scm")

; (primitive-load "results/pln-attractions.scm")
(primitive-load "results/pln-step-6.scm")
(define input-file "dw_pairs.txt")
(define result-file "results/pln-intensional-similarities.scm")

;; ---------- GO pairs ---------- ;;
(pln-load 'empty)
(pln-add-rule-by-name "intensional-similarity-direct-introduction-rule")

(let* ((port (open-input-file input-file))
       (line (read-line port)))
  (while (not (eof-object? line))
    (let* ((pairs (string-split line #\,))
           (n1 (Concept (car pairs)))
           (n2 (Concept (cadr pairs)))
           (target (IntensionalSimilarity n1 n2))
           (results
             (cog-outgoing-set
               (pln-bc
                 target
                 #:maximum-iterations 10
               )
             )
           ))
      (write-atoms-to-file result-file results)
      (format #t "Done for pairs: ~a\nResults: ~a\n" pairs results)
    )
    (set! line (read-line port))
  )
  (close-port port)
)

(define fp (open-output-file "pln_results_all.txt"))
(for-each
  (lambda (x)
    (display
      (format #f "~a,~a,~a\n"
        (cog-name (gar x))
        (cog-name (gdr x))
        (if (equal? (stv 1 0) (cog-tv x))
          0
          (cog-mean x)
        )
      )
      fp
    )
  )
  (cog-get-atoms 'IntensionalSimilarityLink)
)
(close-port fp)
