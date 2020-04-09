;; Simple intensional reasoning test

;; Parameters
(define rs 0)                           ; Random seed
(define ss 1)                         ; Subsampled portion of the KBs
(define mi 10)                      ; Maximum number of iterations
(define cp 1)                           ; Complexity penalty

;; Load modules
(use-modules (opencog))
(use-modules (opencog randgen))
(use-modules (opencog logger))
(use-modules (opencog ure))
(use-modules (opencog pln))
(use-modules (opencog bioscience))
(load "bio-as-utils.scm")

;; Parameters string
(define param-str (string-append
                   "-rs=" (number->string rs)
                   "-ss=" (number->string ss)
                   "-mi=" (number->string mi)
                   "-cp=" (number->string cp)))

(define log-filename
  (string-append "log/manhin-intentional-reasoning.log"))

;; (cog-logger-set-timestamp! #f)
;; (cog-logger-set-sync! #t)
(cog-logger-set-level! "debug")
(cog-logger-set-filename! log-filename)
;; (ure-logger-set-timestamp! #f)
;; (ure-logger-set-sync! #t)
(ure-logger-set-level! "debug")
(ure-logger-set-filename! log-filename)

;; Load kbs
(define db-lst (load-kbs (list
                          "results/manhin-preprocessing.scm"
                         )
                         #:subsmp ss))

(cog-logger-debug "(length db-lst) = ~a" (length db-lst))

;; Load PLN
(pln-load 'empty)

;; For now calculate the intensional inheritance between GO concepts
;; as a test.

; JJJ
(use-modules (ice-9 rdelim))

(define input-file "dw_pairs.txt")
(define result-file "results/intsim.scm")

(define bp (ConceptNode "GO:0042100"))
(define tp (ConceptNode "GO:0042098"))
(define bd (ConceptNode "GO:0030183"))
(define td (ConceptNode "GO:0030217"))

;; ---------- GO pairs ---------- ;;
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
                 #:maximum-iterations mi
                 #:complexity-penalty cp
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

#!
;; ---------- Extensional Similarity Intensional Difference ---------- ;;
(define gos (filter (lambda (x) (string-prefix? "GO:" (cog-name x))) (cog-get-atoms 'ConceptNode)))
(define output-file "new-pairs.txt")
(define output-port (open-output-file output-file))

(pln-add-rule-by-name "intensional-difference-direct-introduction-rule")

(do ((i 0 (+ i 2)))
    ((>= (+ 2 i) (length gos)))
  (let* ((n1 (list-ref gos i))
         (n2 (list-ref gos (+ 1 i)))
         (target (IntensionalDifference n1 n2))
         (results
           (cog-outgoing-set
             (pln-bc
               target
               #:maximum-iterations mi
               #:complexity-penalty cp
             ))))
    (format #t "Done checking: ~a\n" target)
    (if (and (> (cog-mean target) 0) (> (cog-confidence target) 0))
      (display (format #f "~a,~a\n" (cog-name (gar target)) (cog-name (gdr target))) output-port)
    )
  )
)
(close-port output-port)

(let* ((port (open-input-file input-file))
       (line (read-line port)))
  (while (not (eof-object? line))
    (let* ((pairs (string-split line #\,))
           (n1 (Concept (car pairs)))
           (n2 (Concept (cadr pairs)))
           (target (IntensionalDifference n1 n2))
           (results
             (cog-outgoing-set
               (pln-bc
                 target
                 #:maximum-iterations mi
                 #:complexity-penalty cp
               )
             )
           ))
      (write-atoms-to-file result-file results)
      (format #t "Done for pairs: ~a\n" pairs)
    )
    (set! line (read-line port))
  )
  (close-port port)
)

(format #t "Generating members...\n")
(load "../pln/opencog/pln/rules/intensional/intensional-difference-member-introduction.scm")
(define results-2 (cog-outgoing-set (cog-execute! intensional-difference-member-introduction-rule)))
(write-atoms-to-file result-file results-2)

(format #t "Calculating similarities...\n")
(load "../pln/opencog/pln/rules/extensional/extensional-similarity-direct-introduction.scm")
(define results-3 (cog-outgoing-set (cog-execute! extensional-similarity-direct-introduction-rule)))
(write-atoms-to-file result-file results-3)
#!

#!
(define bp-tp (IntensionalDifference bp tp))
(define bd-td (IntensionalDifference bd td))

(define results-1 (pln-bc bp-tp
                          #:maximum-iterations mi
                          #:complexity-penalty cp))
(define results-2 (pln-bc bd-td
                          #:maximum-iterations mi
                          #:complexity-penalty cp))

;(pln-add-rule-by-name "intensional-difference-member-introduction-rule")
;
;(pln-add-rule-by-name "extensional-similarity-direct-introduction-rule")
;
;(define results-3 (pln-bc (ExtensionalSimilarity bp-tp)
;                          #:maximum-iterations mi
;                          #:complexity-penalty cp))
;(define results-4 (pln-bc (ExtensionalSimilarity bd-td)
;                          #:maximum-iterations mi
;                          #:complexity-penalty cp))

;; All results
(define results
  (append
    (cog-outgoing-set results-1)
    (cog-outgoing-set results-2)
    (cog-outgoing-set results-3)
    (cog-outgoing-set results-4)
  )
)

(define scm-filename
  (string-append "results/manhin-intentional-reasoning.scm"))
(write-atoms-to-file scm-filename results)

(define bps  (length (cog-get-root (ConceptNode "GO:0042100"))))
(define tps  (length (cog-get-root (ConceptNode "GO:0042098"))))
(define bds  (length (cog-get-root (ConceptNode "GO:0030183"))))
(define tds  (length (cog-get-root (ConceptNode "GO:0030217"))))

;; Run backward chainer to produce intensional links.
; (define Y (Variable "$Y"))
;; (define target (IntensionalInheritance X Y))
;; (define target (IntensionalSimilarity X Y))
; (define target (IntensionalDifference X Y))
(define results-1 (pln-bc (IntensionalDifference bp tp)
                          #:maximum-iterations mi
                          #:complexity-penalty cp))
(define results-2 (pln-bc (IntensionalDifference bd td)
                          #:maximum-iterations mi
                          #:complexity-penalty cp))
(define results-3 (pln-bc (IntensionalDifference tp bp)
                          #:maximum-iterations mi
                          #:complexity-penalty cp))
(define results-4 (pln-bc (IntensionalDifference td bd)
                          #:maximum-iterations mi
                          #:complexity-penalty cp))
(define rr (list))

(define (cc g1 g2)
;  (define kk (Predicate "count-key"))
  (define result (pln-bc (IntensionalDifference g1 g2)
          #:maximum-iterations mi
          #:complexity-penalty cp))
;  (cog-set-value! g1 kk (cog-incoming-set g1))
;  (cog-set-value! g2 kk (cog-incoming-set g2))
  (format #t "\n===(~a, ~a)\n~a~a=~a\n" (cog-incoming-size g1) (cog-incoming-size g2) g1 g2 (cog-tv (IntensionalDifferenceLink g1 g2)))
  (set! rr (append rr (cog-outgoing-set result)))
)

;; All results
(define results
  (append
    (cog-outgoing-set results-1)
    (cog-outgoing-set results-2)
    (cog-outgoing-set results-3)
    (cog-outgoing-set results-4)
    rr
  )
)

;; Write results in file
(define scm-filename
  (string-append "results/intentional-tb.scm"))

(define gos (filter (lambda (x) (string-prefix? "GO:" (cog-name x))) (cog-get-atoms 'ConceptNode)))
(define threshold 2)

(do ((i 0 (+ i 1)))
;    ((> i 10000))
    ((= i (length gos)))
  (if (> threshold (abs (- bps (length (cog-get-root (list-ref gos i))))))
    (cc bp (list-ref gos i)))

  (if (> threshold (abs (- bds (length (cog-get-root (list-ref gos i))))))
    (cc bd (list-ref gos i)))

  (if (> threshold (abs (- tps (length (cog-get-root (list-ref gos i))))))
    (cc tp (list-ref gos i)))

  (if (> threshold (abs (- tds (length (cog-get-root (list-ref gos i))))))
    (cc td (list-ref gos i)))
)

(write-atoms-to-file scm-filename results)

(format #t "No. of IDs: ~a\nResults:\n~a"
(length (cog-get-atoms 'IntensionalDifferenceLink))
(filter (lambda (x) (not (equal? (stv 1 0) (cog-tv x)))) (cog-get-atoms 'IntensionalDifferenceLink))
)

; JJJ
; (format #t "bp tp: ~a\n" (cog-tv (IntensionalDifferenceLink bp tp)))
; (format #t "bd td: ~a\n" (cog-tv (IntensionalDifferenceLink bd td)))
; (format #t "tp bp: ~a\n" (cog-tv (IntensionalDifferenceLink tp bp)))
; (format #t "td bd: ~a\n" (cog-tv (IntensionalDifferenceLink td bd)))
!#
