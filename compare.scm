(use-modules (opencog) (ice-9 rdelim) (srfi srfi-1))

(define go-pairs (list))
(define pln-results (list))
(define dw-results (list))
(define match-pairs (list))
(define mismatch-pairs (list))

(call-with-input-file "pln_results_all.txt"
  (lambda (fp)
    (let ((line (read-line fp)))
      (while (not (eof-object? line))
        (let* ((contents (string-split line #\,))
               (go1 (list-ref contents 0))
               (go2 (list-ref contents 1))
               (go-pair (string-join (list go1 go2) ","))
               (result (list-ref contents 2)))
          (set! go-pairs (append go-pairs (list (cons go1 go2))))
          (set! pln-results (assoc-set! pln-results go-pair result))
        )
        (set! line (read-line fp))
      )
    )
  )
)

(call-with-input-file "dw_results_all.txt"
  (lambda (fp)
    (let ((line (read-line fp)))
      (while (not (eof-object? line))
        (let* ((contents (string-split line #\,))
               (go1 (list-ref contents 0))
               (go2 (list-ref contents 1))
               (go-pair (string-join (list go1 go2) ","))
               (result (list-ref contents 2)))
;          (set! go-pairs (append go-pairs (list (cons go1 go2))))
          (set! dw-results (assoc-set! dw-results go-pair result))
        )
        (set! line (read-line fp))
      )
    )
  )
)

(set! go-pairs (delete-duplicates go-pairs))
(format #t "Total no. of pairs: ~a\n" (length go-pairs))

(define (get-results margin size-margin)
  (define match-cnt 0)
  (define mismatch-cnt 0)

  (define size-match-cnt 0)
  (define size-mismatch-cnt 0)

  (set! match-pairs (list))
  (set! mismatch-pairs (list))

  (do ((i 0 (1+ i)))
      ((= i (length go-pairs)))
    (let* ((pair (list-ref go-pairs i))
           (go1 (car pair))
           (go2 (cdr pair))
           (k1 (format #f "~a,~a" go1 go2))
           (k2 (format #f "~a,~a" go2 go1))
           (pln-r1 (assoc-ref pln-results k1))
           (pln-r2 (assoc-ref pln-results k2))
           (dw-r1 (assoc-ref dw-results k1))
           (dw-r2 (assoc-ref dw-results k2))
           (pln-result (if pln-r1 pln-r1 pln-r2))
           (dw-result (if dw-r1 dw-r1 dw-r2)))
      (if (<= (abs (- (string->number pln-result) (string->number dw-result))) margin)
        (begin
          (set! match-cnt (1+ match-cnt))
          (set! match-pairs (append match-pairs (list pair))))
        (begin
          (set! mismatch-cnt (1+ mismatch-cnt))
          (set! mismatch-pairs (append mismatch-pairs (list pair))))
      )
    )
  )

  (format #t "Match: ~a (~,2f%)\nMismatch: ~a\n"
    match-cnt
    (exact->inexact (/ match-cnt (length go-pairs)))
    mismatch-cnt
  )

  (do ((i 0 (1+ i)))
      ((= i (length match-pairs)))
    (let* ((pair (list-ref match-pairs i))
           (go1 (Concept (car pair)))
           (go2 (Concept (cdr pair)))
           (go1-size (length (filter (lambda (x) (equal? (gdr x) go1)) (cog-incoming-by-type go1 'MemberLink))))
           (go2-size (length (filter (lambda (x) (equal? (gdr x) go2)) (cog-incoming-by-type go2 'MemberLink))))
           (size-diff (abs (- go1-size go2-size)))
           (size-avg (/ (+ go1-size go2-size) 2)))
;      (if (<= (abs (- go1-size go2-size)) size-margin)
      (if (or (= size-avg 0) (<= (/ size-diff size-avg) size-margin))
        (set! size-match-cnt (1+ size-match-cnt))
        (set! size-mismatch-cnt (1+ size-mismatch-cnt))
      )
    )
  )

  (format #t "For match-pairs:\nSize match: ~a (~,2f%)\nSize mismatch: ~a\n"
    size-match-cnt
    (exact->inexact (* (/ size-match-cnt (length match-pairs)) 100))
    size-mismatch-cnt
  )

  (set! size-match-cnt 0)
  (set! size-mismatch-cnt 0)

  (do ((i 0 (1+ i)))
      ((= i (length mismatch-pairs)))
    (let* ((pair (list-ref mismatch-pairs i))
           (go1 (Concept (car pair)))
           (go2 (Concept (cdr pair)))
           (go1-size (length (cog-incoming-by-type go1 'MemberLink)))
           (go2-size (length (cog-incoming-by-type go2 'MemberLink)))
           (size-diff (abs (- go1-size go2-size)))
           (size-avg (/ (+ go1-size go2-size) 2)))
;      (if (<= (abs (- go1-size go2-size)) size-margin)
      (if (or (= size-avg 0) (<= (/ size-diff size-avg) size-margin))
        (set! size-match-cnt (1+ size-match-cnt))
        (set! size-mismatch-cnt (1+ size-mismatch-cnt))
      )
    )
  )

  (format #t "For mismatch-pairs:\nSize match: ~a (~,2f%)\nSize mismatch: ~a\n"
    size-match-cnt
    (exact->inexact (* (/ size-mismatch-cnt (length mismatch-pairs)) 100))
    size-mismatch-cnt
  )
)

(display (get-results 0.1 0.25)) (newline)
(display (get-results 0.2 0.25)) (newline)
(display (get-results 0.3 0.25)) (newline)
