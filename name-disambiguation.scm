(use-modules (ice-9 rdelim))

(define symbol-alist '())

(define* (download-symbols #:optional
            (symbol-filename
                (string-append
                    "gene_names_"
                    (strftime "%d-%m-%y" (localtime (current-time)))
                    ".txt")))
"
Download all the gene symbols -- Approved symbol, Previous symbols,
and Alias symbols, from genenames.org.
"
    (define base-url "https://www.genenames.org/cgi-bin/download/custom?")
    (define args
        (list "col=gd_app_sym&"
              "col=gd_prev_sym&"
              "col=gd_aliases&"
              "status=Approved&"
              "status=Entry%20Withdrawn&"
              "hgnc_dbtag=on&"
              "order_by=gd_app_sym_sort&"
              "format=text&submit=submit"))

    (system (string-append "curl '" base-url (string-join args "") "' --output " symbol-filename))
)

(define (read-symbols-from-file filename)
    (let ((port (open-input-file filename))
          ; Skip the first line
          (line (read-line port)))
        (set! line (read-line port))
        (while (not (eof-object? line))
            (let ((symbols
                    (delete
                        ""
                        (string-split
                            (string-delete #\space line)
                            (lambda (c) (or (char=? c #\tab) (char=? c #\,)))
                        )
                    )
                 ))
                (for-each
                    (lambda (symbol)
                        (set! symbol-alist
                            (assoc-set! symbol-alist symbol (delete symbol symbols)))
; (format #t "----- Looking at symbol: ~a\nand the values: ~a\n\n" symbol (assoc-ref symbol-alist symbol))
                    )
                    symbols
                )
            )
            (set! line (read-line port))
        )
        (close-port port)
    )
)
