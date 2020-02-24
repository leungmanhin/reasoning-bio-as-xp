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
