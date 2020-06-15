;; Set of utilities for filtering the dataset and the results

;; Load modules
(use-modules (srfi srfi-1))
(use-modules (opencog))
(use-modules (opencog exec))
(use-modules (opencog ure))
(use-modules (opencog logger))
(use-modules (opencog bioscience))

;; For debugging
;; (cog-logger-set-stdout! #t)
;; (cog-logger-set-sync! #t)

;; Helpers
(define (true-predicate x) #t)

(define (false-predicate x) #f)

(define (null-mean? x)
  (and (cog-atom? x) (< (cog-mean x) 1e-16)))

(define (non-null-mean? x)
  (and (cog-atom? x) (< 1e-16 (cog-mean x))))

;; Whether all nodes of the given link have non null mean
(define (all-nodes-non-null-mean? LINK)
  (let* ((nodes (cog-get-all-nodes LINK)))
    (every non-null-mean? nodes)))

(define (scope? x)
  (cog-subtype? 'ScopeLink (cog-type x)))

(define (lst? x)
  (cog-subtype? 'ListLink (cog-type x)))

(define (and? x)
  (cog-subtype? 'AndLink (cog-type x)))

(define (present? x)
  (cog-subtype? 'PresentLink (cog-type x)))

(define (eval? x)
  (cog-subtype? 'EvaluationLink (cog-type x)))

(define (eval-pred-name? name x)
  (and (eval? x)
       (equal? (cog-name (gar x)) name)))

(define (eval-GO_namespace? x)
  (eval-pred-name? "GO_namespace" x))

(define (member? x)
  (cog-subtype? 'MemberLink (cog-type x)))

(define (subset? x)
  (cog-subtype? 'SubsetLink (cog-type x)))

;; TODO: better use scheme basename function
(define (rm-extension fn ext)
  (if (string-suffix? (string-append "." ext) fn)
      (substring fn 0 (- (string-length fn) 4))
      fn))

(define (string-starts-with? str prefix)
  (if (< (string-length str) (string-length prefix))
    #f
    (let* ((start 0)
           (end (string-length prefix))
           (str-prefix (substring str start end)))
      (equal? str-prefix prefix))))

(define (smp? A)
  (and (eq? (cog-type A) 'ConceptNode)
       (string-starts-with? (cog-name A) "SMP0")))

(define (go? A)
  (and (eq? (cog-type A) 'ConceptNode)
       (string-starts-with? (cog-name A) "GO:")))

(define (gene? A)
  (and (eq? (cog-type A) 'GeneNode)))

(define (GO_term? A)
  (and (eq? (cog-type A) 'ConceptNode)
       (equal? (cog-name A) "GO_term")))

(define (inheritance-GO_term? A)
  (and (eq? (cog-type A) 'InheritanceLink)
       (GO_term? (gdr A))))

(define (get-smps)
  (filter smp? (cog-get-atoms 'ConceptNode)))

(define (get-genes)
  (cog-get-atoms 'GeneNode))

(define (get-go-categories)
  (filter go? (cog-get-atoms 'ConceptNode)))

(define (go-subset? S)
  (and (subset? S) (go? (gar S)) (go? (gdr S))))

(define (get-go-subsets)
  (filter go-subset? (cog-get-atoms 'SubsetLink)))

(define (get-members C)
"
  Given a concept node C, return all its members
"
  (let* ((member-links (cog-filter 'MemberLink (cog-incoming-set C)))
         (member-of-C? (lambda (x) (equal? C (gdr x))))
         (members (map gar (filter member-of-C? member-links))))
    members))

(define (get-member-links EL-TYPE C-TYPE)
  (let* ((mbr-links (cog-get-atoms 'MemberLink #f))
         (valid-types? (lambda (x) (and (equal? (cog-type (gar x)) EL-TYPE)
                                        (equal? (cog-type (gdr x)) C-TYPE)))))
    (filter valid-types? mbr-links)))

(define (get-cardinality C)
"
  Giveb a concept node C, return its number of members
"
  (length (get-members C)))

(define (number->hexstr n)
  (format #f "~x" n))

(define (bool->string b)
  (if b "t" "f"))

(define (get-pattern eval-pattern)
  (cog-outgoing-atom (gdr eval-pattern) 0))

(define (get-body pattern)
  (cog-outgoing-atom pattern 1))

(define (get-clauses body)
  (if (or (eq? (cog-type body) 'PresentLink)
          (eq? (cog-type body) 'AndLink))
      (cog-outgoing-set body)
      body))

(define (mk-rand-selector prob)
"
  Return a random selector that is true with the given probability
"
  (lambda (x) (<= (cog-randgen-randfloat) prob)))

(define (load-filter-in pred-in? filename)
"
  1. Load filename in an auxiliaury atomspace
  2. Grab all atoms
  3. Only retain the valid ones according to pred-in?
  4. Copy the valid atoms in the current atomspace
  5. Return the list of the copied atoms

  An auxiliary atomspace is used to avoid deleting atoms, which can be
  quite costly.
"
  (let* (;; Load file in a temporary atomspace
         (base-as (cog-set-atomspace! (cog-new-atomspace)))
         (dummy (primitive-load filename))

         ;; Filter in atoms satisfying pred
         (atoms (filter pred-in? (cog-get-atoms 'Atom #t)))

         ;; Copy admissible atoms in the base atomspace
         (base-atoms (cog-cp base-as atoms))

         ;; Discard the temporary atomspace
         (dummy (cog-set-atomspace! base-as)))
    base-atoms))

(define* (load-kb kb-filename
                  #:key
                  (subsmp 1)
                  (filter-in true-predicate)
		  (filter-out false-predicate))
"
  Load knowledge base, optionally perform some filtering and return
  a list of all atoms (not just root atoms) loaded.

  Note that filters only affect root atoms as far as their
  belonging to the atomspace is concerned, however they affect
  all atoms as far as their belonging to the resulting list is
  concerned. TODO: maybe we want to return only roots, though
  that might be undesirable if the list is passed as db to the
  pattern miner and we want to mine subgraphs.

  Usage: (load-kb kb-filename
                  #:subsmp ssp
                  #:filter-in pred-in?
                  #:filter-out pred-out?)

  kb-filename: Scheme file containing the knowledge base.

  ssp [optional, default=1]: Probability of randomly retaining
  atoms.

  pred-in? [optional, default=true-predicate]: Predicate that
  atoms must satisfy in order to be retained.

  pred-out? [optional, default=false-predicate]: Predicate that
  atoms must contradict in order to be retained.
"
  (let* (;; Define filter for admissible atoms
         (rand-selected? (mk-rand-selector subsmp))
         (admissible? (lambda (x) (and
				    (rand-selected? x)
				    (filter-in x)
				    (not (filter-out x))))))
    (load-filter-in admissible? kb-filename)))

(define* (load-kbs kbs-filenames
                   #:key
                   (subsmp 1)
		   (filter-in true-predicate)
                   (filter-out false-predicate))
"
  Like load-kb but takes a list of filenames.

  Note that filters only affect root atoms as far as their
  belonging to the atomspace is concerned, however they affect
  all atoms as far as their belonging to the resulting list is
  concerned. TODO: maybe we want to return only roots, though
  that might be undesirable if the list is passed as db to the
  pattern miner and we want to mine subgraphs.

  Usage: (load-kbs kb-filenames
                   #:subsmp ssp
                   #:filter-in pred-in?
                   #:filter-out pred-out?)

  kb-filenames: Scheme list of files containing knowledge bases.

  ssp [optional, default=1]: Probability of randomly retaining atoms.

  pred-in? [optional, default=true-predicate]: Predicate that atoms
  must satisfy in order to be retained.

  pred-out? [optional, default=false-predicate]: Predicate that atoms
  must contradict in order to be retained.
"
  (concatenate (map (lambda (x) (load-kb x
					 #:subsmp subsmp
					 #:filter-in filter-in
                                         #:filter-out filter-out))
                    kbs-filenames)))

(define (add-extra-smp-go-terms)
  ;; Small Molecule Pathway concept
  (let* ((smps (get-smps))
         (smp-cpt (Concept "SMP_term"))
         (gos (get-go-categories))
         (go-cpt (Concept "GO_term")))
    (append (map (lambda (x) (Inheritance x smp-cpt)) smps)
            (map (lambda (x) (Inheritance x go-cpt)) gos))))

(define (write-atoms-to-file filename a-lst)
"
  Write a list of atoms into filename.
"
  (let* ((port (open-file filename "w")))
    (for-each (lambda (x) (write x port)) a-lst)
    (close-port port))

  *unspecified*)

(define (load-pattern-evaluations filename)
"
  1. Load the given scheme file.
  2. Only retain pattern evaluations with
    2.1 non null surprisingness
    2.2 pairs of SMP and GO concept and gene variabale member
"
  (let* (;; Define filter for admissible atoms
         (surp-eval? (lambda (x) (or (eval-pred-name? "nisurp" x)
                                     (eval-pred-name? "jsdsurp" x))))
         (member-clauses? (lambda (x)
                            (let* ((body (get-body (get-pattern x)))
                                   (clauses (get-clauses body))
                                   (are-members (map member? clauses)))
                              (not (member #f are-members)))))
         (pos-mean? (lambda (x) (< 0 (cog-mean x))))
         (admissible? (lambda (x) (and
                                   (pos-mean? x)
                                   (surp-eval? x)
                                   (member-clauses? x)))))
    (load-filter-in admissible? filename)))

(define (extract-GO-SMP-pair pattern-eval)
  (let* ((body (get-body (get-pattern pattern-eval)))
         (clauses (get-clauses body))
         (cpts (map gdr clauses)))
    cpts))

(define (true-subset-inverse S)
"
  Given a subset with a true value

  Subset (stv 1 1)
    A <ATV>
    B <BTV>

  Return

  Subset <TV>
    B <BTV>
    A <ATV>

  where TV is calculated as follows

  TV.strength = (ATV.strength * ATV.count) / (BTV.strength * BTV.count)
  TV.count = (BTV.strength * BTV.count)

  Which is technically correct since (Subset A B) is true.
"
(let* ((A (gar S))
       (B (gdr S))
       (ATV (cog-tv A))
       (BTV (cog-tv B))
       (A-positive-count (* (cog-tv-mean ATV) (cog-tv-count ATV)))
       (B-positive-count (* (cog-tv-mean BTV) (cog-tv-count BTV)))
       (TV-strength (if (< 0 B-positive-count)
                        (exact->inexact (/ A-positive-count B-positive-count))
                        1))
       (TV-count B-positive-count)
       (TV-confidence (count->confidence TV-count))
       (TV (stv TV-strength TV-confidence)))
  (Subset TV B A)))

(define (gt-zero-confidence? A)
"
  Return #t iff A's confidence is greater than 0
"
  (> (cog-confidence A) 0))

(define (gt-zero-mean? A)
"
  Return #t iff A's mean is greater than 0
"
  (> (cog-mean A) 0))

(define (gt-zero-mean-and-confidence? A)
  (and (gt-zero-confidence? A) (gt-zero-mean? A)))

;; Trace querying utilities
;;
;; TODO: move to URE

;; Print the inference steps leading to this target

(define (get-direct-steps-to-target target)
"
  Return all inference steps directly inferring the given target, in
  the following format:

  (Set
    (List <rule-1> <source-1> <iteration-1>)
    ...
    (List <rule-n> <source-n> <iteration-n>))
"
  (let* ((pattern (Execution
                    (Variable "$rule")
                    (List
                      (Variable "$source")
                      (Variable "$iteration"))
                    target))
         (vardecl (VariableList
                    (TypedVariable (Variable "$rule") (Type 'DefinedSchemaNode))
                    (Variable "$source")
                    (TypedVariable (Variable "$iteration") (Type 'NumberNode))))
         (gl (Get vardecl pattern)))
    (cog-execute! gl)))

(define (get-direct-steps-from-source source)
"
  Return all inference steps directly inferred from the give source, in
  the following format:

  (Set
    (List <rule-1> <target-1> <iteration-1> )
    ...
    (List <rule-n> <target-n> <iteration-n>))
"
  (let* ((pattern (Execution
                    (Variable "$rule")
                    (List
                      source
                      (Variable "$iteration"))
                    (Variable "$target")))
         (vardecl (VariableList
                    (TypedVariable (Variable "$rule") (Type 'DefinedSchemaNode))
                    (Variable "$target")
                    (TypedVariable (Variable "$iteration") (Type 'NumberNode))))
         (gl (Get vardecl pattern)))
    (cog-execute! gl)))

(define (get-trails-to-target-rec target . inners)
"
  Return all inference trails leading to the given target, in the
  following format:

  (Set
    (List
      (List <rule-11> <inter-11> <iteration-11>)
      ...
      (List <rule-1m> <inter-1m> <iteration-1m>))
    ...
    (List
      (List <rule-n1> <inter-n1> <iteration-n1>)
      ...
      (List <rule-nm> <inter-nm> <iteration-nm>)))
"
  (let* ((get-inner (lambda (s) (gdr s))) ; Get the inner target of a step
         (direct-steps (get-direct-steps-to-target target))
         ;; Remove cycles
         (inners? (lambda (s) (member (get-inner s) inners)))
         (not-inners? (lambda (s) (not (inners? s))))
         (direct-steps-no-cycles (filter not-inners? (cog-outgoing-set direct-steps)))
         ;; Given a direct inference step, find the trails going to
         ;; that inference step, and append the inference step to them
         (get-trails (lambda (s)
                       (let* ((inrs (if (inners? s) inners (cons (get-inner s) inners))))
                         (cog-outgoing-set (apply get-trails-to-target-rec (cons (get-inner s) inrs))))))
         (append-step-to-trail (lambda (t s)
                                 (List (cog-outgoing-set t) s)))
         (append-step-to-trails (lambda (ts s)
                                  (if (null? ts)
                                      (List s)
                                      (map (lambda (t) (append-step-to-trail t s)) ts))))
         (get-trails-with-direct-step (lambda (s)
                                        (let* ((ts (get-trails s)))
                                          (append-step-to-trails ts s)))))
    (Set (map get-trails-with-direct-step direct-steps-no-cycles))))

(define (get-trails-to-target target)
  (get-trails-to-target-rec target target))

(define (get-trails-from-source-rec source . inners)
"
  Return all inference trails coming from the given source, in the
  following format:

  (Set
    (List
      (List <rule-11> <inter-11> <iteration-11>)
      ...
      (List <rule-1m> <inter-1m> <iteration-1m>))
    ...
    (List
      (List <rule-n1> <inter-n1> <iteration-n1>)
      ...
      (List <rule-nm> <inter-nm> <iteration-nm>)))
"
  (let* ((get-inner (lambda (s) (gdr s))) ; Get the inner target of a step
         (direct-steps (get-direct-steps-from-source source))
         ;; Remove cycles
         (inners? (lambda (s) (member (get-inner s) inners)))
         (not-inners? (lambda (s) (not (inners? s))))
         (direct-steps-no-cycles (filter not-inners? (cog-outgoing-set direct-steps)))
         ;; Given a direct inference step, find the trails going to
         ;; that inference step, and append the inference step to them
         (get-trails (lambda (s)
                       (let* ((inrs (if (inners? s) inners (cons (get-inner s) inners))))
                         (cog-outgoing-set (apply get-trails-from-source-rec (cons (get-inner s) inrs))))))
         (prepend-step-to-trail (lambda (t s) (List s (cog-outgoing-set t))))
         (prepend-step-to-trails (lambda (ts s)
                                   (if (null? ts)
                                       (List s)
                                       (map (lambda (t) (prepend-step-to-trail t s)) ts))))
         (get-trails-with-direct-step (lambda (s)
                                        (let* ((ts (get-trails s)))
                                          (prepend-step-to-trails ts s)))))
    (Set (map get-trails-with-direct-step direct-steps-no-cycles))))

(define (get-trails-from-source source)
  (get-trails-from-source-rec source source))
