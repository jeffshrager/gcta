;;; (load (compile-file "vtb.lisp"))

;;; To Do:

;;; Creating a range of treatment option types, as might rationally be
;;; created by a tumor board. Output file usually will be sent to
;;; julia analysis code (wanal,jl)

;;; The drug target feature vector (tfn) is the numerical value coding
;;; a an n-bit vector. 

(defstruct tx n name cat tfn)

(defvar *o* nil)
(defparameter *tfgenes* 30)
(defparameter *txs* nil)
(defparameter *cats* '(:imuno :trgtd :chemo))
(defparameter *ntxs-per-cat* 20)

(defun init-txs ()
  ;; Create the treatments
  (loop for cat in *cats*
	with txn = 0
	do (loop for ntx below *ntxs-per-cat*
		 as tx = (make-tx :n (incf txn)
				   :name (random-drug-name cat)
				   :cat cat
				   :tfn (random (expt 2 *tfgenes*)))
		 do (push tx *txs*))))

(defun random-drug-name (cat)
  (let ((proposed-name (propose-random-drug-name cat)))
    (if (not (find proposed-name *txs* :key #'tx-name))
	proposed-name
      (random-drug-name cat))))

(defparameter *drug-name-parts*
  '(("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z")
    ("1" "2" "3" "4" "5" "6" "7" "8" "9" "9")
    ("alfa" "beta" "cato" "dlta" "etho" "fara" "gama" "htro" "iota" "jelo" "keto" "lama" "meta" "neto" "octo" "peta"
     "qrka" "rtro" "sega" "tata" "umma" "velo" "wako" "xeno" "yeti" "zeta")
    ))

(defun propose-random-drug-name (cat)
  (loop for list in *drug-name-parts*
	with name = ""
	do (setf name (format nil "~a~a" name (nth (random (length list)) list)))
	finally (return (format nil "~a~a" name (cdr (assoc cat '((:chemo . "kem") (:trgtd . "tar") (:imuno . "mim"))))))))

;;; Tumor board decision making. TBs always consider the treatments by
;;; category, selects 1-3 plausibles in each cat ranking them, and
;;; then returns the whole list ranked. There are more algorithmically
;;; efficient ways to do this, but we want to do it in a way that
;;; models the way tbs do it.

;;; Tweaking slightly changes the score. Depending upon how strong the
;;; recommendation is we probablistically either up or down regulate
;;; it so that when there is less overlap, there's more variance of
;;; opinion. (This would be better done in Julia with
;;; distribution. After tweaking we twiddle, that is, randomize within
;;; isoscore sets. These represent two different phenomena: tweaking
;;; is uncertainty of efficacy, whereas twiddling i.e., randoming w/i
;;; eq sets represents process randomness.

(defvar *score->txs* (make-hash-table :test #'equal))

(defun tweak-and-twiddle (scored-txs)
  ;; Fun function!! :-)
  (clrhash *score->txs*)
  ;; First Tweak (and load up the twiddle table)...
  (loop for (score . tx) in 
	(loop for (score . tx) in scored-txs
	      ;; This is a sort of UUUgly but simple hack. The prob. that we
	      ;; twiddle the score is based on the overlap where lower
	      ;; scores are more likely to return a 0...
	      collect (cons (if (zerop (random score)) 
				;; And here we either increase it or decrease by 1
				(if (zerop (random 2)) (1+ score) (1- score))
			      score)
			    tx))
	do (push tx (gethash score *score->txs*)))
  ;; Then Twiddle (from the twiddle table):
  (loop for (score . txs) in 
	(sort 
	 (loop for score being the hash-keys of *score->txs*
	       using (hash-value txs)
	       collect (cons score (shuffle-list txs)))
	 #'> :key #'car)
	append (loop for tx in txs collect (cons score tx))))

;;; by Paul Nielse (From http://computer-programming-forum.com/50-lisp/2ff7fc6d728e6fd3.htm)

(defun shuffle-list (list) 
  (let ((new-list (copy-list list)))    ; Skip the copying if destructive 
    (do ((lst new-list (cdr lst))) 
        ((null lst) new-list) 
      (rotatef (car lst) (nth (random (length lst)) lst))))) 

;;; Here we choose up to 3, but stop if there's a skip by more than
;;; one value

(defun top-n-txs (ranked-score.txs)
  (loop for n below 3
	with last-score = (caar ranked-score.txs)
	as ranked-score.tx in ranked-score.txs
	as (score . tx) in ranked-score.txs
	until (> (abs (- score last-score)) 1)
	do (setf last-score score)
	collect ranked-score.tx))

(defun score-txs (ptfn)
  (let ((ptfv (write-to-string ptfn :base 2)))
    (loop for tx in *txs*
	 as ttfv = (write-to-string (tx-tfn tx) :base 2)
	 collect (cons (tfv-overlap ptfv ttfv) tx))))
  
(defun tfv-overlap (ptfv ttfv)
  (loop for pb across ptfv
	as tb across ttfv
	if (string-equal pb tb)
	sum 1))

;;; There are three different ways to report the results. You always
;;; get 5 options, but don't know which algorithm was used to create
;;; the result. (IRL you might actually be told how it was created.)

(defparameter *n-pros* 3) ;; Should be either the same as number of cats, or MxNCATS
(defparameter *n-cons* 2)

(defun create-recommendation (p i ptfn)
  (let* ((scored-txs (tweak-and-twiddle (score-txs ptfn)))
	 (ranked-txs (sort scored-txs #'> :key #'car))
	 (alg (nth (random 3) '(rec-n-pros
				rec-one-pro-per-cat+cons
				rec-combine-as-pros+cons
				)))
	 )
    (let* ((r (funcall (symbol-function alg) ranked-txs)))
      (loop for p/c in '(:pro :con)
	    do
	    (if (eq p/c :pro) (format t "~%tb_~a [~a]:~%" i alg) (format t "   --~%"))
	    (loop for stx in (cdr (assoc p/c r))
		     as score = (car stx)
		     as tx = (cdr stx)
		     as tn from 1 by 1
		     do (format *o* "~a	~a	~a	~a	~a	~a	~a	~a	~a	~a	~a~%"
				p ptfn i p/c tn score (tx-n tx) (tx-cat tx) (tx-tfn tx) (tx-name tx) alg)
		     (format t "   ~a_~a: score:~a :: ~a~%" p/c tn score tx)))
      (list alg r))))

;;; In each case the tbrecs is (score . tx), and the categorized
;;; tbrecs (cat+tbrecs) is ((:cat tx tx tx) ...). You need to be
;;; careful as to whether you want to keep or strip the cats.
	     
(defun rec-n-pros (ranked-txs)
  `((:pro . ,(first-n (+ *n-pros* *n-cons*) ranked-txs))
    (:con . nil)))

(defun rec-combine-as-pros+cons (ranked-txs)
  `((:pro . ,(first-n *n-pros* ranked-txs))
    (:con . ,(worst-n-in-order *n-cons* ranked-txs))))

(defun first-n (n l) (loop for k below n as i in l collect i))

(defun worst-n-in-order (n ranked-txs)
  ;; resort then reverse back get the the tail back in good->poor order
  (reverse (first-n n (sort ranked-txs #'< :key #'car))))

(defun rec-one-pro-per-cat+cons (ranked-txs)
  ;; And the cats are sorted in accord with the score of the top entry
  `((:pro . ,(sort (loop for cat in *cats*
		   collect (loop for tx in ranked-txs
				 when (eq cat (tx-cat (cdr tx)))
				 do (return tx)))
		   #'> :key #'car))
    (:con . ,(worst-n-in-order *n-cons* ranked-txs))))
  
(defun run ()
  (with-open-file
   (*o* (format nil "results/vtb-~a.xls" (get-universal-time)) :direction :output)
   (format *o* "p	ptfn	i	p/c	tn	score	txn	cat	tfn	name	alg~%")
   (init-txs)
   (loop for p below 100 ;; For 100 patients
	 as genotype = (random (expt 2 *tfgenes*))
	 do 
	 (format t "~%================ ~a ~a ================~%" p genotype)
	 (loop for i below 10  ;; Get 10 opinions per patient
	       collect (create-recommendation p i genotype))
	 )))


;;; Simple footrule, based on
;;; http://theory.stanford.edu/~sergei/slides/www10-metrics.pdf, but
;;; takes into account a weight matrix, and the possibility of missing
;;; elements, and non-unique scores. 

(defun swf-distance (l1 l2 &key (weightfn #'swf-default-weightfn)
			&aux (len1 (length l1)) (len2 (length l2)) (missing-score (min len1 len2)) model tail)
  ;; We always use the longer list as the model
  (if (> len2 len1)
      (setq model l2 tail l1)
    (setq model l1 tail l2))
  (loop for me in model
	as mp from 0 by 1
	as weight = (funcall weightfn mp)
	sum (* weight (abs (- mp (or (position me tail :test #'equal) missing-score))))))

(defun swf-default-weightfn (i) (/ 1.0 (1+ i)))

;;; Kendell's W (from: https://github.com/ugolbck/kendall-w/blob/master/kendall_w/kendall_w.py)

    ;; 0 indicates no agreement and 1 indicates unanimous agreement.
    ;; Parameters
    ;; ---------
    ;; data : list
    ;;     List of lists with shape (n_items * n_annotators)
    ;; Return
    ;; ---------
    ;; W : float
    ;;     Kendall's W [0:1]
    ;; Example
    ;; ---------
    ;; annotations = [
    ;;     [1, 1, 1, 2], # item 1
    ;;     [2, 2, 2, 3], # item 2
    ;;     [3, 3, 3, 1], # item 3
    ;; ]
    ;; # Annotator #4 disagrees with the other annotators
    ;; # Annotators #1, #2, #3 agree
    ;; W = kendall_w(annotations)
    ;; # output: 0.4375

;; (defun kw (data)
;;   (let* ((m (length (car data)))  ;; Number of annotators
;; 	 (n (length data)) ;; Number of items
;; 	 (sums (mapcar #'(lambda (l) (apply #'+ l)) data)) ;; Sum of each item ranks
;; 	 (rbar (float (/ (apply #'+ sums) n))) ;; Mean of ranking sums
;; 	 ;; Sum of squared deviations from the mean
;; 	 (s (loop for sum in data
;; 		  as ssum = (apply #'+ sum)
;; 		  sum (expt (- ssum rbar) 2)))
;; 	 )
;;     (/ (* 12 S) (* (expt m 2) (- (expt n 3) n)))))
    
    ;; # Tests
    ;; if not all(len(i) == m for i in data):
    ;;     raise ValueError("Items must all have the same number of annotators.\
    ;;         At least one sublist of argument 'data' has different length than\
    ;;         the first sublist.")
    ;; if m <= 1:
    ;;     raise ValueError("Kendall's W is irrevelent for only one annotator,\
    ;;         try adding more lists to argument 'data'.")
    ;; if m == 2:
    ;;     warnings.warn("Kendall's W is adapted to measure agreement between\
    ;;         more than two annotators. The results might not be reliable in\
    ;;         this case.", Warning)
    ;; # Tests
    ;; if n <= 1:
    ;;     raise ValueError("Kendall's W is irrevelent for only one item,\
    ;;         try adding more sublists to argument 'data'.")



;(print (kw '((1  1  1  2) (2  2  2  3) (3  3  3  1))))

;; (defun wtest ()
;;   (labels ((create-wtest-matrix
;; 	    ()
;; 	    ;; We start with perfect harmony and then start swapping things
;; 	    (let ((m (loop for i below 6 collect (loop for j below 10 collect (1+ i)))))
;; 	      (loop for i below (random 20)
;; 		    as a = (random 6)
;; 		    as b = (random 6)
;; 		    as p = (random 10)
;; 		    as ****************
;; 		    )
;; 	      m)
;; 	    ))
;; 	  (loop for i below 100
;; 		as data = (create-wtest-matrix)
;; 		do (print (list data (kw data)))))
;;   )
;; (wtest)

;(trace get-tbrecs worst-n-in-order)
(run)
