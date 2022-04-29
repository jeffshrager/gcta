;; (load (compile-file "rg.lisp"))
(load "lhstats.dx32fsl")
(defparameter *nbuckets* 3)

(defvar *trace* nil)

(defvar *bucket->real-pr* (make-hash-table :test #'equal))
(defvar *bucket->observed-rgs* (make-hash-table :test #'equal))

(defvar *memory* "Car is latest experimental result, as: (bucket-number . color)")

(defparameter *stopping-threshold* 0.01)
(defparameter *n-for-stats* 100)

(defvar *meta-run-memory-lengths* nil)

(defparameter *selection-functions* "filled in at end of module after fns are defined")

(defparameter *stopping-thresholds* '(0.1 0.05 0.01))

(defparameter *sd-lengths* '(40 20 10))

(defun meta-run (&key (n 100) (initvals (list 25 50 75)))
  (with-open-file
   (o (format nil "~a-data.xls" (get-universal-time)) :direction :output)
   (format o "stopthresh	selfn	N	Mean	sd	Min	Max	SDn~%")
   (loop for sdn in *sd-lengths*
	 do
	 (loop for selfn in *selection-functions*
	       do
	       (setf *n-for-stats* sdn)
	       (loop for stopthresh in *stopping-thresholds*
			do (setf *stopping-threshold* stopthresh)
			(setf *meta-run-memory-lengths* nil
			      *nbuckets* (length initvals))
			(print (list selfn stopthresh))
			(loop for r below n
			      do (push (run selfn #'stop?low-stdev :init-fixed? initvals) *meta-run-memory-lengths*))
			(format o "~a	~a	~a	~a	~a	~a	~a	~a~%"
				stopthresh (function-name selfn)
				(length *meta-run-memory-lengths*)
				(float (stats::mean *meta-run-memory-lengths*))
				(stats::sd *meta-run-memory-lengths*)
				(reduce #'min *meta-run-memory-lengths*)
				(reduce #'max *meta-run-memory-lengths*)
				sdn
				))))))

(defun run (selection-fn stopping-fn &key (init-fixed? nil))
  (clrhash *bucket->real-pr*)
  (clrhash *bucket->observed-rgs*)
  (loop for bkt below *nbuckets*
	as nr = (if init-fixed?
		    (pop init-fixed?)
		    (1+ (random 100)) ;; avoid /0
		  )
	do
	(setf (gethash bkt *bucket->real-pr*) (/ nr 100.0)) ;; this is pr
	(setf (gethash bkt *bucket->observed-rgs*) (list (cons 0 1))) ;; Initial hypothesis about /0 err!
	)
  ;(dht *bucket->real-pr*)
  (setf *memory* nil)
  (loop for n from 1 by 1
	until (done? stopping-fn) 
	do (update-stats (funcall selection-fn))
	;(when (zerop (mod n 100)) (print n)) 
	)
  '(loop for b below *nbuckets*
	as rgs = (gethash b *bucket->observed-rgs*)
	as clean = (first-n *n-for-stats* (mapcar #'(lambda (rg) (float (/ (car rg) (+ (car rg) (cdr rg)))))
						  (delete t rgs  :test #'(lambda (ignore obs) (= 0 (cdr obs))))))
	do (format t "bucket ~a: n=~a, real=~a, obs: last=~a, mean(~a)=~a, sd(~a)=~a~%"
		   b
		   (length rgs)
		   (gethash b *bucket->real-pr*)
		   (float (/ (caar rgs) (+ (caar rgs) (cdar rgs))))
		   *n-for-stats*
		   (stats::mean clean)		   
		   *n-for-stats*
		   (stats::sd clean)		   
		   ))
  (length *memory*)
  )

;;; This needs to estimate the strength of belief, either
;;; statistically, or perhaps by counting stability, otherwise you can
;;; accidentally hit the value, esp. if the target is a small
;;; rational, like 0.75 = 3/4 -- which you could accidentally hit in 7
;;; tries!

(defun done? (stopping-fn)
  (loop for bucket being the hash-keys of *bucket->observed-rgs*
	using (hash-value rgs)
	unless (funcall stopping-fn bucket rgs)
	do (return nil)
	finally (return t)))

(defun stop?close-enough (bucket rgs)
  "Stop when the difference between predicted and actual is less than the stopping threshold"
  (let ((sd (abs (- (/ (caar rgs) (+ (caar rgs) (cdar rgs))) (gethash bucket *bucket->real-pr*)))))
    (and (> (cdar rgs) 0)
	 (> sd 0)
	 (< sd *stopping-threshold*))))

#|

There's a subtlty re having non-zero observation counts. On the one
hand we want to avoid divides by zero, whereas there really could be
no observations of success in a really bad drug. For the moment we
just assume that successes are the demoninator, and that we assume
fail if that's zero. But this isn't quite right, and has to be thought
through.

|#

(defun stop?low-stdev (bucket rgs)
  "Stop when the standard-deviation of the whole observation set is less than the stopping threshold. 
   We only bother where there are more than 5 observations, and only calculate the sd on the last *n-for-stats* values."
  (setq rgs (first-n *n-for-stats* rgs))
  (and (valid-observations? rgs)
       (< (rgsd rgs) *stopping-threshold*)))

(defun valid-observations? (rgs)
  (and (> (length rgs) 5)
       (> (cdar rgs) 0) ;; at least one observation of success (avoids /0)
       ))

(defun rgsd (rgs)
  "Gets the sd of the first *n-for-stats* of a list of pairs: (a . b) as fractions."
  (stats::sd (mapcar #'(lambda (rg) (float (/ (car rg) (+ (car rg) (cdr rg)))))
		     (delete t (first-n *n-for-stats* rgs)  :test #'(lambda (ignore obs) (= 0 (cdr obs)))))))

(defun first-n (n l)
  (loop for i below n
	as elt in l
	collect elt))

(defun select-randomly ()
  (random *nbuckets*))

(defun select-non-completed ()
  (loop as b = (select-randomly) ;; avoid head bias
	when (not (stop?low-stdev b (gethash b *bucket->observed-rgs*)))
	do (return b)))


(defun select-largest-sd () ;; Implies non-completed!
  (let ((selection 
	 ;; Unless they all have a valid SD, reduce to non-completion
	 (if (loop for b below *nbuckets*
		   unless (valid-observations? (gethash b *bucket->observed-rgs*))
		   do (return nil)
		   finally (return t))
	     ;; Have to exclude completed buckets.
	     (loop for b below *nbuckets*
		   with maxsd = nil
		   with maxbucket = nil
		   as rgs = (gethash b *bucket->observed-rgs*)
		   unless (stop?low-stdev b rgs)
		   do (let ((thissd (rgsd rgs)))
			(if (or (null maxsd) (> thissd maxsd))
			    (setf maxbucket b maxsd thissd)))
		   finally (return maxbucket))
	   (select-non-completed)
	   )))
    (when *trace*
      (loop for b below *nbuckets*
	    when (valid-observations? (gethash b *bucket->observed-rgs*))
	    do (print (cons b (rgsd (gethash b *bucket->observed-rgs*)))))
      (print selection)
      )
    selection))

(defun update-stats (selection)
  (let ((r? (< (/ (random 100) 100.0) (gethash selection *bucket->real-pr*)))
	(obs (copy-list (car (gethash selection *bucket->observed-rgs*))))
	)
    (if r? (incf (car obs)) (incf (cdr obs)))
    (push obs (gethash selection *bucket->observed-rgs*))
    (push (list selection r? obs) *memory*)))

(defun dht (table &optional (n 10))
  (maphash #'(lambda (key value)
	       (when (zerop (decf n)) (return-from dht))
	       (format t "~s: ~s~%" key value)	       
	       )
	   table))
	 
(defparameter *selection-functions*
  (list #'select-non-completed #'select-largest-sd #'select-randomly))

(setf *trace* nil)
(meta-run)
