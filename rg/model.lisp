;; (load (compile-file "model.lisp"))
(load "lhstats.dx32fsl")
(defparameter *ntxs* 3)

(defvar *trace* nil)

(defvar *txs* nil)

(defparameter *stopping-threshold* 0.01)
(defparameter *n-for-stats* 100)

(defvar *meta-run-memory-lengths* nil)

(defparameter *selection-functions* "filled in at end of module after fns are defined")

;;; Status is :run :grad :fail :pause
;;; Effect is the "real" effect from which observations are drawn
(defstruct tx status effect n+ n-) 

(defparameter *stopping-thresholds* '(0.1 0.05 0.01))

(defparameter *stat-lengths* '(50))

(defun meta-run (&key (n 100) (initvals (list 25 50 75)))
  (with-open-file
   (o (format nil "results/~a-data.xls" (get-universal-time)) :direction :output)
   (format o "stopthresh	selfn	N	Mean	sd	Min	Max	SDn~%")
   (loop for sdn in *stat-lengths*
	 do
	 (loop for selfn in *selection-functions*
	       do
	       (setf *n-for-stats* sdn)
	       (loop for stopthresh in *stopping-thresholds*
			do (setf *stopping-threshold* stopthresh)
			(setf *meta-run-memory-lengths* nil
			      *ntxs* (length initvals))
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

(defvar *memory* "Experimental result for a given run, as: (tx . outcome)")

(defun run (selection-fn stopping-fn &key (init-fixed? nil))
  (setf *txs* 
	(loop for txn below *ntxs*
	      as nr = (if init-fixed?
			  (pop init-fixed?)
			(1+ (random 100)) ;; avoid /0
			)
	      collect (make-tx :status :run :effect (/ nr 100.0) :n+ 0 :n- 1)))
  (setf *memory* nil)
  (loop for n from 1 by 1
	until (done? stopping-fn) 
	do (update-stats (funcall selection-fn))
	;(when (zerop (mod n 100)) (print n)) 
	)
  (length *memory*)
  )

;;; This needs to estimate the strength of belief, either
;;; statistically, or perhaps by counting stability, otherwise you can
;;; accidentally hit the value, esp. if the target is a small
;;; rational, like 0.75 = 3/4 -- which you could accidentally hit in 7
;;; tries!

(defun done? (stopping-fn)
  (loop for tx in *txs*
	unless (funcall stopping-fn tx (tx-n+ tx) (tx-n- tx))
	do (return nil)
	finally (return t)))

(defun stop?close-enough (tx n+ n-)
  "Stop when the difference between predicted and actual is less than the stopping threshold"
  (let ((sd (abs (- (/ n+ (+ n+ n-)) (tx-effect tx)))))
    (and (> n- 0)
	 (> sd 0)
	 (< sd *stopping-threshold*))))

#|

This is highly confused! It won't work! The sd isn't defined over a
slowly changing set of counts! (or ratios that result from those bcs
there really just one slowly changing ratio, not a set of real-valued
samples!

There's a subtlty re having non-zero observation counts. On the one
hand we want to avoid divides by zero, whereas there really could be
no observations of success in a really bad drug. For the moment we
just assume that successes are the demoninator, and that we assume
fail if that's zero. But this isn't quite right, and has to be thought
through.

|#

(defun stop?low-stdev (tx n+ n-)
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
  (random *ntxs*))

(defun select-non-completed ()
  (loop as b = (select-randomly) ;; avoid head bias
	when (not (stop?low-stdev b (gethash b *tx->observed-rgs*)))
	do (return b)))


(defun select-largest-sd () ;; Implies non-completed!
  (let ((selection 
	 ;; Unless they all have a valid SD, reduce to non-completion
	 (if (loop for b below *ntxs*
		   unless (valid-observations? (gethash b *tx->observed-rgs*))
		   do (return nil)
		   finally (return t))
	     ;; Have to exclude completed txs.
	     (loop for b below *ntxs*
		   with maxsd = nil
		   with maxtx = nil
		   as rgs = (gethash b *tx->observed-rgs*)
		   unless (stop?low-stdev b rgs)
		   do (let ((thissd (rgsd rgs)))
			(if (or (null maxsd) (> thissd maxsd))
			    (setf maxtx b maxsd thissd)))
		   finally (return maxtx))
	   (select-non-completed)
	   )))
    (when *trace*
      (loop for b below *ntxs*
	    when (valid-observations? (gethash b *tx->observed-rgs*))
	    do (print (cons b (rgsd (gethash b *tx->observed-rgs*)))))
      (print selection)
      )
    selection))

(defun select-smallest-sd () ;; Implies non-completed!
  "This is an anti-pattern, here to test what happens when folks choose the choices with the best evidence.
Note that here we do NOT care whether the tx has reached min threshold, but just keep piling on! (see ***). 
Problem is that then it'll NEVER end, so we have a second order threshold which is x/<delta-factor> of the target."
  (let ((selection 
	 ;; Unless they all have a valid SD, reduce to non-completion
	 (if (loop for b below *ntxs*
		   unless (valid-observations? (gethash b *tx->observed-rgs*))
		   do (return nil)
		   finally (return t))
	     ;; Have to exclude completed txs.
	     (loop for b below *ntxs*
		   with minsd = nil
		   with mintx = nil
		   with *stopping-threshold* = (/ *stopping-threshold* 10.0)
		   as rgs = (gethash b *tx->observed-rgs*)
		   unless (stop?low-stdev b rgs)
		   do (let ((thissd (rgsd rgs)))
			(if (or (null minsd) (< thissd minsd))
			    (setf mintx b minsd thissd)))
		   finally (return mintx))
	   (select-non-completed)
	   )))
    (when *trace*
      (loop for b below *ntxs*
	    when (valid-observations? (gethash b *tx->observed-rgs*))
	    do (print (cons b (rgsd (gethash b *tx->observed-rgs*)))))
      (print selection)
      )
    selection))

(defun update-stats (selection)
  (let ((r? (< (/ (random 100) 100.0) (gethash selection *tx->real-pr*)))
	(obs (copy-list (car (gethash selection *tx->observed-rgs*))))
	)
    (if r? (incf (car obs)) (incf (cdr obs)))
    (push obs (gethash selection *tx->observed-rgs*))
    (push (list selection r? obs) *memory*)))

(defun dht (table &optional (n 10))
  (maphash #'(lambda (key value)
	       (when (zerop (decf n)) (return-from dht))
	       (format t "~s: ~s~%" key value)	       
	       )
	   table))
	 
(defparameter *selection-functions*
  (list #'select-non-completed #'select-largest-sd #'select-smallest-sd #'select-randomly))

(setf *trace* nil)
(meta-run)
