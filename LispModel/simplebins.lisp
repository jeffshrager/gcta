;; (load (compile-file "simplebins.lisp"))

(defparameter *nbins* 3)
(defparameter *nopts/bin* 4)
(defparameter *ncats/opt* 2)
(defparameter *npatients* 100)

(defvar *bins* nil)
(defvar *actual-winner-nth-bin+sum* nil)
;;(defvar *pick-sequences* nil)
(defvar *dead* nil)
(defvar *live* nil)
(defvar *bin-record* nil)

(defun init-bins ()
  (setf *bins*
	(loop for bins below *nbins*
	      collect (loop for n below *nopts/bin*
			    collect (random *ncats/opt*))))
  (setq *actual-winner-nth-bin+sum* 
	(loop for bin in *bins*
	      as n from 1 by 1
	      with max = -1
	      with winner = 0
	      as sum = (reduce #'+ bin)
	      when (> sum max)
	      do (setf winner n max sum)
	      finally (return (cons winner sum))))
  (setf *bin-record* (loop for i below *nbins* collect 0))
  ;;(print (length (setf *pick-sequences* (pick-sequences)))) 
  )

#|

;;; This is old code used to create all possible pick sequences.

(defun pick-sequences ()
  (loop for combo in (all-combinations (* *nbins* *nopts/bin*) *nbins*)
	;; Make sure that you don't over pick from any bin.
	when (not (loop for n-1 below *nbins*
			as n = (1+ n-1)
			if (> (count n combo) *nopts/bin*)
			do (return t)))
	collect combo))

(defun all-combinations (len lim)
  (mapcan #'identity (all-combos2 len lim)))

(defun all-combos2 (len lim)
  (cond ((= len 0) (list (list nil)))
	(t (loop for i from 1 to lim
		 collect (insert-all i (all-combos2 (1- len) lim))))))

(defun insert-all (what in)
  (loop for elt in in
	append (loop for subelt in elt
		     collect (cons what subelt))))

|#

(defun run (nruns)
  (with-open-file 
   (o (format nil "results/~a.xls" (get-universal-time))
      :direction :output)
   (loop for run below nruns
	 do 
	 (setf *live* 0 *dead* 0)
	 (print run)
	 (init-bins)
	 (loop for p below *npatients*
	       as which-bin = (choose-bin)
	       as choice = (nth (random *nopts/bin*) (nth which-bin *bins*))
	       do 
	       (if (zerop choice)
		   (record which-bin :success)
		 (record which-bin :fail)))
	 (format o "~a	~a~%" *live* *dead*)
	 )))

(defun choose-bin ()
  (loop for i below *nbins*
	with max = -1
	with maxn = 0
	as v in *bin-record*
	when (> v max)
	do (setf max v maxn i)
	finally (return maxn)))

(defun record (which-bin s/f)
  (case s/f
	(:success 
	 (incf *live*)
	 (incf (nth which-bin *bin-record*)))
	(:fail
	 (incf *dead*)
	 (decf (nth which-bin *bin-record*)))
	))

(run 1000)
