;;; (load (compile-file "life.lisp"))
(declaim (optimize (debug 3)))
(defparameter *l* 71) ;; Want there to be a central element (for no really good reason)
(defparameter *c* (1- (round (/ *l* 2))))
(defvar *body* nil)

(defvar *results.xls* nil)

(defstruct asite
  var ;; an integer -- 0 indicates that there's no tumor here
  nut-count ;; counts the number of UNfilled cells around the current one
  next-var ;; The type that this will become on update
  )

(defparameter *nseeds* 4)

(defparameter *seeds* '(nil a b c d e f g h i j k l)) ;; Just a convenience
(defparameter *chars* "-abcdefghijk")

(defun init (&key (size 11))
  (setq *l* size)
  (setf *body* (make-array (list *l* *l*)))
  (loop for x below *l*
	do (loop for y below *l*
		 do (setf (aref *body* x y)
			  (make-asite :var nil :nut-count 0))))
  ;; Seed a new tumor in the center, of var type 1
  (loop for n from 1 to *nseeds*
	as seed in *seeds*
	do (setf (asite-var (aref *body* (+ 2 (random (- *l* 4))) (+ 2 (random (- *l* 4))))) seed))
  )

;;; Run forward one time unit

(defun tic ()
  ;; Figure the nut status of every tumor. First clear all the
  ;; nut counts. (This could be done in one pass...but f'it)
  (loop for x below *l*
	do (loop for y below *l*
		 do (setf (asite-nut-count (aref *body* x y)) 0)))
  ;; Now walk through the whole array again and add one to each cell
  ;; around each (protect edge effects) -- but only when there's not other cell there (i.e., it's 0)
  (loop for x from 1 to (- *l* 2)
	do (loop for y from 1 to (- *l* 2)
		 do
		 (when (null (asite-var (aref *body* x y)))
		   (incf (asite-nut-count (aref *body* (1- x) (1- y)))) ;; Above left
		   (incf (asite-nut-count (aref *body* (1- x) y))) ;; Above
		   (incf (asite-nut-count (aref *body* (1- x) (1+ y)))) ;; Above right
		   (incf (asite-nut-count (aref *body* x (1- y)))) ;; left
		   (incf (asite-nut-count (aref *body* x (1+ y)))) ;; right
		   (incf (asite-nut-count (aref *body* (1+ x) (1- y)))) ;; Below left
		   (incf (asite-nut-count (aref *body* (1+ x) y))) ;; Below
		   (incf (asite-nut-count (aref *body* (1+ x) (1+ y)))) ;; Below right
		   )))
  ;; Now spread the tumor to nearby cells when there is a lot of
  ;; nut. For a given cell, the neighboring tumor with the
  ;; highest nut support wins (uniform random on ties).
  (spread)
  )

(defun spread ()
  ;; Compute next vars
  (loop for x from 1 to (- *l* 2)
	do (loop for y from 1 to (- *l* 2)
		 do (setf (asite-next-var (aref *body* x y))
			  (next-var x y))))
  ;; Now update (unless it's going to update to 0 -- this should be filtered above)
  (loop for x from 1 to (- *l* 2)
	do (loop for y from 1 to (- *l* 2)
		 as next = (asite-next-var (aref *body* x y))
		 do (setf (asite-var (aref *body* x y)) next))))

(defun next-var (x y)
  ;; Select a cell to replce based on total nut count, probabalistically.
  (let* ((vars (surrounding-vars x y))
	 (sum (reduce #'+ (mapcar #'cdr vars))))
    ;; Norm the nut counts
    (loop for elt in vars
	  do (setf (cdr elt)
		   (truncate (* 100 (float (/ (cdr elt) sum))))))
    ;; Now select one probabalistically
    (loop with rand = (random 100)
	  with sum = 0
	  for (var . threshold) in vars
	  do (incf sum threshold)
	  when (>= sum rand)
	  do (return var))))
	  
;;; Returns (var . nut-count) for all the surrounding cells. 

(defun surrounding-vars (x y)
  (summarize-vars 
   (macrolet
    ((svget (x y)
	    `(let ((asite (aref *body* ,x ,y)))
	       (cons (asite-var asite) (asite-nut-count asite)))))
    (loop for elt in 
	  (list (svget (1- x) (1- y)) ;; Above left
		(svget (1- x) y)      ;; Above
		(svget (1- x) (1+ y)) ;; Above right
		(svget x (1- y))      ;; left
		(svget x (1+ y))      ;; right
		(svget (1+ x) (1- y)) ;; Below left
		(svget (1+ x) y)      ;; Below
		(svget (1+ x) (1+ y)) ;; Below right
		)
	  unless (null (car elt))
	  collect elt))))

;;; Takes the output from surrounding-vars, as: ((1 . 2) (1 . 1) (1
;;; . 3) (1 . 1) (0 . 1) (1 . 0) (0 . 1) (1 . 0)) and combines the nut
;;; counts for each variant, as: ((5 . 1) (2 . 0)) [Note that we're only
;;; counting the CARs of the above data!]

(defun summarize-vars (vars)
  (loop for elt in (loop with r = nil
			 as (var . nut) in vars
			 as elt = (find var r :key #'car)
			 if elt do (incf (cdr elt) nut)
			 else do (push (cons var nut) r)
			 finally (return r))
	unless (zerop (cdr elt))
	collect elt))

;;; Utils

(defun show (&key verbose &aux r)
  (when verbose (format t "~%(variant/nutcount)~%   "))
  (when verbose (loop for y below *l* do (format t " ~a  " y)))
  (format t "~%")
  (loop for x below *l*
	do
	(when verbose (format t "~%~s: " x))
	(loop for y below *l*
	      as asite = (aref *body* x y)
	      as var = (asite-var asite)
	      as char = (aref *chars* (position var *seeds*))
	      do (if (assoc var r) (incf (cdr (assoc var r))) (push (cons var 1) r))
	      (if verbose (format t "~a/~a " char (asite-nut-count asite))
		(format t "~c " char))
	      )
	(format t "~%")
	)
  (format t "~%Totals: ~a~%" r)
  (loop for (nil . n) in r
	do (format *results.xls* "~a	" n))
  (format *results.xls* "~%")
  )

(defun test (&key verbose tics size)
  (with-open-file
   (*results.xls* (format nil "results/~a.xls" (get-universal-time))
		  :direction :output)
   (init :size size)
   (loop as i from 1 to tics
	 do
	 (tic)
	 (show :verbose verbose))))

(test :verbose nil :tics 100 :size 51)
