; (load "hreader.lisp")
(defvar *hypotheses* nil)

(defun init ()
  (with-open-file 
   (i "hypotheses.lisp")
   (setf *hypotheses* (read i))
   ))

(defun test ()
  (init)
  (pprint *hypotheses*)
  (loop for (m* drug-sets) in *hypotheses*
	do (format t "~%
For this mutation set: ~a, use these drugs:~%" m*)
	(loop for ds in drug-sets
	      do (format t "   ~a~%" ds))))

(test)


