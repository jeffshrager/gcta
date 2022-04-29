;;; (load "reshape.lisp")

;;; Just a hack to reshuffle results

(defun string-split (string &key (delimiter #\space))
  (let ((substrings '())
        (length (length string))
        (last 0))
    (flet ((add-substring (i)
             (push (subseq string last i) substrings)))
      (dotimes (i length)
        (when (eq (char string i) delimiter)
          (add-substring i)
          (setq last (1+ i))))
      (add-substring length)
      (nreverse substrings))))

(defun nthpart (n f)
  (nth n (string-split f :delimiter #\_)))

(defvar *t* (make-hash-table :test #'equal))
(defvar *keys* nil)

(setf *print-pretty* nil *print-length* nil)

(defun reshape (pattern outfile column maxlines)
  (clrhash *t*)
  (setf *keys* nil)
  (loop for fn in (mapcar #'pathname-name (directory pattern))
	as l = (nthpart 5 fn)
	as pws = (nthpart 6 fn)
	as rule = (nthpart 7 fn)
	as key = (list l pws rule)
	do (with-open-file
	    (i (format nil "results/~a.xls" fn))
	    (push key *keys*)
	    (setf (gethash key *t*)
		  (loop for line = (string-split (read-line i nil) :delimiter #\tab)
			until (null (nth (1- column) line))
			collect line))))

(with-open-file
   ;; Can't call it idperm because then it'll get swept by the above directory search
   (o outfile :direction :output :if-exists :supersede)
   (loop for key in *keys*
	 do (format o "~a	" key))
   (format o "~%")
   (loop for i below maxlines ;; UUU
	 do (loop for key in *keys*
		  do (format o "~a	" (second (nth i (gethash key *t*)))))
	 (format o "~%")))
  )

(reshape "results/*idperm*.xls" "results/allidps.xls" 2 1500)
(reshape "results/*mvtbcdistest*.xls" "results/allcdistest.xls" 3 10000)




