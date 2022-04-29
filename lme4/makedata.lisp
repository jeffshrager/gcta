;;; (load (compile-file "makedata.lisp"))

(defun make-data (fn &key (nsubjects 100) (minobs 3) (maxobs 100)
			(pobs 0.1)
			(a_trend 0.314159) ;(/ (rno 3 10) 10.0))
			(noise_range 0.3)
			)
  (with-open-file
   (o (format nil "~a_s~a_n~a.tsv" fn nsubjects (truncate (* 100 noise_range)))
	      :direction :output :if-exists :supersede)
   (format o "SID	Tx	Day	KPS~%")
   (loop for sid below nsubjects
	 as kps = (rno 50 80)
	 as tx = (nth (random 2) '(:a :b))
	 ;with a_trend = (+ (/ (random 100) 33.0) a_trend)
	 with b_trend = (- a_trend)
	 do (loop for day below (+ minobs (random maxobs))
		  do
		  (when (< (/ (random 100) 100.0) pobs)
		    (format o "s~a	~a	~a	~a~%" sid tx day kps))
		  ;; Update kps and die as needed
		  (setf kps
			(+ (* (nth (random 2) '(1.0 -1.0))
			      (* (/ (random 100) 10.0))
			      noise_range)
			   (+ kps (case tx (:a a_trend) (t b_trend)))
			   ))
		  (if (> kps 100.0) (setf kps 100.0))
		  (if (< kps 0.0) (setf kps 0.0))
		  ))))

(defun rno (low high)
  (+ low (/ (random (* 100 (abs (- low high)))) 100.0)))

;(make-data (get-universal-time))
