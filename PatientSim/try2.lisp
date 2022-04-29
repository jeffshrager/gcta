;; (load (compile-file "try2.lisp"))

(defun run (&key (n 1000))
  (with-open-file
   (o "~/Downloads/try2.xls" :direction :output :if-exists :supersede)
   (loop for x below n
	 do (format o "~A~%" (+ (sin (/ 100.0)) (sin (/ x 300)) (sin (/ x 10.0))))
	 )))

(run)

