;;; (load "pwys.lisp")

(setf *model-graph* 
  '((size + cln3)
    (cln3 + mbf)
    (cln3 + sbf)
    (sbf + cln12)
    (cln12 - sic1)
    (cln12 - cdh1)
    (sic1 - clb56)
    (sic1 - clb12)
    (mbf + clb56)
    (clb56 + mcm1/sff) ; weird yellow/green in the paper
    (clb56 - sic1)
    (clb56 + clb12) ; weird yellow/green in the paper
    (clb56 - cdh1)
    (cdh1 - clb12)
    (clb12 - cdh1)
    (clb12 + mcm1/sff)
    (clb12 + cdc20&14)
    (clb12 - sic1)
    (clb12 - swi5)
    (clb12 - mbf)
    (clb12 - sbf)
    (mcm1/sff + cdc20&14)
    (mcm1/sff + swi5)
    (mcm1/sff + clb12)
    (cdc20&14 + swi5)
    (cdc20&14 + cdh1)
    (cdc20&14 - clb12)
    (cdc20&14 + sic1)
    (cdc20&14 - clb56)
    (swi5 + sic1)
    ))

(defun graph-proteins (graph)
  (remove-duplicates 
    (loop for (from +- to) in graph
	  append (list from to))))

(defun proteins-directly-controlling (desired-to-protein graph)
  (remove-duplicates 
    (loop for (from-protein +- to-protein) in graph
           WHEN (EQUAL to-protein desired-to-protein)
           collect from-protein)))

(defun graph-walker (to-prot graph &optional seen)
  (let ((direct-controllers (proteins-directly-controlling to-prot graph)))
    (cond ((null direct-controllers) (list to-prot))
          (t (loop for prot in direct-controllers
                    unless (member prot seen)
                    collect (cons to-prot 
                                  (graph-walker prot graph 
                                               (cons prot seen))))))))

(defvar *pwys* nil)

(defun explode-pathways (gwr) ;; graph-walker-result
  (setq *pwys* nil)
  (mapcar #'(lambda (sub) (ep2 sub nil)) gwr)
  *pwys*)

(defun ep2 (gwr collector)
  (cond ((and (= 2 (length gwr)) (symbolp (first gwr)) (symbolp (second gwr)))
	 (pushnew (append (reverse collector) gwr) *pwys* :test #'tree-equal))
	(t (mapcar #'(lambda (sub) (ep2 sub (cons (car gwr) collector))) (cdr gwr)))))

(defvar *gene->pathways* (make-hash-table :test #'equal))

(defun collect-genes->pathways (pwys)
  (clrhash *gene->pathways*)
  (loop for pwy in pwys
	do (loop for gene in pwy
		 do (push pwy (gethash gene *gene->pathways*)))))


(collect-genes->pathways (explode-pathways (graph-walker 'sic1 *model-graph*)))

(loop for gene being the hash-keys of *gene->pathways*
      using (hash-value pwys)
      do (print (list gene (length pwys))))


