;; (load "mets.lisp")
(load "lhstats.dx32fsl")

(defun p* (&key ma mb mab mn wa wb wab wn)
  (let* (
	 (w+ (+ wa wb wab wn))
	 (m+ (+ ma mb mab mn))
	 (d+ (+ ma wa mab wab))
	 (c+ (+ mb wb mab mab))
	 (b+ (+ mab wab))
	 (++ (+ w+ m+))
	 (pw (float (/ w+ ++)))
	 (pm (float (/ m+ ++)))
	 (pd (float (/ d+ ++)))
	 (pc (float (/ c+ ++)))
	 (pawx (float (/ (+ wa wab) w+)))
	 (pwax (float (/ (+ wa wab) d+)))
	 )
    `(
      (:ma ,ma)
      (:mb ,mb)
      (:mab ,mab)
      (:mn ,mn)
      (:wa ,wa)
      (:wb ,wb)
      (:wab ,wab)
      (:wn ,wn)
      (:w+ ,w+)
      (:m+ ,m+)
      (:d+ ,d+)
      (:c+ ,c+)
      (:b+ ,b+)
      (:++ ,++)
      (:pw ,pw)
      (:pm ,pm)
      (:pd ,pd)
      (:pc ,pc)
      (:pawx ,pawx)
      (:pwax ,pwax)
      (:pwa ,(/ (* pawx pw) pd))
      (:paw ,(/ (* pwax pd) pw))
      (:pwax ,pwax)
      (:pawx ,pawx)
      )))

;(defparameter *args* '(:ma 97 :mb 50 :mab 45 :mn 108 :wa 47 :wb 34 :wab 86 :wn 33))
(defparameter *args* '(:ma 91 :mb 89 :mab 92 :mn 88 :wa 92 :wb 88 :wab 90 :wn 89))

(defun run ()
  (let ((**** (apply #'p* *args*)))
    (pprint ****)
    (pprint (apply #'fisher *args*))
    (loop as (arg n) on *args* by #'cddr
	  with *pwax = (second (assoc :pwax ****))
	  with *pawx = (second (assoc :pawx ****))
	  as args+1 = (rplarg arg (+ 1 n) *args*)
	  as res = (apply #'p* args+1)
	  as pwax = (second (assoc :pwax res))
	  as pawx = (second (assoc :pawx res))
	  do (format t "~%~%=======================")
	  (pprint `(,arg :d-pwax ,(- pwax *pwax) :d-pawx ,(- pawx *pawx)))
	  (pprint args+1)
	  (apply #'fisher args+1))))

(defun rplarg (arg n+1 args)
  (cond ((null args) nil)
	((eq (car args) arg) (append (list arg n+1) (cddr args)))
	(t (append (list (car args) (second args))
		   (rplarg arg n+1 (cddr args))))))

(defun fisher (&key ma mb mab mn wa wb wab wn)
  (let ((contingency-table (make-array '(2 4))))
    (setf (aref contingency-table 0 0) ma 
          (aref contingency-table 0 1) mb
          (aref contingency-table 0 2) mab 
          (aref contingency-table 0 3) mn
	  (aref contingency-table 1 0) wa
          (aref contingency-table 1 1) wb
          (aref contingency-table 1 2) wab
          (aref contingency-table 1 3) wn
	  )
    (pprint (stats::chi-square-test-rxc contingency-table))
    ;;(loop for tails in '(:both :positive :negative)
    ;;	  do (pprint `(,tails ,(stats::chi-square-test-rxc contingency-table :tails tails))))
    ))

(run)

