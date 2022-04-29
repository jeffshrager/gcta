;;; (load (compile-file "weaver.lisp"))

#|

The weaver simulates the range of possible cohort reorganizations that
can happen with each new patient. The first patient is obviously an
n-of-1 experiment, but when the next (second) patient arrives, there
is a decision to make as to whether to create a cohort (i.e., give the
same tx as the first pt), or do something separate. And, these become
more complex as there are more patients in the system. To make matters
even more complex, you don't need to have pts actually arrive; At
every time point, pts can fail at tx, whole drugs can fail or
graduate, and other things. All we're doing here is to count the
option space. The theory is that it get wider (of course) at first,
but then narrows down as you start to learn which treatements work for
which subsets of features space.

|#

(defparameter *npts* 10) ;; Both #pts and #cycles -- one pt arrives on each cycle
(defvar *ptids* 0)
(defparameter *ntxs* 5)

(defvar *weave* nil)

(defvar *pts* nil)
(defstruct pt id health) ;; Health is 0-1 dead-ned
(defvar *txs* nil)
(defstruct tx id)
(defvar *cohorts* nil)
(defstruct cohort tx pts)

(defun run ()
  (loop for ptn below *npts*
	as newpt = (make-pt :id (incf *ptids*) :health 0.7)
	
