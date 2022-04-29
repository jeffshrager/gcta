;; (load (compile-file "pairsci.lisp"))

(defstruct agent id x y memory)

;;; The width and height of the maze. Both must be odd to preserve the
;;; edges.

;;; (Thought: The maze metaphor could be overkill. Why can't we just
;;; be looking for, say, (0,0) starting at a random point in space,
;;; and making random moves? The actual path someone else took is what
;;; gives the search it's MCMC-like nature, but is thisreally
;;; relevant? Probably only if we're going to follow others' actual
;;; paths somehow, but really all we need to do is to chase them to
;;; the front of their search and get to the goal ahead of
;;; them. Soccer Science is played on an open field, not a maze. (See
;;; also: Payette N (2009) Simulating Science: A Multiagent Model of
;;; Scientific Evolution. Paper presented at the conference on
;;; Modeling Science: Understanding, forecasting, and communicating
;;; the science system. N. Payette2009Simulating Science: A Multiagent
;;; Model of Scientific Evolution. Paper presented at the conference
;;; on Modeling Science: Understanding, forecasting, and communicating
;;; the science system.Amsterdam (October 6–9, 2009). Amsterdam
;;; (October 6–9, 2009)))

(defparameter *width* 17)
(defparameter *height* 9)

(defvar maze)
(defvar *goal*)
(defvar *agents* nil)

;;; BAR, Bias Against Replication, is used to decide when to retrace
;;; either one's own (I-BAR), or anyone's (G-BAR) steps. The BAR value
;;; is used in the expression: (zerop (random <bar value>)), and if
;;; this is true (i.e., one in n times it'll be zero), then we allow
;;; replication, otherwise we reject the move. See MOVE-AGENT for the
;;; detailed code.

(defvar *i-bar*)
(defparameter *i-bar-range* '(1 1 1)) ; use n n 1 to just do n 

;; The global memory tells us where others have already tread
(defvar *global-memory* nil)

(defvar *g-bar*)
(defparameter *g-bar-range* '(1 5 1)) ; use n n 1 to just do n

(defvar *nagents*)
(defparameter *nagents-range* '(1 10 2)) ; use n n 1 to just do n

(defparameter *nruns* 1000)

(defvar *runs* nil)

(defparameter *display-maze* nil) ; :final :every-move <number> -- any other value never displays

(defun run ()
  (format t "~%nagents	i-bar	g-bar	mean~%") 
  (loop for *nagents* from (first *nagents-range*) to (second *nagents-range*) by (or (third *nagents-range*) 1)
	do (loop for *i-bar* from (first *i-bar-range*) to (second *i-bar-range*) by (or (third *i-bar-range*) 1)
		 do (loop for *g-bar* from (first *g-bar-range*) to (second *g-bar-range*) by (or (third *g-bar-range*) 1)
			  do (setq *runs* nil)
			  (loop for i below *nruns*
				do (push (run-maze) *runs*)
				(when (equal :final *display-maze*) (display-maze))
				)
			  (format t "~a	~a	~a	~a~%" *nagents* *i-bar* *g-bar* (/ (reduce #'+ *runs*) (float (length *runs*))))
			  ))))

(defparameter *agent-id-base* 64)

(defun generate-agents ()
  (setq *agents* 
	(loop for id from 1 to *nagents*
	      collect (make-agent :id id
				  :x 1 :y 0))))

(defun run-maze ()
  (generate-maze)
  (generate-agents)
  (setf *global-memory* nil)
  (loop with i = 0
	do (loop for a in *agents*
	      do (move-agent a)
	      (incf i)
	      ;; Return from the run anytime any agent reaches the goal
	      (when (equal :every-move *display-maze*) (display-maze))
	      (when (and (numberp *display-maze*) (zerop (mod i *display-maze*))) (display-maze))
	      (if (and (= (agent-x a) (car *goal*))
		       (= (agent-y a) (cdr *goal*)))
		  (return-from run-maze i)))))

(defun move-agent (a)
  (let* ((x1 (agent-x a))
	 (y1 (agent-y a)))
    (loop as x2 = (+ x1 (- (random 3) 1))
	  as y2 = (+ y1 (- (random 3) 1))
	  when (and (and (> x2 0) (< x2 *width*)) ;; Dont run off the maze
		    (and (> y2 0) (< y2 *height*))
		    (= (aref maze x2 y2) 0) ;; Or into walls! BAR= bias
		    ;; against replication, either globally (*g-bar*)
		    ;; or locally (individually: *i-bar*)
		    (if (member (cons x2 y2) (agent-memory a) :test #'equal)
			(zerop (random *i-bar*))
		      (if (member (cons x2 y2) *global-memory* :test #'equal)
			  (zerop (random *g-bar*))
			t))
		    )
	  do (return (progn (setf (agent-x a) x2 (agent-y a) y2)
			    (pushnew (cons x2 y2) (agent-memory a) :test #'equal)
			    (pushnew (cons x2 y2) *global-memory* :test #'equal))))))
	  
; The maze builder is based on one by Joe Wingbermuehle (2003).

(defun carve-maze (x y)
  (let ((d (random 4)))
    (dotimes (c 4)
      (let* ((cd (mod (+ c d) 4))
	     (dv (cond
		  ((= cd 0) (list 1 0))
		  ((= cd 1) (list 0 1))
		  ((= cd 2) (list -1 0))
		  (t        (list 0 -1))))
	     (x1 (+ x (car dv)))
	     (y1 (+ y (cadr dv)))
	     (x2 (+ x1 (car dv)))
	     (y2 (+ y1 (cadr dv)))
	     )
	(if (and (and (> x2 0) (< x2 *width*))
		 (and (> y2 0) (< y2 *height*)))
	    (if (and (= (aref maze x1 y1) 1)
		     (= (aref maze x2 y2) 1))
		(let ()
		  (setf (aref maze x1 y1) 0)
		  (setf (aref maze x2 y2) 0)
		  (carve-maze x2 y2)
                  )))))))

(defun generate-maze ()
   (setq *random-state* (make-random-state t))
   (setf (aref maze 1 1) 0)
   (carve-maze 1 1)
   (setf (aref maze 1 0) 0)
   (setf (aref maze (- *width* 1) (- *height* 2)) 0)
   (setf *goal* (cons (- *width* 1) (- *height* 2)))
   )

(defun display-maze ()
  (dotimes (y *height*)
    (dotimes (x *width*)
      (if (= (aref maze x y) 1)
	  (princ "[]")
	(if (equal (cons x y) *goal*)
	    (princ "**")	       
	  (loop for a in *agents*
		if (and (= (agent-x a) x)
			(= (agent-y a) y))
		do (progn (pac a) (pac a) (return nil))
		else if (member (cons x y) (agent-memory a) :test #'equal)
		do (progn (pac a) (princ #\!) (return nil))
		finally (princ "  ")))))
    (terpri)
    ))

(defun pac (a)
  (princ (code-char (+ *agent-id-base* (agent-id a)))))

(run)