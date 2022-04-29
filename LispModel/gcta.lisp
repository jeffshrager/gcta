; (load (compile-file "gcta.lisp"))

;;; III Variants have less and less effect when they are drugged for longer?

(eval-when 
 (compile load)
 (unless (probe-file "lhstats.dx32fsl")
   (compile-file "lhstats.lisp"))
 (unless (find-package 'STATISTICS)
   (load "lhstats.dx32fsl"))
 (unless (probe-file "normal.dx32fsl")
   (compile-file "normal.lisp"))
 (unless (ignore-errors (symbol-function 'normal-random-integer))
   (load "normal.dx32fsl"))
 )

(defvar *verbose?* t) ;; verbose flag

(defun formatv? (format &rest args)
  (when *verbose?*
    (apply #'format `(t ,format ,@args))))

;;; =========================================================
;;; The literature (aka. the log)

(defvar *log* nil)

(defstruct log-entry
  year event pid did health)

;;; =========================================================
;;; PATIENTS AND DISEASES (VARIANTS)

;;; Disease's are caused by non-normal values in the patient's variant
;;; vector. Patietns also have a "health" (effectively "tumor load" or
;;; "disease state"), which is incremented each year in accord with
;;; the variants and their weights. 

(defstruct patient
  id
  variants
  age
  regimen
  health ;; This gets recomputed every cycle (year)
  )

(defparameter *variant-vector-length* 100) 
(defparameter *variant-normal-value* 0)
(defparameter *variant-abnormal-values* '(1))
(defparameter *n-abvals* (length *variant-abnormal-values*))
(defparameter *starting-age* 50)
(defparameter *%-abnormal-variant* 10) ;; 0-100 this is on a PER VECTOR ELEMENT BASIS, *NOT* over the whole vector

;;; Variant weights: Most of the variants as non-drivers, but some are
;;; drivers. On a given run, which variants are drivers is precomputed
;;; and fixed. So the main goal of the discovery problem is to
;;; efficiently find out which are the drivers, and then efficiently
;;; find out which drugs are best to use, usually the ones that
;;; preferentially target the drivers. Weight is between 0 and 1.

(defvar *variant-weights* nil)

(defparameter *mean-variant-weight* -0.6) ;; Should be negative in order for variants to decrease health
(defparameter *variant-weight-sdev* 0.3)

(defun set-variant-weights ()
  (setq *variant-weights*
	(scaled-normal :n *variant-vector-length* :m *mean-variant-weight* :d *variant-weight-sdev*)))

(defvar *latest-assigned-patient-id* 0)
(defvar *latest-assigned-drug-id* 0)

(defun new-patient ()
  (make-patient :id (incf *latest-assigned-patient-id*)
		:health 0.0 ;; ??? Maybe make this a random variable?
		:age *starting-age*
		:variants
		(loop for i below *variant-vector-length*
		      if (< (random 100) *%-abnormal-variant*)
		      collect (nth (random *n-abvals*) *variant-abnormal-values*)
		      else collect *variant-normal-value*)))

;;; =========================================================
;;; DRUGS AND TARGETS

;;; Drugs have one or more targets, which are specific locations in
;;; the variant space, and both a positive and negative effect
;;; size. The positive effect size improves the patient's state of
;;; health in accord with the number of matches between drug and
;;; actual variant. The negative effect size make the patient worse,
;;; regardless of biomarker matches. These model effective treatment
;;; v. side effects. A drug cocktail simply adds these (!!! FFF adding
;;; is probably not really the correct function)

(defstruct drug 
  id 
  targets 
  +effect-per-match
  -effect-per-target
  )

;;; With this and the variant set so very low, the probabiliy of a
;;; drug being helpful is going to be very low (specifically, it
;;; should be the product of this number and the % variant density)

(defparameter *%-targets-per-drug* 5) 
(defparameter *general-cancer-effect-per-variant* -0.5)
(defparameter *mean-+effect-per-match* 2.0)
(defparameter *mean--effect-per-target* -2.0)
(defparameter *effect-size-devs* 0.1) ;; Applies to both (??? separate these?)

(defun new-drug ()
  (make-drug 
   :id (incf *latest-assigned-drug-id*)
   :targets (loop for i below *variant-vector-length*
		  if (< (random 100) *%-targets-per-drug*)
		  collect (nth (random *n-abvals*) *variant-abnormal-values*)
		  else collect *variant-normal-value*)
   :+effect-per-match (first (scaled-normal :n 1 :m *mean-+effect-per-match* :d *effect-size-devs*))
   :-effect-per-target (first (scaled-normal :n 1 :m *mean--effect-per-target* :d *effect-size-devs*))
   ))

(defun new-personalized-drug (p)
  (make-drug 
   :id (incf *latest-assigned-drug-id*)
   :targets (patient-variants p)
   :+effect-per-match (first (scaled-normal :n 1 :m *mean-+effect-per-match* :d *effect-size-devs*))
   :-effect-per-target (first (scaled-normal :n 1 :m *mean--effect-per-target* :d *effect-size-devs*))
   ))

;;; Health is a value between -1 and +1. Everyone starts at 0 and then
;;; their health improves or declines in accord with a scaled sigmoid
;;; function, as pushed around the number of properly targeted
;;; variants (CTVs), missed variants (MVs), and drug targets
;;; (DTs). The CTVs improve health whereas the MVs and DTs make health
;;; worse.

(defun cCTVs=contrib-from-correctly-targeted-variants (d p)
  (* (drug-+effect-per-match d)
     (loop for v in (patient-variants p)
	   as vw in *variant-weights*
	   as g in (drug-targets d)
	   when (and (not (equal v *variant-normal-value*))
		     (not (equal g *variant-normal-value*)))
	   sum vw)))
(defun cMVs=contrib-from-missed-variants (d p)
  (* *general-cancer-effect-per-variant*
     (loop for v in (patient-variants p)
	   as g in (drug-targets d)
	   when (and (not (equal v *variant-normal-value*))
		     (equal g *variant-normal-value*))
	   sum 1)))
(defun cDTs=contrib-from-drug-targets (d p) ;; to figure side-effects, which are per-target
  (* (drug--effect-per-target d)
     (loop for g in (drug-targets d)
	   when (not (equal g *variant-normal-value*))
	   sum 1)))

;;; The general health modulator either adds or subtracts points
;;; depeneding on your age. Before this age you get improvement points
;;; (i.e., negative counts) depending on how young you are, but after
;;; it you start declining. So younger folks actually get slowly
;;; better, all things being equal, but older one get slowly worse.
(defparameter *ghm-zero-age* 70)

(defun cGHM-contrib-from-General-Health-Modulator (p)
  (- (/ (- (patient-age p) *ghm-zero-age*) 10.0)))

(defun total-health-delta-per-drug (d p)
  ;; Another way to think of the first two terms is just the number of
  ;; variants - the number correctly targeted, and the last term (side
  ;; effects) is the cost of targeting overall.
  (let* ((cCTVs (cCTVs=contrib-from-correctly-targeted-variants d p))
	 (cMVs (cMVs=contrib-from-missed-variants d p))
	 (cDTs (cDTs=contrib-from-drug-targets d p))
	 (result (+ cCTVs cMVs cDTs))
	 )
    (formatv? "~%For patient ~a (cGHM=~a) & drug ~a: cCTVs = ~a, cMVs = ~a, cDTs = ~a --> ~a~%" 
	    (patient-id p) (cGHM-contrib-from-General-Health-Modulator p) (drug-id d) cCTVs cMVs cDTs result)
    result))

;;; And now we use that total number to change the patient's health
;;; through this master equation, and then renormalize to -1/+1. The
;;; way this works is that for every one "hit", you get x% closer to
;;; whichever end. So, zero would leave you where you are (except for
;;; the general health modulator (see above), negative counts would
;;; put you closer to -1 (unhealthy) and positive would put you closer
;;; to +1 (healthy), but you can't every reach the end points. There
;;; is a death point near -1 (but obviously no super health point). 

(defun update-patient-health (p)
  (let ((h (patient-health p))
	(delta (/ (min 99 ;; farthest we're allowed to go!
		       (truncate 
			(loop for d in (patient-regimen p)
			      with h = (patient-health p)
			      with delta = (cGHM-contrib-from-General-Health-Modulator p)
			      ;; WWW ??? Variants are going to be
			      ;; multiply counted * n drugs !!!
			      ;; Should end up with the summary by
			      ;; either creating a pseudo drug, or
			      ;; gathgering the unmatch variants and
			      ;; adding them once.
			      do (when (< h *health-threshold*) 
				   (incf delta (total-health-delta-per-drug d p)))
			      finally (return delta))))
		  100.0)))
    (setf (patient-health p)
	  (if (not (= 0 delta))
	      (if (> delta 0)
		  (+ h (abs (* delta (- 1.0 h))))
		(- h (abs (* delta (- -1.0 h)))))
	    h))
    (formatv? "For patient ~a, pre-update health = ~a, delta = ~a, new health = ~a~%" 
	    (patient-id p) h delta (patient-health p))))

;;; Okay, so with all the above, we should be able to run a basic
;;; simulation.

(defparameter *n-years-to-run* 50)
(defparameter *n-population* 100)
(defparameter *n-drugs* 100)

;; When a person reaches this health value, he gets kicked out of the population
(defparameter *death-threshold* -0.9) 

;; When a person reaches this health value, he you stop drugging him,
;; or if he goes below, you start again.
(defparameter *health-threshold* -0.25) 

(defvar *p* nil) ;; This is the main global population store
(defvar *d* nil) ;; All drugs
(defvar *h* nil) ;; heavan or hell -- you decide!

(defun sim ()
  ;; Init
  (setf *latest-assigned-patient-id* 0)
  (setq *h* nil) ;; clear out the dead
  (set-variant-weights)
  (formatv? "~%*variant-weights* = ~s~%" *variant-weights*)
  ;; Drugs
  (setq *d* (loop for i below *n-drugs* collect (new-drug)))
  ;; Patients
  (setq *p* 
	(loop for i below *n-population* 
	      as d in *d*
	      collect (let ((p (new-patient)))
			(setf (patient-regimen p) 
			      ;;(list (new-personalized-drug p))
			      (list (new-drug))
			      )
			p)))
  (when *verbose?* (pprint *p*))
  ;; EXECUTION
  (loop for y below *n-years-to-run*
	do 
	(formatv? "~%~%************************* ~a ************************~%~%" y)
	(loop for p in *p*
		 do (update-patient-health p)
		 #+nil (dp p)
		 (incf (patient-age p)))
	(formatv? "****** Health summary: ")
	(mapcar #'(lambda (p) (formatv? "~a " (truncate (* 100 (patient-health p))))) *p*)
	(setf *p* ;; Bring out your dead!
	      (loop for p in *p*
		    as ph = (patient-health p)
		    as le = (car (push (make-log-entry :year y :pid (patient-id p) :health ph :event nil 
											      :did (mapcar #'drug-id (patient-regimen p))) *log*))
		    if (> ph *death-threshold*)
		    collect (progn (setf (log-entry-event le) :live) p)
		    ;; No collection in the else case so patient is
		    ;; dropped from *p* and pushed onto *h*
		    else do (progn (setf (log-entry-event le) :dead) (push p *h*))
		    ))
	(formatv? "~%~%******* Population is now: ~a e-souls.~%~%" (length *p*))
	))

(defun dp (&optional p) ;; describe-patient (or all patients)
  (cond ((null p) (setq p *p*))
	((numberp p) (setq p `(,(find p *p* :key #'patient-id) ,(find p *h* :key #'patient-id))))
	((not (listp p)) (setq p (list p))))
  (loop for p in p
	when p
	do 
	(let ((*print-length* nil)
	      (*print-width* nil)
	      (*PRINT-PRETTY* nil))
	  (pprint p))))

(defun dh () (dp *h*))

;;; Experimental jig

(defun run (&key (reps 5))
  (setq *log* nil)
  ;;  (loop for *mean-variant-weight* from 0.0 to 1.0 by 0.1
  ;;	do (loop for *variant-weight-sdev* from 0.2 to 1.0 by 0.1
  (let ((p* (loop for i below reps do (sim) collect (length *p*))))
    (format t "@ mvw=~a, vwsdev=~a, p*=~a, mean(p*)=~a~%" 
	    *mean-variant-weight* *variant-weight-sdev* p* (float (STATISTICS:MEAN-SD-N p*)))))

(untrace)
;(trace total-health-delta-count cCTVs=contrib-from-correctly-targeted-variants 
; cMVs=contrib-from-missed-variants cDTs=contrib-from-drug-targets nGHM-n-from-General-Health-Modulator)
(setq *verbose?* nil)
(run)
