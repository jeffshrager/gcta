;;; +=========================================================================+
;;; | Copyright (c) 2009-2010 Jeff Shrager (jshrager@stanford.edu)            |
;;; |                                                                         |
;;; | Permission is hereby granted, free of charge, to any person obtaining   |
;;; | a copy of this software and associated documentation files (the         |
;;; | "Software"), to deal in the Software without restriction, including     |
;;; | without limitation the rights to use, copy, modify, merge, publish,     |
;;; | distribute, sublicense, and/or sell copies of the Software, and to      |
;;; | permit persons to whom the Software is furnished to do so, subject to   |
;;; | the following conditions:                                               |
;;; |                                                                         |
;;; | The above copyright notice and this permission notice shall be included |
;;; | in all copies or substantial portions of the Software.                  |
;;; |                                                                         |
;;; | THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         |
;;; | EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      |
;;; | MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  |
;;; | IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    |
;;; | CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    |
;;; | TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       |
;;; | SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  |
;;; +=========================================================================+

;;; To run: 
;;; 1. Create a subdirectory called "results"
;;; 2. Start any Common Lisp. (Tested on Allegro Common Lisp (v 8.1) by Franz Inc., 
;;;    and CCL by Clozure Associates.)
;;; 3. Enter: (load (compile-file "WHATEVER_YOU_CALLED_THIS_FILE")) to the lisp reader.
;;;    (Obviously, change WHATEVER_YOU_CALLED_THIS_FILE to whatever you called the file
;;;     when you downloaded it. The quotes ("...") are required in the above command, 
;;;     so it should look something like: (load (compile-file "WHATEVER_YOU_CALLED_THIS_FILE"))
;;; The program should self-execute upon load, producing a number of outputs in 
;;;    the results subdirectory. If it does not self-execute, you can start it manually
;;;    by entering: (plos) to the lisp reader.

(defparameter *run-label* nil)

(defparameter *params* 
  '(
    *disease-name*
    *disease-p-of-diagnosis-per-person-per-year*
    *disease-initial-level-at-diagnosis*
    *disease-progression-multiplier-per-year*
    *disease-trials*
    *disease-death-threshold*
    *most-effective-drug-power*
    *treatment-effectiveness-distribution-model*
    *initial-compounds-that-address-a-given-disease*
    *discovery-model*
    *additional-non-functional-drugs*
    *treatment-effectiveness-mean*
    *stochastic-treatment-effectiveness-sd*
    *do-not-repeat-treatments* 
    *years-to-run*
    *n-people*
    *run-replicates*
    *treatment-model* 
    *initial-treatment-model-if-need-be*
    *disease-state-measurement-error-sd*
    *measures-are-precise*
    *publication-regimen*
    *how-many-sequential-declines-counts-as-lengthy-experience*
    *how-many-back-pubs-to-look-at*
    *n-clinical-trials-at-once*
    *n-patients-per-trial*
    *trial-length-in-years*
    ))

;;; General parameters:
(defparameter *years-to-run* 100)
(defparameter *n-people* 1000)
(defparameter *run-replicates* 10)
    
;;; Disease specs
(defparameter *disease-name* "cancer")
(defparameter *disease-p-of-diagnosis-per-person-per-year* 0.01)
(defparameter *disease-initial-level-at-diagnosis* 0.1)
(defparameter *disease-progression-multiplier-per-year* 1.2)
(defparameter *disease-trials* nil)
(defparameter *disease-death-threshold* 0.5)

;;; :even :stochastic :golf-course
(defparameter *treatment-effectiveness-distribution-model* :even)
(defparameter *initial-compounds-that-address-a-given-disease* 10)
;;; low for :even model
(defparameter *most-effective-drug-power* 0.6) 
;;; Values < 1 are helpful, otherwise harmful! -- so most drugs have no primary effect
(defparameter *treatment-effectiveness-mean* 0.9)
(defparameter *additional-non-functional-drugs* 0)
;;; Only relevant when randomly distributed drug effect.
(defparameter *stochastic-treatment-effectiveness-sd* 0.2) 

(defparameter *discovery-model* :no-discovery)
;;; [[Code removed for PLoS publication as not relevant.]]

;;; [[Code removed for PLoS publication as not relevant.]]

; Implies do not treat after you run out of possibilities
(defparameter *do-not-repeat-treatments* nil)

;;; The treatment model says how to switch drugs when when there is
;;; failure to improve.  Options are: :random :placebo :no-switch
;;; :switch-to-best-published-treatment :trials :hill-climb (This is
;;; like random by always switch to the next best treatment) The
;;; initial treatment model is only relevant in the no-switch case.
(defparameter *treatment-model* :random) 
(defparameter *initial-treatment-model-if-need-be* nil) 
(defparameter *disease-state-measurement-error-sd* 2.0)
(defparameter *measures-are-precise* nil)
(defparameter *publication-regimen* :publish-only-after-lengthy-experience)
(defparameter *how-many-sequential-declines-counts-as-lengthy-experience* 1)
(defparameter *how-many-back-pubs-to-look-at* 100)
(defparameter *n-clinical-trials-at-once* 3)
(defparameter *n-patients-per-trial* 100)
(defparameter *trial-length-in-years* 5)

;;; The labels here get appended to those given as the first arg to
;;; sim, with "+" between them. This is a conveniece for labeling.
(defvar *augment-labels-stack* nil)

(defun augmented-label (label)
  (loop for aug-label in *augment-labels-stack*
	do (setf label
		 (format nil "~a+~a" label aug-label))
	finally (return label)))

(defvar *master-average-tbl* (make-hash-table :test #'equal))

(defvar *summary-stream* nil) 

(defvar *the-literature* nil)
(defvar *trials* nil)
(defvar *next-birth-id* nil)

;;; [[Code removed for PLoS publication as not relevant.]]

(defstruct disease 
  name 
  p-of-diagnosis-per-person-per-year 
  active-trials 
  completed-trials
  initial-level-at-diagnosis 
  progression-multiplier-per-year 
  death-threshold 
  drugs 
  trials)

(defstruct person 
  id 
  scientist? 
  disease 
  disease-state 
  drug 
  disease-state-history
  ;; [[Code removed for PLoS publication as not relevant.]]
  drugs-tried ;; This is slightly redundant for efficiency (could get it from the clinical-record)
  treatment-history ;; inverse chron order list of clinical-record
  out-of-options-do-not-treat ;; Gets set when, for example, you have no-repeat and are out of drugs.
  )

(defstruct clinical-record
  treatment
  date-started
  date-ended)

(defstruct drug 
  id 
  effect ;; 1.0 is no effect, smaller numbers are better
  next-best-treatment ;; used for hill-climbing treatment model
  )

;;; [[Code removed for PLoS publication as not relevant.]]

(defmacro sim (label &rest letargs)
  `(let ,letargs
     (setq *run-label* (augmented-label ,label)) ; Global so that summarize can get it above
     (print ',letargs)
     (run)
     ))

(defmacro test-set (auglabel args &body body)
  `(let ((*augment-labels-stack*
	  (cons ,auglabel *augment-labels-stack*))
	 ,@args)
     ,@body))

(defvar *run-name* "")
(defmacro run-test-set  (name &body body)
  `(with-open-file 
    (*summary-stream* ,(format nil "results/~a-~a-SUMMARY.xls" (get-universal-time) name)
		      :direction :output)
    (setq *run-name* ,name)
    (clrhash *master-average-tbl*)
   ;; ------------ Test code begins here -----------------
    ,@body
    (final-averages-report)
    ))

(defun plos ()
  (setq *disease*
	(make-disease 
	 :name *disease-name*
	 :p-of-diagnosis-per-person-per-year *disease-p-of-diagnosis-per-person-per-year*
	 :initial-level-at-diagnosis *disease-initial-level-at-diagnosis*
	 :progression-multiplier-per-year *disease-progression-multiplier-per-year*
	 ;; [[Code removed for PLoS publication as not relevant.]]
	 :death-threshold *disease-death-threshold*
	 ))
  (run-test-set 
   "plos" ;; different length of sequential decline
   (let ((*run-replicates* 250))
     (loop for *how-many-sequential-declines-counts-as-lengthy-experience*
	   from 1 to 12 by 1
	   do
	   (sim 
	    (format nil "bestpub+pubafterexperience+seqdec~a" 
		    *how-many-sequential-declines-counts-as-lengthy-experience*)
	    (*treatment-model* :best-published-treatment)
	    (*publication-regimen* :publish-only-after-lengthy-experience)
	    ))
     ))
  )

(defparameter *results-vars* 
  '(*final-population*
    *people-with-disease*
    *percent-with-diease*
    ))

(defvar *final-population* ())
(defvar *people-with-disease* ())
(defvar *percent-with-diease* ())

(defun clear-results-vars ()
  (loop for var in *results-vars* 
	do (set var nil)))

(defvar *logging-level* 0)
(defvar *log-file* nil)
(defvar *log-file-name* nil)

(defmacro log! (level format &rest vars)
  ;; The level is optional. If it's a string, assume that it's actually
  ;; the format; regroups and change a ten to ten ones...
  (if (stringp level)
      (setq vars (cons format vars)
	    format level
	    level 0))
  `(when (>= *logging-level* ,level)
     (progn
       (format t ,format ,@vars)
       (format *log-file* ,format ,@vars)
       )))

(defparameter *disease* nil)

(defvar *people* nil)

(defun init-population ()
  (setq *next-birth-id* 1)
  (setq *people*
	(loop for p below *n-people*
	      collect (make-new-person 0.0))))

(defun make-new-person (p-scienist)
  ;; Note that we include the placebo here because for some people,
  ;; the placebo might be the personal best.
  (let* ((drugs (cons *placebo* (disease-drugs *disease*)))
	 (ndrugs (length drugs))
	 (person (make-person 
		  :id (prog1 *next-birth-id* (incf *next-birth-id*))
		  :scientist? (<= (/ (random 100) 100.0) p-scienist)
		  ;; [[Code removed for PLoS publication as not relevant.]]
		  ))
	 )
    ;; [[Code removed for PLoS publication as not relevant.]]
    person
    ))

(defun report-params ()
  (loop for param in *params*
	do (log! 0 (format nil "~a = ~~a~~%" param)
		 (symbol-value param)))
  (log! 2 "*diease* = ~s~%" *disease*)
  )

;;; Collects the run data (in reverse of replicates -- usually not relevant)
;;; and then at the end dumps them as a spreadsheet.

(defvar *trace-data-table* (make-hash-table :test #'equal))

(defun run ()
  (clear-results-vars)
  (sleep 1) ; ensure that previous run file is closed!
  (with-open-file 
   (*log-file* 
    (format nil "~a.log" (setq *log-file-name* 
			       (format nil "results/~a-~a" (get-universal-time) *run-label*)))
	       :direction :output)
   ;; [[Code removed for PLoS publication as not relevant.]]
   (report-params)
   (clrhash *trace-data-table*)
   (loop for n below *run-replicates*
	 do 
	 (format t "------------ ~a ------------~% " n)
	 (run1)
	 )
   (report-results)
   (record-in-results-spreadsheet)
   (record-trace-results)
   ))

(defun record-trace-results ()
  (let ((trace-file-name (format nil "~a.xls" *log-file-name*)))
      (log! 0 "Dumping trace spreadsheet in ~a" trace-file-name)
      (with-open-file 
       (tfo trace-file-name :direction :output)
       (loop for year below *years-to-run*
	     as sum = 0.0
	     with means = nil
	     as data = (gethash year *trace-data-table*)
	     do 
	     (format tfo "~a~c" year #\tab)
	     (loop for (pop nil) in data 
		   do 
		   (incf sum pop)
		   (format tfo "~a~c" pop #\tab))
	     (let ((mean (/ sum (length data))))
	       (format tfo "<-Pop,~c~a~c%Sick->~c" #\tab mean #\tab #\tab)
	       (push mean means)
	       )
	     (format tfo "~a~c" year #\tab)
	     (loop for (nil sick) in data
		      do (format tfo "~a~c" sick #\tab))
	     (format tfo "~%")
	     finally (setf (gethash *run-label* *master-average-tbl*)
			   means)
	     )
       )))

(defun final-averages-report ()
  (with-open-file 
   (o (format nil "results/~a-~a-MEANS.xls" (get-universal-time) *run-name*)
      :direction :output :if-exists :supersede)
   (loop for label being the hash-keys of *master-average-tbl*
	 using (hash-value values)
	 do 
	 (format o "~a~c" label #\tab)
	 (loop for value in (reverse values)
	       do (format o "~a~c" value #\tab))
	 (format o "~%")
	 )))

(defun record-in-results-spreadsheet ()
  (let ((*print-pretty* nil))
    (with-open-file 
     (rs "results.xls" :direction :output :if-exists :append :if-does-not-exist :create)
     (format rs "logfile~c"  #\tab)
     (loop for param in *params*
	   do (format rs "~a~c" param #\tab))
     (loop for var in *results-vars* 
	   do (format rs "~a-mean~c~a-sd~c~a-se~c~a-high~c~a-low~c" 
		      var #\tab var #\tab var #\tab var #\tab var #\tab))
     (format rs "~%")
     (format rs "~a~c" *log-file-name* #\tab)
     (loop for param in *params*
	   do (format rs "~a~c" (symbol-value param) #\tab))
     (loop for var in *results-vars* 
	   as vals = (sort (symbol-value var) #'>)
	   do (format rs "~a~c~a~c~a~c~a~c~a~c" 
		      (mean vals) #\tab 
		      (standard-deviation vals) #\tab
		      (standard-error vals) #\tab
		      (first vals) #\tab 
		      (first (last vals)) #\tab)
	   )
     (format rs "~%")
     )))

(defvar *year* nil)

(defun run1 (&aux n-with-disease n-people)
  ;; Temporary sanity test
  (if (and *do-not-repeat-treatments* 
	   (not (eq *treatment-model* :random)))
      (error "*do-not-repeat-treatments* only applies to the :random *treatment-model*"))
  (setup-disease-drug-relationships)
  (init-population)
  (clear-the-literature)
  (loop for *year* below *years-to-run*
	as n-people = (length *people*)
	as n-with-disease = (loop for person in *people*
				when (person-disease-state person)
				sum 1)
	as p-with-disease = (if *people* 
				(/ (float n-with-disease) n-people)
			      0.0)
	do 
	(log! 1 "===== Year ~a; population = ~a, ~a (~a%) have disease =====~%" 
	      *year* n-people n-with-disease p-with-disease
	      )
	(push (list n-people p-with-disease)
	      (gethash *year* *trace-data-table*))
	(loop for person in *people*
	      do (visit-the-doctor! person))
	(remove-dead-people)
	(discover-new-treatments)
	finally 
	(push n-people *final-population*)
	(push n-with-disease *people-with-disease*)
	(push p-with-disease *percent-with-diease*)
	))

(defun discover-new-treatments ()
  (case *discovery-model* 
	(:no-discovery nil)
	;; [[Code removed for PLoS publication as not relevant.]]
	;; [[Code removed for PLoS publication as not relevant.]]
	;; [[Code removed for PLoS publication as not relevant.]]
	(t (error "Unkown discovery model: ~a" *discovery-model*))))

;;; In the :add-hill-climbing-links-occassionally discovery model, we
;;; according to *discovery-p-of-adding-hill-climbing-link* add one
;;; hill climbing link to the drug information. This isn't really a
;;; very good model of discovery -- needs to be better thought
;;; through. FFF

;;; [[Code removed for PLoS publication as not relevant.]]

(defun report-results ()
  (loop for var in *results-vars* 
	as vals = (sort (symbol-value var) #'>)
	do 
	(log! 0 "~a m ~a sd ~a se ~a h ~a l ~a~%" var (mean vals) 
	      (standard-deviation vals) (standard-error vals)
	      (first vals) (first (last vals)))
	(when (and *summary-stream* (eq var '*final-population*))
	  (format *summary-stream* "~a~c~a~c~a~%" 
		  *run-label* #\tab (mean vals) #\tab (standard-error vals)))))

(defun visit-the-doctor! (person)
  (let* ((current-disease (person-disease person)))
    ;; First off we update the person's disease state if they are diagnosed.
    (when current-disease
      (let* ((true-disease-state (person-disease-state person))
	     (measured-disease-state 
	      (measure true-disease-state *disease-state-measurement-error-sd*)))
	;; Save the current state in their medical record. THIS IS MEASURED, NOT TRUE!
	(push (list (person-drug person) measured-disease-state)
	      (person-disease-state-history person))
	;; Now update their disease state. This is the TRUE diease state, NO MEASUREMENT ERROR!
	;; Note that this will make the diease ALWAYS get worse, and then the drugs will
	;; (one hopes) treat it! 
	(setf (person-disease-state person)
	      (* true-disease-state
		 (disease-progression-multiplier-per-year current-disease))))
      )
    ;; Next, if they are on a treatment regimen, update according to
    ;; that drug effectiveness. 
    (let ((on-drug (person-drug person)))
      (when on-drug
	(setf (person-disease-state person)
	      (max 0.00001 ;; make sure that it never bottoms out; 0.0 will blow up later.
		   (* (person-disease-state person)
		      (drug-effect (person-drug person)))))))
    ;; Consider publishing!
    (maybe-publish? person)
    ;; If they aren't getting better try a different treatment, unless they're out of options.
    (when (not (improving? person)) 
      (unless (person-out-of-options-do-not-treat person)
	(treat person current-disease)))
    ;; Finally report the state of a diagnosed individual.
    (when current-disease
      (log! 3 "~a (on ~a) new ~a state = ~a~%"
	    (person-id person)
	    (drug-id (person-drug person))
	    (disease-name current-disease)
	    (person-disease-state person))
      )
    ;; Give the person a disease with some probabiliy if they
    ;; don't already have one.
    (unless current-disease
      (when (<= (/ (random 100) 100.0)
		(disease-p-of-diagnosis-per-person-per-year *disease*))
	(setf (person-disease person) *disease*
	      (person-disease-state person) (disease-initial-level-at-diagnosis *disease*))
	(setf current-disease (person-disease person))
	(log! 2 "Person #~a diagnosed with ~a~%" 
	      (person-id person)
	      (disease-name current-disease))
	;; Now put them on a treatment regimen.
	(treat person current-disease) 
	))
    ))

(defun improving? (person)
  (let ((history (mapcar #'second (person-disease-state-history person))))
    (or (null (cdr history)) ; not enough records to make a determination
	(< (first history) (second history)))))

(defparameter *placebo* (make-drug :id "placebo" :effect 1.0))

(defun treat (person disease)
  ;; If they're already on a drug, just change them (unless no-switch),
  ;; otherwise, force them on a drug (first treatment!)
  ;; according to the *initial-treatment-model-if-need-be* parameter.
  (if (person-drug person) ;; Already on a drug...
      ;; Change unless it's a no-switch treatment model.
      (unless (eq :no-switch *treatment-model*)
	(set-drug person disease))
    ;; Not on anything yet, need to choose and initial drug.
    (let ((*treatment-model* 
	   ;; Only consider the initial treatment model if the
	   ;; main model is :no-switch or :hill-climb
	   (case *treatment-model*
		 ;; Temp bind the initial treatment model.
		 ((:no-switch :hill-climb) 
		  (or *initial-treatment-model-if-need-be*
		      (error "In :no-switch and :hill-climb, *initial-treatment-model-if-need-be* must be set!")
		      ))
		 (t *treatment-model*))))
      (set-drug person disease)))
  (log! 2 "Treating ~a with ~a~%" (person-id person) (drug-id (person-drug person)))
  )

(defun set-drug (person disease)
  (let* ((new-drug
	  (case *treatment-model* 
		;; [[Code removed for PLoS publication as not relevant.]]
		;; [[Code removed for PLoS publication as not relevant.]]
		;; [[Code removed for PLoS publication as not relevant.]]
		(:best-published-treatment 
		 (best-published-treatment disease))
		(t (error "In (TREAT ~a ~a) treatment model was ~a, which is unknown."
			  person disease *treatment-model*))
		)))
    ;; Now update and change the clinical history appropriately.
    (cond ((null new-drug) ;; can't find anything to switch to (do-not-repeat is on!)
	   ;; [Note that you'll never be here if the no-treat flag was already set!]
	   ;; Close the record on the last treatment.
	   (setf (clinical-record-date-ended 
		  (first (person-treatment-history person)))
		 *year*)
	   ;; Set their "do not treat" flag...
	   (setf (person-out-of-options-do-not-treat person) t)
	   ;; And wave them a fond fair well....
	   )
	  ((eq new-drug (person-drug person))) ;; Same drug, do nothing
	  ((null (person-drug person)) ;; first drug
	   ;; Create a clinical record entry for the new drug...
	   (push (make-clinical-record
		  :treatment (list new-drug)
		  :date-started *year*)
		 (person-treatment-history person))
	   ;; And treat!
	   (pushnew new-drug (person-drugs-tried person))
	   (setf (person-drug person) new-drug)
	   )
	  (t ;; Already on a drug...
	   ;; Close the record on the previous drug...
	   (setf (clinical-record-date-ended 
		  (first (person-treatment-history person)))
		 *year*)
	   ;; Now treat, as above, and add record entries.
	   ;; Create a clinical record entry for the new drug...
	   (push (make-clinical-record
		  :treatment (list new-drug)
		  :date-started *year*)
		 (person-treatment-history person))
	   ;; And treat!
	   (pushnew new-drug (person-drugs-tried person))
	   (setf (person-drug person) new-drug))
	  ) ;; Close cond
    ))

(defun best-published-treatment (disease)
  (loop for pub in *the-literature*
	as p below *how-many-back-pubs-to-look-at* 
	with best-treatment = nil
	as this-delta = (publication-delta pub) ;; WWW BBB This is a MEAN, NOT A DELTA!!! (I think)
	as best-delta = 10000.0	;; don't know what this has to be to start?
	when (< this-delta best-delta)
	do (setf best-treatment (publication-drug pub))
	finally (return (or best-treatment
			    (nth (random (length (disease-drugs disease))) (disease-drugs disease))))
	))

(defun remove-dead-people ()
  (setq *people*
	(loop for person in *people*
	      when (or (null (person-disease-state person))
		       (< (person-disease-state person)
			  (disease-death-threshold (person-disease person))))
	      collect person))
  )
				   
;;; Only a small number of drug treat a given disease, and with
;;; varying degrees of success. The diease model 'knows' these
;;; bindings, but WE don't! 

;;; The drug effect is a down-regulator on the disease.

(defun setup-disease-drug-relationships ()
  (log! 2 "Setting up disease-drug relationships.")
  ;; We'll assume that there are so many drugs relative to the
  ;; number selected, that we don't have to test for
  ;; duplication.
  (setf (disease-drugs *disease*) nil) ;; clear from previous runs
  (case *treatment-effectiveness-distribution-model*
	;; [[Code removed for PLoS publication as not relevant.]]
	(:even 
	 ;; Use the value given value as a low end, and create n even spaces
	 (loop for effect from *most-effective-drug-power*
	       to (+ *treatment-effectiveness-mean* 
		     (- *treatment-effectiveness-mean* 
			*most-effective-drug-power*))
	       by (/ (* 2.0 
			(- *treatment-effectiveness-mean* 
			   *most-effective-drug-power*))
		     *initial-compounds-that-address-a-given-disease*)
	       as n from 1 by 1
	       as ddr = (make-drug :id n :effect effect)
	       do (log! 2 "Even DDR: ~s~%" ddr)
	       (push ddr (disease-drugs *disease*))))
	;; [[Code removed for PLoS publication as not relevant.]]
	(t (error "Unknown *treatment-effectiveness-distribution-model* = ~s"
		  *treatment-effectiveness-distribution-model*))
	)
  ;; Now add additional drugs that do nothing (if desired)
  (loop for i below *additional-non-functional-drugs*
	as n from (1+ (length (disease-drugs *disease*))) by 1
	as ddr = (make-drug :id n :effect 1.0)
	do (log! 2 "Ineffective DDR: ~s~%" ddr)
	(push ddr (disease-drugs *disease*)))
  ;; Finally, create hill-climbing relationships by filling in the next-best-treatment slots.
  ;; (Unless the discovery model prevents, in which case these will be created dynamically.)
  (case *discovery-model*
	(:no-discovery (loop for (t1 . +) on (sort (copy-list (disease-drugs *disease*))
						   #'> :key #'drug-effect)
			     do (setf (drug-next-best-treatment t1)
				      (car +))))
	;; [[Code removed for PLoS publication as not relevant.]]
	(t (error "Unknown discovery model: ~a" *discovery-model*)))
  )

;;; Gaussian random number generation.  From:
;;; http://www.bearcave.com/misl/misl_tech/wavelets/hurst/random.html

(defun measure (mean sigma)
  (let (x y)
    (declare (double x y r2))
    (flet 
     ((r () (/ (random 1000) 1000.0)))
     (let ((val
	    (loop with r2 = 0.0 with r1 = 0.0 
		  while (or (> r1 1.0) (= r2 0.0))
		  do 
		  ;; choose x,y in uniform square (-1,-1) to (+1,+1)
		  (setf x (+ -1 (* 2 (r))))
		  (setf y (+ -1 (* 2 (r))))
		  ;; see if it is in the unit circle
		  (setf r2 (+ (* x x) (* y y)))
		  finally (return (* sigma y (sqrt (* -2.0 (/ (log r2) r2)))))
		  )))
       (+ mean (if (complexp val) (imagpart val) val))))))

;;; All publications go on in "the literature", which is just a stack
;;; of published case reports, etc. There's no distinction made here
;;; between blogs and papers, case reports and reviews. FFF maybe
;;; there should be?

(defun clear-the-literature ()
  (setq *the-literature* nil))

(defstruct publication drug delta case)

(defun maybe-publish? (person)
  (let* ((history (person-disease-state-history person))
	 (state-history (mapcar #'second (person-disease-state-history person)))
	 (drug-history (mapcar #'first (person-disease-state-history person)))
	 (current-drug (person-drug person))
	 (delta (when (second state-history) (/ (first state-history) (second state-history))))
	 )
    (when delta
      (case *publication-regimen*
	    ;; [[Code removed for PLoS publication as not relevant.]]
	    (:publish-only-after-lengthy-experience
	     (when (nthcdr *how-many-sequential-declines-counts-as-lengthy-experience* history)
	       (let ((experiences 
		      (loop for n below *how-many-sequential-declines-counts-as-lengthy-experience*
			    as entry in history
			    collect entry))
		     (experiences-back-one ;; Used to calculate the delta later.
		      (loop for n below *how-many-sequential-declines-counts-as-lengthy-experience*
			    as entry in (cdr history)
			    collect entry)))
		 (if (and 
		      ;; Has to all be the same drug
		      (loop with d1 = (caar experiences)
			    as (d2) in (cdr experiences)
			    unless (eq d1 d2)
			    do (return nil)
			    finally (return t))
		      ;; And always decreasing disease state.
		      (loop as (nil v1) in experiences
			    as (nil v2) in (cdr experiences)
			    unless (< v1 v2)
			    do (return nil)
			    finally (return t)))
		     (push (make-publication :drug (caar experiences)
					     ;; !!! This was WRONG. The paper claims it's a delta, but it was a mean!
					     ;; ??? Which way should this go? Should it be ABS?
;; 					     :delta (- (/ (reduce #'+ (mapcar #'second experiences)) ;; Mean of now
;; 							  *how-many-sequential-declines-counts-as-lengthy-experience*)
;; 						       (/ (reduce #'+ (mapcar #'second experiences-back-one)) ;; Mean of back one
;; 							  *how-many-sequential-declines-counts-as-lengthy-experience*)
;;						       )
					     ;; !!! This is the (published!) WRONG version:
					     :delta (/ (reduce #'+ (mapcar #'second experiences)) ;; Mean of now
						       *how-many-sequential-declines-counts-as-lengthy-experience*)
					     ;; :case person
					     )
			   *the-literature*)))))
	    (t (error "in MAYBE-PUBLISH?: Unknown *publication-regimen*: ~s" *publication-regimen*))
	    ))))

;;; Imported from Jeff's stats for PLoS publication completeness:

(defun sum (l &aux (sum 0))
  (dolist (n l) (incf sum n)) sum)
(defun mean (l)
  (/ (sum l) (float (length l))))
(defun standard-deviation (l)
  (if (null (cdr l))
      0.0
    (sqrt
     (let ((m (mean l)))
       (* (/ 1.0 (1- (length l)))
	  (sum (mapcar #'(lambda (x) (expt (- x m) 2)) l))
	  )))))
(defun standard-error (l)
  (/ (standard-deviation l) (sqrt (length l))))

;;; Self-execution on load:

(plos)

