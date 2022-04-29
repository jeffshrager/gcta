#|

Complilation sequence:

/Users/jeffshrager/Desktop/cc/GCTA/sim/RL/sim: lisp
(load (compile-file "compile-simpx.lisp"))
(save-application "compile-simpx.exe" :toplevel-function #'main :prepend-kernel t) 
rm -rf experiments/201707081100-new-choice-algorithm/
tf
python master-script.py
: lisp
? (load (compile-file "compile-simpx.lisp"))
? (save-application "compile-simpx.exe" :toplevel-function #'main :prepend-kernel t) 
: cd ..
: tf
: python master-script.py

|#

;; To-Dos
;; For Jeff: Why does it seem to randomly pick 1 or 2 drugs in the first round?
;; Add in costs (add drug costs as an adversity to subtract from our k)
;; Reading in and using the hypotheses output by nathaniel 
;; Currently we change on decline for drugs, could we create a more elaborate condition?
;; Choose-drugs needs to be updated for multiple matches - getting best option - and also no matches
;; Fix untried-mono, binary-1-at-time, and untried-binary to adapt to changing tumors - or fix the way patient drug histories are stored
;; ?? Use untried treatment pairs, with replacement (i.e. all combos) ?? - modelling drug synergy? Include a synergy score for combos?
;; Fix the hypotheses functions so that matrix only read once each main loop or something and remove other ones that are now unused
;; Fix the defparameters to be scoped globally at the top --> clean up the code 
;; fix the *markers* for the new set-up style

;;; Simulating patients for Musella "Virtual Trials" application. When
;;; someone is diagnosed they start out with a karnofsky (k) score of
;;; 70-90, and then, dependings on treatment and tumor aggressiveness,
;;; either reach 0, and are dead, or reach 100, and just stay there.

;; TO RUN - We use the Clozure Common Lisp implementation of Common Lisp 
;; Navigate to the folder containing the relevant files and start the lisp listener (REPL)
;; Enter:
;; (load (compile-file "simpx.lisp"))
;; Hit enter to run

(eval-when (:compile-toplevel :load-toplevel :execute)
	   (unless (find-package :statistics) (load "lhstats"))
           (ql:quickload "cl-json")
           )

;; ==============================================================================================
;; Treatment Rules - on decline and stops when limit reached
;;		random-monotherapy: randomly selects 1 drug to give the patient, with replacement.
;;		untried-monotherapy: randomly selects 1 drug to give the patient, without replacement (i.e. can't recieve a drug you have already received).
;; 		random-binary-cocktail: randomly selects 2 drugs to give the patient, with replacement.
;; 		binary-cocktail-changing-1-at-a-time: maintains two drugs at all times, randomly changing one out without replacement. 
;;											  If no other options, stick with current cocktail.
;;		untried-binary-cocktail: randomly selects 2 drugs to give the patient, without replacement. 
;;		targeted-therapy: ultimately, it will use a set of hypotheses supplied by experts/classifier to assign drugs to patients based on their mutations. 
;;						  Drugs are still randomly assigned, but skewed towards the theoretically more effective drugs. 

;;; ==============================================================================================
;;; Object definitions

;;; The patient history is in rev. chron order so that the first thing is the most recent treatment 
(defstruct (patient (:print-function
		     (lambda (p s k)
		       (format s "((stype . patient) (id . ~a) (birthyear . ~a) (gender . ~a) (mdx . ~a) (k . ~a) (tumor . ~a) (diagnosisdate . ~a) (drugs . ~a) (history  ~a))"
			       (patient-id p) (patient-birth-year p) (patient-gender p) 
			       (patient-mbdx p) (patient-k p) (patient-tumor p)
			       (patient-dx-utime p)
			       (patient-drugs p) (patient-history p)))))
  id
  birth-year
  gender
  mbdx ;; mbdx = months beyond dx (starts at 0)
  k ;; k = karnofsky score
  tumor 
  dx-utime
  drugs
  history
  hypotheses 
  previous-drugs
  )

(defstruct (drug (:print-function
		  (lambda (d s k)
		    (format s "((stype . drug) (name . ~a) (markers . ~b) (adversity . ~a) (mrh . ~a))"
			    (drug-name d) (drug-markers d) (drug-adversity d) (drug-mrh d)))))
  name markers adversity mrh) ;; mrh ==  molecular-response-history

(defun name->drug (name)
	(loop for drug in *drugs* 
		when (equal (drug-name drug) name)
		do (return drug)))

(defun param-reader(location version)
	(with-open-file (input location)
      	(let* ((data (json:decode-json input)))
      		(case version
      			(:simulator
	      		   	(setf params (loop for param in (rest (assoc :simulator data))
	      		   		collect param))
	      			(return-from param-reader params))
	      		(:run-label
	      			(return-from param-reader (rest (assoc :run-label data))))
	      		(:sim-output
	      			(return-from param-reader (rest (assoc :sim-output data))))
	      		(:classifier-output
	      			(return-from param-reader (rest (assoc :classifier-output data))))))))

;; Turns a string into a keyword
(defun make-keyword (name) (intern (string-upcase name) "KEYWORD")) 

;; Kills the cl-json method for converting camel-case to lisp (since it naturally inserts "*" or "+" into text - and we don't want that)
(defun cl-json:camel-case-to-lisp (s) (string-upcase (string-trim ":" s)))(defun cl-json:camel-case-to-lisp (s) (string-upcase (string-trim ":" s)))

;; Retrieved from: http://faculty.washington.edu/dbp/SAPACLISP-1.x/matrix.lisp
(defun list-to-2d-array (list)
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))

;; Retrieved from: http://faculty.washington.edu/dbp/SAPACLISP-1.x/matrix.lisp
(defun multiply-matrix-and-vector (a-matrix b-vector &key (result (make-array (nth 0 (array-dimensions a-matrix)))))
  (let ((m (nth 0 (array-dimensions a-matrix)))
        (n (length b-vector)))
    (dotimes (i m result)
      (setf (aref result i) 0.0)
      (dotimes (j n)
	(incf (aref result i)
	      (* (aref a-matrix i j) (aref b-vector j)))))))

;; Retrieved from:
;; http://faculty.washington.edu/dbp/SAPACLISP-1.x/matrix.lisp

(defun multiply-two-matrices
  (a-matrix
   b-matrix
   &key
   (result
    (make-array
     (list (nth 0 (array-dimensions a-matrix))
	   (nth 1 (array-dimensions b-matrix))))))
  (let ((m (nth 0 (array-dimensions a-matrix)))
        (n (nth 1 (array-dimensions b-matrix)))
        (common (nth 0 (array-dimensions b-matrix))))
    (dotimes (i m result)
      (dotimes (j n)
        (setf (aref result i j) 0.0)
        (dotimes (k common)
          (incf (aref result i j)
                (* (aref a-matrix i k) (aref b-matrix k j))))))))

(defvar *hypotheses-matrix* nil)
(defvar *drug-dict* nil)

(defun read-hyps-2 (location)
  "Pulls the LAST hypotheses and drug dictionary."
  (with-open-file
   (i location)
   (loop for j = (handler-case (json:decode-json-from-source i) (end-of-file nil))
	 as n from 1 by 1
         until (null j)
         do
	 (setf *hypotheses-matrix* (list-to-2d-array (rest (assoc :BEST-MATRIX j)))
	       *drug-dict* (print (rest (assoc :DRUG-DICT j)))))))

;;; ==========================================================================================
;;; Globals - Note that many parameters are now contained in the "inputs/parameters.json" file

(defparameter *parameters* (param-reader "../parameters.json" :simulator))

(defparameter *sim-version* (get-universal-time))
(defparameter *patient-id-counter* 0) ;; Super global reset on each load
(defparameter *fcounter* 1000) ;; File differentiator -- start here and incf'ed for each file in a run
(defparameter *run-label* (param-reader "../parameters.json" :run-label))
(defparameter *sim-output* (param-reader "../parameters.json" :sim-output))
(defparameter *classifier-output* 
  (format nil "../experiments/~a/~a/" *run-label* (param-reader "../parameters.json" :classifier-output)))

(defvar *mrhlog* nil)
(defvar *mrhlog-as-list* nil)
(defvar *dump-pxs-in-json?* nil)

;; The genetic dictionary assigns a gensym to each combination in
;; order that they have names in the mrh log file. III ***
;; Importantly, the keys are UNORDERED! In order (so to speak) to make
;; that happen, we do a numberical encoding trick on them and sort
;; it. UUU) The only reason we need this is because not all
;; combinations are actually seen. In fact, very few of the actual
;; combinations are seen since most genetics are combos of just 2 or 3
;; mutations. So this saves us from having to list all gazillion
;; combinations of up to (length *markers*) markers, when we only see
;; a fraction of those.

(defvar *genetics* (make-hash-table :test #'equal))

;; Probability for each mutation that it will appear in the tumor
;; (Note: I use 'marker' and 'mutation' interchangably.)
(defparameter *p-mutation* (rest (assoc :p-mutation *parameters*))) ;; @0.25 you tend to get two markers
(defvar *param-header* nil)
(defvar *log* t)

(defparameter *txrules* (loop for rule in (rest (assoc :txrules *parameters*)) collect (make-keyword rule))) 
(defparameter *txrule* (make-keyword (rest (assoc :txrule *parameters*)))) ;; Has the current rule (Careful of the difference btwn the singular and plural here!)

;;; Tumor aggresiveness is correlated with the number of
;;; mutations. You get *k-score-decr-per-mutation-per-month* points
;;; for each mutation substracted from your k score per month.
(defparameter *k-score-decr-per-mutation* (rest (assoc :k-score-decr-per-mutation *parameters*)))
(defparameter *k-score-incr-per-drug-mutation-overlap* (rest (assoc :k-score-incr-per-drug-mutation-overlap *parameters*)))
(defparameter *p-mutational-escape* (rest (assoc :p-mutational-escape *parameters*)))
(defparameter *hypotheses-threshold* 0.3)

(defparameter *n-patients* (rest (assoc :n-patients *parameters*)))
(defvar *pxs* nil)
(defparameter *nmonths* (rest (assoc :nmonths *parameters*)))
(defvar *p* nil) ;; global for debugging, will contain the last patient record

(defvar *date* 0) ;; Everyone starts at 0 and every cycle is a virtual month.

(defvar *hypotheses* nil) ;; Holds the initial hypotheses ("priors") generated by classifier or experts

;;; =======================================================
;;; Drugs, trials, and other contextual parameters

;; Mutation is really simple; if you have the key in the markers list,
;; then the tumor has that mutation. WWW The order of these is used to
;; create genetic codes later in the ... code ... this usually won't
;; matter, unless you're trying to make comparisons between runs.

(defparameter *markers* 
  (loop for name in (rest (assoc :markers *parameters*)) 
	collect (make-keyword (string-trim "+" name)))) 

;;; Treatments are each assigned pairs of markers that they work well
;;; with. If the patient's tumor has both markers, then you get
;;; maximum effect from the drug, otherwise, you get progressively
;;; less benefit. Each drug also has an adversity score, which gets
;;; SUBTRACTED each time the drug is used.

(defparameter *less-confounded-drugs*
  (list
   ;; All combinations of markers except allomeric pairs (e.g., no IHD1 and IDH2)
   (make-drug :name :Bevacizumab :markers '(:EGFR :IDH1) :adversity 1)
   (make-drug :name :Temozolomide :markers '(:IDH1 :IDH2) :adversity 1)
   (make-drug :name :Cabazitaxel :markers '(:IDH2 :HLA-A1) :adversity 1)
   (make-drug :name :TCAR :markers '(:HLA-A1 :HLA-A2) :adversity 1)
   (make-drug :name :Disatinib :markers '(:HLA-A2 :MGMT-Promoter) :adversity 1)
   (make-drug :name :Nivolumab :markers '(:EGFR) :adversity 1)
   (make-drug :name :Doxorubicin :markers '(:IDH1) :adversity 1)
   (make-drug :name :varlilumab :markers '(:IDH2) :adversity 1) 
   (make-drug :name :Durvalumab :markers '(:MGMT-Promoter) :adversity 1)
   (make-drug :name :Sorafenib :markers '(:HLA-A1) :adversity 1)
   (make-drug :name :Pembrolizumab :markers '(:HLA-A2) :adversity 1)
   ))

(defparameter *not-confounded-drugs*
  (list
   ;; All combinations of markers except allomeric pairs (e.g., no IHD1 and IDH2)
   (make-drug :name :Bevacizumab :markers '(:EGFR) :adversity 1)
   (make-drug :name :Temozolomide :markers '(:IDH1) :adversity 1)
   (make-drug :name :Cabazitaxel :markers '(:IDH2) :adversity 1)
   (make-drug :name :TCAR :markers '(:HLA-A1) :adversity 1)
   (make-drug :name :Disatinib :markers '(:HLA-A2) :adversity 1)
   (make-drug :name :Nivolumab :markers '(:MGMT-Promoter) :adversity 1)
   ))

(defparameter *trials*
  '(("NCT01109095" "TCAR")
    ("NCT02758366" "Doxorubicin")
    ("NCT01588769" "good immune markers" "ALECSAT")
    ("NCT02458508" "Melanocortin Receptor 4" "Radio-chemotherapy with temozolomide")
    ("NCT01280552" "HLA-A1 or HLA-A2 positive" "ICT-107 Immunotherapy")
    ("NCT02546102" "HLA-A2 positive" "ICT-107 Immunotherapy")
    ("NCT02664363" "EGFRvIII (del exons2-7) EGFRvIII" "CAR T cells")
    ("NCT01269424" "good immune markers" "O6-Benzylguanine and Temozolomide With MGMTP140K Genetically Modified Blood Stem Cells")
    ("NCT02977780" "IDH1 R132H negative,  MGMT promoter is unmethylated" "Temozolomide, Neratinib,CC-115,Abemaciclib")
    ("NCT03018288" "IDH wild type and MGMT promoter is unmethylated" "pembrolizumab + vaccine HSPPC-96")
    ("NCT02667587" "MGMT methylated or indeterminate tumor subtype" "Nivolumab")
    ("NCT02617589" "Unmethylated MGMT" "Nivolumab + TMZ")
    ("NCT02658981" "MGMT methylation status" "Anti-LAG-3  or Anti-CD137 Alone or with Anti-PD-1")
    ("NCT02327078" "EGFR and ALK status" "Nivolumab + Epacadostat")
    ("NCT02152982" "MGMT hypermethylation status" "Veliparib + TMZ")
    ("NCT01480479" "EGFRvIII positive" "Rindopepimut/GM-CSF (virus)")
    ("NCT02573324" "EGFR amplification" "ABT-414")
    ("NCT02414165" "IDH mutation status" "Toca 511, a Retroviral Replicating Vector +  Toca FC")
    ("NCT01491893" "MGMT status" "PVSRIPO")
    ))

;;; =======================================================
;;; Create patients and their tumors

(defun genpx ()
  (make-patient 
   :id (incf *patient-id-counter*)
   :birth-year (- 2017 (if (zerop (random 2)) 
			   (+ (random 10))
			 (+ (random 20) 50)))
   :gender (if (zerop (random 2)) :m :f)
   :tumor (gentumor)
   :k (+ 70 (random 20)) ;; Patients always present between 70 and 90
   :mbdx 0
   ;; After Jan 1 2015 == 3629091661
   :dx-utime (+ 3629091661 (random 76414218)) ;; Should go two years or so into the future
   :drugs nil ;; Keeps track of what drugs the px is on NOW
   :previous-drugs nil
   :history `(((:action . :dx) (:date . 0))))) ;; Every history entry starts with an action key.

;; Tumor load translates to karnofsky score.

(defun gentumor ()
  (if (eq :sanity-check *txrule*)
      ;; Under sanity-testing always generates single-marker tumors
      (list (nth (random 6) *markers*))
    ;; Normal generation
    (loop as ms = (loop for marker in *markers*
			if (< (/ (random 100) 100.0) *p-mutation*)
			collect marker) ;; Can't be nil
	  until ms
	  finally (progn (gene-code ms)
			 (return ms)))
    ))

;;; =======================================================
;;; Updating disease status

;;; Tumor aggresiveness is correlated with the number of
;;; mutations. You get *k-score-decr-per-mutation-per-month* points
;;; for each mutation substracted from your k score per month.

(defun update-k (px)
  ;; Decr k for disease load, then incr for treatments
  (loop for mutation in (patient-tumor px)
	do (decf (patient-k px) *k-score-decr-per-mutation*))
  ;; Now incr for overlapping treatments
  (loop for tx in (patient-drugs px)
	as dms = (drug-markers tx)
	do (loop for m in (patient-tumor px)
		 when (member m dms)
		 do (incf (patient-k px) *k-score-incr-per-drug-mutation-overlap*))
	;; While were here, do adversity
	(decf (patient-k px) (drug-adversity tx)))
  ;; Now chop the ceiling and floor
  (setf (patient-k px) (max (min (patient-k px) 100) 0))
  ;; Finally, stochastically mutate the tumor
  (if (> *p-mutational-escape* (/ (random 100) 100.0))
      (setf (patient-tumor px) (gentumor)))
  )

;;; =======================================================
;;; Drug decision making

;;; This extremely unhygenic macro is used to compute drug
;;; updates. You get a bunch of vars that you can (unhygenically) use,
;;; and you have to leave the new drugs list in DRUGS. Whatever DRUGS
;;; is at the end will get reset back into the px record. However, if
;;; you don't change DRUGS (that is, if it ends up being the same as
;;; the current px drugs) then no change is recorded.  Note that the
;;; "only change on decline" is built in to this; perhaps some day it
;;; should be conditionalized.

#|
(make-drug :name :Bevacizumab :markers '(:EGFR :IDH1) :adversity 1)
   (make-drug :name :Temozolomide :markers '(:IDH1 :IDH2) :adversity 1)
   (make-drug :name :Cabazitaxel :markers '(:IDH2 :HLA-A1) :adversity 1)
   (make-drug :name :TCAR :markers '(:HLA-A1 :HLA-A2) :adversity 1)
   (make-drug :name :Disatinib :markers '(:HLA-A2 :MGMT-Promoter) :adversity 1)
   (make-drug :name :Nivolumab :markers '(:EGFR) :adversity 1)
   (make-drug :name :Doxorubicin :markers '(:IDH1) :adversity 1)
   (make-drug :name :varlilumab :markers '(:IDH2) :adversity 1) 
   (make-drug :name :Durvalumab :markers '(:MGMT-Promoter) :adversity 1)
   (make-drug :name :Sorafenib :markers '(:HLA-A1) :adversity 1)
   (make-drug :name :Pembrolizumab :markers '(:HLA-A2) :adversity 1)
   {"EGFR": 0, "IDH1": 1, "IDH2": 2, "HLA-A1": 3, "HLA-A2": 4, "MGMT-PROMOTER": 5}
   |#

(defparameter *fml* ;; trival (level 0)
  '(
    (1 :Nivolumab 5)
    (1 :Cabazitaxel -5)
    (2 :Nivolumab 5)
    (2 :Cabazitaxel -5)
    (4 :Nivolumab 5)
    (4 :Cabazitaxel -5)
    (8 :Nivolumab 5)
    (8 :Cabazitaxel -5)
    (16 :Nivolumab 5)
    (16 :Cabazitaxel -5)
    (32 :Nivolumab 5)
    (32 :Cabazitaxel -5)
    )
  )

(defparameter *fml* ;; unconfounded (level 1)
    '(
     (1 :Nivolumab 5)
     (1 :Cabazitaxel -5)
     (2 :Cabazitaxel 5)
     (2 :Doxorubicin 5)
     (4 :Doxorubicin -5)
     (4 :Varlilumab 5)
     (8 :Varlilumab -5)
     (8 :Sorafenib 5)
     (16 :Sorafenib -5)
     (16 :Pembrolizumab 5)
     (32 :Pembrolizumab -5)
     (32 :Durvalumab 5)
  ))

(defparameter *fml* ;; confounded (level 2)
    '(
     (1 :Nivolumab 5)
     (1 :Cabazitaxel -5)
     (2 :Cabazitaxel 5)
     (2 :Nivolumab -5)
     (4 :Nivolumab 5)
     (4 :Doxorubicin -5)
     (8 :Doxorubicin 5)
     (8 :Nivolumab -5)
     (16 :Doxorubicin 5)
     (16 :Disatinib -5)
     (32 :Disatinib 5)
     (32 :Doxorubicin -5)
  ))

(eval-when 
 (:compile-toplevel :load-toplevel :execute)
 (defmacro redrug (px initexpr &body new-drugs-body)
   ;; Sets up some unhygenic vars that can be used in the body: ks, ds, cds, well as drugs (which should be changed)
   `(let* ((ks-5 (histscan-fast-5 :scan ,px #'cdar)) ;; the reverse chron of k scores
	   ;;(ds (reduce #'append (histscan :changed-tx ,px #'(lambda (l) (cdr (assoc :result l)))))) ;; all previously used drugs
	   (cds (patient-drugs ,px))  ;; and cds -- current drugs
	   (drugs (patient-drugs ,px)) ;; If this isn't reset by the body, we don't record a change.
	   (histlen (length ks-5))
	   (previous-drugs (patient-previous-drugs ,px))
	   )
      (setf drugs
	    (if (= 1 histlen) ,initexpr
	      ;; Change when the patient is declining Pretty simple test for the
	      ;; moment -- just looks at the first two!
	      (progn
		;; In the molecular response history of each drug in the
		;; current set record how the px's k changed for these
		;; drugs. Confusingly, we do this twice, once with the
		;; cocktail elements separated, and once together. This is
		;; somewhat confusing and needs to be sorted out on the back
		;; end. First the separate drugs in the cocktail.
		(if (eq *txrule* :sanity-check)
		    (let* ((line (nth (random (length *fml*)) *fml*))
			   (tumor (first line))
			   (drug (second line))
			   (dk (third line))
			   (external `(0 0 ,tumor 1 ,drug 0 0 ,dk))
			   )
		      (format *mrhlog* "~{~a~^  ~}~%" external)
		      (setf (fourth external) 2)
		      (format *mrhlog* "~{~a~^  ~}~%" external)
		      )
		  (loop for drug in drugs
			as dn from 1 by 1 ;; Cocktail counter
			do 
			;; Record for learning (and general tracing)
			(push (list (patient-tumor ,px) (first ks-5) (second ks-5)) (drug-mrh drug)) ;; [1]
			(push (list (patient-tumor ,px) (first ks-5) (second ks-5) drug) *mrhlog-as-list*)
			(format *mrhlog* "~a	~a	~a	~a	~a	~a	~a	~a	~a~%" 
				(patient-id ,px)
				histlen
				(gene-code (patient-tumor ,px))
				dn
				(drug-name drug)
				(gene-code (drug-markers drug))
				(first ks-5)
				(second ks-5)
				(- (first ks-5) (second ks-5))
				)
			))
		;; Now redrug if things are going south!
		(if (< (first ks-5) (second ks-5))
		    (progn ,@new-drugs-body)))
	      ))
      ;; Here DRUGS should be either changed by the body, or not. If
      ;; it turns out to be the same as what's already there, then we
      ;; record a :NO-CHANGE
      (unless (set-equal drugs (patient-drugs ,px)) ;; Wasn't changed!
	(loop for drug in drugs do (pushnew drug (patient-previous-drugs ,px)))
	(record ,px :changed-tx (cons :rule *txrule*) (cons :result drugs))
	(setf (patient-drugs ,px) drugs))
      )))

(defun assess-and-treat (px)
  ;; Scan (so to speak) the patient ... basically just read off their k score
  (record px :scan `(:k= . ,(patient-k px)))
  ;; Treat according to the rule
  (case *txrule*
	(:random-binary-cocktail (redrug px (n-random 2 *drugs*) 
					 :random-binary-cocktail (n-random 2 *drugs*)))
	(:random-monotherapy
	 (redrug px (n-random 1 *drugs*) 
		 :random-monotherapy (list (nth (random (length *drugs*)) *drugs*))))
	(:untried-monotherapy
	 ;; Change to a drug that isn't in the ds list.
	 (redrug px (n-random 1 *drugs*) 
		 :untried-monotherapy
		 (let* ((possible-drugs (set-difference *drugs* previous-drugs))
			(new-drug (when possible-drugs (nth (random (length possible-drugs)) possible-drugs))))
		   (if new-drug (list new-drug) drugs))))
	(:binary-cocktail-changing-1-at-a-time
	 ;; Maintain 2 at all times, randomly changing out one of the pair (on decline).
	 (redrug px (n-random 2 *drugs*) 
		 :binary-cocktail-changing-1-at-a-time 
		 (binary-cocktail-changing-1-at-a-time drugs previous-drugs cds)))
	(:untried-binary-cocktail
	 (redrug px (n-random 2 *drugs*) 
		 :untried-binary-cocktail (untried-binary-cocktail drugs previous-drugs cds)))
	(:targeted-therapy ;; Currently starting on binary cocktail
	 (redrug px (n-random 2 *drugs*) 
		 :targeted-therapy (choose-drugs (patient-tumor px)))) 
	(:hypotheses-therapy
	 (redrug px (hypotheses-therapy (patient-tumor px))
		 :hypotheses-therapy (hypotheses-therapy (patient-tumor px))))
	(:sanity-check
	 (redrug px (sanity-check (patient-tumor px))
		 :sanity-check (sanity-check (patient-tumor px))))
	(t (error "Bad txrule in assess-and-treat ~a" *txrule*))
	)
  ;; Stop if they're dead
  (when (or (zerop (patient-k px)) (= 100 (patient-k px))) (throw 'done px))
  )

(defun binary-cocktail-changing-1-at-a-time (drugs ds cds)
  (let* ((possible-drugs (set-difference *drugs* ds)))
    (if possible-drugs ;; If there's nothing else we can do, just leave things as is.
	(let ((random-new-drug (nth (random (length possible-drugs)) possible-drugs)))
	  (if (< (length cds) 2) ;; Early on there won't be 2 in the binary-cocktail -- just add one randomly
	      (setf drugs (cons random-new-drug cds))
	  ;; Otherwise (=2) kick one out randomly, and add a new one randomly
	  (list random-new-drug (nth (random 2) cds))))
      ;; Need to return the current list if no change.
      drugs)))

(defun untried-binary-cocktail (drugs ds cds)
	(let* ((possible-drugs (set-difference *drugs* ds)))
		(if (and possible-drugs (>= (length possible-drugs) 2)) ;; if there are no untried drugs, leave it as is. currently we miss one drug (since there are 11)
			(let* ((random-new-drug-1 (nth (random (length possible-drugs)) possible-drugs))
			       (foo (delete random-new-drug-1 possible-drugs))
			       (random-new-drug-2 (nth (random (length possible-drugs)) possible-drugs)))
			(setf drugs (list random-new-drug-1 random-new-drug-2)))))
	drugs)

(defun choose-drugs (muts)
 	;; what to do with no-match - try another treatment model or assume matches all of them
 	;; what happens if matches mutliple - answers: best match criteria (can have multiple best matches -> combines suggestions)
 	;; The above needs to be implemented
 	(loop for drug in (choose-most-likley-drugs (matching-drug-lists muts))
 		collect (name->drug drug)))

 (defun matching-drug-lists (tmuts)
 	(loop for (dmuts drugs) in *hypotheses*
 		when (mut-match? tmuts dmuts)
 		append drugs)) ;; if you have multiple matches -> append drug lists 

(defun mut-match? (tmuts dmuts &aux (matches 0))
  (loop for (dm mk) in dmuts
	do (case mk 
		 (1 (loop for tm in tmuts
			  when (equalp dm tm) do (return (incf matches))))
		 (-1 (loop for tm in tmuts
			  when (equalp dm tm) do (return-from mut-match? nil))) ;; to exit the function immediately for contra-indications
		 (0 :continue)))
  (if (> matches 0) t nil))
 	
(defun choose-most-likley-drugs (drugs &aux (probs 0) (total 0))
  (let* ((stow 
	  (if (= (length drugs) 0)
	      (loop for (dmuts previous-drugs) in *hypotheses*
		    append previous-drugs)
	    drugs)))
    (loop for (score previous-drugs) in stow 
	  do (incf total score))
    (loop for (score previous-drugs) in stow
	  with target = (/ (random 10000) 10000)
	  do (incf probs (/ score total))
	  (when (<= target probs) 
	    (return-from choose-most-likley-drugs previous-drugs)))))

(defvar *mut-mat* nil)

(defun hypotheses-therapy (tmuts &aux ds-index)
  (if *hypotheses-matrix*
      (progn
	;; Setup weight matrix:
	(setf *mut-mat* (make-array '(6 6) :initial-element 0))
	(loop for i from 0 to 5 do (setf (aref *mut-mat* i i) 1))
	(setf *hypotheses-matrix* (multiply-two-matrices *mut-mat* *hypotheses-matrix*))
	;; Setup tumor vector:
	(setf pmuts (make-array '(1 6) :initial-element 0))
	;; Sim and Learner have reversed endein for the bit sequence
	(loop for mut in (reverse tmuts) do (setf (aref pmuts 0 (position mut *markers*)) 1))
	;; Compute the scores by the mat mult:
	(setf drug-p-vector (1xN-array->list (multiply-two-matrices pmuts *hypotheses-matrix*)))
	;; Choose 2 according to the probabilities (%%% Inefficient in various ways FFF)
	;; This is tricky if it can't pick 2 bcs they're all zeros except one (or none)!
	;; In that case we either choose randomly (if there are none) or duplicate (for 1).
	(loop with r = nil
	      with palist = (mapcar #'cons drug-p-vector *drug-dict*)
	      as k below 100 ;; Keep from looping out!
	      until (= 2 (length r)) ;; Early stopping rule if we find 2
	      do (pushnew (palist-choose palist) r)
	      finally (return (mapcar #'name->drug 
				      (case (length r)
					    (2 r)
					    (0 (n-random 2 *drugs*))
					    (1 (cons (first r) r)) ;; Duplicate
					    ))))
	)
    ;; No hyps: Just throw the dice!
    (n-random 2 *drugs*)
    )
  )

(defun 1xN-array->list (a)
  (loop for x below (second (array-dimensions a))
	collect (aref a 0 x)))

(defun first-n (n l)
  (loop for m below n as elt in l collect elt))

(defun sanity-check (tmuts)
  ;; We need to make a small number of errors here otherwise there's
  ;; no variance and the learner has nothing to learn from!
  (if (zerop (random 10)) (n-random 1 *drugs*)
    (loop for tm in tmuts
	  do (loop for drug in *drugs* 
		   do (if (equalp (car (drug-markers drug)) tm) 
			  (return-from sanity-check (list drug)))))))

(defun histscan (target-key px fn)
  (loop for (nil key . rest) in (patient-history px)
	when (eq target-key (cdr key))
	collect (funcall fn rest)))

(defun histscan-fast-5 (target-key px fn)
  (loop for (nil key . rest) in (patient-history px)
  	as limit below 5
	when (eq target-key (cdr key))
	collect (funcall fn rest)))
  
;;; =======================================================
;;; Utils

;;; Take a p(alist) like: ((0.5 . a) (0.25 . b) ...). First it
;;; normalizes it, a then it chooses a, b, ... etc in accord with the
;;; probabilities.

(defun palist-choose (palist)
  ; (print palist)
  (let* ((norm (float (reduce #'+ (mapcar #'car palist))))
	 (target (/ (random 1000) 1000.0)))
    (loop for (p . label) in palist
	  with sum = 0
	  do (incf sum (/ p norm))
	  (if (< target sum)
	      (return label)))))

(defun utime->standard (utime)
  (multiple-value-bind
   (second minute hour date month year day daylight-p zone)
   (decode-universal-time utime)
   (format nil "~4,'0d-~2,'0d-~2,'0dT00:00:00.000Z" year month date)))

(defun remdups (l &key (test #'eq))
  (loop for l+ on l
	unless (member (car l+) (cdr l+))
	collect (car l+)))

(defun dht (table &optional (n 10))
  (maphash #'(lambda (key value)
	       (when (zerop (decf n)) (return-from dht))
	       (format t "~s: ~s~%" key value)	       
	       )
	   table))

(defun record (px key &rest entry)
  (push `((:date . ,*date*) (:action . ,key) ,@entry) (patient-history px)))

(defun all-list-combinations (plists)
  (cond ((null plists) (list nil))
	(t (mapcan #'(lambda (elt) (insert-all-with-copies elt (all-list-combinations (cdr plists)))) (car plists)))))

(defun insert-all-with-copies (what intos)
  (loop for into in intos
	collect (cons what (copy-list into))))

(defun set-equal (a b)
  (and (equal (length a) (length b))
       (null (set-difference a b))
       (null (set-difference b a)) ;; I don't think you need both of these conditions since the lengths are the same
       ))

(defun n-random (n l) ;; With repeats!
  (loop for m below n collect (nth (random (length l)) l)))

;; Hypothesis reader: reads in hypotheses generated by ML algorithm as
;; "priors" for mutation-drug combos

(defun hreader ()
  (with-open-file 
   (i "inputs/hypotheses.lisp")
   (setf *hypotheses* (read i))
   ))

;;; Q&D Graphic of px's k dynamics -- mostly for debugging
;;;
;;; (defun plot-k-history (px)
;;;   (format *log* "~%")
;;;   (loop for (nil (nil . key) . rest) in (patient-history px)
;;; 	do (case key 
;;; 		 (:scan 
;;; 		  (loop for i below (cdar rest) do (format t "."))
;;; 			(format *log* "~a~%" (cdar rest)))
;;; 		 (:changed-tx (unless (member '(:result . :no-change) rest :test #'equal)
;;; 				(format *log* "Changed treatments: ~a~%" (cdr (assoc :result rest)))))
;;; 		 )))

;;; =======================================================
;;; Learning which drugs go with which molecular signatures.

(defvar *tms->des* (make-hash-table :test #'equal)) ;; *tumor-markers->drug-effects*
(defun molecular-signature (drug)
  (clrhash *tms->des*)
  (loop for (markers knew kold) in (drug-mrh drug)
	do (push (- knew kold) (gethash markers *tms->des*)))
  (sort 
   (loop for markers being the hash-keys of *tms->des*
	 using (hash-value kdels)
	 collect (list (if (cdr kdels) (float (STATISTICS:MEAN kdels)) (car kdels))
		       (if (cdr kdels) (STATISTICS:STANDARD-ERROR-OF-THE-MEAN kdels) 0.0)
		       markers))
   #'> :key #'car))

;;; =======================================================
;;; Genetic symbolification

;; Take the *markers* list as a bit position, and assign a binary key
;; to each combination of markers. nil -> 0, and anything else is a
;; unique order-independent number for a given marker combination.

(defun gene-code (markers)
  (setf markers (remdups markers)) ;; Make sure we're not double counting -- very confusing!
  (let ((code ( / (reduce #'+ (mapcar #'(lambda (marker) (expt 2 (length (member marker *markers*)))) markers)) 2)))
    (unless (gethash code *genetics*) (push markers (gethash code *genetics*)))
    code ))
	    
;;; =======================================================
;;; Main fns

(defun one-px (ncycles)
  ;;(format *mrhlog* "---	---	---	---	---	---	---	---~%")
  (catch 'done ;; "I see dead people."
    (loop with px = (genpx)
	  as *date* below ncycles
	  do 
	  (assess-and-treat px)
	  (update-k px)
	  finally (return px))))

(defun run-trial ()
  (format t "~%--- Running Trial ---~%")
  ;; Clear drug mrhs
  (loop for drug in *drugs* do (setf (drug-mrh drug) nil))
  (setf *pxs* nil)
  (with-open-file
   (o (format nil "../experiments/~a/pxs.json" *run-label*) :direction :output :if-exists :append :if-does-not-exist :create)
   (setq *foo*
	 (loop for np below *n-patients*
	       with final-n = (- *n-patients* 1)
	       do 
	       (setf *p* (one-px *nmonths*))
	       (push *p* *pxs*)
	       (case *dump-pxs-in-json?*
		     (:simple 
		      ;; This ugliness is because *p* is a thing of type
		      ;; string, but it wants to be a list for json
		      ;; encoding. UUU
		      (format o "~%~a~%" (json:encode-json-to-string
					  (read-from-string (format nil "~s" *p*))))
				)
		     (:rbstyle 
		      (if (= np 0) (format o "[~%"))
		      (dump-patient-rb-style *p* o)
		      (if (= np final-n) 
			  (format o "~%]~%")
			(format o ",~%")))
		     (t :no-output)
		     ))))
  )

(defun dump-patient-rb-style (p o)
  ;; This is a highly specific style that Robert wanted
  (format o 
	  "{caseId: \"~a\",
 source: \"Sim_version_~a\",
 cancer_group: \"Brain\",
 cancer_type: \"Glioblastoma\",
 birthdate: ~s,
 gender: \"~a\",
 diagnosis_date: ~s,
 case_history: 
   [
"
	  (patient-id p)
	  *sim-version*
	  (format nil "~a-01-01T00:00:00.000Z" (patient-birth-year p))
	  (patient-gender p)
	  (utime->standard (patient-dx-utime p))
	  )
  (let ((history (reverse (patient-history p))))
    (loop for elt in history
	  as future on (cdr history) ;; This runs forward for computing tx durations UUU
	  as k from 0 by 1
	  with end-k = (- (length (patient-history p)) 2) ;; Don't put a comma on this one
	  as type = (cdr (assoc :action elt))
	  as now-date = (cdr (assoc :date elt))
	  as event-date-utime = (+ (patient-dx-utime p) (* now-date 2678400)) 
	  as event-standard-date = (utime->standard event-date-utime)
	  do 
	  (format o "   {event_type: \"~a\", event_date: \"~a\"" type event-standard-date)
	  (case type
		(:dx :no-further-output)
		(:scan (format o ", kscore: \"~a\"" (cdr (assoc :k= elt))))
		(:CHANGED-TX 
		 (format o ", 
     drugs:[")
		 (loop for d in (cdr (assoc :result elt))
		       as k from 0 by 1
		       with end-k = (- (length (cdr (assoc :result elt))) 1)
		       do (format o "{name: \"~a\",start_date: \"~a\", end_date: \"~a\", adverse_severity: \"~a\"}" 
				  (drug-name d) event-standard-date 
				  (utime->standard ;; OMG! This is so hairy!
				   (+ event-date-utime 
				      (* 2678400 (- (redruging-future-date d future) now-date))))
				  (drug-adversity d))
		       (if (= k end-k) (format o "") (format o ",
            "))
		       )
		 (format o "]")
		 )
		)
	  (format o "}")
	  (if (= k end-k) (format o "~%") (format o ",~%"))
	  ))
  (format o "  ]
}"))

;;; Hairily figure the end date of a drug by walking forward in the
;;; history until the drug isn't in the list anylonger (or the history
;;; ends). UUU
(defun redruging-future-date (drug future)
  (loop for felt in future
	with last-date = (cdr (assoc :date felt)) ;; This starts here and gets updated in progress
	as date = (cdr (assoc :date felt))
	until (and (eq :changed-tx (cdr (assoc :action felt)))
		   (not (member drug (cdr (assoc :result felt))))) ;; I hope that this works by eq -- should!
	do (setf last-date date)
	finally (return last-date)))

(defvar *tms->effs* (make-hash-table :test #'equal)) ;; *tumor-markers->drug-effectivenesses* (or whatever) 

(defvar *tmr->drugs* (make-hash-table :test #'equal))
(defun report (&aux (ts+ (incf *fcounter*)))
  (format t "~%--- Reporting ---~%")
  ;; p of using each tx on each tumor type
  ;; Each k-summary (ks) in kss is: ((m1 m2 ...) klatest kprevious) -- set up here: [1]
  (clrhash *tmr->drugs*)
  (with-open-file 
   (o (format nil "../experiments/~a/~a-p-drugs-per-tumor.xls" *run-label* (get-universal-time))
      :direction :output :if-exists :append :if-does-not-exist :create)
   ;; (patient-tumor ,px) (first ks-5) (second ks-5) drug)
   (format o "Tumor	")
   (loop for drug in *drugs* do (format o "~a	" (drug-name drug)))
   (format o "~%")
   (loop for (tm k1 k2 d) in *mrhlog-as-list*
	 do (push (drug-name d) (gethash tm *tmr->drugs*)))
   (loop for tm being the hash-keys of *tmr->drugs*
	 using (hash-value drugs)
	 as maxp = -1
	 as maxd = nil
	 as l = (float (length drugs))
	 do 
	 (format o "~a	" (gene-code tm))
	 (loop for drug in *drugs* 
	       as drug-name = (drug-name drug)
	       as p = (/ (count drug-name drugs) l)
	       do (format o "	~a" p)
	       (if (> p maxp) (setf maxp p maxd drug-name)))
	 (format o "	~a~%" maxd))
   )
  ;; How patients' k-scores fared
  #+nil ;; !!!!!!!!!!! Temporarily out !!!!!!!!!!!!!
  (let* ((kss (mapcar #'(lambda (px) (histscan-fast-5 :scan px #'cdar)) *pxs*))
	 (kss! (reduce #'append kss)))
    (format *log* "~a	~a~%" (float (STATISTICS:MEAN kss!)) (STATISTICS:STANDARD-ERROR-OF-THE-MEAN kss!))
    (with-open-file 
     (o (format nil "../experiments/~a/~a-trajectories.xls" *run-label* ts+) :direction :output :if-exists :append :if-does-not-exist :create)
     (print *param-header* o) (terpri o)
     (loop for ks in kss
	   do (loop for i in (reverse ks)
		    do (format o "~a	" i)) ;; add in the ts here !!!!!!!!!!!!!!!!!!!!!
	   (format o "~%")))
    )
  ;; Drug targeting
  #+nil ;; !!!!!!!!!!! Temporarily out !!!!!!!!!!!!!
  (with-open-file 
   (o (format nil "../experiments/~a/~a-drug-targets.xls" *run-label* ts+) :direction :output :if-exists :append :if-does-not-exist :create)
   (print *param-header* o) (terpri o)
   (loop for drug in *drugs*
	 as molsig = (molecular-signature drug)
	 do 
	 ;; Note that this possibly confusingly replaces the mrh with its new value.
	 (setf (drug-mrh drug) molsig)
	 (pprint drug o)
	 ;; This is a slightly hairy inversion, putting the drugs to the tumor markers
	 (loop for (m se tmks) in molsig
	       do (push (list (drug-name drug) m se)
			(gethash tmks *tms->effs*)))
	 )
   (format o "~%=======================~%")
   (loop for tmks being the hash-keys of *tms->effs*
	 using (hash-value drug-effs)
	 do (pprint (list tmks (sort drug-effs #'> :key #'second)) o))
   ))

(defun run ()
  (format t "~%--- Starting Run ---~%")
  ;; ** If there are hypotheses, overrule the sanity check (which only matters for the first round anyway!)
  (if (eq *txrule* :sanity-check)
      (if *hypotheses-matrix*
	  (progn 
	  (setq *txrule* :hypotheses-therapy)
	  (format t "~% ************* OVER-RULING :SANITY-CHECK WITH :HYPOTHESES-THERAPY *************** ~%"))))
  (setf *drugs* 
	(case (print *txrule*)
	      (:sanity-check *not-confounded-drugs*)
	      (t (loop for drug in (rest (assoc :drugs *parameters*))
		       collect (make-drug :name (make-keyword (string-trim "+" (nth 0 drug))) 
					  :markers (loop for mk in (cdr (nth 1 drug)) collect (make-keyword mk)) 
					  :adversity (cdr (nth 2 drug)))))))
  (format t "Using *drugs*: ~s~%" *drugs*)
  (clrhash *genetics*)
  ;; Put the drugs in the genetics table
  (loop for drug in *drugs* do (gene-code (drug-markers drug)))
  ;; Core execution loop, open output files...
  (format t "~%--- Starting Core Execution Loops ---~%")
  (with-open-file
   (*mrhlog* (format nil "../experiments/~a/~a" *run-label* *sim-output*) 
	     :direction :output :if-exists :append :if-does-not-exist :create)
   ;;(format *mrhlog* "Patient	TimeStep	TumorCode	ComPos	Drug	DrugCode	KNow	KPrev	KDiff~%")
   (setf *mrhlog-as-list* nil)
   (with-open-file 
    (*log* (format nil "../experiments/~a/stats.xls" *run-label*) 
	   :direction :output :if-exists :append :if-does-not-exist :create)
    (loop for pn in *parameters* do (format *log* "~a	" pn))
    (format *log* "mean-k	k-stderr~%")
    ;; And go for it! ...
    (run-trial)
    (report)))
  ;; Dump genetics
  (with-open-file 
   (o (format nil "../experiments/~a/genetics.tsv" *run-label*) :direction :output :if-exists :append :if-does-not-exist :create)
   (loop for code being the hash-keys of *genetics*
	 using (hash-value markerss)
	 do (format o "~a	~{~a~^ ~}~%" code (car markerss))))
  )
  
;;; =======================================================
;;; Startup calls
(defun main () 
  (format t "~%--- Starting Sim! ---~%")
  (setf *random-state* (make-random-state t)) ;; randomize
  (untrace)
  (setf *dump-pxs-in-json?* nil) ;; :simple :rbstyle :none (or nil)
  (read-hyps-2 *classifier-output*) ;; holds the weight matrix produced by the classifier in *hypotheses-matrix*
  (setf *drug-dict* (loop for x below 11 collect (loop for drug in *drug-dict* when (= x (cdr drug)) do (return (car drug)))))
  (format t "*hypotheses-matrix* is: ~s~%" *hypotheses-matrix*)
  (run)
  (format t "~%--- Sim Done! ---~%")
  )
