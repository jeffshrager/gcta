;; (load (compile-file "simple.lisp"))

;;; As run6 ... 

;;; (The runcode MUST be manually entered at all time -- resist the urge to default it out!!)
(defun run7 (runcode &key (ntypes 3) (mvols 10) (cycles 1000) (detection 10000)
		     (undetection 10) (metastasis (* detection 0.8)) (p-mutation-per-100 40) (tailoff 0.99))
  (with-open-file
   (o (format nil "results/run~a_~a_ty~a_vo~a_cy~a_de~a_ud~a_me~a_pmp~a_tail~a.xls" runcode (get-universal-time) ntypes mvols cycles detection undetection metastasis p-mutation-per-100 tailoff)
      :direction :output :if-exists :supersede)
   (format o "PS	NMets	LargestMet	TotVol	")
   (loop for m from 1 to mvols do (format o "~avol	~atype	" m m))
   (loop for n from 1 to ntypes do (format o "tx~a	" n)) ;; Types and treatments are the same
   (format o "~%")
   (loop for i below cycles
	 ;; Each volume is a pair: (ncells . type)
	 with volumes = (let ((vols (make-array mvols :initial-element nil)))
			  (setf (aref vols 0) (cons 1 1))
			  vols)
	 ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	 ;; Treatments and types are 1:1 at the moment
	 with txs = (make-array ntypes :initial-element 1.0)
	 do
	 ;; Termination check
	 (when (zerop (loop for (cells . type) across volumes when cells sum cells))
	   (format t "~%************* RUN TERMINATED @ CYCLE ~a -- NO EVIDENCE OF DISEASE! *************~%" i)
	   (return :done))
	 ;; Check for metastatsis, and metastasize (spread) as necessary/possible -- and maybe mutate
	 (loop for vol across volumes
	       as (cells . type) = vol
	       if (and cells (> cells metastasis))
	       do
	       (let* ((possible-locations (loop for entry across volumes as p from 0 by 1 unless entry collect p))
		      (selected-location (when possible-locations (nth (random (length possible-locations)) possible-locations))))
		 (when selected-location
		   (setf (aref volumes selected-location)
			 (cons 1 (if (< (random 100) p-mutation-per-100)
				     (random ntypes)
				   type))))))
	 ;; Detection and treatment -- increase treatment for the
	 ;; relevant type; if the number of cells is under the
	 ;; UNdetection limit, reduce treatment (but not above 1.0 -- no treatment)
	 (loop for (cells . type) across volumes
	       as p from 0 by 1
	       do
	       (when cells
		 (if (> cells detection) (setf (aref txs type) (* (aref txs type) 0.95))
		   (if (< cells undetection) (setf (aref txs type) (min (/ (aref txs type) tailoff) 1.0))))
		 (unless (<= detection metastasis) ;;Don't let detection go below metastatis
		   (setf detection (* 0.9 detection))) ;; Also increase vigalance!
		 ))
	 ;; Finally update volumes, less tx multipliers
	 (loop for vol across volumes
	       as (cells . type) = vol
	       when cells
	       do (let ((newvol (* (+ cells cells) (aref txs type))))
		    ;; Threshold at the bottom to zero
		    (setf (car vol) (if (< newvol 1.0) 0.0 newvol)))) 
	 ;; And repoert!
	 (let* ((allvols (loop for (cells . nil) across volumes collect cells))
	       (justmets (delete nil allvols)))
	   (format o "~a	" (performance-score allvols detection))
	   (format o "~a	" (length justmets))     ; NMets
	   (format o "~a	" (apply #'max justmets)) ;; Largest met
	   (format o "~a	" (reduce #'+ justmets)) ;; Total volume
	   )
	 (loop for (cells . type) across volumes do (format o "~a	~a	" (or cells "") (or type "")))
	 (loop for x across txs do (format o "~a	" (or x "")))
	 (format o "~%")
	 )))

(defun performance-score (mets detection)
  (loop for met in mets
	with ps = 100
	as m from 1 by 1
	as important? = (oddp m)
	if (and met important? (> met detection))
	do (setf ps (* ps 0.8))
	finally (return ps)))

(run7 7 :ntypes 5 :mvols 25 :cycles 1000 :detection 10000 :metastasis (* 10000 0.8) :p-mutation-per-100 40 :tailoff 0.998)

#|

;;; As run5, but treatments tail off when tumors become
;;; undetectable. The tailoff is a number that, if the cells are under
;;; detection threhsold, we DIVIDE by the tail constant (to a max of
;;; 1.0 -- recall that 1.0 is no treatment bcs it's a multiplier!) All
;;; this is sort of upside down and so hard to understand. (And, I
;;; predict, eventually hard to maintain.)

;;; Nb. this isn't really the way treatment works...I guess that sort
;;; of goes w/o saying, but I mean that treatment isn't really based
;;; on the tumor volume, but usually on the first derivative of the
;;; volume -- that is, if the tumor is increasing in size (or number
;;; of mets) then you do aggressive treatment. But that would require
;;; keeping met size histories...maybe in a future version! FFF

;;; (The runcode MUST be manually entered at all time -- resist the urge to default it out!!)
(defun run6 (runcode &key (ntypes 3) (mvols 10) (cycles 1000) (detection 10000)
		     (undetection 10) (metastasis (* detection 0.8)) (p-mutation-per-100 40) (tailoff 0.99))
  (with-open-file
   (o (format nil "results/run~a_~a_ty~a_vo~a_cy~a_de~a_ud~a_me~a_pmp~a_tail~a.xls" runcode (get-universal-time) ntypes mvols cycles detection undetection metastasis p-mutation-per-100 tailoff)
      :direction :output :if-exists :supersede)
   (format o "PS	NMets	LargestMet	TotVol	")
   (loop for m from 1 to mvols do (format o "~avol	~atype	" m m))
   (loop for n from 1 to ntypes do (format o "tx~a	" n)) ;; Types and treatments are the same
   (format o "~%")
   (loop for i below cycles
	 ;; Each volume is a pair: (ncells . type)
	 with volumes = (let ((vols (make-array mvols :initial-element nil)))
			  (setf (aref vols 0) (cons 1 1))
			  vols)
	 ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	 ;; Treatments and types are 1:1 at the moment
	 with txs = (make-array ntypes :initial-element 1.0)
	 do
	 ;; Termination check
	 (when (zerop (loop for (cells . type) across volumes when cells sum cells))
	   (format t "~%************* RUN TERMINATED @ CYCLE ~a -- NO EVIDENCE OF DISEASE! *************~%" i)
	   (return :done))
	 ;; Check for metastatsis, and metastasize (spread) as necessary/possible -- and maybe mutate
	 (loop for vol across volumes
	       as (cells . type) = vol
	       if (and cells (> cells metastasis))
	       do
	       (let* ((possible-locations (loop for entry across volumes as p from 0 by 1 unless entry collect p))
		      (selected-location (when possible-locations (nth (random (length possible-locations)) possible-locations))))
		 (when selected-location
		   (setf (aref volumes selected-location)
			 (cons 1 (if (< (random 100) p-mutation-per-100)
				     (random ntypes)
				   type))))))
	 ;; Detection and treatment -- increase treatment for the
	 ;; relevant type; if the number of cells is under the
	 ;; UNdetection limit, reduce treatment (but not above 1.0 -- no treatment)
	 (loop for (cells . type) across volumes
	       as p from 0 by 1
	       do
	       (when cells
		 (if (> cells detection) (setf (aref txs type) (* (aref txs type) 0.95))
		   (if (< cells undetection) (setf (aref txs type) (min (/ (aref txs type) tailoff) 1.0))))
		 (unless (<= detection metastasis) ;;Don't let detection go below metastatis
		   (setf detection (* 0.9 detection))) ;; Also increase vigalance!
		 ))
	 ;; Finally update volumes, less tx multipliers
	 (loop for vol across volumes
	       as (cells . type) = vol
	       when cells
	       do (let ((newvol (* (+ cells cells) (aref txs type))))
		    ;; Threshold at the bottom to zero
		    (setf (car vol) (if (< newvol 1.0) 0.0 newvol)))) 
	 ;; And repoert!
	 (let* ((allvols (loop for (cells . nil) across volumes collect cells))
	       (justmets (delete nil allvols)))
	   (format o "~a	" (performance-score allvols detection))
	   (format o "~a	" (length justmets))     ; NMets
	   (format o "~a	" (apply #'max justmets)) ;; Largest met
	   (format o "~a	" (reduce #'+ justmets)) ;; Total volume
	   )
	 (loop for (cells . type) across volumes do (format o "~a	~a	" (or cells "") (or type "")))
	 (loop for x across txs do (format o "~a	" (or x "")))
	 (format o "~%")
	 )))

(defun performance-score (mets detection)
  (loop for met in mets
	with ps = 100
	as m from 1 by 1
	as important? = (oddp m)
	if (and met important? (> met detection))
	do (setf ps (* ps 0.8))
	finally (return ps)))

(run6 6 :ntypes 5 :mvols 25 :cycles 1000 :detection 10000 :metastasis (* 10000 0.8) :p-mutation-per-100 40 :tailoff 0.998)

;;; As run4, but includes performance score computation, where simply
;;; even nth volumes hurt a little, whereas odd ones hurt a lot more.

(defun run5 (&key (ntypes 3) (mvols 10) (cycles 1000) (detection 10000) (metastasis (* detection 0.8)) (p-mutation-per-100 40))
  (with-open-file
   (o (format nil "results/~a_ty~a_vo~a_cy~a_de~a_me~a_pmp~a.xls" (get-universal-time) ntypes mvols cycles detection metastasis p-mutation-per-100)
      :direction :output :if-exists :supersede)
   (loop for m from 1 to mvols do (format o "~avol	~atype	" m m))
   (format o "Largest Met	Tot Vol	N Mets	PS	")
   (loop for n from 1 to ntypes do (format o "tx~a	" n)) ;; Types and treatments are the same
   (format o "~%")
   (loop for i below cycles
	 ;; Each volume is a pair: (ncells . type)
	 with volumes = (let ((vols (make-array mvols :initial-element nil)))
			  (setf (aref vols 0) (cons 1 1))
			  vols)
	 ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	 ;; Treatments and types are 1:1 at the moment
	 with txs = (make-array ntypes :initial-element 1.0)
	 do
	 ;; Termination check
	 (when (zerop (loop for (cells . type) across volumes when cells sum cells))
	   (format t "~%************* RUN TERMINATED @ CYCLE ~a -- NO EVIDENCE OF DISEASE! *************~%" i)
	   (return :done))
	 ;; Check for metastatsis, and metastasize (spread) as necessary/possible -- and maybe mutate
	 (loop for vol across volumes
	       as (cells . type) = vol
	       if (and cells (> cells metastasis))
	       do
	       (let* ((possible-locations (loop for entry across volumes as p from 0 by 1 unless entry collect p))
		      (selected-location (when possible-locations (nth (random (length possible-locations)) possible-locations))))
		 (when selected-location
		   (setf (aref volumes selected-location)
			 (cons 1 (if (< (random 100) p-mutation-per-100)
				     (random ntypes)
				   type))))))
	 ;; Detection and treatment -- increase treatment for the relevant type
	 (loop for (cells . type) across volumes
	       as p from 0 by 1
	       when (and cells (> cells detection))
	       do
	       (setf (aref txs type) (* (aref txs type) 0.95))
	       (unless (<= detection metastasis) ;;Don't let detection go below metastatis
		 (setf detection (* 0.9 detection))) ;; Also increase vigalance!
	       )
	 ;; Finally update volumes, less tx multipliers
	 (loop for vol across volumes
	       as (cells . type) = vol
	       when cells
	       do (setf (car vol) (* (+ cells cells) (aref txs type))))
	 ;; And repoert!
	 (loop for (cells . type) across volumes do (format o "~a	~a	" (or cells "") (or type "")))
	 (let* ((allvols (loop for (cells . nil) across volumes collect cells))
	       (justmets (delete nil allvols)))
	   (format o "~a	" (apply #'max justmets)) ;; Largest met
	   (format o "~a	" (reduce #'+ justmets)) ;; Total volume
	   (format o "~a	" (length justmets))     ; NMets
	   (format o "~a	" (performance-score allvols detection))
	   )
	 (loop for x across txs do (format o "~a	" (or x "")))
	 (format o "~%")
	 )))

(defun performance-score (mets detection)
  (loop for met in mets
	with ps = 100
	as m from 1 by 1
	as important? = (oddp m)
	if (and met important? (> met detection))
	do (setf ps (* ps 0.8))
	finally (return ps)))

(run5 :ntypes 5 :mvols 25 :cycles 1000 :detection 10000 :metastasis (* 10000 0.8) :p-mutation-per-100 40)

;;; Same as run3, but goes to specific results files.

(defun run4 (&key (ntypes 3) (mvols 10) (cycles 1000) (detection 10000) (metastasis (* detection 0.8)) (p-mutation-per-100 40))
  (with-open-file
   (o (format nil "results/ty~a_vo~a_cy~a_de~a_me~a_pmp~a.xls" ntypes mvols cycles detection metastasis p-mutation-per-100)
      :direction :output :if-exists :supersede)
   (loop for m from 1 to mvols do (format o "v~a	t~a	" m m))
   (format o "totvol	")
   (loop for n from 1 to ntypes do (format o "tx~a	" n)) ;; Types and treatments are the same
   (format o "~%")
   (loop for i below cycles
	 ;; Each volume is a pair: (ncells . type)
	 with volumes = (let ((vols (make-array mvols :initial-element nil)))
			  (setf (aref vols 0) (cons 1 1))
			  vols)
	 ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	 ;; Treatments and types are 1:1 at the moment
	 with txs = (make-array ntypes :initial-element 1.0)
	 do
	 ;; Check for metastatsis, and metastasize (spread) as necessary/possible -- and maybe mutate
	 (loop for vol across volumes
	       as (cells . type) = vol
	       if (and cells (> cells metastasis))
	       do
	       (let* ((possible-locations (loop for entry across volumes as p from 0 by 1 unless entry collect p))
		      (selected-location (when possible-locations (nth (random (length possible-locations)) possible-locations))))
		 (when selected-location
		   (setf (aref volumes selected-location)
			 (cons 1 (if (< (random 100) p-mutation-per-100)
				     (random ntypes)
				   type))))))
	 ;; Detection and treatment -- increase treatment for the relevant type
	 (loop for (cells . type) across volumes
	       as p from 0 by 1
	       when (and cells (> cells detection))
	       do
	       (setf (aref txs type) (* (aref txs type) 0.95))
	       (setf detection (* 0.9 detection)) ;; Also increase vigalance!
	       )
	 ;; Finally update volumes, less tx multipliers
	 (loop for vol across volumes
	       as (cells . type) = vol
	       when cells
	       do (setf (car vol) (* (+ cells cells) (aref txs type))))
	 ;; And repoert!
	 (loop for (cells . type) across volumes do (format o "~a	~a	" (or cells "") (or type "")))
	 (format o "~a	" (loop for (cells . type) across volumes when cells sum cells))
	 (loop for x across txs do (format o "~a	" (or x "")))
	 (format o "~%")
	 )))

;(run4 :ntypes 5 :mvols 25 :cycles 1000 :detection 10000 :metastasis (* 10000 0.8) :p-mutation-per-100 40)

;;; Similar to run2, but with n types of tumor, and m volumes, and a
;;; random volume selection and "implantation" probability. We assume
;;; that the correct tx for the tumor type is adopted, at least
;;; initially, for each new tumor implantation. (As tumors do not, in
;;; this version, evolve in place, this is a good assumption. Types
;;; should be "variants", but "v" here means vol.)

(defun run3 (&key (ntypes 3) (mvols 10) (cycles 1000) (detection 10000) (metastasis (* detection 0.8)) (p-mutation-per-100 40))
  (with-open-file
   (o "out.xls" :direction :output :if-exists :supersede)
   (loop for m from 1 to mvols do (format o "v~a	t~a	" m m))
   (format o "totvol	")
   (loop for n from 1 to ntypes do (format o "tx~a	" n)) ;; Types and treatments are the same
   (loop for i below cycles
	 ;; Each volume is a pair: (ncells . type)
	 with volumes = (let ((vols (make-array mvols :initial-element nil)))
			  (setf (aref vols 0) (cons 1 1))
			  vols)
	 ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	 ;; Treatments and types are 1:1 at the moment
	 with txs = (make-array ntypes :initial-element 1.0)
	 do
	 ;; Check for metastatsis, and metastasize (spread) as necessary/possible -- and maybe mutate
	 (loop for vol across volumes
	       as (cells . type) = vol
	       if (and cells (> cells metastasis))
	       do
	       (let* ((possible-locations (loop for entry across volumes as p from 0 by 1 unless entry collect p))
		      (selected-location (when possible-locations (nth (random (length possible-locations)) possible-locations))))
		 (when selected-location
		   (setf (aref volumes selected-location)
			 (cons 1 (if (< (random 100) p-mutation-per-100)
				     (random ntypes)
				   type))))))
	 ;; Detection and treatment -- increase treatment for the relevant type
	 (loop for (cells . type) across volumes
	       as p from 0 by 1
	       when (and cells (> cells detection))
	       do
	       (setf (aref txs type) (* (aref txs type) 0.95))
	       (setf detection (* 0.9 detection)) ;; Also increase vigalance!
	       )
	 ;; Finally update volumes, less tx multipliers
	 (loop for vol across volumes
	       as (cells . type) = vol
	       when cells
	       do (setf (car vol) (* (+ cells cells) (aref txs type))))
	 ;; And repoert!
	 (loop for (cells . type) across volumes do (format o "~a	~a	" (or cells "") (or type "")))
	 (format o "~a	" (loop for (cells . type) across volumes when cells sum cells))
	 (loop for x across txs do (format o "~a	" (or x "")))
	 (format o "~%")
	 )))

;;; Similar to run1, but here the volumes are considered as parts of
;;; the body, and you don't get any tumor at all in the empty parts
;;; until one of the tumors reaches detection threshold (equally
;;; interpretable as metastisis threshold). This creates a somewhat
;;; interesting and surprising "M"-shaped pattern of total tumor load
;;; because of the delayed emergence of the distant mets, and at the
;;; same time the greatly reduced (by that time) detection threshold.

(defun run2 ()
    (with-open-file
     (o "simple2.xls" :direction :output :if-exists :supersede)
     (format o "v1	v2	v3	vtot	tx1	tx2	tx3~%")
     (loop for i below 1000
	   with volumes = (list 1 0 0)
	   with txs = (list 1.0 1.0 1.0)  ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	   with detection = 10000
	   as high = (apply #'max volumes)
	   do (loop for vol in volumes
		    as p from 0 by 1
		    if (> vol detection) ;; when this volume is above detection
		    do
		    ;; Metastisize ... assuming there's anywhere to go!
		    (let ((p (position 0.0 volumes))) (when p (incf (nth p volumes))))
		    (setf (nth p txs) (* (nth p txs) 0.95)) ;; add treatment
		    (setf detection (* detection 0.95)) ;; and reduce the detection threshold
		    )
	   (loop for vol in volumes
		 for tx in txs
		 as p from 0 by 1
		 do (setf (nth p volumes) (* tx (+ vol vol)))) ;; Vol doubles times treatment reducer
	   (loop for v in volumes do (format o "~a	" v))
	   (format o "~a	" (reduce #'+ volumes))
	   (loop for x in txs do (format o "~a	" x))
	   (format o "~%")
	   )))

;;; Just keep doubling each volume independently, and when each hits
;;; the detection threshold, again independently, treat it,
;;; independently, and keep pressing harder and harder on the
;;; treatments (independently!) until they get below detection
;;; threshold. Also, the detection theshold is reduced each time tx is
;;; increased, simulating increased surveillance vigilance.

(defun run1 ()
    (with-open-file
     (o "simple1.xls" :direction :output :if-exists :supersede)
     (format o "v1	v2	v3	vtot	tx1	tx2	tx3~%")
     (loop for i below 1000
	   with volumes = (list 10 2 1)
	   with txs = (list 1.0 1.0 1.0) ;; Note NO EFFECT = 1.0 (this is multiplied in for treatment)
	   with detection = 10000
	   as high = (apply #'max volumes)
	   do (loop for vol in volumes
		    as p from 0 by 1
		    if (> vol detection) ;; when this volume is above detection
		    do
		    (setf (nth p txs) (* (nth p txs) 0.95)) ;; add treatment
		    (setf detection (* detection 0.95)) ;; and reduce the detection threshold
		    )
	   (loop for vol in volumes
		 for tx in txs
		 as p from 0 by 1
		 do (setf (nth p volumes) (* tx (+ vol vol)))) ;; Vol doubles times treatment reducer
	   (loop for v in volumes do (format o "~a	" v))
	   (format o "~a	" (reduce #'+ volumes))
	   (loop for x in txs do (format o "~a	" x))
	   (format o "~%")
	   )))
		 
		 
|#
