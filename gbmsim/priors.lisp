(defun priors (drugs ds cds px)
	(let* ((possible-drugs (set-difference *drugs* ds)))
		(if (and possible_drugs (= (gene-code (nth 0 patient-hypotheses px)) (gene-code patient-tumor px)))
			(loop for (m* drug-sets) in *hypotheses* with match? nil
				if (= (m* gene-code) (px-tumor gene-code)) ;; test all mutations listed in *hypotheses* for px mutation
					do (setf drugs (nth 1 (pop drug-sets))) ;; set px drugs to the highest prior and remove it from drug-sets
					(return drugs) ;; breaks the loop
			)

;; the below works 
;; note that *hypotheses* is preserved 
;; perhaps we can assign each px a set of hypotheses for their particular mutation but this would perhaps be too messy
;; current desired tests:
;;		simply picking the best
;; 		randomly selecting according to "prior" distribution
(defun test (markers)
  ;;(init)
  (pprint *hypotheses*)
  (loop for (m* drug-sets) in *hypotheses*
    do (format t "~%")
    (format t "this mutation set: ~a, use these drugs:~%" m*)
    (loop for ds in drug-sets
      do (format t "   ~a~%" ds)
    ))
  (format t "BELOW MEEEEEEEEEEEEEEEEEEE ~%~%")
  (loop for (muts drug-sets) in *hypotheses* 
    ;;initially (let* ((mutations '())))
    do (loop for (m* indicator) in muts 
    	;;initially (let ((mutations '()) (contra-mutations '())))
    	with mutations = ()
    	with contra-mutations = ()
    	do (format t "This ~a mutation set can be broken down as:~%" muts)
    	   (format t "~a with an indicator of ~a~%" m* indicator)
    	   (case indicator
    	   		(1 (push m* mutations))
    	   		(-1 (push m* contra-mutations))
    	   	)


    	   (format t "~a~%" mutations)
    	   (format t "contra-mutations ~a~%" contra-mutations)
    	) 
  ))

;;do (let* ((drogas (loop for setn in drug-sets maximizing (nth
;;	)))))

;;;;;;; currently nonfunctional
;; what if I looped through drug-sets to produce the separate list of accuracies and drug combos 

;;if (= (gene-code m*) (gene-code markers)) ;; if we fail, we could move on to some other thing, such as random binary cocktail or whatever (untried or something)
  ;;    do (loop for setn in drug-sets 
    ;;  		maximizing (nth 0 setn) into drogas
      ;;		do (print (nth 0 setn)) 
      	;; )
      ;;finally (format t "These are the drogas: ~a ~%" drogas)
      ;;do (print (nth 0 drug-sets))
      ;;   (setf drug (nth 1 (pop drug-sets))) ;; this returns the first drug set, i.e. the one equal to 0.5! tight 
      ;;   (print drug)
      ;;   (print drug-sets)
      ;;   (let* ((match? t)))
      ;;(return drug)
      ;; (format t "I should not appear") ;; if uncommented it will not appear... for demo purposes only 
    ;;finally  (unless match? (print m*));; if we fail every single one we need to come up with some sort of option, this currently just executes for all the misses