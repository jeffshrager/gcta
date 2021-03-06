;;; (load (compile-file "reshape.lisp"))

;;; (loop for k being the hash-keys of *t* using (hash-value v) do (print (list k v)))
;;; Just a hack to reshuffle results

(defparameter
  *vtbs*
  '(
    ;; PtNo Concordance AlgDistMeanVSConsensus AlgDistStDev Consensus Judges
    (41 1.343 1.102 0.588 "[11, 24, 26]" ( "[11, 26, 55]" "[11, 26, 55]" "[11, 24, 55]" "[55, 39, 24]" "[11, 55, 39]" "[24, 28, 26]" "[26, 11, 39]" "[41, 26, 24]" "[24, 41, 55]" "[24, 41, 11]"))
    (76 1.280 0.931 0.663 "[28, 56, 43]" ( "[56, 47, 38]" "[38, 43, 45]" "[28, 56, 26]" "[56, 28, 26]" "[26, 47, 6]" "[28, 26, 45]" "[28, 43, 47]" "[28, 56, 43]" "[28, 43, 6]" "[38, 43, 9]"))
    (72 1.253 1.031 0.389 "[26, 40, 22]" ( "[40, 22, 26]" "[26, 2, 52]" "[40, 26, 22]" "[2, 22, 52]" "[40, 26, 22]" "[52, 26, 22]" "[22, 40, 26]" "[40, 4, 52]" "[26, 2, 52]" "[26, 2, 52]"))
    (63 1.224 1.095 0.415 "[55, 26, 35]" ( "[35, 26, 40]" "[35, 55, 12]" "[55, 35, 40]" "[26, 55, 35]" "[26, 35, 40]" "[40, 55, 7]" "[55, 35, 26]" "[40, 55, 26]" "[26, 40, 55]" "[26, 55, 14]"))
    (46 1.189 0.834 0.597 "[19, 29, 5]" ( "[41, 29, 5]" "[19, 29, 60]" "[19, 60, 5]" "[5, 60, 41]" "[29, 19, 60]" "[19, 29, 60]" "[19, 29, 41]" "[19, 60, 5]" "[41, 19, 5]" "[29, 5, 19]"))
    (70 1.123 0.822 0.549 "[35, 30, 38]" ( "[35, 19, 53]" "[30, 19, 38]" "[30, 35, 19]" "[35, 30, 38]" "[38, 19, 53]" "[35, 38, 30]" "[35, 19, 53]" "[38, 30, 35]" "[35, 30, 53]" "[30, 38, 35]"))
    (24 1.116 0.772 0.562 "[11, 51, 40]" ( "[11, 40, 55]" "[11, 51, 40]" "[11, 40, 51]" "[55, 11, 40]" "[11, 51, 55]" "[51, 11, 40]" "[40, 51, 11]" "[11, 40, 51]" "[55, 11, 40]" "[51, 55, 11]"))
    (23 1.094 0.875 0.506 "[24, 18, 28]" ( "[24, 18, 26]" "[18, 24, 55]" "[24, 26, 18]" "[18, 28, 41]" "[26, 24, 28]" "[18, 28, 55]" "[24, 26, 28]" "[18, 24, 41]" "[24, 18, 55]" "[28, 24, 18]"))
    (50 1.091 0.880 0.412 "[19, 6, 20]" ( "[6, 19, 20]" "[20, 19, 6]" "[19, 45, 22]" "[19, 45, 22]" "[6, 20, 19]" "[19, 20, 6]" "[6, 45, 22]" "[20, 45, 22]" "[19, 20, 6]" "[6, 19, 20]"))
    (54 0.956 0.797 0.529 "[29, 20, 40]" ( "[40, 29, 20]" "[29, 20, 40]" "[20, 40, 29]" "[29, 20, 40]" "[29, 20, 45]" "[40, 29, 20]" "[20, 29, 45]" "[20, 29, 45]" "[40, 29, 20]" "[40, 20, 29]"))
    (1 0.889 0.573 0.510 "[10, 32, 9]" ( "[10, 9, 48]" "[10, 43, 50]" "[32, 4, 10]" "[10, 48, 9]" "[10, 43, 4]" "[10, 32, 48]" "[9, 32, 48]" "[10, 32, 43]" "[10, 32, 9]" "[10, 32, 48]"))
    (40 0.865 0.575 0.559 "[1, 24, 27]" ( "[24, 1, 27]" "[1, 24, 27]" "[27, 24, 3]" "[1, 24, 52]" "[1, 24, 52]" "[1, 27, 24]" "[1, 27, 24]" "[1, 24, 27]" "[27, 24, 1]" "[24, 27, 1]"))
    (34 0.835 0.597 0.604 "[3, 47, 59]" ( "[3, 59, 47]" "[3, 47, 29]" "[59, 3, 47]" "[59, 47, 29]" "[3, 47, 59]" "[59, 3, 47]" "[3, 47, 32]" "[3, 47, 32]" "[3, 47, 59]" "[59, 47, 3]"))
    (12 0.828 0.569 0.555 "[17, 9, 46]" ( "[17, 9, 24]" "[17, 9, 46]" "[17, 9, 46]" "[9, 17, 46]" "[17, 24, 46]" "[17, 9, 46]" "[9, 24, 46]" "[9, 24, 46]" "[9, 17, 46]" "[17, 24, 46]"))
    (56 0.828 0.594 0.541 "[20, 15, 39]" ( "[15, 39, 19]" "[20, 15, 21]" "[15, 39, 44]" "[20, 39, 41]" "[20, 15, 19]" "[15, 20, 39]" "[15, 20, 19]" "[20, 39, 53]" "[20, 15, 39]" "[20, 15, 39]"))
    (31 0.822 0.600 0.582 "[20, 40, 53]" ( "[53, 40, 20]" "[20, 40, 53]" "[53, 40, 20]" "[20, 40, 53]" "[20, 40, 53]" "[20, 53, 40]" "[40, 20, 53]" "[40, 53, 20]" "[20, 40, 53]" "[53, 20, 40]"))
    (79 0.818 0.663 0.431 "[48, 23, 41]" ( "[48, 41, 22]" "[48, 22, 7]" "[48, 41, 23]" "[48, 22, 23]" "[23, 41, 7]" "[48, 22, 23]" "[23, 41, 7]" "[48, 22, 41]" "[48, 23, 7]" "[48, 41, 23]"))
    (4 0.809 0.544 0.514 "[33, 36, 20]" ( "[33, 36, 20]" "[36, 20, 45]" "[33, 36, 6]" "[33, 36, 6]" "[33, 20, 45]" "[36, 33, 6]" "[33, 36, 29]" "[36, 7, 45]" "[33, 36, 29]" "[33, 7, 45]"))
    (91 0.809 0.625 0.506 "[45, 51, 40]" ( "[45, 40, 6]" "[45, 40, 6]" "[51, 45, 40]" "[51, 40, 6]" "[51, 45, 44]" "[51, 40, 6]" "[45, 51, 44]" "[45, 51, 44]" "[45, 51, 44]" "[45, 51, 44]"))
    (19 0.801 0.472 0.616 "[16, 35, 58]" ( "[16, 35, 58]" "[16, 58, 35]" "[16, 35, 58]" "[7, 35, 58]" "[35, 58, 16]" "[16, 35, 58]" "[16, 35, 58]" "[58, 16, 35]" "[16, 35, 58]" "[16, 58, 35]"))
    (14 0.787 0.530 0.417 "[1, 8, 54]" ( "[1, 8, 39]" "[8, 1, 3]" "[8, 1, 3]" "[1, 45, 31]" "[1, 8, 3]" "[1, 8, 54]" "[8, 45, 30]" "[1, 54, 31]" "[1, 8, 3]" "[1, 54, 39]"))
    (60 0.772 0.631 0.442 "[26, 28, 55]" ( "[26, 28, 55]" "[28, 26, 37]" "[28, 26, 37]" "[26, 28, 37]" "[26, 55, 2]" "[26, 55, 2]" "[28, 55, 7]" "[28, 26, 55]" "[28, 26, 37]" "[26, 28, 39]"))
    (39 0.753 0.691 0.440 "[4, 7, 38]" ( "[7, 4, 2]" "[7, 4, 2]" "[7, 38, 55]" "[4, 7, 38]" "[7, 4, 2]" "[4, 38, 55]" "[4, 7, 38]" "[7, 4, 2]" "[4, 38, 55]" "[7, 4, 38]"))
    (82 0.753 0.471 0.467 "[9, 6, 54]" ( "[9, 54, 37]" "[6, 9, 29]" "[9, 54, 6]" "[9, 6, 54]" "[9, 6, 56]" "[6, 54, 29]" "[9, 6, 54]" "[9, 56, 39]" "[9, 6, 54]" "[9, 56, 36]"))
    (3 0.750 0.581 0.469 "[9, 45, 37]" ( "[9, 45, 37]" "[9, 37, 45]" "[9, 37, 45]" "[9, 37, 45]" "[45, 37, 9]" "[45, 9, 37]" "[9, 45, 37]" "[45, 37, 9]" "[9, 37, 45]" "[37, 9, 45]"))
    (78 0.720 0.620 0.489 "[53, 15, 9]" ( "[15, 53, 38]" "[9, 53, 38]" "[53, 15, 46]" "[53, 15, 9]" "[15, 53, 37]" "[53, 15, 46]" "[15, 53, 37]" "[15, 53, 46]" "[53, 15, 9]" "[15, 53, 38]"))
    (22 0.688 0.426 0.279 "[8, 3, 31]" ( "[8, 3, 1]" "[3, 8, 31]" "[8, 31, 59]" "[8, 3, 1]" "[8, 31, 41]" "[8, 3, 1]" "[8, 27, 59]" "[8, 31, 51]" "[8, 3, 31]" "[3, 8, 31]"))
    (92 0.685 0.517 0.471 "[9, 26, 22]" ( "[9, 26, 56]" "[26, 9, 10]" "[26, 9, 22]" "[26, 9, 10]" "[9, 26, 22]" "[9, 26, 56]" "[26, 22, 10]" "[9, 26, 56]" "[9, 22, 56]" "[9, 26, 56]"))
    (71 0.681 0.433 0.298 "[48, 9, 50]" ( "[48, 50, 38]" "[48, 9, 50]" "[48, 7, 50]" "[48, 50, 38]" "[9, 48, 38]" "[9, 48, 50]" "[48, 50, 38]" "[48, 9, 23]" "[48, 9, 38]" "[48, 9, 23]"))
    (53 0.678 0.457 0.281 "[27, 35, 15]" ( "[27, 15, 56]" "[35, 27, 16]" "[27, 35, 56]" "[27, 15, 56]" "[27, 35, 16]" "[27, 43, 17]" "[27, 15, 56]" "[27, 35, 16]" "[27, 35, 17]" "[35, 27, 15]"))
    (66 0.666 0.441 0.492 "[13, 5, 50]" ( "[5, 28, 50]" "[13, 50, 36]" "[13, 5, 50]" "[13, 5, 47]" "[5, 13, 50]" "[13, 5, 50]" "[13, 5, 50]" "[13, 5, 50]" "[5, 13, 50]" "[5, 13, 50]"))
    (28 0.633 0.463 0.298 "[44, 45, 1]" ( "[44, 45, 40]" "[45, 44, 40]" "[45, 44, 53]" "[44, 1, 36]" "[44, 1, 36]" "[44, 45, 53]" "[44, 45, 53]" "[44, 45, 40]" "[44, 1, 40]" "[44, 1, 36]"))
    (83 0.625 0.452 0.419 "[2, 23, 14]" ( "[2, 4, 14]" "[2, 14, 4]" "[2, 4, 14]" "[2, 23, 52]" "[2, 4, 14]" "[2, 23, 52]" "[2, 23, 52]" "[14, 23, 52]" "[2, 23, 52]" "[2, 23, 52]"))
    (73 0.618 0.450 0.393 "[56, 31, 30]" ( "[56, 31, 30]" "[56, 31, 30]" "[31, 56, 35]" "[31, 56, 16]" "[56, 31, 35]" "[31, 56, 30]" "[56, 27, 30]" "[56, 31, 27]" "[31, 56, 16]" "[56, 31, 35]"))
    (67 0.601 0.491 0.355 "[30, 19, 23]" ( "[19, 30, 23]" "[19, 30, 56]" "[30, 19, 56]" "[19, 30, 23]" "[30, 23, 19]" "[30, 19, 56]" "[19, 30, 23]" "[19, 30, 56]" "[30, 19, 23]" "[30, 19, 5]"))
    (89 0.592 0.403 0.400 "[1, 30, 56]" ( "[1, 30, 56]" "[1, 56, 30]" "[1, 60, 54]" "[1, 30, 60]" "[1, 56, 30]" "[1, 30, 54]" "[56, 1, 54]" "[1, 60, 3]" "[1, 30, 54]" "[1, 30, 3]"))
    (64 0.592 0.361 0.338 "[39, 32, 4]" ( "[39, 4, 59]" "[32, 39, 4]" "[32, 39, 25]" "[39, 4, 48]" "[39, 32, 25]" "[39, 32, 4]" "[39, 4, 48]" "[39, 32, 4]" "[39, 32, 4]" "[39, 32, 25]"))
    (84 0.572 0.445 0.358 "[44, 59, 7]" ( "[44, 59, 7]" "[59, 44, 7]" "[44, 59, 7]" "[44, 59, 52]" "[44, 59, 7]" "[59, 44, 7]" "[59, 44, 7]" "[44, 7, 37]" "[59, 44, 7]" "[59, 44, 7]"))
    (33 0.569 0.438 0.367 "[44, 33, 45]" ( "[44, 48, 33]" "[44, 33, 3]" "[44, 45, 22]" "[44, 33, 3]" "[44, 33, 3]" "[45, 44, 48]" "[44, 45, 33]" "[44, 45, 33]" "[44, 3, 33]" "[44, 33, 3]"))
    (97 0.555 0.342 0.357 "[8, 42, 60]" ( "[42, 8, 55]" "[42, 8, 60]" "[8, 42, 60]" "[8, 42, 60]" "[8, 42, 22]" "[8, 60, 31]" "[8, 42, 31]" "[8, 42, 31]" "[42, 8, 60]" "[8, 42, 60]"))
    (35 0.551 0.334 0.545 "[2, 55, 8]" ( "[2, 55, 8]" "[2, 55, 8]" "[2, 55, 8]" "[55, 4, 2]" "[2, 55, 37]" "[2, 55, 37]" "[2, 55, 37]" "[2, 55, 37]" "[2, 55, 8]" "[55, 17, 2]"))
    (9 0.548 0.463 0.215 "[53, 20, 54]" ( "[53, 20, 32]" "[53, 45, 14]" "[53, 46, 45]" "[53, 10, 24]" "[53, 54, 46]" "[53, 45, 54]" "[53, 20, 24]" "[53, 46, 54]" "[53, 54, 46]" "[53, 20, 23]"))
    (42 0.519 0.316 0.423 "[39, 50, 32]" ( "[39, 50, 16]" "[50, 39, 32]" "[39, 50, 32]" "[39, 50, 32]" "[39, 50, 32]" "[39, 50, 32]" "[50, 39, 32]" "[39, 50, 15]" "[39, 50, 5]" "[50, 32, 39]"))
    (21 0.518 0.434 0.192 "[54, 36, 39]" ( "[54, 36, 45]" "[54, 45, 39]" "[54, 39, 51]" "[54, 51, 39]" "[54, 39, 1]" "[54, 36, 45]" "[54, 36, 13]" "[54, 51, 36]" "[54, 51, 36]" "[54, 40, 39]"))
    (55 0.514 0.309 0.324 "[7, 27, 23]" ( "[7, 27, 5]" "[7, 23, 27]" "[7, 27, 23]" "[27, 7, 23]" "[7, 27, 23]" "[27, 7, 23]" "[27, 7, 23]" "[7, 27, 47]" "[7, 27, 41]" "[7, 27, 23]"))
    (88 0.496 0.404 0.215 "[26, 60, 43]" ( "[26, 51, 43]" "[26, 43, 1]" "[26, 60, 51]" "[26, 43, 1]" "[26, 46, 43]" "[26, 60, 1]" "[26, 60, 46]" "[26, 60, 1]" "[26, 46, 1]" "[26, 51, 60]"))
    (58 0.496 0.386 0.238 "[58, 23, 60]" ( "[58, 48, 60]" "[58, 23, 1]" "[58, 23, 1]" "[58, 48, 23]" "[58, 60, 48]" "[58, 30, 23]" "[58, 23, 60]" "[58, 60, 30]" "[58, 23, 30]" "[58, 30, 14]"))
    (0 0.495 0.400 0.210 "[2, 3, 41]" ( "[2, 3, 4]" "[2, 21, 41]" "[2, 3, 4]" "[2, 28, 41]" "[2, 3, 21]" "[2, 3, 1]" "[2, 21, 41]" "[2, 31, 41]" "[2, 31, 41]" "[2, 31, 41]"))
    (44 0.466 0.284 0.441 "[44, 38, 7]" ( "[44, 38, 7]" "[44, 38, 35]" "[44, 38, 7]" "[44, 38, 7]" "[38, 44, 7]" "[38, 44, 7]" "[38, 7, 44]" "[44, 38, 7]" "[44, 38, 7]" "[44, 38, 7]"))
    (68 0.464 0.339 0.230 "[30, 46, 36]" ( "[30, 46, 10]" "[30, 36, 29]" "[30, 36, 29]" "[30, 29, 46]" "[30, 46, 36]" "[30, 26, 36]" "[30, 26, 36]" "[30, 46, 26]" "[30, 46, 13]" "[30, 46, 26]"))
    (61 0.454 0.310 0.300 "[28, 1, 51]" ( "[28, 1, 30]" "[28, 1, 51]" "[28, 1, 51]" "[28, 51, 30]" "[28, 31, 51]" "[28, 1, 51]" "[28, 31, 30]" "[28, 32, 30]" "[28, 30, 1]" "[28, 1, 51]"))
    (90 0.448 0.391 0.414 "[55, 17, 14]" ( "[17, 55, 14]" "[55, 17, 14]" "[55, 17, 14]" "[55, 17, 14]" "[17, 55, 14]" "[17, 55, 24]" "[17, 55, 14]" "[55, 17, 14]" "[55, 17, 14]" "[17, 55, 14]"))
    (30 0.416 0.304 0.292 "[54, 20, 28]" ( "[54, 20, 28]" "[54, 19, 28]" "[54, 19, 28]" "[54, 20, 28]" "[54, 43, 20]" "[54, 43, 19]" "[54, 19, 28]" "[54, 20, 28]" "[54, 20, 28]" "[54, 20, 43]"))
    (93 0.413 0.346 0.239 "[5, 6, 38]" ( "[5, 38, 47]" "[5, 38, 47]" "[5, 38, 47]" "[5, 6, 16]" "[5, 6, 16]" "[5, 6, 16]" "[5, 6, 47]" "[5, 38, 47]" "[5, 47, 27]" "[5, 6, 38]"))
    (36 0.409 0.256 0.382 "[11, 45, 37]" ( "[11, 45, 37]" "[11, 45, 37]" "[11, 45, 37]" "[11, 45, 37]" "[45, 11, 43]" "[11, 45, 37]" "[45, 11, 37]" "[45, 11, 37]" "[11, 45, 42]" "[11, 45, 37]"))
    (57 0.401 0.306 0.286 "[26, 38, 35]" ( "[26, 47, 2]" "[26, 47, 2]" "[26, 38, 35]" "[26, 35, 38]" "[26, 47, 2]" "[26, 38, 35]" "[26, 38, 35]" "[26, 38, 35]" "[26, 35, 38]" "[26, 35, 38]"))
    (29 0.399 0.280 0.222 "[24, 55, 32]" ( "[24, 55, 32]" "[24, 32, 55]" "[24, 50, 32]" "[24, 50, 6]" "[24, 55, 6]" "[24, 55, 50]" "[24, 32, 55]" "[24, 55, 50]" "[24, 55, 32]" "[24, 32, 55]"))
    (11 0.385 0.230 0.299 "[22, 8, 41]" ( "[22, 8, 45]" "[22, 8, 41]" "[8, 22, 45]" "[22, 8, 41]" "[22, 8, 41]" "[22, 45, 8]" "[22, 41, 8]" "[22, 8, 45]" "[22, 8, 45]" "[22, 8, 41]"))
    (37 0.383 0.295 0.237 "[3, 27, 36]" ( "[3, 45, 27]" "[3, 36, 27]" "[3, 36, 27]" "[3, 27, 36]" "[3, 27, 45]" "[3, 45, 36]" "[3, 27, 36]" "[3, 36, 27]" "[3, 45, 36]" "[3, 27, 36]"))
    (6 0.380 0.294 0.235 "[22, 43, 17]" ( "[22, 17, 43]" "[22, 43, 17]" "[22, 43, 9]" "[22, 9, 17]" "[22, 9, 43]" "[22, 9, 43]" "[22, 43, 17]" "[22, 17, 43]" "[22, 17, 43]" "[22, 43, 17]"))
    (59 0.374 0.370 0.212 "[16, 58, 23]" ( "[16, 23, 58]" "[16, 58, 23]" "[16, 30, 58]" "[16, 30, 58]" "[16, 23, 58]" "[16, 58, 23]" "[16, 30, 23]" "[16, 23, 58]" "[16, 30, 58]" "[16, 23, 58]"))
    (86 0.368 0.270 0.235 "[22, 9, 21]" ( "[22, 9, 21]" "[22, 45, 9]" "[22, 21, 9]" "[22, 9, 45]" "[22, 9, 45]" "[22, 45, 9]" "[22, 21, 45]" "[22, 9, 21]" "[22, 9, 21]" "[22, 21, 9]"))
    (95 0.366 0.320 0.246 "[48, 12, 25]" ( "[48, 25, 38]" "[48, 25, 38]" "[48, 12, 38]" "[48, 12, 25]" "[48, 25, 38]" "[48, 25, 38]" "[48, 12, 38]" "[48, 38, 12]" "[48, 12, 38]" "[48, 12, 25]"))
    (17 0.354 0.242 0.187 "[7, 33, 19]" ( "[7, 20, 19]" "[7, 19, 33]" "[7, 33, 55]" "[7, 33, 55]" "[7, 33, 18]" "[7, 33, 55]" "[7, 33, 20]" "[7, 33, 19]" "[7, 33, 55]" "[7, 19, 18]"))
    (81 0.349 0.266 0.250 "[38, 27, 3]" ( "[38, 27, 3]" "[38, 27, 26]" "[38, 3, 44]" "[38, 27, 26]" "[38, 27, 3]" "[38, 27, 26]" "[38, 3, 44]" "[38, 3, 44]" "[38, 3, 44]" "[38, 27, 3]"))
    (25 0.331 0.210 0.254 "[25, 17, 54]" ( "[25, 17, 54]" "[25, 17, 31]" "[25, 39, 31]" "[25, 17, 54]" "[25, 17, 54]" "[25, 54, 17]" "[25, 54, 31]" "[25, 17, 54]" "[25, 17, 54]" "[25, 54, 17]"))
    (47 0.325 0.253 0.218 "[16, 50, 47]" ( "[16, 50, 35]" "[16, 47, 50]" "[16, 47, 27]" "[16, 47, 36]" "[16, 50, 27]" "[16, 47, 50]" "[16, 47, 50]" "[16, 50, 47]" "[16, 50, 47]" "[16, 50, 47]"))
    (20 0.323 0.208 0.173 "[19, 40, 20]" ( "[19, 3, 40]" "[19, 40, 48]" "[19, 40, 3]" "[19, 20, 40]" "[19, 40, 5]" "[19, 40, 3]" "[19, 20, 40]" "[19, 40, 20]" "[19, 40, 20]" "[19, 40, 55]"))
    (65 0.298 0.203 0.192 "[50, 3, 41]" ( "[50, 3, 19]" "[50, 3, 13]" "[50, 3, 19]" "[50, 3, 29]" "[50, 41, 13]" "[50, 3, 19]" "[50, 3, 13]" "[50, 3, 41]" "[50, 54, 3]" "[50, 3, 41]"))
    (10 0.294 0.196 0.189 "[3, 19, 8]" ( "[3, 48, 23]" "[3, 19, 8]" "[3, 19, 48]" "[3, 19, 43]" "[3, 19, 5]" "[3, 19, 8]" "[3, 8, 19]" "[3, 19, 48]" "[3, 19, 43]" "[3, 19, 43]"))
    (26 0.267 0.150 0.316 "[23, 44, 6]" ( "[23, 44, 6]" "[23, 44, 6]" "[23, 44, 6]" "[44, 23, 6]" "[23, 44, 6]" "[44, 23, 6]" "[23, 44, 6]" "[23, 44, 6]" "[23, 44, 6]" "[23, 44, 6]"))
    (2 0.261 0.169 0.271 "[28, 4, 55]" ( "[4, 28, 12]" "[28, 4, 55]" "[28, 4, 12]" "[28, 4, 12]" "[28, 4, 55]" "[28, 4, 55]" "[28, 4, 2]" "[28, 4, 55]" "[28, 4, 2]" "[28, 4, 12]"))
    (18 0.239 0.145 0.209 "[43, 19, 17]" ( "[43, 19, 17]" "[43, 19, 17]" "[43, 19, 17]" "[43, 17, 19]" "[43, 19, 17]" "[43, 17, 32]" "[43, 19, 17]" "[43, 17, 19]" "[43, 19, 18]" "[43, 19, 17]"))
    (32 0.233 0.144 0.171 "[54, 29, 53]" ( "[54, 53, 29]" "[54, 53, 29]" "[54, 29, 53]" "[54, 29, 5]" "[54, 29, 53]" "[54, 29, 5]" "[54, 29, 53]" "[54, 29, 53]" "[54, 53, 29]" "[54, 29, 53]"))
    (45 0.220 0.143 0.193 "[31, 6, 52]" ( "[31, 6, 52]" "[31, 6, 20]" "[31, 6, 20]" "[31, 6, 22]" "[31, 23, 20]" "[31, 6, 52]" "[31, 6, 23]" "[31, 6, 52]" "[31, 6, 22]" "[31, 6, 52]"))
    (13 0.219 0.138 0.145 "[17, 39, 26]" ( "[17, 39, 60]" "[17, 39, 26]" "[17, 39, 60]" "[17, 39, 60]" "[17, 26, 39]" "[17, 39, 26]" "[17, 39, 60]" "[17, 39, 26]" "[17, 39, 26]" "[17, 26, 39]"))
    (52 0.214 0.127 0.198 "[31, 20, 47]" ( "[31, 9, 56]" "[31, 20, 47]" "[31, 20, 45]" "[31, 20, 45]" "[31, 20, 56]" "[31, 20, 47]" "[31, 20, 26]" "[31, 20, 47]" "[31, 20, 47]" "[31, 20, 47]"))
    (48 0.213 0.148 0.158 "[6, 16, 56]" ( "[6, 16, 56]" "[6, 16, 44]" "[6, 56, 35]" "[6, 16, 56]" "[6, 16, 56]" "[6, 16, 44]" "[6, 16, 12]" "[6, 16, 12]" "[6, 16, 44]" "[6, 16, 44]"))
    (75 0.200 0.150 0.194 "[22, 43, 9]" ( "[22, 9, 43]" "[22, 43, 9]" "[22, 9, 43]" "[22, 9, 43]" "[22, 43, 9]" "[22, 43, 9]" "[22, 9, 43]" "[22, 43, 9]" "[22, 43, 9]" "[22, 43, 9]"))
    (5 0.199 0.131 0.113 "[51, 14, 13]" ( "[51, 14, 36]" "[51, 14, 36]" "[51, 14, 13]" "[51, 14, 13]" "[51, 14, 36]" "[51, 13, 14]" "[51, 14, 45]" "[51, 14, 27]" "[51, 14, 13]" "[51, 14, 45]"))
    (49 0.198 0.109 0.231 "[47, 41, 23]" ( "[47, 41, 23]" "[47, 41, 23]" "[47, 41, 23]" "[47, 23, 5]" "[47, 41, 23]" "[47, 23, 1]" "[47, 41, 23]" "[47, 41, 23]" "[47, 41, 23]" "[47, 41, 23]"))
    (51 0.192 0.133 0.165 "[23, 48, 19]" ( "[23, 48, 43]" "[23, 48, 19]" "[23, 48, 19]" "[23, 48, 19]" "[23, 48, 43]" "[23, 48, 19]" "[23, 43, 48]" "[23, 48, 52]" "[23, 48, 52]" "[23, 48, 52]"))
    (7 0.170 0.100 0.124 "[20, 33, 12]" ( "[20, 33, 57]" "[20, 33, 12]" "[20, 33, 12]" "[20, 33, 12]" "[20, 33, 12]" "[20, 33, 12]" "[20, 33, 55]" "[20, 33, 57]" "[20, 33, 57]" "[20, 12, 33]"))
    (98 0.135 0.109 0.075 "[9, 40, 7]" ( "[9, 40, 7]" "[9, 40, 6]" "[9, 40, 55]" "[9, 40, 36]" "[9, 40, 36]" "[9, 40, 7]" "[9, 40, 7]" "[9, 40, 46]" "[9, 40, 49]" "[9, 40, 36]"))
    (69 0.131 0.069 0.126 "[31, 54, 47]" ( "[31, 47, 54]" "[31, 54, 47]" "[31, 54, 47]" "[31, 54, 47]" "[31, 54, 19]" "[31, 54, 47]" "[31, 54, 47]" "[31, 54, 47]" "[31, 54, 47]" "[31, 54, 19]"))
    (27 0.122 0.094 0.081 "[17, 22, 41]" ( "[17, 22, 52]" "[17, 22, 8]" "[17, 22, 41]" "[17, 22, 41]" "[17, 22, 52]" "[17, 22, 55]" "[17, 22, 41]" "[17, 22, 41]" "[17, 22, 52]" "[17, 22, 55]"))
    (85 0.122 0.094 0.081 "[30, 50, 1]" ( "[30, 50, 46]" "[30, 50, 1]" "[30, 50, 1]" "[30, 50, 16]" "[30, 50, 16]" "[30, 50, 46]" "[30, 50, 29]" "[30, 50, 1]" "[30, 50, 16]" "[30, 50, 1]"))
    (43 0.111 0.078 0.082 "[36, 6, 10]" ( "[36, 6, 43]" "[36, 6, 10]" "[36, 6, 9]" "[36, 6, 10]" "[36, 6, 10]" "[36, 6, 9]" "[36, 6, 54]" "[36, 6, 10]" "[36, 6, 9]" "[36, 6, 10]"))
    (16 0.108 0.078 0.082 "[31, 2, 43]" ( "[31, 2, 27]" "[31, 2, 43]" "[31, 2, 51]" "[31, 2, 51]" "[31, 2, 43]" "[31, 2, 27]" "[31, 2, 43]" "[31, 2, 43]" "[31, 2, 51]" "[31, 2, 43]"))
    (80 0.107 0.053 0.123 "[48, 38, 31]" ( "[48, 38, 31]" "[48, 38, 31]" "[48, 38, 31]" "[48, 38, 31]" "[48, 38, 31]" "[48, 38, 31]" "[48, 38, 31]" "[48, 31, 38]" "[48, 38, 4]" "[48, 38, 31]"))
    (96 0.107 0.053 0.123 "[60, 1, 26]" ( "[60, 1, 26]" "[60, 1, 26]" "[60, 1, 26]" "[60, 1, 26]" "[60, 1, 54]" "[60, 1, 26]" "[60, 26, 1]" "[60, 1, 26]" "[60, 1, 26]" "[60, 1, 26]"))
    (8 0.101 0.062 0.081 "[23, 41, 5]" ( "[23, 41, 5]" "[23, 41, 36]" "[23, 41, 5]" "[23, 41, 52]" "[23, 41, 5]" "[23, 41, 36]" "[23, 41, 5]" "[23, 41, 5]" "[23, 41, 17]" "[23, 41, 5]"))
    (38 0.101 0.078 0.082 "[33, 44, 12]" ( "[33, 44, 57]" "[33, 44, 12]" "[33, 44, 57]" "[33, 44, 12]" "[33, 44, 27]" "[33, 44, 57]" "[33, 44, 57]" "[33, 44, 12]" "[33, 44, 12]" "[33, 44, 12]"))
    (99 0.101 0.078 0.082 "[26, 9, 28]" ( "[26, 9, 43]" "[26, 9, 28]" "[26, 9, 28]" "[26, 9, 43]" "[26, 9, 28]" "[26, 9, 28]" "[26, 9, 28]" "[26, 9, 43]" "[26, 9, 46]" "[26, 9, 43]"))
    (62 0.097 0.062 0.081 "[38, 1, 50]" ( "[38, 1, 50]" "[38, 1, 50]" "[38, 1, 50]" "[38, 1, 47]" "[38, 1, 50]" "[38, 1, 23]" "[38, 1, 23]" "[38, 1, 50]" "[38, 1, 47]" "[38, 1, 50]"))
    (74 0.056 0.031 0.066 "[23, 41, 19]" ( "[23, 41, 29]" "[23, 41, 29]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]" "[23, 41, 19]"))
    (87 0.056 0.031 0.066 "[22, 9, 20]" ( "[22, 9, 20]" "[22, 9, 19]" "[22, 9, 20]" "[22, 9, 20]" "[22, 9, 20]" "[22, 9, 19]" "[22, 9, 20]" "[22, 9, 20]" "[22, 9, 20]" "[22, 9, 20]"))
    (15 0.031 0.016 0.049 "[41, 39, 19]" ( "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 19]" "[41, 39, 18]" "[41, 39, 19]" "[41, 39, 19]"))
    (94 0.031 0.016 0.049 "[6, 52, 23]" ( "[6, 52, 30]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]" "[6, 52, 23]"))
    (77 0.000 0.000 0.000 "[31, 51, 10]" ( "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]" "[31, 51, 10]"))
    ))

(defparameter *scorekeys*
  '(("tailharmpws" . "ssfr") ("all1pws" . "ssfr") ("randpws" . "ssfr")
    ("tailharmpws" . "ltgt") ("all1pws" . "ltgt") ("randpws" . "ltgt")
    ))

(defun string-split (string &key (delimiter #\space))
  (let ((substrings '())
        (length (length string))
        (last 0))
    (flet ((add-substring (i)
             (push (subseq string last i) substrings)))
      (dotimes (i length)
        (when (eq (char string i) delimiter)
          (add-substring i)
          (setq last (1+ i))))
      (add-substring length)
      (nreverse substrings))))

(defun nthpart (n f)
  (nth n (string-split f :delimiter #\_)))

(defvar *t* (make-hash-table :test #'equal))
(defvar *keys* nil)

(setf *print-pretty* nil *print-length* nil)

(defun reshape-idperms (pattern outfile column maxlines)
  (clrhash *t*)
  (setf *keys* nil)
  (loop for fn in (mapcar #'pathname-name (directory pattern))
	as l = (nthpart 5 fn)
	as pws = (nthpart 6 fn)
	as rule = (nthpart 7 fn)
	as key = (list l pws rule)
	do (with-open-file
	    (i (format nil "results/~a.xls" fn))
	    (push key *keys*)
	    (setf (gethash key *t*)
		  (loop for line = (string-split (read-line i nil) :delimiter #\tab)
			until (null (nth (1- column) line))
			collect line))))
  ;; Combine into a master file...
  (with-open-file
   (o outfile :direction :output :if-exists :supersede)
   (loop for key in *keys*
	 do (format o "~a:~a	" (car key) (reformat-key (cons (cadr key) (caddr key)))))
   (format o "~%")
   (loop for i below maxlines ;; UUU
	 do (loop for key in *keys*
		  do (format o "~2$	" (let ((v (nth (1- column) (nth i (gethash key *t*)))))
					    (if v (read-from-string v) "NaN"))))
	 (format o "~%")))
  )

;;; Depends on having *t* pre-loaded by reshape-idperms!

(defun texify-3-wide-idperms ()
  (reshape-idperms "results/*idperms_3*.xls" "results/idps3.xls" 2 1500) ;; Sets up *t*
  (with-open-file
   (o "results/idptables.tex" :direction :output :if-exists :supersede)
   (format o "~%\\begin{longtable}{| l | l || l | l || l | l || l | l || l | l || l | l |}~%\\hline~%")
   (loop for key in *scorekeys*
	 with p = nil
	 do
	 (if p (format o " & ") (setf p t))
	 (format o "\\multicolumn{2}{|c|}{\\textbf{~a}}" (reformat-key key))
	 finally (format o "\\\\~%\\hline~%"))
   (loop for n below (length (gethash '("3" "randpws" "ssfr") *t*)) ;; Use this as a proxy for all (???)
	 do 
	 (loop for key in *scorekeys*
	       with p = nil
	       as (id score) = (nth n (gethash (list "3" (car key) (cdr key)) *t*))
	       do
	       (if p (format o " & ") (setf p t))
	       (format o "~a & ~2$ " id (read-from-string score))
	       finally (format o "\\\\~%\\hline~%")))
   (format o "\\end{longtable}~%~%")
   )
  )

(defun texify-multicons (expid)
  (clrhash *t*)
  (loop for file in (mapcar #'pathname-name (directory (format nil "results/*~a_multicon*.lisp" expid)))
	as (a b c d e f g h i j) = (string-split file :delimiter #\_)
	do (setf (gethash (cons g h) *t*)
		 (with-open-file
		  (i (format nil "results/~a.lisp" file))
		  (loop for line = (read i nil nil)
			until (null line)
			collect line))))
  (with-open-file
   (o "results/vtbtables.tex" :direction :output :if-exists :supersede)
   (format o "
\\begin{longtable}{| l || l | l | l | l | l | l | }
\\hline
\\multicolumn{7}{|c|}{\\textbf{Legend for concordance results subtable}} \\\\
\\hline
\\hline
")
   (loop for key in *scorekeys*
	 as i from 1 by 1
	 do (format o "& ~a " (reformat-key key))
	 (when (= i 3) (format o " & Consensus & CxMean & CxDev\\\\~%\\hline~%")))
   (format o "&  & \\multicolumn{2}{r|}{(above: tailharmpws.ssfr)}\\\\~%\\hline~%\\end{longtable}~%~%")
   (format o "~%\\begin{longtable}{| l || l | l | l | l | l | l | }~%")
   (loop for (PtNo Concordance AlgDistMeanVSConsensus AlgDistStDev Consensus js) in (sort *vtbs* #'< :key #'car)
	 do (format o "~%\\hline~%~a " ptno)
	 (loop for j in js
	       as i from 1 by 1
	       do (format o " & ~a" j)
	       (when (= i 5) (format o " \\\\~%"))
	       )
	 (format o "\\\\~%\\hline~%")
	 (loop for key in *scorekeys*
	       as i from 1 by 1
	       as vs  = (second (assoc ptno (gethash key *t*)))
	       do (format o "& ~a " vs)
	       ;; Nb. we get these from the hard-coded data above -- comes from the tailharmpws.ssfr run
	       (when (= i 3) (format o "& ~a & ~s & ~s\\\\~%\\hline~%" Consensus AlgDistMeanVSConsensus AlgDistStDev))
	       finally (format o "\\\\~%\\hline~%")
	       ))
   (format o "\\end{longtable}~%")
   (format o "\\begin{longtable}{l l l}")
   (loop for k being the hash-keys of *t*
	 using (hash-value v)
	 do
	 (setq v (reverse v))
	 (format o "~a\\\\ ~%" (reformat-key k))
	 (format o " & Near & ~{~a, ~} ...\\\\~%" (reformat-scores (first-n 6 v)))
	 (format o " & Far & ... ~{~a, ~} \\\\~%\\hline~%\\hline~%\\smallskip~%" (reformat-scores (reverse (first-n 6 (reverse v)))))
	 )
   (format o "\\end{longtable}~%")
   ))

(defun first-n (n l) (loop for i below n as j in l collect j))

(defun reformat-key (key) (format nil "~a.~a" (car key) (cdr key)))

(defun reformat-scores (ss)
  (loop for (a b) in ss collect
	(format nil "~a=~a" a b)))

(defun fancy-score (s)
  (let ((s (read-from-string s)))
    (if (= 0.0 s)
	;;(format nil "{\\color{red}~2$}" s)
	(format nil "~2$" s) ;; temporarily dyked
      (format nil "~2$" s))))

(defun extract-comparator ()
  (with-open-file
   (o "results/compare.xls" :direction :output :if-exists :supersede)
   (format o "PtNo	Concordance	AlgDistMeanVSConsensus	AlgDistStDev~%")
   (loop for (PtNo Concordance AlgDistMeanVSConsensus AlgDistStDev Consensus js) in (sort *vtbs* #'< :key #'car)
	 do (format o "~a	~a	~a	~a~%"
		    PtNo Concordance AlgDistMeanVSConsensus AlgDistStDev))))


;;; (Don't call the file "...idperm..." or the dir will sweep it in!)
;(reshape-idperms "results/*idperm*.xls" "results/allidps.xls" 2 1500)
;(reshape-idperms "results/*mvtbcdistest*.xls" "results/allcdistest.xls" 3 10000)

;(texify-3-wide-idperms)

;;; (reshape-multicons reuses *t* for its own purposes)
;(texify-multicons "Apr29")
(extract-comparator)
