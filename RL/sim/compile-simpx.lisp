(defun main()
    (load (compile-file "simpx.lisp"))
    (save-application "simpx.exe" :toplevel-function #'main :prepend-kernel t)
)