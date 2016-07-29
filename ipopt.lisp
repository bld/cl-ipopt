(in-package :cl-ipopt)

(defcallback eval_f :int ((n :int) (x :pointer) (new_x :int) (obj_value :pointer))
  (assert (= n 4))
  (let ((x0 (mem-aref x :double 0))
	(x1 (mem-aref x :double 1))
	(x2 (mem-aref x :double 2))
	(x3 (mem-aref x :double 3)))
    (setf (mem-ref obj_value :double) (+ (* x0 x3 (+ x0 x1 x2)) x2)))
  1)

(defcallback eval_grad_f :int ((n :int) (x :pointer) (new_x :int) (grad_f :pointer))
  (assert (= n 4))
  (let ((x0 (mem-aref x :double 0))
	(x1 (mem-aref x :double 1))
	(x2 (mem-aref x :double 2))
	(x3 (mem-aref x :double 3)))
    (setf (mem-aref grad_f :double 0) (+ (* x0 x3) (* x3 (+ x0 x1 x2)))
	  (mem-aref grad_f :double 1) (* x0 x3)
	  (mem-aref grad_f :double 2) (+ (* x0 x3) 1)
	  (mem-aref grad_f :double 3) (* x0 (+ x0 x1 x2))))
  1)

(defcallback eval_g :int ((n :int) (x :pointer) (new_x :int) (m :int) (g :pointer))
  (assert (= n 4))
  (assert (= m 2))
  (let ((x0 (mem-aref x :double 0))
	(x1 (mem-aref x :double 1))
	(x2 (mem-aref x :double 2))
	(x3 (mem-aref x :double 3)))
    (setf (mem-aref g :double 0) (* x0 x1 x2 x3)
	  (mem-aref g :double 1) (+ (* x0 x0) (* x1 x1) (* x2 x2) (* x3 x3))))
  1)

(defcallback eval_jac_g :int ((n :int) (x :pointer) (new_x :int) (m :int) (nele_jac :int) (irow :pointer) (jcol :pointer) (values :pointer))
  (let ((x0 (mem-aref x :double 0))
	(x1 (mem-aref x :double 1))
	(x2 (mem-aref x :double 2))
	(x3 (mem-aref x :double 3)))
    (if (null-pointer-p values) ;; Return Jacobian structure
	(setf (mem-aref irow :int 0) 0 (mem-aref jcol :int 0) 0
	      (mem-aref irow :int 1) 0 (mem-aref jcol :int 1) 1
	      (mem-aref irow :int 2) 0 (mem-aref jcol :int 2) 2
	      (mem-aref irow :int 3) 0 (mem-aref jcol :int 3) 3
	      (mem-aref irow :int 4) 1 (mem-aref jcol :int 4) 0
	      (mem-aref irow :int 5) 1 (mem-aref jcol :int 5) 1
	      (mem-aref irow :int 6) 1 (mem-aref jcol :int 6) 2
	      (mem-aref irow :int 7) 1 (mem-aref jcol :int 7) 3)
	(setf
	 (mem-aref values :double 0) (* x1 x2 x3) ; 0,0
	 (mem-aref values :double 1) (* x0 x2 x3) ; 0,1
	 (mem-aref values :double 2) (* x0 x1 x3) ; 0,2
	 (mem-aref values :double 3) (* x0 x1 x2) ; 0,3
	 
	 (mem-aref values :double 4) (* 2 x0) ; 1,0
	 (mem-aref values :double 5) (* 2 x1) ; 1,1
	 (mem-aref values :double 6) (* 2 x2) ; 1,2
	 (mem-aref values :double 7) (* 2 x3)))) ; 1,3
  1)

(defcallback eval_h :int ((n :int) (x :pointer) (new_x :int) (obj_factor :double) (m :int) (lmb :pointer) (new_lmb :int) (nele_hess :int) (irow :int) (jcol :int) (values :pointer))
  (let ((x0 (mem-aref x :double 0))
	(x1 (mem-aref x :double 1))
	(x2 (mem-aref x :double 2))
	(x3 (mem-aref x :double 3))
	(lmb0 (mem-aref lmb :double 0))
	(lmb1 (mem-aref lmb :double 1)))
    (if (null-pointer-p values) ;; Return Hessian structure
	(let ((i 0))
	  (dotimes (row 4)
	    (dotimes (col (1+ row))
	      (setf (mem-aref irow :int i) row
		    (mem-aref jcol :int i) col)
	      (incf i)))
	  (assert (= i nele_hess)))
	;; Return the values. Symmetric, fill lower left
	;; Fill objective portion
	(progn
	  (setf (mem-aref values :double 0) (* obj_factor 2 x3) ; 0,0
		(mem-aref values :double 1) (* obj_factor x3) ; 1,0
		(mem-aref values :double 2) 0 ; 1,1
		
		(mem-aref values :double 3) (* obj_factor x3) ; 2,0
		(mem-aref values :double 4) 0 ; 2,1
		(mem-aref values :double 5) 0 ; 2,2
		
		(mem-aref values :double 6) (* obj_factor (+ (* 2 x0) x1 x2)) ; 3,0
		(mem-aref values :double 7) (* obj_factor x0) ; 3,1
		(mem-aref values :double 8) (* obj_factor x0) ; 3,2
		(mem-aref values :double 9) 0) ; 3,3
	  ;; Add portion for the 1st constraint
	  (incf (mem-aref values :double 1) (* lmb0 x2 x3)) ; 1,0
	  
	  (incf (mem-aref values :double 3) (* lmb0 x1 x3)) ; 2,0
	  (incf (mem-aref values :double 4) (* lmb0 x0 x3)) ; 2,1
	  
	  (incf (mem-aref values :double 6) (* lmb0 x1 x2)) ; 3,0
	  (incf (mem-aref values :double 7) (* lmb0 x0 x2)) ; 3,1
	  (incf (mem-aref values :double 8) (* lmb0 x0 x1)) ; 3,2
	  
	  ;; Add portion for the 2nd constraint
	  (incf (mem-aref values :double 0) (* lmb1 2)) ; 0,0
	  (incf (mem-aref values :double 2) (* lmb1 2)) ; 1,1
	  (incf (mem-aref values :double 5) (* lmb1 2)) ; 2,2
	  (incf (mem-aref values :double 9) (* lmb1 2))))) ; 3,3
  1)

(defun test ()
  (let ((n 4)
	(m 2))
    (with-foreign-objects
	((x_l :double n)
	 (x_u :double n)
	 (g_l :double m)
	 (g_u :pointer m)
	 (x :pointer n)
	 (mult_x_l :double n)
	 (mult_x_u :double n)
	 (obj :double))

      (setf (mem-aref g_l :double 0) 25d0 (mem-aref g_u :double 0) 2d19
	    (mem-aref g_l :double 1) 40d0 (mem-aref g_u :double 1) 40d0)

      (dotimes (i n)
	(setf (mem-aref x_l :double i) 1d0
	      (mem-aref x_u :double i) 5d0))

      (let ((nlp (createipoptproblem n x_l x_u m g_l g_u 8 10 0 (callback eval_f) (callback eval_g) (callback eval_grad_f) (callback eval_jac_g) (callback eval_h))))
	;; Set some options
	(addipoptnumoption nlp "tol" 1d-9)
	(addipoptstroption nlp "mu_strategy" "adaptive")
	;; Allocate space for initial point and set values
	(setf (mem-aref x :double 0) 1d0
	      (mem-aref x :double 1) 5d0
	      (mem-aref x :double 2) 5d0
	      (mem-aref x :double 3) 1d0)
	;; Solve problem
	(ipoptsolve nlp x (null-pointer) obj (null-pointer) mult_x_l mult_x_u (null-pointer))
	))))
