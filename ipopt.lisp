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
	(x3 (mem-aref x :double 3)))
    (if (null-pointer-p values) ;; Return Hessian structure
	(let ((i 0))
	  (dotimes (row 4)
	    (dotimes (col (1+ row))
	      (setf (mem-aref irow :int i) row
		    (mem-aref jcol :int i) col)
	      (incf i)))
	  (assert (= i nele_hess)))
	))
  1)

(defun test ()
  (let ((n 4)
	(m 2))
    (with-foreign-objects
	((x_l :double n)
	 (x_u :double n)
	 (g_l :double m)
	 (g_u :pointer m)
	 (status :int)
	 (x :pointer n)
	 (mult_x_l :double n)
	 (mult_x_u :double n)
	 (obj :double))
      )))
