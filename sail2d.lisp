(in-package :cl-ipopt)

(defcallback sail2d_f :int ((n :int) (x :pointer) (new_x :int) (obj_value :pointer))
  (setf (mem-aref obj_value :double) (mem-aref x :double (1- n)))
  1)

(defcallback sail2d_grad_f :int ((n :int) (x :pointer) (new_x :int) (grad_f :pointer))
  (dotimes (i (1- n))
    (setf (mem-aref grad_f :double i) 0d0))
  (setf (mem-aref grad_f :double (1- n)) 1d0)
  1)

(defclass sail-data ()
  ((mu :initarg :mu)
   (lightness :initarg :lightness)
   (sia :initarg :sia)))

(defun sail2d-eom (tm y p)
  (with-slots (mu lightness sia) p
    (let ((r (aref y 0))
	  (th (aref y 1))
	  (vr (aref y 2))
	  (vt (aref y 3))
	  (co (cos sia))
	  (si (sin sia)))
      (vector
       vr
       (/ vt r)
       (+ (/ (expt vt 2) r)
	  (* mu (/ (1- (* lightness (expt co 2) (abs co))) (expt r 2))))
       (- (/ (* mu lightness (expt co 2) si) (expt r 2))
	  (/ (* vr vt) r))))))

(defun sail2d-integrate (x)
  (let* ((mu 1d0)
	 (lightness 0.05d0)
	 (n (1- (length x)))
	 (x0 #(1d0 0d0 0d0 1d0))
	 (t0 0d0)
	 (tf (aref x n))
	 (tsegment (/ (- tf t0) n)))
    (loop for i below n
       for tfi = (+ t0 tsegment) then (+ tfi tsegment)
       for t0i = t0 then (+ t0i tsegment)
       for x0i = x0 then xfi
       for sia across x
       for res = (rka #'sail2d-eom t0i tfi x0i :param (make-instance 'sail-data :mu mu :lightness lightness :sia sia))
       for xfi = (second (car (last res)))
       append res)))

(defun sail2d-g (x)
  (let ((rf 1.5d0)
	(vrf 0d0)
	(vtf 0d0)
	(xf (second (car (last (sail2d-integrate x))))))
    (vector (- rf (aref xf 0))
	    (- vrf (aref xf 2))
	    (- vtf (aref xf 3)))))

(defcallback sail2d_g :int ((n :int) (x :pointer) (new_x :int) (m :int) (g :pointer))
  (let ((gv (sail2d-g (from-foreign-vector x :double n))))
    (set-foreign-vector gv g :double m))
  1)

(defcallback sail2d_jac_g :int ((n :int) (x :pointer) (new_x :int) (m :int) (nele_jac :int) (irow :pointer) (jcol :pointer) (values :pointer))
  (assert (null-pointer-p values)) ; always estimate Jacobian
  ;; Return Jacobian structure: m x n: 3 x n
  (set-structure (jacobian-structure m n) irow jcol)
  1)

(defcallback sail2d_h :int ((n :int) (x :pointer) (new_x :int) (obj_factor :double) (m :int) (lmb :pointer) (new_lmb :int) (nele_hess :int) (irow :pointer) (jcol :pointer) (values :pointer))
  (if (null-pointer-p values)
      ;; Return Hessian structure: n x n
      (set-structure (hessian-structure n) irow jcol)
      ;; Otherwise return Hessian
      (dotimes (i nele_hess)
	(setf (mem-aref values :double i) 0d0)))
  1)

(defun sail2d-test ()
  (let ((n 6)
	(m 3))
    (with-foreign-objects
	((x_l :double n)
	 (x_u :double n)
	 (g_l :double m)
	 (g_u :double m)
	 (x :pointer n)
	 (mult_x_l :double n)
	 (mult_x_u :double n)
	 (obj :double))

      ;; Constraint functions are all equality: set upper/lower to 0
      (setf (mem-aref g_l :double 0) 0d0 (mem-aref g_u :double 0) 0d0
	    (mem-aref g_l :double 1) 0d0 (mem-aref g_u :double 1) 0d0
	    (mem-aref g_l :double 2) 0d0 (mem-aref g_u :double 2) 0d0)

      ;; Bound on sun incidence angles -pi/2 to pi/2
      (dotimes (i (1- n))
	(setf (mem-aref x_l :double i) (- (/ pi 2))
	      (mem-aref x_u :double i) (/ pi 2)))

      ;; Bounds on tf 5 to 100
      (setf (mem-aref x_l :double (1- n)) 5d0
	    (mem-aref x_u :double (1- n)) 100d0)

      ;; Create NLP problem
      (let ((nlp (createipoptproblem n x_l x_u m g_l g_u (length (jacobian-structure m n)) (length (hessian-structure n)) 0 (callback sail2d_f) (callback sail2d_g) (callback sail2d_grad_f) (callback sail2d_jac_g) (callback sail2d_h))))

	;; Set options
	;;(addipoptnumoption nlp "tol" 1d-9)
	(addipoptintoption nlp "max_iter" 1000)
	;;(addipoptstroption nlp "mu_strategy" "adaptive")
	;;(addipoptstroption nlp "hessian_approximation" "limited-memory")
	(addipoptstroption nlp "jacobian_approximation" "finite-difference-values")
	

	;; Allocate space for initial point and set values
	;; Initial angles: (atan (/ (sqrt 2)))
	(dotimes (i (1- n))
	  (setf (mem-aref x :double i) (atan (/ (sqrt 2d0)))))
	;; Initial time of flight: 10 TU
	(setf (mem-aref x :double (1- n)) 10d0)

	;; Solve problem
	(ipoptsolve nlp x (null-pointer) obj (null-pointer) mult_x_l mult_x_u (null-pointer))
	(let* ((xv (from-foreign-vector x :double n))
	       (f (mem-ref obj :double))
	       (g (sail2d-g xv))
	       (traj (sail2d-integrate xv)))
	  (export-trajectory traj "sail2d_traj.dat")
	  (values xv f g))
	))))

(defun export-trajectory (traj file)
  (with-open-file (s file :direction :output :if-exists :supersede)
    (loop for (tm x) in traj
       do (format s "~a ~a ~a ~a ~a~%" tm (aref x 0) (aref x 1) (aref x 2) (aref x 3)))))
