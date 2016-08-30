(in-package :cl-ipopt)

(defun jacobian-structure (m n)
  (loop for i below m
     append (loop for j below n
	       collect (list i j))))

(defun hessian-structure (n)
  (let ((i 0))
    (loop for row below n
       append (loop for col below (1+ row)
		 collect (list row col)
		 do (incf i)))))

(defun grad (f x)
  (let* ((eps (sqrt double-float-epsilon))
	 (n (length x))
	 (gf (make-array n)))
    (dotimes (i n)
      (let ((x+1 (copy-seq x))
	    (x-1 (copy-seq x)))
	(incf (aref x+1 i) eps)
	(decf (aref x-1 i) eps)
	(setf (aref gf i) (/ (- (funcall f x+1) (funcall f x-1)) 2 eps))))
    gf))

(defun from-foreign-vector (fv type size)
  "Make list vector from foreign vector given pointer, type, and size"
  (let ((v (make-array size)))
    (dotimes (i size)
      (setf (aref v i) (mem-aref fv type i)))
    v))

(defun set-foreign-vector (v fv type size)
  "Assign Lisp vector elements to foreign pointer of given type and size"
  (dotimes (i size)
    (setf (mem-aref fv type i) (aref v i))))

(defun set-structure (str irow jcol)
  (loop for strk in str
     for k below (length str)
     do (setf (mem-aref irow :int k) (first strk)
	      (mem-aref jcol :int k) (second strk))))

