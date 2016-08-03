(defpackage :cl-ipopt
  (:use :cl :cffi :bld-ode))

(in-package :cl-ipopt)

(define-foreign-library libipopt
  (:unix "libipopt.so.1")
  (t (:default "IpOpt-vc10")))

(use-foreign-library libipopt)
