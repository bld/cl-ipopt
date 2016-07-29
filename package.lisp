(defpackage :cl-ipopt
  (:use :cl :cffi))

(in-package :cl-ipopt)

(define-foreign-library libipopt
  (:unix "libipopt.so.1")
  (t (:default "IpOptFSS")))

(use-foreign-library libipopt)
