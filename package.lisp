(defpackage :cl-ipopt
  (:use :cl :cffi))

(in-package :cl-ipopt)

(define-foreign-library libipopt
  (t (:default "IpOptFSS")))

(use-foreign-library libipopt)
