(asdf:defsystem :cl-ipopt
  :author "Ben Diedrich"
  :license "MIT"
  :description "Common Lisp CFFI interface to IPOPT"
  :depends-on ("cffi" "bld-ode")
  :serial t
  :components
  ((:file "package")
   (:file "swig")
   (:file "ipopt")
   (:file "hs071")
   (:file "sail2d")))
