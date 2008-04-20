
(asdf:defsystem #:cl-bio-rucksack
  :name "cl-bio-rucksack"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version #.(with-open-file
                 (vers (merge-pathnames "version.lisp-expr" *load-truename*))
               (read vers))
  :licence "BSD"
  :depends-on (cl-bio ch-asdf rucksack)
  :components
  ((:module :rucksack
            :components
            ((:module :rucksack-data)
             (:cl-source-file "defpackage")
             (:cl-source-file "rucksack" :depends-on ("defpackage"))))))