;;; PDB (Protein Data Bank) File Parser
;;;
;;; Copyright (c) 2007 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;;
;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.
;;;
;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;

;;; PDB Files consist of a number of "records" which are, in turn, one
;;; or more 80-column lines of ASCII text.
;;;
;;; The reference for the file format can be found at:
;;;   http://www.wwpdb.org/documentation/format23/v2.3.html
;;;
;;; There are a number classes of records:
;;;
;;;  * single (non-continuable) records. Records of these types may
;;;    only appear once (per type) and are not "continuable" (see below).
;;;
;;;  * single continuable records. Records of these types may only
;;;    appear once (per type), but are continuable across mutiple lines.
;;;
;;;  * multiple records. multiple records (per type) are permitted and
;;;    these are continuable across multiple lines.
;;;
;;; * extended line multiple records
;;;
;;; * grouping records
;;; 

(in-package :bio)

;;; since we need to look ahead to see if we are continuing lines or
;;; not, establish a special variable to hold the value of the current
;;; line so that we only need to read it once.
(defparameter *current-line* nil)
(defparameter *current-entry* nil)

(defvar *pdb-string-buffer* nil)

(defun fast-pdb-symbol (string start end)
  "Extracts a string from the line and returns it as a shared Keyword Symbol.
   Never allocates a new string if the symbol already exists!"
  (declare (type string string) (type fixnum start end)
           (optimize (speed 3) (safety 0)))
  (let ((buf-idx 0))
    (loop for i from start below end
          for char = (char string i)
          do (unless (char= char #\Space)
               (setf (char *pdb-string-buffer* buf-idx) char)
               (incf buf-idx)))
    (if (= buf-idx 0)
        nil
        (intern (subseq *pdb-string-buffer* 0 buf-idx) :keyword))))

(declaim (inline strict-pdb-float))
(defun strict-pdb-float (string start end &optional default-value)
  "Parses a fixed-width PDB float, strictly rejecting malformed characters.
   Returns DEFAULT-VALUE if the field is entirely blank."
  (declare (type string string) (type fixnum start end)
           (optimize (speed 3) (safety 0)))
  (let ((sign 1.0d0) 
        (val 0.0d0) 
        (dec-val 0.0d0) 
        (dec-div 1.0d0) 
        (in-dec nil)
        (has-digits nil))
    (loop for i from start below end
          for char = (char string i)
          do (cond
               ((char= char #\Space) nil)
               ((char= char #\-) (setf sign -1.0d0))
               ((char= char #\+) (setf sign 1.0d0))
               ((char= char #\.) (setf in-dec t))
               (t 
                (let ((code (char-code char)))
                  (if (and (>= code 48) (<= code 57))
                      (let ((digit (- code 48)))
                        (setf has-digits t)
                        (if in-dec
                            (setf dec-val (+ (* dec-val 10.0d0) digit)
                                  dec-div (* dec-div 10.0d0))
                            (setf val (+ (* val 10.0d0) digit))))
                      (error "Malformed PDB float at column ~D: Invalid character '~A' in string ~S" 
                             i char (fast-pdb-symbol string start end)))))))
    ;; Return the calculated float, or the fallback if it was completely blank
    (if has-digits
        (* sign (+ val (/ dec-val dec-div)))
        default-value)))

(declaim (inline strict-pdb-coord))
(defun strict-pdb-coord (string start end)
  "Parses a PDB float, but explicitly rejects missing/blank values with an error.
   Used for X, Y, and Z coordinates which MUST exist."
  (declare (type string string) (type fixnum start end)
           (optimize (speed 3) (safety 0)))
  (let ((val (strict-pdb-float string start end)))
    (if val
        val
        (error "CRITICAL PARSE ERROR: Missing coordinate between columns ~D and ~D in string ~S" 
               start end string))))

(declaim (inline strict-pdb-int))
(defun strict-pdb-int (string start end)
  "Parses a fixed-width PDB integer, strictly returning a fixnum.
   Rejects decimals and invalid characters."
  (declare (type string string) (type fixnum start end)
           (optimize (speed 3) (safety 0)))
  (let ((sign 1)
        (val 0)
        (has-digits nil))
    (declare (type fixnum sign val))
    (loop for i from start below end
          for char = (char string i)
          do (cond
               ((char= char #\Space) nil)
               ((char= char #\-) (setf sign -1))
               ((char= char #\+) (setf sign 1))
               (t
                (let ((code (char-code char)))
                  (if (and (>= code 48) (<= code 57))
                      (progn
                        (setf has-digits t)
                        (setf val (+ (* val 10) (- code 48))))
                      (error "Malformed PDB integer at column ~D: Invalid character '~A' in string ~S" 
                             i char (fast-pdb-symbol string start end)))))))
    (when has-digits
      (* sign val))))

(defun read-pdb-record-name (string)
  "Extracts the first 6 characters as a keyword, with zero consing and no regex."
  (let ((len (length string)))
    (when (> len 0)
      (fast-pdb-symbol string 0 (min 6 len)))))

(declaim (inline make-resi-id))
(defun make-resi-id (chain-id residue-seq-number)
  "Creates a standard Lisp cons cell representing a residue (Chain . SeqNumber)."
  (cons chain-id residue-seq-number))

(defun generate-resi-range (chain start-seq end-seq)
  "Generates a list of residue cons cells from start to end."
  (loop for seq from start-seq to end-seq
        collect (make-resi-id chain seq)))

(defclass pdb-entry ()
  ((classification :initarg :classification :accessor classification :initform nil)
   (dep-date :initarg :dep-date :accessor dep-date :initform nil)
   (id-code :initarg :id-code :accessor id-code :initform nil)
   (obsolete :initarg :obsolete :accessor obsolete :initform nil)
   (title :initarg :title :accessor title :initform nil)
   (molecules :initarg :molecules :accessor molecules :initform nil)
   (chains :initarg :chains :accessor chains :initform nil)
   (atom-hash :initarg :atom-hash :accessor atom-hash :initform (make-hash-table))
   (sites :initarg :sites :accessor sites :initform nil)
   (helices :initarg :helices :accessor helices :initform nil)
   (sheets :initarg :sheets :accessor sheets :initform nil)
   (cispeps :initarg :cispeps :accessor cispeps :initform nil)
   (ssbonds :initarg :ssbonds :accessor ssbonds :initform nil)
   (links :initarg :links :accessor links :initform nil)))

(defclass pdb-molecule ()
  ((id :initarg :id :accessor :id)
   (name :initarg :name :accessor name)
   (chains :initarg :chains :accessor chains)
   (fragments :initarg :fragments :accessor fragments)
   (synonym :initarg :synonym :accessor synonym)
   (ec :initarg :ec :accessor ec)
   (engineered :initarg :engineered :accessor engineered)
   (mutation :initarg :mutation :accessor mutation)
   (other-details :initarg :other-details :accessor other-details)))

(defclass pdb-chain ()
  ((name :initarg :name :accessor name)
   (sequence :initarg :sequence :accessor chain-sequence)))

(defclass pdb-record ()
  ((record-name :reader record-name :allocation class)
   (continuable :initarg :continuable :accessor continuable
                :initform nil))
  (:documentation "class for holding information about pdb-records
  while parsing PDB files."))

(defclass continuable-pdb-record (pdb-record)
  ((continuable :initarg :continuable :accessor continuable
                :initform t)
   (continuation-columns :initarg :continuation-columns
                         :accessor continuation-columns
                         :initform '(8 10))
   (field-columns :initarg :field-columns
                  :accessor field-columns
                  :initform '(10 80))
   (lines :initarg :lines :accessor lines :initform nil)
   (data :initarg :data :accessor data :initform nil))
  (:documentation "subclass of pdb-record for holding information
  about contiuable pdb-records for use while parsing PDB files."))

(defgeneric start-pdb-record (record-name line &key entry)
  (:documentation "Reads the first line of a PDB record"))

(defgeneric continue-pdb-record (record line cont &key entry))

(defgeneric finish-pdb-record (record &key entry))

(defgeneric read-continuation (record line)
  (:documentation "generic function for reading the second and later
  lines of a continued record of a PDB file."))

(defmethod read-continuation ((record continuable-pdb-record) line)
  (destructuring-bind (start end)
      (continuation-columns record)
    (strict-pdb-int line start end)))

(defmethod continue-pdb-record ((record continuable-pdb-record) line cont
                                &key (entry *current-entry*))
  (declare (ignore entry))
  (destructuring-bind (start end) (field-columns record)
    (let* ((len (length line))
           ;; Safely bound both start and end
           (actual-start (min (1+ start) len))
           (actual-end (min end len))
           (continuation-lines (subseq line actual-start actual-end)))
      (push continuation-lines (lines record)))))

(defmethod finish-pdb-record (record &key (entry *current-entry*))
  (declare (ignore entry))
  record)

(defmethod finish-pdb-record ((record continuable-pdb-record)
                              &key (entry *current-entry*))
  (declare (ignore entry))
  (setf (data record)
        (apply #'concatenate 'string
               (mapcar (lambda (s) (string-trim '(#\Space #\Tab) s)) 
                       (nreverse (lines record)))))
  record)

(defmethod start-pdb-record (record-name line &key (entry *current-entry*))
  (declare (ignore entry))
  (progn
    (format t "Ignoring ~A~%" line)))

(defmethod start-pdb-record ((record-name (eql :header)) line
                           &key (entry *current-entry*))
  (let ((classification (fast-pdb-symbol line 10 50))
        (dep-date (fast-pdb-symbol line 50 59))
        (id-code (fast-pdb-symbol line 62 66)))
    (setf (classification entry)
          (remove-trailing-spaces classification))
    (setf (dep-date entry)
          (remove-trailing-spaces dep-date))
    (setf (id-code entry)
          (remove-trailing-spaces id-code)))
  nil)

(defmethod start-pdb-record ((record-name (eql :obslte)) line
                           &key (entry *current-entry*))
  (let ((record-name (fast-pdb-symbol line 0 6))
        (continuation (fast-pdb-symbol line 8 10))
        (date (fast-pdb-symbol line 11 20))
        (id-code (fast-pdb-symbol line 21 25))
        (rid-code-1 (fast-pdb-symbol line 31 35))
        (rid-code-2 (fast-pdb-symbol line 36 40))
        (rid-code-3 (fast-pdb-symbol line 41 45))
        (rid-code-4 (fast-pdb-symbol line 46 50))
        (rid-code-5 (fast-pdb-symbol line 51 55))
        (rid-code-6 (fast-pdb-symbol line 56 60))
        (rid-code-7 (fast-pdb-symbol line 61 65))
        (rid-code-8 (fast-pdb-symbol line 66 70)))))

(defclass pdb-title (continuable-pdb-record)
  ((record-name :initform :TITLE :allocation :class)))

(defmethod start-pdb-record ((record-name (eql :title)) line
                          &key (entry *current-entry*))
  (declare (ignore entry))
  (let* ((record (make-instance 'pdb-title))
         (start (first (field-columns record)))
         (end (second (field-columns record)))
         (len (length line)))
    ;; Bound the subseq to the actual length of the line
    (setf (lines record) (list (subseq line (min start len) (min end len))))
    record))

(defmethod finish-pdb-record ((record pdb-title) &key (entry *current-entry*))
  (call-next-method)
  (setf (title entry)
        (data record)))

(defstruct (pdb-atom (:conc-name nil) 
                     (:constructor make-pdb-atom))
  "A lightweight, statically-typed struct replacing the heavy CLOS atom classes."
  rec-name        
  (atom-number 0 :type fixnum) 
  atom-name
  alt-loc 
  residue-name 
  chain-id
  (residue-seq-number 0 :type fixnum) 
  insertion-code
  (x-coord 0.0d0 :type double-float) 
  (y-coord 0.0d0 :type double-float) 
  (z-coord 0.0d0 :type double-float)
  ;; Use (or type null) for optional fields that might return NIL
  (occupancy 0.0d0 :type (or double-float null)) 
  (temp-factor 0.0d0 :type (or double-float null))
  element-symbol 
  (charge 0 :type (or fixnum null)))

(defclass pdb-compound (continuable-pdb-record)
  ((record-name :initform :COMPND :allocation :class)))

(defmethod start-pdb-record ((record-name (eql :compnd)) line
                          &key (entry *current-entry*))
  (declare (ignore entry))
  (let* ((record (make-instance 'pdb-compound))
         (start (first (field-columns record)))
         (end (second (field-columns record)))
         (len (length line)))
    ;; Bound the subseq to the actual length of the line
    (setf (lines record) (list (subseq line (min start len) (min end len))))
    record))

(defclass pdb-remark (pdb-record)
  ((remark-num :initarg :remark-num :accessor remark-num)
   (text :initarg :text :accessor text)))

(defmethod start-pdb-record ((record-name (eql :remark)) line
                             &key (entry *current-entry*))
  (declare (ignore entry))
  (let* ((record (make-instance 'pdb-remark))
         (len (length line)))
    ;; Safely grab the remark number (columns 7-10)
    (when (>= len 10)
      (setf (remark-num record) (strict-pdb-int line 7 10)))
    ;; Safely grab the text (columns 11-79)
    (when (> len 11)
      (setf (text record) (subseq line 11 (min 79 len))))
    record))

(defstruct pdb-site
  (seq-num 0 :type fixnum) site-id (num-res 0 :type fixnum) residues)

(defstruct pdb-helix
  (ser-num 0 :type fixnum) helix-id 
  init-res-name init-chain (init-seq-num 0 :type fixnum) init-icode
  end-res-name end-chain (end-seq-num 0 :type fixnum) end-icode
  (helix-class 0 :type fixnum) comment (length 0 :type (or fixnum null)))

(defstruct pdb-sheet
  (strand 0 :type fixnum) sheet-id (num-strands 0 :type fixnum)
  init-res-name init-chain (init-seq-num 0 :type fixnum) init-icode
  end-res-name end-chain (end-seq-num 0 :type fixnum) end-icode
  (sense 0 :type fixnum)
  cur-atom cur-res-name cur-chain (cur-res-seq 0 :type (or fixnum null)) cur-icode
  prev-atom prev-res-name prev-chain (prev-res-seq 0 :type (or fixnum null)) prev-icode)

(defstruct pdb-cispep
  (ser-num 0 :type fixnum)
  pep1 chain1 (seq-num1 0 :type fixnum) icode1
  pep2 chain2 (seq-num2 0 :type fixnum) icode2
  (mod-num 0 :type (or fixnum null))
  (measure 0.0d0 :type (or double-float null)))

(defstruct pdb-ssbond
  (ser-num 0 :type fixnum)
  res-name1 chain1 (seq-num1 0 :type fixnum) icode1
  res-name2 chain2 (seq-num2 0 :type fixnum) icode2
  sym1 sym2 (length 0.0d0 :type (or double-float null)))

(defstruct pdb-link
  atom1 alt1 res-name1 chain1 (seq-num1 0 :type fixnum) icode1
  atom2 alt2 res-name2 chain2 (seq-num2 0 :type fixnum) icode2
  sym1 sym2 (length 0.0d0 :type (or double-float null)))

(defun parse-pdb-atom-record (line entry)
  (let ((rec-name           (fast-pdb-symbol line 0 6))
        (atom-number        (strict-pdb-int line 6 11))
        (atom-name          (fast-pdb-symbol line 12 16))
        (alt-loc            (fast-pdb-symbol line 16 17))
        (residue-name       (fast-pdb-symbol line 17 20))
        (chain-id           (fast-pdb-symbol line 21 22))
        (residue-seq-number (strict-pdb-int line 22 26))
        (insertion-code     (fast-pdb-symbol line 26 27))
        (x-coord            (strict-pdb-coord line 30 38))
        (y-coord            (strict-pdb-coord line 38 46))
        (z-coord            (strict-pdb-coord line 46 54))
        (occupancy          (strict-pdb-float line 54 60))
        (temp-factor        (strict-pdb-float line 60 66))
        (element-symbol     (fast-pdb-symbol line 76 78))
        (charge             (strict-pdb-int line 78 80)))
        
    (let ((atom (make-pdb-atom :rec-name rec-name 
                               :atom-number atom-number 
                               :atom-name atom-name
                               :alt-loc alt-loc 
                               :residue-name residue-name 
                               :chain-id chain-id
                               :residue-seq-number residue-seq-number 
                               :insertion-code insertion-code
                               :x-coord x-coord 
                               :y-coord y-coord 
                               :z-coord z-coord
                               :occupancy occupancy 
                               :temp-factor temp-factor
                               :element-symbol element-symbol 
                               :charge charge)))
      (setf (gethash atom-number (atom-hash entry)) atom))))

(defmethod start-pdb-record ((record-name (eql :atom)) line
                             &key (entry *current-entry*))
  (parse-pdb-atom-record line entry))

(defmethod start-pdb-record ((record-name (eql :hetatm)) line
                             &key (entry *current-entry*))
  (parse-pdb-atom-record line entry))

;; --- SITE ---
(defmethod start-pdb-record ((record-name (eql :site)) line &key (entry *current-entry*))
  (let ((seq-num (strict-pdb-int line 7 10))
        (site-id (fast-pdb-symbol line 11 14))
        (num-res (strict-pdb-int line 15 17))
        (residues nil))
    ;; Extract up to 4 residues per line [cite: 97, 102]
    (loop for (res chain seq icode) in '((18 22 23 27) (29 33 34 38) (40 44 45 49) (51 55 56 60))
          do (let ((res-name (fast-pdb-symbol line res (+ res 3))))
               (when res-name
                 (push (list :res-name res-name
                             :chain (fast-pdb-symbol line chain (+ chain 1))
                             :seq-num (strict-pdb-int line seq (+ seq 4))
                             :icode (fast-pdb-symbol line icode (+ icode 1)))
                       residues))))
    (push (make-pdb-site :seq-num seq-num :site-id site-id :num-res num-res 
                         :residues (nreverse residues))
          (sites entry)))
  nil)

;; --- HELIX ---
(defmethod start-pdb-record ((record-name (eql :helix)) line &key (entry *current-entry*))
  (push (make-pdb-helix
         :ser-num      (strict-pdb-int line 7 10)     :helix-id     (fast-pdb-symbol line 11 14)
         :init-res-name (fast-pdb-symbol line 15 18)  :init-chain   (fast-pdb-symbol line 19 20)
         :init-seq-num (strict-pdb-int line 21 25)    :init-icode   (fast-pdb-symbol line 25 26)
         :end-res-name (fast-pdb-symbol line 27 30)   :end-chain    (fast-pdb-symbol line 31 32)
         :end-seq-num  (strict-pdb-int line 33 37)    :end-icode    (fast-pdb-symbol line 37 38)
         :helix-class  (strict-pdb-int line 38 40)
         :comment      (string-trim '(#\Space) (subseq line 40 (min 70 (length line))))
         :length       (strict-pdb-int line 71 76))
        (helices entry))
  nil)

;; --- SHEET ---
(defmethod start-pdb-record ((record-name (eql :sheet)) line &key (entry *current-entry*))
  (push (make-pdb-sheet
         :strand       (strict-pdb-int line 7 10)     :sheet-id     (fast-pdb-symbol line 11 14)
         :num-strands  (strict-pdb-int line 14 16)
         :init-res-name (fast-pdb-symbol line 17 20)  :init-chain   (fast-pdb-symbol line 21 22)
         :init-seq-num (strict-pdb-int line 22 26)    :init-icode   (fast-pdb-symbol line 26 27)
         :end-res-name (fast-pdb-symbol line 28 31)   :end-chain    (fast-pdb-symbol line 32 33)
         :end-seq-num  (strict-pdb-int line 33 37)    :end-icode    (fast-pdb-symbol line 37 38)
         :sense        (strict-pdb-int line 38 40)
         :cur-atom     (fast-pdb-symbol line 41 45)   :cur-res-name (fast-pdb-symbol line 45 48)
         :cur-chain    (fast-pdb-symbol line 49 50)   :cur-res-seq  (strict-pdb-int line 50 54)
         :cur-icode    (fast-pdb-symbol line 54 55)
         :prev-atom    (fast-pdb-symbol line 56 60)   :prev-res-name (fast-pdb-symbol line 60 63)
         :prev-chain   (fast-pdb-symbol line 64 65)   :prev-res-seq (strict-pdb-int line 65 69)
         :prev-icode   (fast-pdb-symbol line 69 70))
        (sheets entry))
  nil)

;; --- CISPEP ---
(defmethod start-pdb-record ((record-name (eql :cispep)) line &key (entry *current-entry*))
  (push (make-pdb-cispep
         :ser-num  (strict-pdb-int line 7 10)
         :pep1     (fast-pdb-symbol line 11 14) :chain1   (fast-pdb-symbol line 15 16)
         :seq-num1 (strict-pdb-int line 17 21)  :icode1   (fast-pdb-symbol line 21 22)
         :pep2     (fast-pdb-symbol line 25 28) :chain2   (fast-pdb-symbol line 29 30)
         :seq-num2 (strict-pdb-int line 31 35)  :icode2   (fast-pdb-symbol line 35 36)
         :mod-num  (strict-pdb-int line 43 46)  :measure  (strict-pdb-float line 53 59))
        (cispeps entry))
  nil)

;; --- SSBOND ---
(defmethod start-pdb-record ((record-name (eql :ssbond)) line &key (entry *current-entry*))
  (push (make-pdb-ssbond
         :ser-num   (strict-pdb-int line 7 10)
         :res-name1 (fast-pdb-symbol line 11 14) :chain1   (fast-pdb-symbol line 15 16)
         :seq-num1  (strict-pdb-int line 17 21)  :icode1   (fast-pdb-symbol line 21 22)
         :res-name2 (fast-pdb-symbol line 25 28) :chain2   (fast-pdb-symbol line 29 30)
         :seq-num2  (strict-pdb-int line 31 35)  :icode2   (fast-pdb-symbol line 35 36)
         :sym1      (string-trim '(#\Space) (subseq line 59 65))
         :sym2      (string-trim '(#\Space) (subseq line 66 72))
         :length    (strict-pdb-float line 73 78))
        (ssbonds entry))
  nil)

;; --- LINK ---
(defmethod start-pdb-record ((record-name (eql :link)) line &key (entry *current-entry*))
  (push (make-pdb-link
         :atom1     (fast-pdb-symbol line 12 16) :alt1      (fast-pdb-symbol line 16 17)
         :res-name1 (fast-pdb-symbol line 17 20) :chain1    (fast-pdb-symbol line 21 22)
         :seq-num1  (strict-pdb-int line 22 26)  :icode1    (fast-pdb-symbol line 26 27)
         :atom2     (fast-pdb-symbol line 42 46) :alt2      (fast-pdb-symbol line 46 47)
         :res-name2 (fast-pdb-symbol line 47 50) :chain2    (fast-pdb-symbol line 51 52)
         :seq-num2  (strict-pdb-int line 52 56)  :icode2    (fast-pdb-symbol line 56 57)
         :sym1      (string-trim '(#\Space) (subseq line 59 65))
         :sym2      (string-trim '(#\Space) (subseq line 66 72))
         :length    (strict-pdb-float line 73 78))
        (links entry))
  nil)

(defun find-first-char (char string &key (start 0) (escape-char #\\))
  (let ((pos (position char string :start start)))
    (when pos (if (and (plusp pos)
                       (char= (char string (1- pos)) escape-char))
                  (find-first-char char string :start (1+ pos) :escape-char escape-char)
                  pos))))

(defun parse-specification (string &optional (start 0))
  (let ((colon (position #\: string :start start)))
    (when colon
      (cons 
       (intern (string-trim '(#\Space #\Tab) (subseq string start colon)) :keyword)
       (string-trim '(#\Space #\Tab) (subseq string (1+ colon)))))))

(defun parse-specification-list (string)
  (let ((results nil))
    (dolist (spec (cl-ppcre:split ";" string))
      (let* ((clean-spec (string-trim '(#\Space #\Tab) spec))
             (colon (position #\: clean-spec)))
        (when colon
          (let ((key (intern (string-trim '(#\Space #\Tab) (subseq clean-spec 0 colon)) :keyword))
                (val (string-trim '(#\Space #\Tab) (subseq clean-spec (1+ colon)))))
            (push (list key val) results)))))
    (nreverse results)))

(defmethod finish-pdb-record ((record pdb-compound) &key (entry *current-entry*))
  (declare (ignore entry))
  (call-next-method)
  (print (parse-specification-list (print (data record)))))

(defgeneric write-pdb-record (record-type data stream))

(defun write-pdb-header (entry stream)
  (format stream "~A~10,0T~A~50,0T~A~62,0T~A~80,0T~%"
          "HEADER"
          (classification entry)
          (dep-date entry)
          (id-code entry)))

(defun write-continuable-pdb-record (record-type data stream
                                     &key
                                     (continuation-start 8)
                                     (field-start 10)
                                     (line-length 60))
  (loop with i = 0
     for line from 1
     while (< i (length data))
     do (let ((end (min (length data)
                        (+ i line-length))))
          (format stream "~?"
                  (format nil "~~A~~~D,0T~~@[~~2D ~~]~~~D,0T~~A~~80,0T~~%"
                          continuation-start
                          field-start)
                  (list (typecase record-type
                          (string (string-upcase record-type))
                          (symbol (symbol-name record-type)))
                        (when (> line 1) line)
                        (fast-pdb-symbol data i end))))
     (incf i (if (zerop i)
                 line-length
                 (1- 59)))))

(defmethod write-pdb-record ((record-type (eql :title)) title stream)
  (write-continuable-pdb-record record-type title stream))

(defparameter *continuable-records*
  '(:author :caveat :compnd :expdata :keywds :obslte :source :sprsde
    :title :anisou :atom :cispep :conect :dbref :helix :het :hetsyn
    :link :modres :mtrix1 :mtrix2 :mtrix3 :revdat :eqadv :seqres
    :sheet :sigatm :siguij :site :ssbond :tvect :formul :hetatm :hetnam))

(defun read-pdb-record (stream)
  (let* ((record-key (read-pdb-record-name *current-line*))
         (record (when record-key 
                   (start-pdb-record record-key *current-line*))))
    (if record
        (progn
          (if (continuable record)
              (loop for str = (setf *current-line* (read-line stream nil nil))
                 while str
                 do
                 (if (eq (record-name record)
                         (read-pdb-record-name str))
                     (let ((cont (read-continuation record str)))
                       (if cont
                           (continue-pdb-record record str cont)
                           (return)))
                     (return)))
              (setf *current-line* (read-line stream nil nil)))
          (finish-pdb-record record))
        (setf *current-line* (read-line stream nil nil)))))

(defun read-pdb-stream (stream)
  ;; bind *current-line* per thread to make this thread-safe
  (let ((*current-line* (read-line stream nil nil))
        (*current-entry* (make-instance 'pdb-entry))
        (*pdb-string-buffer* (make-string 80 :initial-element #\Space)))
    (loop while *current-line*
          do (read-pdb-record stream))
    (resolve-links *current-entry* *deferred-links*)
    *current-entry*))

(defun read-pdb-file (file)
  (with-open-file (stream file)
    (read-pdb-stream stream)))

(defun write-pdb-stream (stream entry)
  (write-pdb-header entry stream)
  (let ((title (title entry)))
    (when title (write-pdb-record :title title stream))))
