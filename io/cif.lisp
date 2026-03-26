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

(in-package :bio)

(defparameter *cif-trim-bag* '(#\' #\" #\Space #\Tab #\Return #\Newline)
  "The mathematically safe list of characters to trim from mmCIF values.")

(defparameter *cif-loop-dispatch-table*
  '(("_atom_site."           . process-atom-site-row)            ;; Atoms
    ("_struct_conf."         . process-struct-conf-row)          ;; Helices
    ("_struct_sheet_range."  . process-struct-sheet-row)         ;; Sheets
    ("_struct_site_gen."     . process-struct-site-row)          ;; Sites
    ("_struct_conn."         . process-struct-conn-row)          ;; Disulfides & Links
    ("_struct_mon_prot_cis." . process-struct-mon-prot-cis-row)) ;; Cis-Peptides
  "Maps mmCIF loop prefixes to the functions that process their rows.")

(defvar *cif-buffer* nil
  "Thread-local adjustable string buffer for the lexer. Must be bound dynamically.")

(defvar *cif-string-pool* nil
  "An ephemeral hash table to deduplicate strings and prevent keyword bloat.")

(defun pool-string (str)
  "Returns a deduplicated string pointer. Fast as EQ, safe for the GC."
  (if *cif-string-pool*
      (or (gethash str *cif-string-pool*)
          (setf (gethash str *cif-string-pool*) str))
      str))

(declaim (inline reset-cif-buffer))
(defun reset-cif-buffer ()
  (setf (fill-pointer *cif-buffer*) 0))

(declaim (inline push-cif-char))
(defun push-cif-char (char)
  (vector-push-extend char *cif-buffer*))

;; ----------------------------------------------------------------------
;; 0. TYPE PARSERS
;; ----------------------------------------------------------------------

(defun cif-int (alist key &optional (default 0))
  "Extracts a CIF value and mathematically guarantees a FIXNUM."
  (let ((val (cdr (assoc key alist))))
    (cond ((integerp val) val)
          ((stringp val) (or (parse-integer val :junk-allowed t) default))
          ((symbolp val) (or (parse-integer (string val) :junk-allowed t) default))
          (t default))))

(defun cif-str (alist key &optional (default ""))
  "Extracts a CIF value and guarantees a standard STRING."
  (let ((val (cdr (assoc key alist))))
    (if val (princ-to-string val) default)))

(defun cif-sym (alist key &optional default)
  "Extracts a CIF value and guarantees a KEYWORD SYMBOL."
  (let ((val (cdr (assoc key alist))))
    (values ;; Explicitly strip the secondary status value returned by 'intern'
     (cond ((keywordp val) val)
           ((symbolp val) (intern (string val) "KEYWORD"))
           ((stringp val) (intern val "KEYWORD"))
           (t default)))))

(defun cif-float (alist key &optional (default 0.0d0))
  "Extracts a CIF value and mathematically guarantees a DOUBLE-FLOAT."
  (let ((val (cdr (assoc key alist))))
    (cond
      ((floatp val) (coerce val 'double-float))
      ((consp val) (coerce (car val) 'double-float))
      ((integerp val) (coerce val 'double-float))
      ((stringp val) 
       (let* ((*read-default-float-format* 'double-float)
              (parsed (ignore-errors (read-from-string val))))
         (if (numberp parsed) (coerce parsed 'double-float) default)))
      (t default))))

(defun parse-cif-number (str)
  "Parses a CIF number, handling standard uncertainty and scientific notation safely.
   Returns (values float-value float-uncertainty)."
  (when (or (string= str "?") (string= str "."))
    (return-from parse-cif-number (values nil nil)))
  
  (let* ((e-pos (position-if (lambda (c) (member c '(#\e #\E))) str))
         (exp-val (if e-pos (or (parse-integer (subseq str (1+ e-pos)) :junk-allowed t) 0) 0))
         (base-str (if e-pos (subseq str 0 e-pos) str))
         (paren-pos (position #\( base-str))
         (val-str (if paren-pos (subseq base-str 0 paren-pos) base-str))
         (unc-str (when paren-pos (subseq base-str (1+ paren-pos) (position #\) base-str))))
         ;; SAFETY ADDED HERE
         (parsed-val (let ((*read-default-float-format* 'double-float))
                       (ignore-errors (read-from-string val-str)))))
    
    (unless (numberp parsed-val)
      (return-from parse-cif-number (values nil nil)))
      
    (let ((val-float (coerce parsed-val 'double-float)))
      (values
       (* val-float (expt 10.0d0 exp-val))
       (if unc-str
           (let* ((unc-int (or (parse-integer unc-str :junk-allowed t) 0))
                  (dot-pos (position #\. val-str))
                  (decimals (if dot-pos (- (length val-str) dot-pos 1) 0)))
             (* unc-int (expt 10.0d0 (- exp-val decimals))))
           0.0d0)))))

(defun parse-cif-float (str)
  "Wraps a CIF number into a (value . uncertainty) cons cell."
  (multiple-value-bind (v u) (parse-cif-number str)
    (if v (cons v u) nil)))

(defun split-whitespace (str)
  "Helper to safely tokenize a string by spaces and newlines."
  (let ((tokens nil)
        (start 0)
        (len (length str)))
    (loop
      ;; Find the start of the next word
      (let ((s (position-if-not (lambda (c) (member c '(#\Space #\Tab #\Newline #\Return))) str :start start)))
        (unless s (return (nreverse tokens))) ;; No more words? We're done.
        
        ;; Find the end of the word
        (let ((e (position-if (lambda (c) (member c '(#\Space #\Tab #\Newline #\Return))) str :start s)))
          (push (subseq str s (or e len)) tokens)
          (unless e (return (nreverse tokens))) ;; Reached the end of the string? We're done.
          
          ;; Advance the pointer
          (setf start e))))))

(defun parse-range (str)
  "Parses a CIF range into a (min . max) cons cell, perfectly handling negative numbers."
  (when (or (string= str "?") (string= str ".")) (return-from parse-range nil))
  (let* ((clean-str (remove-if (lambda (c) (member c '(#\Space #\Tab #\Newline))) str))
         ;; Find a comma, colon, or a hyphen that is NOT the first character
         (split-pos (or (position #\, clean-str)
                        (position #\: clean-str)
                        (position #\- clean-str :start 1))))
    (let ((*read-default-float-format* 'double-float))
      (if split-pos
          (cons (read-from-string (subseq clean-str 0 split-pos))
                (read-from-string (subseq clean-str (1+ split-pos))))
          (let ((val (read-from-string clean-str)))
            (cons val val))))))

(defun parse-3x4-matrix (str)
  "Parses a 3x4 matrix into a cons cell: (values-array . uncertainties-array)."
  (when (or (string= str "?") (string= str ".")) (return-from parse-3x4-matrix nil))
  (let ((vals (make-array 12 :element-type 'double-float :initial-element 0.0d0))
        (uncs (make-array 12 :element-type 'double-float :initial-element 0.0d0))
        (idx 0))
    (loop for token in (split-whitespace str)
          while (< idx 12)
          do (multiple-value-bind (v u) (parse-cif-number token)
               (setf (aref vals idx) v)
               (setf (aref uncs idx) u)
               (incf idx)))
    (cons vals uncs)))

(defun parse-3x4-matrices (str)
  "Parses a string of multiple 3x4 matrices by calling parse-3x4-matrix on 12-token chunks."
  (when (or (string= str "?") (string= str ".")) (return-from parse-3x4-matrices nil))
  (let ((tokens (split-whitespace str))
        (result nil))
    (loop while (>= (length tokens) 12)
          do (let ((chunk-str (format nil "~{~A~^ ~}" (subseq tokens 0 12))))
               (push (parse-3x4-matrix chunk-str) result))
             (setf tokens (nthcdr 12 tokens)))
    (nreverse result)))

(defun parse-cif-boolean (str)
  "Parses YES/NO to Lisp T/NIL."
  (string-equal str "YES"))

;; ----------------------------------------------------------------------
;; 1. THE LEXER
;; ----------------------------------------------------------------------

(defun get-cif-token (stream &optional (eof-error-p t) eof-value)
  "Reads the next mmCIF token from the stream. 
   Returns a cons: (TYPE . VALUE), where TYPE is :TAG, :VALUE, :LOOP, :DATA, or :EOF."
  (reset-cif-buffer)
  (let ((char (read-char stream eof-error-p eof-value))
        (at-line-start t))
    
    (loop
      (cond
        ;; 1. Handle EOF
        ((eql char eof-value)
         (return (cons :EOF nil)))
        
        ;; 2. Skip Comments (everything after '#' until a newline)
        ((char= char #\#)
         (loop for c = (read-char stream nil nil)
               until (or (null c) (char= c #\Newline)))
         (setf char (read-char stream eof-error-p eof-value)
               at-line-start t))
        
        ;; 3. Skip Whitespace
        ((member char '(#\Space #\Tab #\Newline #\Return))
         (when (char= char #\Newline)
           (setf at-line-start t))
         (setf char (read-char stream eof-error-p eof-value)))

        ;; 4. Multi-line Text Blocks (Semicolon at start of line)
        ((and at-line-start (char= char #\;))
         (return (cons :VALUE (read-cif-text-block stream))))
        
        ;; 5. Quoted Strings (Single or Double)
        ((or (char= char #\') (char= char #\"))
         (let ((quote-char char))
           (loop for c = (read-char stream eof-error-p eof-value)
                 until (or (eql c eof-value) (char= c quote-char))
                 do (push-cif-char c))
           (return (cons :VALUE (copy-seq *cif-buffer*)))))
        
        ;; 6. Standard Unquoted Tokens (Tags, Values, loop_, data_)
        (t
         (push-cif-char char)
         (loop for c = (peek-char nil stream nil nil)
               while (and c (not (member c '(#\Space #\Tab #\Newline #\Return))))
               do (push-cif-char (read-char stream nil nil)))
         
         (let ((token-str (copy-seq *cif-buffer*)))
           (return 
             (cond
               ;; ADD THIS NEW BLOCK HERE:
               ((and (>= (length token-str) 5) 
                     (string-equal (subseq token-str 0 5) "save_"))
                (cons :SAVE token-str))
               ((and (>= (length token-str) 5) 
                     (string-equal (subseq token-str 0 5) "data_"))
                (cons :DATA token-str))
               ((string-equal token-str "loop_")
                (cons :LOOP nil))
               ((char= (char token-str 0) #\_)
                (cons :TAG token-str))
               (t (cons :VALUE token-str))))))))))

(defun read-cif-text-block (stream)
  "Reads a multi-line text block bounded by ';' at the start of a line."
  (reset-cif-buffer)
  (let ((prev-char #\Newline))
    (loop for c = (read-char stream nil nil)
          while c
          do (if (and (char= prev-char #\Newline) (char= c #\;))
                 ;; We hit the closing semicolon, finish up
                 (return (copy-seq *cif-buffer*))
                 (progn
                   (push-cif-char c)
                   (setf prev-char c))))))

;; ----------------------------------------------------------------------
;; 2. THE PARSER STATE MACHINE
;; ----------------------------------------------------------------------

(defun read-generic-cif (stream)
  "Reads an entire mmCIF file into a generic hash table.
   Keys are tags (strings), values are lists of string values."
  (let ((cif-data (make-hash-table :test 'equal))
        (current-token (get-cif-token stream nil nil)))

    (loop while (and current-token (not (eq (car current-token) :EOF)))
          do (case (car current-token)

               ;; 1. Flat Tag-Value Pair
               (:TAG
                (let ((tag (cdr current-token)))
                  (setf current-token (get-cif-token stream nil nil))
                  (when (eq (car current-token) :VALUE)
                    (push (cdr current-token) (gethash tag cif-data))
                    (setf current-token (get-cif-token stream nil nil)))))

               ;; 2. Loop Block
               (:LOOP
                (let ((headers nil))
                  (setf current-token (get-cif-token stream nil nil))
                  ;; Gather all headers
                  (loop while (eq (car current-token) :TAG)
                        do (push (cdr current-token) headers)
                           (setf current-token (get-cif-token stream nil nil)))
                  (setf headers (nreverse headers))

                  ;; Read values and assign them to the correct column (header)
                  (loop while (eq (car current-token) :VALUE)
                        do (loop for header in headers
                                 do (when (eq (car current-token) :VALUE)
                                      (push (cdr current-token) (gethash header cif-data))
                                      (setf current-token (get-cif-token stream nil nil)))))))

               ;; 3. Ignore DATA blocks or unhandled tokens
               (t (setf current-token (get-cif-token stream nil nil)))))

    ;; Reverse the lists in the hash table so they match the original file order
    (maphash (lambda (k v)
               (setf (gethash k cif-data) (nreverse v)))
             cif-data)
    cif-data))

(defun read-generic-cif-file (file)
  (with-open-file (stream file :direction :input)
    (read-generic-cif stream)))

(defun generate-type-dispatch-table (dictionary-file)
  "A structurally sound DDL2/PDBx compiler that correctly isolates loop types from frame types."
  (let ((explicit-types (make-hash-table :test 'equalp))
        (parent-links   (make-hash-table :test 'equalp))
        (current-names nil) (current-type nil)
        (current-child nil) (current-parent nil))
    
    (with-open-file (stream dictionary-file :direction :input)
      (let ((token (get-cif-token stream nil nil)))
        (loop while (and token (not (eq (car token) :EOF)))
              do (case (car token)
                   
                   (:SAVE
                    (when current-type
                      (dolist (name current-names)
                        (setf (gethash name explicit-types) current-type)))
                    (when (and current-child current-parent)
                      (setf (gethash current-child parent-links) current-parent))
                    (setf current-names nil current-type nil current-child nil current-parent nil)
                    (setf token (get-cif-token stream nil nil)))

                   (:TAG
                    (let ((tag (cdr token)))
                      (setf token (get-cif-token stream nil nil))
                      (when (eq (car token) :VALUE)
                        (let ((val (fast-cif-trim *cif-trim-bag* (cdr token))))
                          (cond
                            ((string-equal tag "_item.name") (push val current-names))
                            ((string-equal tag "_item_type.code") (setf current-type val))
                            ((or (string-equal tag "_item_linked.child_name")
                                 (string-equal tag "_pdbx_item_linked_group_list.child_name"))
                             (setf current-child val))
                            ((or (string-equal tag "_item_linked.parent_name")
                                 (string-equal tag "_pdbx_item_linked_group_list.parent_name"))
                             (setf current-parent val))))
                        (setf token (get-cif-token stream nil nil)))))
                   
                   (:LOOP
                    (let ((headers nil))
                      (setf token (get-cif-token stream nil nil))
                      (loop while (eq (car token) :TAG)
                            do (push (cdr token) headers)
                               (setf token (get-cif-token stream nil nil)))
                      (setf headers (nreverse headers))
                      
                      (let ((child-idx (or (position "_item_linked.child_name" headers :test #'string-equal)
                                           (position "_pdbx_item_linked_group_list.child_name" headers :test #'string-equal)))
                            (parent-idx (or (position "_item_linked.parent_name" headers :test #'string-equal)
                                            (position "_pdbx_item_linked_group_list.parent_name" headers :test #'string-equal)))
                            (name-idx (position "_item.name" headers :test #'string-equal))
                            (type-idx (position "_item_type.code" headers :test #'string-equal)))
                        
                        (loop while (eq (car token) :VALUE)
                              do (let ((row (make-array (length headers))))
                                   (loop for i from 0 below (length headers)
                                         do (setf (aref row i) (fast-cif-trim *cif-trim-bag* (cdr token)))
                                            (setf token (get-cif-token stream nil nil)))
                                   
                                   ;; THE FIX: Isolate row types from frame types
                                   (if type-idx
                                       (when name-idx
                                         (setf (gethash (aref row name-idx) explicit-types) (aref row type-idx)))
                                       (when name-idx
                                         (push (aref row name-idx) current-names)))
                                   
                                   (when (and child-idx parent-idx)
                                     (setf (gethash (aref row child-idx) parent-links) (aref row parent-idx))))))))
                   
                   (t (setf token (get-cif-token stream nil nil)))))))

    (let ((final-dispatch nil))
      (maphash (lambda (item-name type-code)
                 (push (cons item-name (map-code-to-parser type-code)) final-dispatch))
               explicit-types)
      (maphash (lambda (child parent)
                 (unless (gethash child explicit-types)
                   (let ((resolved-type (find-inherited-type parent explicit-types parent-links)))
                     (when resolved-type
                       (push (cons child (map-code-to-parser resolved-type)) final-dispatch)))))
               parent-links)
      final-dispatch)))

(defun find-inherited-type (parent-name types-hash links-hash &optional visited)
  (cond
    ((member parent-name visited :test #'string=) nil) 
    ((gethash parent-name types-hash) (gethash parent-name types-hash))
    (t (let ((next-parent (gethash parent-name links-hash)))
         (if next-parent
             (find-inherited-type next-parent types-hash links-hash (cons parent-name visited))
             nil)))))

(defun map-code-to-parser (code)
  "Maps the complete mmcif_pdbx.dic type codes to strict cl-bio parsers."
  (cond
    ;; --- BOOLEANS ---
    ((string-equal code "boolean") 'parse-cif-boolean)

    ;; --- INTEGERS ---
    ((member code '("int" "positive_int") :test #'string-equal) 
     'strict-pdb-int)

    ;; --- FLOATS (Cons cell of value . uncertainty) ---
    ((string-equal code "float") 'parse-cif-float)

    ;; --- RANGES (Cons cells of min . max) ---
    ((member code '("int-range" "float-range") :test #'string-equal) 
     'parse-range)

    ;; --- MATRICES (Cons cell of array . array) ---
    ((string-equal code "3x4_matrix") 'parse-3x4-matrix)
    ((string-equal code "3x4_matrices") 'parse-3x4-matrices)

    ;; --- SYMBOLS (Fast EQ lookups for identifiers/codes) ---
    ;; We use symbols here because these fields have a strictly limited vocabulary,
    ;; and symbol comparison is O(1) inside struct logic.
    ((member code '("asym_id" "atcode" "code" "code30" "ec-type" "emd_id" 
                    "idname" "name" "pdb_id" "pdb_id_u" "point_group" 
                    "point_group_helical" "point_symmetry" "symop" 
                    "uchar1" "uchar3" "uchar5" "ucode" "uniprot_ptm_id") 
             :test #'string-equal) 
     'fast-pdb-symbol)

    ;; --- STRINGS (Identity fallback to prevent memory/symbol bloat) ---
    ;; Used for unbounded text, lists, dates, and sequences.
    ((member code '("any" "author" "binary" "citation_doi" "date_dep" 
                    "deposition_email" "email" "entity_id_list" "exp_data_doi" 
                    "fax" "id_list" "id_list_spc" "int_list" "line" 
                    "operation_expression" "orcid_id" "pdbx_PDB_obsoleted_db_id" 
                    "pdbx_related_db_id" "pdbx_wavelength_list" "phone" 
                    "seq-one-letter-code" "sequence_dep" "symmetry_operation" 
                    "text" "ucode-alphanum-csv" "uline" "yyyy-mm-dd" 
                    "yyyy-mm-dd:hh:mm" "yyyy-mm-dd:hh:mm-flex") 
             :test #'string-equal)
     'identity)

    ;; --- FALLBACK TRIPWIRE ---
    (t 
     (format t "[WARNING] Unknown dictionary type code '~A' defaulting to IDENTITY.~%" code)
     'identity)))

(defun process-dynamic-row (headers row-values)
  "Dynamically parses a row of mmCIF values using the compiled type dictionary.
   Returns an alist of typed keyword-value pairs."
  (let ((parsed-row nil))
    (loop for header in headers
          for val in row-values
          ;; Skip mmCIF missing value indicators
          unless (or (string= val ".") (string= val "?"))
          do (let* ((parser-symbol (gethash header *cif-type-hash*))
                    ;; Dispatch to the correct high-speed Lisp parser
                    (parsed-val (cond
                                  ((eq parser-symbol 'strict-pdb-int)
                                   (strict-pdb-int val 0 (length val)))
                                  ((eq parser-symbol 'strict-pdb-float)
                                   (strict-pdb-float val 0 (length val)))
                                  ((eq parser-symbol 'fast-pdb-symbol)
                                   (fast-pdb-symbol val 0 (length val)))
                                  ;; Fallback to the raw string if type is 'identity or unknown
                                  (t val))))
               
               ;; Intern the header as a keyword (e.g., :|_struct_title.title|) and pair it
               (push (cons (intern (string-upcase header) :keyword) parsed-val) parsed-row)))
    
    ;; Reverse to maintain the original column order
    (nreverse parsed-row)))

(defun process-atom-site-row (headers row-values entry)
  "Maps a row of mmCIF _atom_site values to our legacy pdb-atom struct."
  (let ((atom-number 0)
        (atom-name nil)
        (alt-loc nil)
        (residue-name nil)
        (chain-id nil)
        (residue-seq-number 0)
        (insertion-code nil)
        (x-coord 0.0d0)
        (y-coord 0.0d0)
        (z-coord 0.0d0)
        (occupancy 0.0d0)
        (temp-factor 0.0d0)
        (element-symbol nil)
        (charge nil)
        (valid-p nil)) ; Flag to ensure we actually parsed an atom ID
        
    (loop for header in headers
          for val in row-values
          ;; Skip mmCIF null indicators
          unless (or (string= val ".") (string= val "?"))
          do (cond
               ((string-equal header "_atom_site.id")
                (setf atom-number (parse-integer val) valid-p t))
               ((string-equal header "_atom_site.label_atom_id")
                (setf atom-name (intern val :keyword)))
               ((string-equal header "_atom_site.label_alt_id")
                (setf alt-loc (intern val :keyword)))
               ((string-equal header "_atom_site.label_comp_id")
                (setf residue-name (intern val :keyword)))
               ((string-equal header "_atom_site.label_asym_id")
                (setf chain-id (intern val :keyword)))
               ((string-equal header "_atom_site.label_seq_id")
                (setf residue-seq-number (parse-integer val)))
               ((string-equal header "_atom_site.pdbx_PDB_ins_code")
                (setf insertion-code (intern val :keyword)))
               ((string-equal header "_atom_site.Cartn_x")
                (setf x-coord (strict-pdb-float val 0 (length val))))
               ((string-equal header "_atom_site.Cartn_y")
                (setf y-coord (strict-pdb-float val 0 (length val))))
               ((string-equal header "_atom_site.Cartn_z")
                (setf z-coord (strict-pdb-float val 0 (length val))))
               ((string-equal header "_atom_site.occupancy")
                (setf occupancy (strict-pdb-float val 0 (length val))))
               ((string-equal header "_atom_site.B_iso_or_equiv")
                (setf temp-factor (strict-pdb-float val 0 (length val))))
               ((string-equal header "_atom_site.type_symbol")
                (setf element-symbol (intern val :keyword)))
               ((string-equal header "_atom_site.pdbx_formal_charge")
                (setf charge (parse-integer val)))))
    
    ;; Only store if it successfully parsed a valid ID
    (when valid-p
      (let ((atom (make-pdb-atom :rec-name :ATOM
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
        (setf (gethash atom-number (atom-hash entry)) atom)))))

(defun read-cif-stream (stream &key (populate-legacy t))
  "The streamlined state machine: fast-paths atoms, dynamically catches everything else,
   and optionally populates legacy structures."
  (let* ((*cif-buffer* (make-array 1024 :element-type 'character :adjustable t :fill-pointer 0))
         (entry (make-instance 'pdb-entry)) 
         (current-token (get-cif-token stream nil nil)))
    
    (loop while (and current-token (not (eq (car current-token) :EOF)))
          do (case (car current-token)
               
               ;; 1. Flat Metadata Tags
               (:TAG
                (let ((tag (cdr current-token)))
                  (setf current-token (get-cif-token stream nil nil))
                  (when (eq (car current-token) :VALUE)
                    (let ((val (cdr current-token)))
                      (unless (or (string= val ".") (string= val "?"))
                        (let* ((parser-symbol (gethash tag *cif-type-hash*))
                               (parsed-val (cond
                                             ((eq parser-symbol 'strict-pdb-int) (strict-pdb-int val 0 (length val)))
                                             ((eq parser-symbol 'strict-pdb-float) (strict-pdb-float val 0 (length val)))
                                             ((eq parser-symbol 'fast-pdb-symbol) (fast-pdb-symbol val 0 (length val)))
                                             (t val))))
                          (setf (gethash (intern (string-upcase tag) :keyword) (cif-metadata entry)) parsed-val))))
                    (setf current-token (get-cif-token stream nil nil)))))
               
               ;; 2. Loops (The Magic Split)
               (:LOOP
                (let ((headers nil))
                  (setf current-token (get-cif-token stream nil nil))
                  (loop while (eq (car current-token) :TAG)
                        do (push (cdr current-token) headers)
                           (setf current-token (get-cif-token stream nil nil)))
                  (setf headers (nreverse headers))
                  
                  (let* ((num-headers (length headers))
                         (first-header (first headers))
                         ;; IS IT AN ATOM?
                         (is-atom (and first-header (>= (length first-header) 11)
                                       (string-equal first-header "_atom_site." :end1 11))))
                    
                    (let ((dynamic-rows nil))
                      (loop while (eq (car current-token) :VALUE)
                            do (let ((row-values nil))
                                 (loop repeat num-headers
                                       do (push (cdr current-token) row-values)
                                          (setf current-token (get-cif-token stream nil nil)))
                                 (setf row-values (nreverse row-values))
                                 
                                 (if is-atom
                                     ;; THE FAST PATH
                                     (process-atom-site-row headers row-values entry)
                                     ;; THE DYNAMIC PATH
                                     (push (process-dynamic-row headers row-values) dynamic-rows))))
                      
                      ;; Group dynamically parsed loops by their prefix
                      (when dynamic-rows
                        (let* ((dot-pos (position #\. first-header))
                               (prefix (if dot-pos 
                                           (intern (string-upcase (subseq first-header 0 dot-pos)) :keyword)
                                           :UNKNOWN_LOOP)))
                          (setf (gethash prefix (cif-metadata entry)) (nreverse dynamic-rows))))))))
               
               (t (setf current-token (get-cif-token stream nil nil)))))
    
    ;; 3. THE OPTIONAL INTERPRETER STEP
    (when populate-legacy
      (populate-legacy-structs entry))
    
    entry))

(defun read-cif-file (file &key (populate-legacy t))
  "Reads an mmCIF file and returns a populated pdb-entry.
   If :populate-legacy is nil, skips building backward-compatible lists (helices, sheets, etc.) for extra speed."
  (with-open-file (stream file :direction :input)
    (read-cif-stream stream :populate-legacy populate-legacy)))

(defun parse-cif-atoms-fast (file-path)
  "A high-speed, zero-metadata mmCIF parser that strictly extracts PDB-ATOM structs 
   into a hash table for rapid data crunching."
  (let ((*cif-buffer* (make-array 1024 :element-type 'character :adjustable t :fill-pointer 0))
        (*cif-string-pool* (make-hash-table :test 'equal))
        (atoms (make-hash-table :test 'eql))
        (atom-count 0))
    
    (with-open-file (stream file-path :direction :input)
      (let ((token (get-cif-token stream nil nil)))
        
        (loop while (and token (not (eq (car token) :EOF)))
              do (case (car token)
                   
                   (:LOOP
                    (let ((headers nil))
                      (setf token (get-cif-token stream nil nil))
                      (loop while (eq (car token) :TAG)
                            do (push (cdr token) headers)
                               (setf token (get-cif-token stream nil nil)))
                      (setf headers (nreverse headers))
                      
                      ;; TARGET ACQUIRED: Is this the atom_site loop?
                      (if (and headers (>= (length (first headers)) 11)
                               (string-equal (subseq (first headers) 0 11) "_atom_site."))
                          
                          ;; --- THE FAST PATH ---
                          (let* ((h-len (length headers))
                                 ;; Prefer 'auth' (author/legacy PDB) over 'label' (internal mmCIF)
                                 (idx-id     (position "_atom_site.id" headers :test #'string-equal))
                                 (idx-name   (or (position "_atom_site.auth_atom_id" headers :test #'string-equal)
                                                 (position "_atom_site.label_atom_id" headers :test #'string-equal)))
                                 (idx-alt    (position "_atom_site.label_alt_id" headers :test #'string-equal))
                                 (idx-res    (or (position "_atom_site.auth_comp_id" headers :test #'string-equal)
                                                 (position "_atom_site.label_comp_id" headers :test #'string-equal)))
                                 (idx-chain  (or (position "_atom_site.auth_asym_id" headers :test #'string-equal)
                                                 (position "_atom_site.label_asym_id" headers :test #'string-equal)))
                                 (idx-seq    (or (position "_atom_site.auth_seq_id" headers :test #'string-equal)
                                                 (position "_atom_site.label_seq_id" headers :test #'string-equal)))
                                 (idx-ins    (position "_atom_site.pdbx_PDB_ins_code" headers :test #'string-equal))
                                 (idx-x      (position "_atom_site.Cartn_x" headers :test #'string-equal))
                                 (idx-y      (position "_atom_site.Cartn_y" headers :test #'string-equal))
                                 (idx-z      (position "_atom_site.Cartn_z" headers :test #'string-equal))
                                 (idx-occ    (position "_atom_site.occupancy" headers :test #'string-equal))
                                 (idx-b      (position "_atom_site.B_iso_or_equiv" headers :test #'string-equal))
                                 (idx-elem   (position "_atom_site.type_symbol" headers :test #'string-equal))
                                 (idx-charge (position "_atom_site.pdbx_formal_charge" headers :test #'string-equal)))
                            
                            (loop while (eq (car token) :VALUE)
                                  do (let ((row (make-array h-len)))
                                       (loop for i from 0 below h-len
                                             do (setf (aref row i) (fast-cif-trim *cif-trim-bag* (cdr token)))
                                                (setf token (get-cif-token stream nil nil)))
                                       
                                       ;; High-speed inline data validation and struct initialization
                                       (macrolet ((get-val (idx)
                                                    `(and ,idx (let ((v (aref row ,idx)))
                                                                 (unless (or (string= v "?") (string= v ".")) v))))
                                                  (fast-float (idx)
                                                    `(let ((v (get-val ,idx)))
                                                       (if v (let ((*read-default-float-format* 'double-float))
                                                               (coerce (read-from-string v) 'double-float)) 
                                                           0.0d0))))
                                         
                                         (let* ((atom-num (let ((v (get-val idx-id))) (if v (or (parse-integer v :junk-allowed t) 0) 0)))
                                                (atom (make-pdb-atom
                                                       :atom-number        atom-num
                                                       :atom-name          (let ((v (get-val idx-name))) (if v (pool-string v)))
                                                       :alt-loc            (let ((v (get-val idx-alt))) (if v (pool-string v)))
                                                       :residue-name       (let ((v (get-val idx-res))) (if v (pool-string v)))
                                                       :chain-id           (let ((v (get-val idx-chain))) (if v (pool-string v)))
                                                       :residue-seq-number (let ((v (get-val idx-seq))) (if v (or (parse-integer v :junk-allowed t) 0) 0))
                                                       :insertion-code     (let ((v (get-val idx-ins))) (if v (pool-string v)))
                                                       :x-coord            (fast-float idx-x)
                                                       :y-coord            (fast-float idx-y)
                                                       :z-coord            (fast-float idx-z)
                                                       :occupancy          (fast-float idx-occ)
                                                       :temp-factor        (fast-float idx-b)
                                                       :element-symbol     (let ((v (get-val idx-elem))) (if v (pool-string v)))
                                                       :charge             (let ((v (get-val idx-charge))) (if v (parse-integer v :junk-allowed t) 0)))))
                                           
                                           (setf (gethash atom-num atoms) atom)
                                           (incf atom-count))))))
                          
                          ;; --- NOT ATOMS: FAST FORWARD ---
                          (loop while (eq (car token) :VALUE)
                                do (setf token (get-cif-token stream nil nil))))))
                   
                   ;; Ignore flat tags and saves
                   (t (setf token (get-cif-token stream nil nil)))))))
    
    (format t "-> Fast-tracked ~D atoms directly to hash!~%" atom-count)
    atoms))

(defun populate-legacy-connections (entry meta)
  "Extracts SSBOND and LINK records from mmCIF _struct_conn."
  (loop for row in (gethash :|_STRUCT_CONN| meta)
        for conn-type = (cdr (assoc :|_STRUCT_CONN.CONN_TYPE_ID| row))
        do (let ((res-1   (cdr (assoc :|_STRUCT_CONN.PTNR1_LABEL_COMP_ID| row)))
                 (chain-1 (cdr (assoc :|_STRUCT_CONN.PTNR1_LABEL_ASYM_ID| row)))
                 (seq-1   (or (cdr (assoc :|_STRUCT_CONN.PTNR1_LABEL_SEQ_ID| row)) 0))
                 (icode-1 (cdr (assoc :|_STRUCT_CONN.PDBX_PTNR1_PDB_INS_CODE| row)))
                 (atom-1  (cdr (assoc :|_STRUCT_CONN.PTNR1_LABEL_ATOM_ID| row)))
                 (sym-1   (cdr (assoc :|_STRUCT_CONN.PTNR1_SYMMETRY| row)))
                 
                 (res-2   (cdr (assoc :|_STRUCT_CONN.PTNR2_LABEL_COMP_ID| row)))
                 (chain-2 (cdr (assoc :|_STRUCT_CONN.PTNR2_LABEL_ASYM_ID| row)))
                 (seq-2   (or (cdr (assoc :|_STRUCT_CONN.PTNR2_LABEL_SEQ_ID| row)) 0))
                 (icode-2 (cdr (assoc :|_STRUCT_CONN.PDBX_PTNR2_PDB_INS_CODE| row)))
                 (atom-2  (cdr (assoc :|_STRUCT_CONN.PTNR2_LABEL_ATOM_ID| row)))
                 (sym-2   (cdr (assoc :|_STRUCT_CONN.PTNR2_SYMMETRY| row)))
                 
                 (length  (or (cdr (assoc :|_STRUCT_CONN.PDBX_DIST_VALUE| row)) 0.0d0)))
             (cond
               ((eq conn-type :|disulf|)
                (push (make-pdb-ssbond
                       :ser-num   (cif-int row :|_STRUCT_CONN.ID|)
                       :res-name1 (cif-sym row :|_STRUCT_CONN.PTNR1_AUTH_COMP_ID|)
                       :chain-id1 (cif-sym row :|_STRUCT_CONN.PTNR1_AUTH_ASYM_ID|)
                       :seq-num1  (cif-int row :|_STRUCT_CONN.PTNR1_AUTH_SEQ_ID|)
                       :icode1    (cif-sym row :|_STRUCT_CONN.PDBX_PTNR1_PDB_INS_CODE|)
                       :res-name2 (cif-sym row :|_STRUCT_CONN.PTNR2_AUTH_COMP_ID|)
                       :chain-id2 (cif-sym row :|_STRUCT_CONN.PTNR2_AUTH_ASYM_ID|)
                       :seq-num2  (cif-int row :|_STRUCT_CONN.PTNR2_AUTH_SEQ_ID|)
                       :icode2    (cif-sym row :|_STRUCT_CONN.PDBX_PTNR2_PDB_INS_CODE|)
                       :sym1      (cif-str row :|_STRUCT_CONN.PTNR1_SYMMETRY|)
                       :sym2      (cif-str row :|_STRUCT_CONN.PTNR2_SYMMETRY|)
                       :length    (cif-float row :|_STRUCT_CONN.CONN_TYPE_ID_DIST|))
                      (ssbonds entry)))
               (conn-type
                (push (make-pdb-link :res-name-1 res-1 :chain-1 chain-1 :seq-num-1 seq-1 :icode-1 icode-1 :atom-1 atom-1
                                     :res-name-2 res-2 :chain-2 chain-2 :seq-num-2 seq-2 :icode-2 icode-2 :atom-2 atom-2
                                     :sym-1 sym-1 :sym-2 sym-2 :length length)
                      (links entry))))))
  (setf (ssbonds entry) (nreverse (ssbonds entry))
        (links entry)   (nreverse (links entry))))

(defun populate-legacy-helices (entry meta)
  "Extracts HELIX records from mmCIF _struct_conf."
  (loop for row in (gethash :|_STRUCT_CONF| meta)
        for conf-type = (string (cdr (assoc :|_STRUCT_CONF.CONF_TYPE_ID| row)))
        when (and (>= (length conf-type) 4) (string-equal (subseq conf-type 0 4) "HELX"))
          do (push (make-pdb-helix
                    :ser-num       (cif-int row :|_STRUCT_CONF.ID|)
                    :helix-id      (cif-sym row :|_STRUCT_CONF.PDBX_PDB_HELIX_ID|)
                    :init-res-name (cif-sym row :|_STRUCT_CONF.BEG_AUTH_COMP_ID|)
                    :init-chain-id (cif-sym row :|_STRUCT_CONF.BEG_AUTH_ASYM_ID|)
                    :init-seq-num  (cif-int row :|_STRUCT_CONF.BEG_AUTH_SEQ_ID|)
                    :init-icode    (cif-sym row :|_STRUCT_CONF.PDBX_BEG_PDB_INS_CODE|)
                    :end-res-name  (cif-sym row :|_STRUCT_CONF.END_AUTH_COMP_ID|)
                    :end-chain-id  (cif-sym row :|_STRUCT_CONF.END_AUTH_ASYM_ID|)
                    :end-seq-num   (cif-int row :|_STRUCT_CONF.END_AUTH_SEQ_ID|)
                    :end-icode     (cif-sym row :|_STRUCT_CONF.PDBX_END_PDB_INS_CODE|)
                    :helix-class   (cif-int row :|_STRUCT_CONF.PDBX_PDB_HELIX_CLASS|)
                    :comment       (cif-str row :|_STRUCT_CONF.DETAILS|)
                    :length        (cif-int row :|_STRUCT_CONF.PDBX_PDB_HELIX_LENGTH|))
                   (helices entry)))
  (setf (helices entry) (nreverse (helices entry))))

(defun populate-legacy-sheets (entry meta)
  "Extracts SHEET records by joining _struct_sheet_range with _struct_sheet_hbond."
  (let ((ranges (gethash :|_STRUCT_SHEET_RANGE| meta))
        (hbonds (gethash :|_STRUCT_SHEET_HBOND| meta)))
    
    (loop for row in ranges
          for sheet-id = (cif-sym row :|_STRUCT_SHEET_RANGE.SHEET_ID|)
          for strand-id = (cif-sym row :|_STRUCT_SHEET_RANGE.ID|)
          
          ;; Perform a relational join: Find the H-bond row where this strand is "range 2"
          do (let ((hb-row (find-if (lambda (hb)
                                      (and (eq (cif-sym hb :|_STRUCT_SHEET_HBOND.SHEET_ID|) sheet-id)
                                           (eq (cif-sym hb :|_STRUCT_SHEET_HBOND.RANGE_ID_2|) strand-id)))
                                    hbonds)))
               
               (push (make-pdb-sheet
                      :strand         (cif-int row :|_STRUCT_SHEET_RANGE.ID|)
                      :sheet-id       sheet-id
                      :num-strands    (cif-int row :|_STRUCT_SHEET.NUMBER_STRANDS|)
                      :init-res-name  (cif-sym row :|_STRUCT_SHEET_RANGE.BEG_AUTH_COMP_ID|)
                      :init-chain-id  (cif-sym row :|_STRUCT_SHEET_RANGE.BEG_AUTH_ASYM_ID|)
                      :init-seq-num   (cif-int row :|_STRUCT_SHEET_RANGE.BEG_AUTH_SEQ_ID|)
                      :init-icode     (cif-sym row :|_STRUCT_SHEET_RANGE.PDBX_BEG_PDB_INS_CODE|)
                      :end-res-name   (cif-sym row :|_STRUCT_SHEET_RANGE.END_AUTH_COMP_ID|)
                      :end-chain-id   (cif-sym row :|_STRUCT_SHEET_RANGE.END_AUTH_ASYM_ID|)
                      :end-seq-num    (cif-int row :|_STRUCT_SHEET_RANGE.END_AUTH_SEQ_ID|)
                      :end-icode      (cif-sym row :|_STRUCT_SHEET_RANGE.PDBX_END_PDB_INS_CODE|)
                      
                      ;; Extract H-bond data if the join was successful
                      :cur-atom       (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_1_AUTH_ATOM_ID|)
                      :cur-res-name   (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_1_AUTH_COMP_ID|)
                      :cur-chain-id   (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_1_AUTH_ASYM_ID|)
                      :cur-res-seq    (cif-int hb-row :|_STRUCT_SHEET_HBOND.RANGE_1_AUTH_SEQ_ID|)
                      :cur-icode      (cif-sym hb-row :|_STRUCT_SHEET_HBOND.PDBX_RANGE_1_PDB_INS_CODE|)
                      :prev-atom      (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_2_AUTH_ATOM_ID|)
                      :prev-res-name  (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_2_AUTH_COMP_ID|)
                      :prev-chain-id  (cif-sym hb-row :|_STRUCT_SHEET_HBOND.RANGE_2_AUTH_ASYM_ID|)
                      :prev-res-seq   (cif-int hb-row :|_STRUCT_SHEET_HBOND.RANGE_2_AUTH_SEQ_ID|)
                      :prev-icode     (cif-sym hb-row :|_STRUCT_SHEET_HBOND.PDBX_RANGE_2_PDB_INS_CODE|))
                     (sheets entry))))
    (setf (sheets entry) (nreverse (sheets entry)))))

(defun populate-legacy-sites (entry meta)
  "Extracts SITE records from mmCIF _struct_site_gen."
  (loop for row in (gethash :|_STRUCT_SITE_GEN| meta)
        for site-id = (cdr (assoc :|_STRUCT_SITE_GEN.SITE_ID| row))
        when site-id
        do (let* ((res-name (cdr (assoc :|_STRUCT_SITE_GEN.LABEL_COMP_ID| row)))
                  (chain    (cdr (assoc :|_STRUCT_SITE_GEN.LABEL_ASYM_ID| row)))
                  (seq-num  (or (cdr (assoc :|_STRUCT_SITE_GEN.LABEL_SEQ_ID| row)) 0))
                  (icode    (cdr (assoc :|_STRUCT_SITE_GEN.PDBX_PDB_INS_CODE| row)))
                  (residue-plist (list :res-name res-name :chain chain :seq-num seq-num :icode icode))
                  (existing-site (find site-id (sites entry) :key #'site-site-id)))
             (if existing-site
                 (progn
                   (incf (site-num-res existing-site))
                   (push residue-plist (site-residues existing-site)))
                 (push (make-pdb-site :site-id site-id :num-res 1 :residues (list residue-plist))
                       (sites entry)))))
  (setf (sites entry) (nreverse (sites entry)))
  (loop for site in (sites entry)
        do (setf (site-residues site) (nreverse (site-residues site)))))

(defun populate-legacy-cispeps (entry meta)
  "Extracts CISPEP records from mmCIF _struct_mon_prot_cis."
  (loop for row in (gethash :|_STRUCT_MON_PROT_CIS| meta)
        for pep1 = (cdr (assoc :|_STRUCT_MON_PROT_CIS.LABEL_COMP_ID| row))
        when pep1
        do (push (make-pdb-cispep 
                  :ser-num   (cif-int row :|_STRUCT_MON_PROT_CIS.PDBX_ID|)
                  :pep1      (cif-sym row :|_STRUCT_MON_PROT_CIS.AUTH_COMP_ID|)
                  :chain-id1 (cif-sym row :|_STRUCT_MON_PROT_CIS.AUTH_ASYM_ID|)
                  :seq-num1  (cif-int row :|_STRUCT_MON_PROT_CIS.AUTH_SEQ_ID|)
                  :icode1    (cif-sym row :|_STRUCT_MON_PROT_CIS.PDBX_PDB_INS_CODE|)
                  ;; Note: CISPEP partner 2 data might be handled differently depending on if you are 
                  ;; pulling from struct_mon_prot_cis or struct_conn. Update these keys as needed!
                  :pep2      (or (cif-sym row :|_STRUCT_MON_PROT_CIS.AUTH_COMP_ID_2|)
                                 (cif-sym row :|_PARTNER_2.AUTH_COMP_ID|))
                  :chain-id2 (or (cif-sym row :|_STRUCT_MON_PROT_CIS.AUTH_ASYM_ID_2|)
                                 (cif-sym row :|_PARTNER_2.AUTH_ASYM_ID|))
                  :seq-num2  (or (cif-int row :|_STRUCT_MON_PROT_CIS.AUTH_SEQ_ID_2|)
                                 (cif-int row :|_PARTNER_2.AUTH_SEQ_ID|))
                  :icode2    (or (cif-sym row :|_STRUCT_MON_PROT_CIS.PDBX_PDB_INS_CODE_2|)
                                 (cif-sym row :|_PARTNER_2.PDBX_PDB_INS_CODE|))
                  :mod-num   (cif-int row :|_STRUCT_MON_PROT_CIS.PDBX_PDB_MODEL_NUM|)
                  :measure   (cif-float row :|_STRUCT_MON_PROT_CIS.PDBX_OMEGA_ANGLE|))
                 (cispeps entry)))
  (setf (cispeps entry) (nreverse (cispeps entry))))

(defun populate-legacy-structs (entry)
  "Master interpreter: reads dynamically parsed mmCIF metadata and delegates 
   to specific helpers to populate legacy PDB lists."
  (let ((meta (cif-metadata entry)))
    (populate-legacy-connections entry meta)
    (populate-legacy-helices entry meta)
    (populate-legacy-sheets entry meta)
    (populate-legacy-sites entry meta)
    (populate-legacy-cispeps entry meta))
  entry)

(defun update-local-dictionary (dictionary-file &key (output-file "io/cif-dictionary.lisp"))
  "Reads the official mmCIF dictionary and generates a permanent Lisp source file 
   containing the compiled type mappings."
  (let ((alist (generate-type-dispatch-table dictionary-file)))
    (with-open-file (stream output-file 
                            :direction :output 
                            :if-exists :supersede
                            :if-does-not-exist :create)
      
      ;; 1. Write the File Header
      (format stream ";;; -------------------------------------------------------------------------~%")
      (format stream ";;; THIS FILE IS AUTO-GENERATED BY CL-BIO.~%")
      (format stream ";;; DO NOT EDIT DIRECTLY.~%")
      (format stream ";;; -------------------------------------------------------------------------~%~%")
      (format stream "(in-package :bio)~%~%")
      
      ;; 2. Write the Alist Variable
      (format stream "(defparameter *cif-type-alist*~%")
      (format stream "  '(~%")
      (loop for (key . val) in alist
            do (format stream "    (\"~A\" . ~S)~%" key val))
      (format stream "   )~%")
      (format stream "  \"Raw association list mapping mmCIF keys to Lisp parser functions.\")~%~%")
      
      ;; 3. Write the Hash Table Generator for fast runtime lookups
      (format stream "(defvar *cif-type-hash* (make-hash-table :test 'equalp)~%")
      (format stream "  \"O(1) lookup table for mmCIF data types.\")~%~%")
      (format stream "(loop for (key . val) in *cif-type-alist*~%")
      (format stream "      do (setf (gethash key *cif-type-hash*) val))~%~%")
      
      (format t "Successfully generated ~A with ~D type definitions!~%" output-file (length alist)))))
