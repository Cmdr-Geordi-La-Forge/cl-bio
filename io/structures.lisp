;;; io/structures.lisp
(in-package :bio)

(defstruct (cif-struct-conn (:type list))
  "Native mmCIF representation of any connection between two entities."
  conn-type-id
  res-name-1 chain-1 seq-num-1 icode-1 atom-1 sym-1
  res-name-2 chain-2 seq-num-2 icode-2 atom-2 sym-2
  length)

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
   (links :initarg :links :accessor links :initform nil)
   (cif-metadata :initarg :cif-metadata :accessor cif-metadata 
                 :initform (make-hash-table :test 'eq))
   (struct-conns :initarg :struct-conns :accessor struct-conns :initform nil)))

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

(defstruct (pdb-site (:conc-name site-))
  (seq-num 0 :type fixnum) site-id (num-res 0 :type fixnum) residues)

(defstruct (pdb-helix (:conc-name helix-))
  (ser-num 0 :type fixnum) helix-id 
  init-res-name init-chain-id (init-seq-num 0 :type fixnum) init-icode
  end-res-name end-chain-id (end-seq-num 0 :type fixnum) end-icode
  (helix-class 0 :type fixnum) comment (length 0 :type (or fixnum null)))

(defstruct (pdb-sheet (:conc-name sheet-))
  (strand 0 :type fixnum) sheet-id (num-strands 0 :type fixnum)
  init-res-name init-chain-id (init-seq-num 0 :type fixnum) init-icode
  end-res-name end-chain-id (end-seq-num 0 :type fixnum) end-icode
  (sense 0 :type fixnum)
  cur-atom cur-res-name cur-chain-id (cur-res-seq 0 :type (or fixnum null)) cur-icode
  prev-atom prev-res-name prev-chain-id (prev-res-seq 0 :type (or fixnum null)) prev-icode)

(defstruct (pdb-cispep (:conc-name cispep-))
  (ser-num 0 :type fixnum)
  pep1 chain-id1 (seq-num1 0 :type fixnum) icode1
  pep2 chain-id2 (seq-num2 0 :type fixnum) icode2
  (mod-num 0 :type (or fixnum null))
  (measure 0.0d0 :type (or double-float null)))

(defstruct (pdb-ssbond (:conc-name ssbond-))
  (ser-num 0 :type fixnum)
  res-name1 chain-id1 (seq-num1 0 :type fixnum) icode1
  res-name2 chain-id2 (seq-num2 0 :type fixnum) icode2
  sym1 sym2 (length 0.0d0 :type (or double-float null)))

(defstruct (pdb-link (:conc-name link-))
  atom1 alt-loc1 res-name1 chain-id1 (seq-num1 0 :type fixnum) icode1
  atom2 alt-loc2 res-name2 chain-id2 (seq-num2 0 :type fixnum) icode2
  sym1 sym2 (length 0.0d0 :type (or double-float null)))
