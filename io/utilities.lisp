;;; Parsing utilities
;;;
;;; Copyright (c) 2006 Cyrus Harmon (ch-lisp@bobobeach.com)
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

(in-package :bio)

;; Shared string buffer for fast symbol interning
(defvar *pdb-string-buffer* (make-string 80 :initial-element #\Space))
(declaim (special *pdb-string-buffer*))

(defun check-end-of-sequence (stream &key (end-char #\>))
  (let ((char (peek-char t stream nil :eof)))
    (when (or (equal char end-char)
              (equal char :eof))
      :eof)))

(defun read-until-char (stream end-char)
  (flet ((read-char-if-not (stream end-char)
           (unless (eq (check-end-of-sequence stream :end-char end-char) :eof)
             (read-char stream))))
    (let ((l))
      (do ((char (read-char-if-not stream end-char) (read-char-if-not stream end-char)))
          ((null char))
        (push char l))
      (coerce (nreverse l) 'string))))

(defun check-end-of-file (stream)
  (equal (peek-char t stream nil :eof) :eof))

(defparameter *line-ending-chars* '(#\Newline #\Return))

(defun read-dos-line (stream)
  (let ((l))
    (do ((char (read-char stream) (read-char stream)))
        ((or (null char)
             (member char *line-ending-chars*)))
      (push char l))
    (let ((p (peek-char nil stream)))
      (cond ((char-equal p #\Newline)
             (read-char stream))))
    (coerce (nreverse l) 'string)))

(defun check-end-of-line (stream)
  (member (peek-char nil stream nil #\Newline) *line-ending-chars*))

(defun remove-initial-spaces (string)
  (string-left-trim '(#\Space #\Tab #\Newline) string))

(defun remove-trailing-spaces (string)
  (string-right-trim '(#\Space #\Tab #\Newline) string))

(defun fast-cif-trim (str)
  "A high-performance trim that avoids allocating memory if the string is already clean."
  (let ((len (length str)))
    (if (or (zerop len)
            (not (or (member (char str 0) *cif-trim-bag*)
                     (member (char str (1- len)) *cif-trim-bag*))))
        str
        (string-trim *cif-trim-bag* str))))

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
