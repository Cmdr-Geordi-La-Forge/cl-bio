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


(in-package #:cl-user)

(defpackage #:bio (:use #:cl)
              (:export

               ;; objects
               #:bio-object
               #:bio-object-set
               #:members
	       
	       #:pdb-entry
	       #:classification
	       #:dep-date
	       #:id-code
	       #:obsolete
	       #:title
	       #:molecules
	       #:chains
	       #:atom-hash
	       
	       #:pdb-atom
	       #:record-name
	       #:atom-number
	       #:atom-name
	       #:alt-loc
	       #:residue-name
	       #:chain-id
	       #:residue-seq-number
	       #:insertation-code
	       #:x-coord
	       #:y-coord
	       #:z-coord
	       #:occupancy
	       #:temp-factor
	       #:element-symbol
	       #:charge
               
               ;; ranges
               #:range
               #:ds-range
             
               #:range-equal
               #:range-contains
               #:range-start
               #:range-end
               #:range-min
               #:range-max
             
               #:+plus-strand+
               #:+unknown-strand+
               #:+minus-strand+
               #:+both-strands+

               ;; descritptors
               #:descriptor
               #:get-descriptors
               #:add-descriptor
               #:remove-descriptor
               
               ;; identifiers
               #:identifier
               #:identifiers
               #:id
               #:type
               #:version
               #:authority

               #:ncbi-gi
               #:get-ncbi-gis

               #:refseq-id
               #:get-refseq-ids
               #:genbank-accession
               #:affymetrix-probe-set-id
               #:flybase-identifier
               #:flybase-gene-identifier
               
               ;; bio-sequences
               #:bio-sequence
               #:residue
               #:residues-string
               #:residues-string-range
               #:seq-length
               #:copy-sequence-range

               ;; sequences with residues
               #:sequence-with-residues
               #:seq-reverse

               ;; sequences with residue codes
               #:sequence-with-residue-codes
               #:residue-code
               #:sequence-encoding-error

               ;; annotated sequences
               #:annotated-sequence

               ;; simple sequences
               #:simple-sequence

               ;; adjustable sequences
               #:adjustable-sequence
               #:insert-residue
               #:insert-residues
               #:insert-residue-codes
               #:append-residue
               #:append-residues
               #:append-residue-codes
               #:delete-residue
               #:delete-residues

               ;; nucleic acid sequences
               #:na-sequence
               #:na-sequence-with-residues

               ;; DNA sequences
               #:dna-sequence
               #:dna-sequence-with-residues
               #:simple-dna-sequence
               #:adjustable-dna-sequence
               #:reverse-complement
               #:make-simple-dna-sequence
               #:make-adjustable-dna-sequence
               #:make-dna-sequence-from-string
               #:make-random-dna-sequence

               ;; RNA sequences
               #:rna-sequence
               #:simple-rna-sequence
               #:adjustable-rna-sequence
               #:rna-sequence-with-residues
               #:make-simple-rna-sequence
               #:make-adjustable-rna-sequence
               #:make-rna-sequence-from-string
               #:make-random-rna-sequence

               ;; amino acid sequences
               #:aa-sequence
               #:simple-aa-sequence
               #:adjustable-aa-sequence
               #:aa-sequence-with-residues
               #:make-simple-aa-sequence
               #:make-adjustable-aa-sequence
               #:make-aa-sequence-from-string
               #:make-random-aa-sequence

               ;; genes
               #:gene
               #:gene-set
               #:genes
               #:gene-type
               #:gene-source
               #:gene-product
               #:gene-products
               #:gene-summary
               
               #:annotation
               #:annotations
               #:get-annotations
               #:exon
               #:cds
               #:sts
               #:repeat-region
               #:source
               #:region

               #:get-genbank-accessions

               #:simple-pairwise-alignment
               #:alpha-sequence
               #:alpha-range
               #:beta-sequence
               #:beta-range
               #:filter-alignments

               ;; artciles
               #:author
               #:author-last-name
               #:author-forenames
               #:author-initials
               
               #:article
               #:article-pmid
               #:article-title
               #:article-authors
               #:article-journal
               #:article-short-journal
               #:article-date
               #:article-volume
               #:article-issue
               #:article-pages
               #:article-abstract
               #:article-mesh-headings
               #:article-doi
               #:article-affiliation

               #:article-set
               #:article-set-articles
               #:article-set-query
               #:article-set-total
               #:article-set-count
               #:article-set-start
               
               ;; utilities
               #:split-string-into-lines-list
               #:split-string-into-lines

               ;; dictionary
               #:lookup
               #:fetch
               
               #:dna->rna
               #:rna->dna
               #:translate

               #:3-letter-residues-string
               #:3-letter-residues-list

               #:read-fasta-sequences
               #:read-fasta-file
               #:write-fasta-file

	       #:read-pdb-file))

(defpackage #:bio-user (:use #:cl #:bio))
