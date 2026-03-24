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
           #:resolution
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
	       #:insertion-code
	       #:x-coord
	       #:y-coord
	       #:z-coord
	       #:occupancy
	       #:temp-factor
	       #:element-symbol
	       #:charge
           
           ;; -----------------------------------------------------------------
           ;; PDB-ENTRY TOPOLOGY LISTS
           ;; -----------------------------------------------------------------
           #:sites #:helices #:sheets #:cispeps #:ssbonds #:links

           ;; -----------------------------------------------------------------
           ;; SITE ACCESSORS
           ;; -----------------------------------------------------------------
           #:site-seq-num #:site-site-id #:site-num-res #:site-residues

           ;; -----------------------------------------------------------------
           ;; HELIX ACCESSORS
           ;; -----------------------------------------------------------------
           #:helix-ser-num #:helix-helix-id 
           #:helix-init-res-name #:helix-init-chain-id #:helix-init-seq-num #:helix-init-icode
           #:helix-end-res-name #:helix-end-chain-id #:helix-end-seq-num #:helix-end-icode
           #:helix-helix-class #:helix-comment #:helix-length

           ;; -----------------------------------------------------------------
           ;; SHEET ACCESSORS
           ;; -----------------------------------------------------------------
           #:sheet-strand #:sheet-sheet-id #:sheet-num-strands
           #:sheet-init-res-name #:sheet-init-chain-id #:sheet-init-seq-num #:sheet-init-icode
           #:sheet-end-res-name #:sheet-end-chain-id #:sheet-end-seq-num #:sheet-end-icode
           #:sheet-sense
           #:sheet-cur-atom #:sheet-cur-res-name #:sheet-cur-chain-id #:sheet-cur-res-seq #:sheet-cur-icode
           #:sheet-prev-atom #:sheet-prev-res-name #:sheet-prev-chain-id #:sheet-prev-res-seq #:sheet-prev-icode

           ;; -----------------------------------------------------------------
           ;; SSBOND ACCESSORS
           ;; -----------------------------------------------------------------
           #:ssbond-ser-num
           #:ssbond-res-name1 #:ssbond-chain-id1 #:ssbond-seq-num1 #:ssbond-icode1
           #:ssbond-res-name2 #:ssbond-chain-id2 #:ssbond-seq-num2 #:ssbond-icode2
           #:ssbond-sym1 #:ssbond-sym2 #:ssbond-length

           ;; -----------------------------------------------------------------
           ;; LINK ACCESSORS
           ;; -----------------------------------------------------------------
           #:link-atom1 #:link-alt-loc1 #:link-res-name1 #:link-chain-id1 #:link-seq-num1 #:link-icode1
           #:link-atom2 #:link-alt-loc2 #:link-res-name2 #:link-chain-id2 #:link-seq-num2 #:link-icode2
           #:link-sym1 #:link-sym2 #:link-length

           ;; -----------------------------------------------------------------
           ;; CISPEP ACCESSORS
           ;; -----------------------------------------------------------------
           #:cispep-ser-num
           #:cispep-pep1 #:cispep-chain-id1 #:cispep-seq-num1 #:cispep-icode1
           #:cispep-pep2 #:cispep-chain-id2 #:cispep-seq-num2 #:cispep-icode2
           #:cispep-mod-num #:cispep-measure
               
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

	       #:read-pdb-file
           #:read-pdb-stream))

(defpackage #:bio-user (:use #:cl #:bio))
