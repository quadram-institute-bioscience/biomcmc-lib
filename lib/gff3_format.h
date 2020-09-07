/* SPDX-License-Identifier: GPL-3.0-or-later 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.

 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file 
 *  \brief GFF3 format 
 *  This library focuses on two units: the chromosome/contig/genome, and the fields (rows) with features. 
 *  
 *  The feature field is the usual, mandatory, 9 columns table. Each feature has a type (CDS, gene, mRNA), the
 *  chromosome/contig/genome it belongs to, its start and end positions, etc.
 *
 *  The chromosome/contig/genome/etc can be described in three places:
 *  1. as a pragma called ##genome-region, before the feature table (optional)
 *  2. as the first column, called "seqid", in the table (mandatory)
 *  3. as a fasta header (sequence name, after the ">") in the ##fasta pragma after the table (optional)
 */

#ifndef _biomcmc_gff3_format_h_
#define _biomcmc_gff3_format_h_

#include "alignment.h"
#include "fortune_cookies.h" 
#include "hashfunctions.h"

typedef struct gff3_file_struct* gff3_t;

typedef struct
{
  char *str;
  uint64_t hash;
  int id;
} gff3_string;

typedef struct 
{
  int start, end;
  uint8_t pos_strand:3, phase:5; 
  gff3_string seqid, source, type, attr_id, attr_parent;
} gff3_fields;

struct gff3_file_struct
{
  gff3_fields *f0, **cds, **gene; // cds and gene are pointers
  int n_f0, n_cds, n_gene;
  char *file_basename;  /*! \brief filename without suffix (file extension) */
  char_vector sequence; /*! \brief from fasta info at end of file; not mandatory */
  char_vector seqname;  /*! \brief from fasta info at end of file; not mandatory */
  hashtable seqname_hash; /*! \brief index from seqname, seq_region, or f0.seqid, in order of preference */
  int *seq_length, *seq_f0_idx; /* \brief same vector size as seqname, has length of each contig and first position in f0[] */
  int ref_counter;
};

gff3_t read_gff3_from_file (char *gff3filename);
void del_gff3_t (gff3_t g3);
gff3_fields find_gff3_field_from_genomic_location (gff3_t g3, const char *seqid, int location); // main function for tatajuba

#endif
