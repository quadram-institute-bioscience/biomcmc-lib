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
 */

#include "gff3_format.h"

gff3_fields gff3_fields_from_char_line (const char *line);
void free_strings_gff3_fields (gff3_fields gff);
gff3_string gff3_string_from_field (const char *start, const char *end);

// static char *feature_type_list[] = {"gene", "cds", "mrna", "region"}; // many more, but I dont use them
// TODO: prokka generates many of these, one per contig: ##sequence-region seqid start end 
// and fasta sequences are one seqid each

compare_gff3_fields_increasing (const void *a, const void *b)
{
  // arbitrary order of distinct genomes 
  if ((*(gff3_fields*)a)->seqid.hash > (*(gff3_fields*)b)->seqid.hash) return 1;
  if ((*(gff3_fields*)a)->seqid.hash < (*(gff3_fields*)b)->seqid.hash) return -1;
  // arbitrary order of feature types (all genes first, then all CDS, etc.)
  if ((*(gff3_fields*)a)->type.hash > (*(gff3_fields*)b)->type.hash) return 1;
  if ((*(gff3_fields*)a)->type.hash < (*(gff3_fields*)b)->type.hash) return -1;
  int result = (*(gff3_fields*)a)->start - (*(gff3_fields*)b)->start;
  if (result) return result; // this is main sorting, by genome location
  return (*(gff3_fields*)a)->end - (*(gff3_fields*)b)->end; // if same type and start, try to sort by end
}

gff3_fields
gff3_fields_from_char_line (const char *line)
{
  gff3_fields gff;
  char *f_start, *f_end;
  int i;

  f_end = line; // will become f_start as son as enters loop 
  for (i = 0; (i < 9) && (f_end != NULL); i++) {
    f_start = f_end; 
    if (f_start[0] == '\t') f_start++;
    f_end = strstr (f_start, '\t');

    switch (i) {
      case 0: // col 1 = SEQID (genome id)
        gff.seqid  = gff3_string_from_field (f_start, f_end); break;
      case 1: // col 2 = source (RefSeq, genbank)
        gff.source = gff3_string_from_field (f_start, f_end); break;
      case 2: // gene, mRNA, CDS
        gff.type   = gff3_string_from_field (f_start, f_end); break;
      case 3: 
        sscanf (f_start, "%d\t", &gff.start); gff.start--; break; // gff3 is one-based but we are zero-based
      case 4: 
        sscanf (f_start, "%d\t", &gff.end);   gff.end--;   break; // gff3 is one-based but we are zero-based
      case 6: // + or - strand (relative to landmark in column 1); can be "." for irrelevant or "?" for unknown 
        if (f_start == '+') gff.pos_strand = 1;
        else if (f_start == '-') gff.pos_strand = 0;
        else gff.pos_strand = 2; break;
      case 7:  // where codon starts, in CDS. It can be 0,1,2 (relative to start if + or to end if -) 
        sscanf (f_start, "%d\t", &gff.phase); break;
      case 8: get_gff3_attributes_from_field (f_start, f_end, &gff.attr_id, &gff.attr_parent); break;
      default: break; // skip 'score' field  (6 of 9)
    };
  }
  if (i < 9) biomcmc_error ("Malformed GFF3 file, found only %d (tab-separated) fields instead of 9 on line \n%s", i, line);
  if (f_end != NULL) biomcmc_error ("Malformed GFF3, more than 9 tab-separating fields found on line\n%s", line);
  return gff;
}

void
free_strings_gff3_fields (gff3_fields gff)
{
  if (gff.seqid.str)   free (gff.seqid.str);
  if (gff.source.str)  free (gff.source.str);
  if (gff.type.str)    free (gff.type.str);
  if (gff.attr_id.str) free (gff.id.str);
  if (gff.attr_parent.str) free (gff.parent.str);
}

gff3_string
gff3_string_from_field (const char *start, const char *end)
{
  gff3_string string;
  type_t len = end - start; // len actually points to final "\t"
  string.str = (char*) biomcmc_malloc (sizeof (char) * len);
  strncpy (string.str, start, len - 1); // do not copy trailing tab
  string.str[len-1] = '\0';
  /* hash is 64 bits, formed by concatenating 2 uint32_t hash values */
  string.hash  = (biomcmc_hashbyte_salted (string.str, len-1, 4) << 32);        // 4 (salt) = CRC algo
  string.hash |= (biomcmc_hashbyte_salted (string.str, len-1, 2) & 0xffffffff); // 2 (salt) = dbj2 algo
  string.id = -1; // defined elsewhere, if relevant
  return string;
}


void /* not actually, I need a struct for this */
read_gff3_from_file (char *gff3filename)
{
  FILE *seqfile;
  char *line = NULL, *line_read = NULL;
  size_t linelength = 0;
  int stage = 0;

  seqfile = biomcmc_fopen (gff3filename, "r");

  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_gff3_string (line)) {
      if ((stage == 0) && strcasestr (line, "##gff-version")) stage = 1;
    }
  }

  fclose (seqfile);
  if (line_read) free (line_read);
