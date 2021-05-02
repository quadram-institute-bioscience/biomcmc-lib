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


int compare_gff3_fields_by_hash (const void *a, const void *b);
int compare_gff3_fields_by_id (const void *a, const void *b);
void free_strings_gff3_fields (gff3_fields gff);

gff3_fields gff3_fields_from_char_line (char *line);
void get_gff3_string_from_field (const char *start, const char *end, gff3_string *string);
uint64_t return_gff3_hashed_string (const char *str, size_t len);
void get_gff3_attributes_from_field (char *start, char *end, gff3_string *attr_id, gff3_string *attr_parent);

gff3_t new_gff3_t (const char *filename);
void add_fields_to_gff3_t (gff3_t g3, gff3_fields gfield);
void gff3_finalise (gff3_t g3, char_vector seq_region);
void merge_seqid_from_fields_and_pragma (gff3_t g3, char_vector *seq_region);
void give_feature_type_id (gff3_t g3);
void generate_feature_type_pointers (gff3_t g3);
void reorder_seqid_charvector_from_pragma (char_vector sid, char_vector *seq_region);
void gff3_generate_seq_vectors (gff3_t g3, hashtable hgs);
int lookup_bruteforce (char_vector haystack, const char *needle);

// static char *feature_type_list[] = {"gene", "cds", "mrna", "region"}; // many more, but I dont use them
// prokka generates many of these, one per contig: ##sequence-region seqid start end 
// and fasta sequences are one seqid each

int
compare_gff3_fields_by_hash (const void *a, const void *b)
{
  // arbitrary order of distinct genomes 
  if (((gff3_fields*)a)->seqid.hash > ((gff3_fields*)b)->seqid.hash) return 1;
  if (((gff3_fields*)a)->seqid.hash < ((gff3_fields*)b)->seqid.hash) return -1;
  // arbitrary order of feature types (all genes then CDS, etc. but based on hashed strings)
  if (((gff3_fields*)a)->type.hash > ((gff3_fields*)b)->type.hash) return 1;
  if (((gff3_fields*)a)->type.hash < ((gff3_fields*)b)->type.hash) return -1;
  int result = ((gff3_fields*)a)->start - ((gff3_fields*)b)->start;
  if (result) return result; // this is main sorting, by genome location
  return ((gff3_fields*)a)->end - ((gff3_fields*)b)->end; // if same type and start, try to sort by end
}

int
compare_gff3_fields_by_id (const void *a, const void *b)
{
  int result = ((gff3_fields*)a)->seqid.id - ((gff3_fields*)b)->seqid.id; 
  if (result) return result; // genomes names in increasing order ("id" comes from char_vector)
  result = ((gff3_fields*)a)->start - ((gff3_fields*)b)->start;
  if (result) return result; // main sorting, by genome location
  result = ((gff3_fields*)a)->end - ((gff3_fields*)b)->end; // if same type and start, try to sort by end
  if (result) return result; 
  result = ((gff3_fields*)a)->type.id -((gff3_fields*)b)->type.id;
  return result; // order of feature types (genes, CDS, etc.) have lowest priority (only if features have same start and end) 
}

static gff3_fields NULL_gff3_fields = {.start = -0xfdfd}; // equivalent to NULL; structs can't be compared (only elements) 

bool
gff3_fields_is_valid (gff3_fields gff)
{
  if (gff.start == NULL_gff3_fields.start) return false;
  return true;
}

gff3_fields
return_null_gff3_field (void)
{
  return NULL_gff3_fields;
}

gff3_fields
gff3_fields_from_char_line (char *line)
{
  gff3_fields gff = NULL_gff3_fields;
  char *f_start, *f_end;
  int i=0, i2;
  // check if proper gff3 fields line, otherwise return NULL before doing anything else
  f_start = line;
  f_end = line + strlen (line); // last column
  for (i = 0; (f_start != NULL) && (f_start < f_end); i++) {
    f_start = strchr (f_start, '\t'); 
    if (f_start) f_start++;
  }
  if (i != 9) return NULL_gff3_fields;

  f_end = line; // will become f_start as son as enters loop 
  for (i = 0; (i < 9) && (f_end != NULL); i++) {
    f_start = f_end; 
    if (f_start[0] == '\t') f_start++;
    f_end = strchr (f_start, '\t');
    if (f_end == NULL) f_end = line + strlen (line); // last column

    switch (i) {
      case 0: // col 1 = SEQID (genome id)
        get_gff3_string_from_field (f_start, f_end, &(gff.seqid)); break;
      case 1: // col 2 = source (RefSeq, genbank)
        get_gff3_string_from_field (f_start, f_end, &(gff.source)); break;
      case 2: // gene, mRNA, CDS
        get_gff3_string_from_field (f_start, f_end, &(gff.type)); break;
      case 3: 
        sscanf (f_start, "%d\t", &gff.start); gff.start--; break; // gff3 is one-based but we are zero-based
      case 4: 
        sscanf (f_start, "%d\t", &gff.end);   gff.end--;   break; // gff3 is one-based but we are zero-based
      case 6: // + or - strand (relative to landmark in column 1); can be "." for irrelevant or "?" for unknown 
        if (*f_start == '+') gff.pos_strand = 1; // positive
        else if (*f_start == '-') gff.pos_strand = 0; // negative
        else if (*f_start == '.') gff.pos_strand = 2; // irrelevant (non-coding)
        else gff.pos_strand = 3;   // unknown ('?')
        break;
      case 7:  // where codon starts, in CDS. It can be 0,1,2 (relative to start if + or to end if -) 
        sscanf (f_start, "%d\t", &i2); gff.phase = i2; break; // bit-fields dont have address
      case 8: get_gff3_attributes_from_field (f_start, f_end, &gff.attr_id, &gff.attr_parent); break;
      default: break; // skip 'score' field  (6 of 9)
    };
  }
  //  if (i < 9) biomcmc_error ("Malformed GFF3 file, found only %d (tab-separated) fields instead of 9 on line \n%s", i, line);
  //  if (f_end != NULL) biomcmc_error ("Malformed GFF3, more than 9 tab-separating fields found on line\n%s", line);
  if (i < 9) {free_strings_gff3_fields (gff); return NULL_gff3_fields; }
  return gff;
}

void
free_strings_gff3_fields (gff3_fields gff)
{
  if (gff.seqid.str)   free (gff.seqid.str);
  if (gff.source.str)  free (gff.source.str);
  if (gff.type.str)    free (gff.type.str);
  if (gff.attr_id.str) free (gff.attr_id.str);
  if (gff.attr_parent.str) free (gff.attr_parent.str);
}

void
get_gff3_string_from_field (const char *start, const char *end, gff3_string *string)
{
  size_t len = end - start + 1; // end actually points to final "\t"; spaces _are_part_of_string_ 
  if (len == 1) { string->str = NULL; string->hash = 0ULL; string->id = -1; return; }

  string->str = (char*) biomcmc_malloc (sizeof (char) * len);
  strncpy (string->str, start, len - 1); // do not copy trailing tab
  string->str[len-1] = '\0';
  /* hash is 64 bits, formed by concatenating 2 uint32_t hash values */
  string->hash = return_gff3_hashed_string (string->str, len-1);
  string->id = -1; // defined when creating char_vector 
  return;
}

uint64_t
return_gff3_hashed_string (const char *str, size_t len)
{
  uint64_t hash = ((uint64_t)(biomcmc_hashbyte_salted (str, len, 4)) << 32); // 4 (salt) = CRC algo
  hash |= (biomcmc_hashbyte_salted (str, len, 2) & 0xffffffffULL);  // 2 (salt) = dbj2 algo
  return hash;
}

void      
get_gff3_attributes_from_field (char *start, char *end, gff3_string *attr_id, gff3_string *attr_parent)
{
  char *s1, *e1;
  s1 = strstr (start, "ID="); // protected attributes start with uppercase 
  if (s1 == NULL) { attr_id->str = NULL; attr_id->hash = 0ULL; attr_id->id = -1; }
  else {
    s1 += 3; // everything after '=', even spaces
    e1 = strchr (s1, ';');
    if (e1 == NULL) e1 = end;
    if (s1 + 1 == e1) { attr_id->str = NULL; attr_id->hash = 0ULL; attr_id->id = -1; }
    else get_gff3_string_from_field (s1, e1, attr_id);
  } 

  s1 = strstr (start, "Parent="); // protected attributes start with uppercase 
  if (s1 == NULL) { attr_parent->str = NULL; attr_parent->hash = 0ULL; attr_parent->id = -1; }
  else { // Parent can be like "Parent=mRNA00001,mRNA00002" (i.e. several parents; here I don't care)
    s1 += 7; // everything after '=', even spaces
    e1 = strchr (s1, ';');
    if (e1 == NULL) e1 = end;
    if (s1 + 1 == e1) { attr_parent->str = NULL; attr_parent->hash = 0ULL; attr_parent->id = -1; }
    else get_gff3_string_from_field (s1, e1, attr_parent);
  } 
  return;
}

gff3_t
read_gff3_from_file (const char *gff3filename)
{
  char *line = NULL, *line_read = NULL, *delim = NULL, *tmpc = NULL;
  size_t linelength = 0;
  int stage = 0, i1, i2, *reg_size = NULL, n_reg_size=0;
  gff3_fields gfield;
  gff3_t g3;
  char_vector seq_region = new_char_vector (1);

  g3 = new_gff3_t (gff3filename);
  tmpc = (char*) biomcmc_malloc (1024 * sizeof (char));

#ifdef HAVE_ZLIB
  gzFile seqfile;
  seqfile = biomcmc_gzopen (gff3filename, "r");
  while (biomcmc_getline_gz (&line_read, &linelength, seqfile) != -1) {
#else
  FILE *seqfile;
  seqfile = biomcmc_fopen (gff3filename, "r");
  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
#endif
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_gff3_line (line)) {
      if ((stage == 0) && (strcasestr (line, "##gff-version") != NULL)) {stage = 1; continue; } // obligatory first line to keep going on

      //else if (stage == 1) { NCBI/refseq does not honour obligation above; some GFF3 do not have gff-version...
      else if (stage < 2) { /* initial pragmas */
        if ((delim = strcasestr (line, "##sequence-region")) != NULL) {
          sscanf (delim, "##sequence-region %s %d %d", tmpc, &i1, &i2); i2--; 
          char_vector_add_string (seq_region, tmpc); // add chromosome/contig name;
          reg_size = (int*) biomcmc_realloc ((int*)reg_size, (n_reg_size + 1) * sizeof (int));
          reg_size[n_reg_size++] = i2; // will be associated to contig name through hashtable 
          continue;
        }
        if (strcasestr (line, "##")) continue; // do nothing atm
        gfield = gff3_fields_from_char_line (line);
        if (gff3_fields_is_valid (gfield)) {
          add_fields_to_gff3_t (g3, gfield);
          stage = 2;
          continue;
        }
      }

      else if (stage == 2) { /* regular fields */
        if (strcasestr (line, "##FASTA")) { stage = 3; continue; }
        if (strcasestr (line, "##")) continue; // do nothing atm
        gfield = gff3_fields_from_char_line (line);
        if (gff3_fields_is_valid (gfield)) { add_fields_to_gff3_t (g3, gfield); continue; }
      }

      else if (stage == 3) {  /* FASTA file */
        if (nonempty_fasta_line (line)) { /* each line can be either the sequence or its name, on a strict order */
          if ((delim = strchr (line, '>'))) char_vector_add_string (g3->seqname, ++delim); /* sequence name (description, in FASTA jargon) */
          else {/* the sequence itself, which may span several lines */
            line = remove_space_from_string (line);
            line = uppercase_string (line);
            char_vector_append_string_big_at_position (g3->sequence, line, g3->seqname->next_avail-1); // counter from name NOT sequence
          }
        } // if non_empty
      } // if stage 3 
    } // if gff3 line not empty
  } // while readline ()
#ifdef HAVE_ZLIB
  gzclose (seqfile);
#else
  fclose (seqfile);
#endif
  char_vector_finalise_big (g3->sequence);

  /* hashtable with genome sizes for all seq_regions (contig/genome/chromosome), but only added to gff3_t after finalisation */
  hashtable hgs = new_hashtable (n_reg_size);
  for (i1 = 0; i1 < n_reg_size; i1++)  insert_hashtable (hgs, seq_region->string[i1], reg_size[i1]); // seq name + genome size
  /* sort, remove unused contig names; and then generate vectors (with genome size etc.) */
  gff3_finalise (g3, seq_region);
  gff3_generate_seq_vectors (g3, hgs);

  del_char_vector (seq_region);
  del_hashtable (hgs);
  if (line_read) free (line_read);
  if (tmpc) free (tmpc);
  if (reg_size) free (reg_size);
  return g3;
}

gff3_t
new_gff3_t (const char *filename)
{
  gff3_t g3 = (gff3_t) biomcmc_malloc (sizeof (struct gff3_file_struct));
  g3->sequence = new_char_vector_big (1);
  g3->seqname = new_char_vector (1);
  g3->f0 = NULL;
  g3->cds = g3->gene = NULL;
  g3->n_f0 = g3->n_cds = g3->n_gene = 0;
  g3->seqname_hash = NULL; // created when finalising 
  g3->seq_f0_idx = g3->seq_length = NULL;
  g3->file_basename = NULL;

  char *last = biomcmc_strrstr (filename, ".gff");
  if (last) {
    size_t len = last - filename;
    g3->file_basename = (char*) biomcmc_malloc ((len+1) * sizeof (char));
    strncpy (g3->file_basename, filename, len);
    g3->file_basename[len] = '\0';
  }

  g3->ref_counter = 1;
  return g3;
}

void
del_gff3_t (gff3_t g3)
{
  if (!g3) return;
  if (--g3->ref_counter) return;
  del_char_vector (g3->sequence);
  del_char_vector (g3->seqname);
  del_hashtable   (g3->seqname_hash);
  if (g3->f0) { // only place where strings are dealloced 
    for (int i = 0; i < g3->n_f0; i++) free_strings_gff3_fields (g3->f0[i]);
    free   (g3->f0);
  }
  if (g3->cds) free  (g3->cds);
  if (g3->gene) free (g3->gene);
  if (g3->seq_length) free (g3->seq_length);
  if (g3->seq_f0_idx) free (g3->seq_f0_idx);
  if (g3->file_basename) free (g3->file_basename);
  free (g3);
}

void
add_fields_to_gff3_t (gff3_t g3, gff3_fields gfield)
{
  g3->f0 = (gff3_fields*) biomcmc_realloc ((gff3_fields*) g3->f0, (g3->n_f0 + 1) * sizeof (gff3_fields));
  g3->f0[g3->n_f0++] = gfield;
}

void
gff3_finalise (gff3_t g3, char_vector seq_region)
{
  /* 1. sort fields, map seqids to char_vector, and point to specific features (cds, gene) */
  give_feature_type_id (g3); // transform cds, genes, etc to 0, 1, etc so that can be properly sorted below
  merge_seqid_from_fields_and_pragma (g3, &seq_region); // updates seq_region to match fields->seqid
  generate_feature_type_pointers (g3); // vectors of pointers (CDS, gene)

  /* 2. now seq_region has updated, authoritative sequence names; */
  char_vector fasta_name = g3->seqname; // just a pointer
  g3->seqname = seq_region; // must increase pointer counter since 
  g3->seqname->ref_counter++; // seq_region will be deleted later
  
  /* 2. if fasta is missing, just copy seqids to seqname */
  if (!g3->sequence->next_avail) {
    del_char_vector (fasta_name); // delete original seqname of length zero
    del_char_vector (g3->sequence); g3->sequence = NULL;
    return;
  }

  int hid, i, *idx;
  /* 3. prepare to reorder fasta, with name order in gff3 */
  idx = (int*) biomcmc_malloc (fasta_name->nstrings * sizeof (int));
  for (i = 0; i < fasta_name->nstrings; i++) {
    hid = lookup_hashtable (g3->seqname_hash, fasta_name->string[i]);
    idx[i] = hid; // negative values if fasta sequence not present in GFF3 
  }

  /* 4. Assume fasta file can have info about only a few of the seqs */
  char_vector seq = new_char_vector (g3->seqname->nstrings);
  for (i = 0; i < fasta_name->nstrings; i++) if (idx[i] >= 0) char_vector_add_string_at_position (seq, g3->sequence->string[i], idx[i]); 
  del_char_vector (g3->sequence);
  g3->sequence = seq;

  del_char_vector (fasta_name); // original fasta names not needed anymore
  if (idx) free (idx);

  return;
}

void
gff3_finalise_OLD (gff3_t g3, char_vector seq_region)
{
  /* 1. sort fields, map seqids to char_vector, and point to specific features (cds, gene) */
  give_feature_type_id (g3); // transform cds, genes, etc to 0, 1, etc so that can be properly sorted below
  merge_seqid_from_fields_and_pragma (g3, &seq_region); // updates seq_region to match fields->seqid
  generate_feature_type_pointers (g3); // vectors of pointers (CDS, gene)

  /* 2. if fasta is missing, just copy seqids to seqname */
  //if ((!g3->sequence->next_avail) || (g3->seqname->next_avail < seq_region->next_avail)) { 
  if (!g3->sequence->next_avail) { 
    //printf ("DBG::no fasta %6d , %6d %6d\n", g3->seqname->nstrings, seq_region->nstrings, seq_region->next_avail);
    del_char_vector (g3->sequence); g3->sequence = NULL;
    del_char_vector (g3->seqname);
    g3->seqname = seq_region;
    g3->seqname->ref_counter++; // seq_region will be deleted later
    //if (g3->seqname && (g3->seqname->next_avail < seq_region->next_avail)) biomcmc_warning ("incomplete fasta pragma in GFF3; ignoring DNA sequences from file\n");
    return;
  }

  /* 3. map seqnames from fasta pragma to seqid from fields; assume fasta can have spurious (extra) seqs */
  int i, n_extra, hid, *order;
  order = (int*) biomcmc_malloc (g3->seqname->nstrings * sizeof (int));
  n_extra = seq_region->nstrings; // OK for us for fasta to have more sequences than needed
  for (i = 0; i < g3->seqname->nstrings; i++) {
    hid = lookup_hashtable (g3->seqname_hash, g3->seqname->string[i]);
    if ((hid < 0) && (n_extra == g3->seqname->nstrings)) {
      biomcmc_warning ("fasta pragma in GFF3 doesn't correspond to field seqids; ignoring DNA sequences from file\n");
      del_char_vector (g3->sequence); g3->sequence = NULL;
      del_char_vector (g3->seqname);
      g3->seqname = seq_region;
      g3->seqname->ref_counter++; // seq_region will be deleted later
      if (order) free (order);
      return;
    }
    else if (hid < 0) order[n_extra++] = i; // last elements, hopefully just extra fasta sequences
    else order[hid] = i;
  }
  /* 4. use order from fields seqids (sequence-region sorted) on fasta */  
  if (seq_region->nstrings < g3->seqname->nstrings) {
    char_vector_reorder_strings_from_external_order (g3->seqname,  order);
    char_vector_reorder_strings_from_external_order (g3->sequence, order);
    char_vector_reduce_to_trimmed_size (g3->seqname,  seq_region->nstrings);	
    char_vector_reduce_to_trimmed_size (g3->sequence, seq_region->nstrings);	
  }

  if (order) free (order);
  return;
}

void
merge_seqid_from_fields_and_pragma (gff3_t g3, char_vector *seq_region)
{ // seq_region will indicate order of seqids; this function updates it to match feature table
  int i, hid;
  char_vector sid = new_char_vector (1);
  qsort (g3->f0, g3->n_f0, sizeof (gff3_fields), compare_gff3_fields_by_hash);

  /* 1. create char_vector with seqid names from feature table */
  char_vector_add_string (sid, g3->f0[0].seqid.str);
  for (i = 1; i < g3->n_f0; i++) 
    if (strcmp(g3->f0[i].seqid.str, sid->string[sid->next_avail-1])) char_vector_add_string (sid, g3->f0[i].seqid.str);

  /* 2. sid has all seqids in feature table, but order is given by seq_region. UPDATED seqid will be seq_region */
  reorder_seqid_charvector_from_pragma (sid, seq_region);
  del_char_vector (sid); // we don't need sid anymore; seq_region has all names, in best possible order 

  /* 3. seq_region has all seqids and in right order; now fields must be updated */
  g3->seqname_hash = new_hashtable ((*seq_region)->nstrings);
  for (i = 0; i < (*seq_region)->nstrings; i++) insert_hashtable (g3->seqname_hash, (*seq_region)->string[i], i);
  for (i = 0; i < g3->n_f0; i++) {
    hid = lookup_hashtable (g3->seqname_hash, g3->f0[i].seqid.str);
    g3->f0[i].seqid.id = hid; // hid should be valid (i.e. not negative) but downstream functions might want to check
  }

  /* 4. now that feature seqid fields have id, sort whole table by genome/contig id and then by location (if you need just one type (CDs, genes) pls 
   * use pointers set by feature_type_pointer() since here they have lowest priority) */
  qsort (g3->f0, g3->n_f0, sizeof (gff3_fields), compare_gff3_fields_by_id);
  return;
}

void
give_feature_type_id (gff3_t g3)
{ // http://www.sequenceontology.org/browser/obob.cgi?rm=term_list&release=current_svn&obo_query=$$
  for (int i = 0; i < g3->n_f0; i++) {
    if ((!strcasecmp (g3->f0[i].type.str, "cds")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000316"))) {
      g3->f0[i].type.id = GFF3_TYPE_cds;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "gene")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000704"))) {
      g3->f0[i].type.id = GFF3_TYPE_gene;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "mRNA")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000234"))) {
      g3->f0[i].type.id = GFF3_TYPE_mRNA;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "exon")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000147"))) {
      g3->f0[i].type.id = GFF3_TYPE_exon;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "polyA_sequence")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000610"))) {
      g3->f0[i].type.id = GFF3_TYPE_polyA_sequence;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "polyA_site")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000553"))) {
      g3->f0[i].type.id = GFF3_TYPE_polyA_site;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "intron")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000188"))) {
      g3->f0[i].type.id = GFF3_TYPE_intron;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "five_prime_UTR")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000204"))) {
      g3->f0[i].type.id = GFF3_TYPE_five_prime_UTR;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "three_prime_UTR")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000205"))) {
      g3->f0[i].type.id = GFF3_TYPE_three_prime_UTR;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "tRNA")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000253"))) {
      g3->f0[i].type.id = GFF3_TYPE_tRNA;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "rRNA")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000252"))) {
      g3->f0[i].type.id = GFF3_TYPE_rRNA;
    }
    else if ((!strcasecmp (g3->f0[i].type.str, "tmRNA")) || (!strcasecmp (g3->f0[i].type.str, "SO:0000584"))) {
      g3->f0[i].type.id = GFF3_TYPE_tmRNA;
    }
    else if (!strcasecmp (g3->f0[i].type.str, "region")) { // simply means 'whole genome/chromosome/contig', not useful generally 
      g3->f0[i].type.id = GFF3_TYPE_region;
    }
    else g3->f0[i].type.id = 0xffff;
  }
  return;
}

void
generate_feature_type_pointers (gff3_t g3)
{
  int i;
  g3->cds  = (gff3_fields**) biomcmc_malloc (g3->n_f0 * sizeof (gff3_fields*));
  g3->gene = (gff3_fields**) biomcmc_malloc (g3->n_f0 * sizeof (gff3_fields*));
  g3->n_gene = g3->n_cds = 0;

  for (i = 0; i < g3->n_f0; i++) {
    if (g3->f0[i].type.id == 0) { // for mapping see give_feature_type_id() above 
      g3->f0[i].attr_id.id = g3->n_cds; // experimental (not used); several CDS (rows) can have same ID in gff3
      g3->cds[g3->n_cds++] = &(g3->f0[i]);
    }
    else if (g3->f0[i].type.id == 1) {
      g3->f0[i].attr_id.id = g3->n_gene; // experimental (not used yet)
      g3->cds[g3->n_gene++] = &(g3->f0[i]);
    }
  }

  if (g3->n_cds > 0)  g3->cds = (gff3_fields**) biomcmc_realloc ((gff3_fields**) g3->cds,  g3->n_cds  * sizeof (gff3_fields*));
  else if (g3->cds)  { free (g3->cds);  g3->cds = NULL; }

  if (g3->n_gene > 0) g3->gene = (gff3_fields**) biomcmc_realloc ((gff3_fields**) g3->gene, g3->n_gene * sizeof (gff3_fields*));
  else if (g3->gene) { free (g3->gene); g3->gene = NULL; }
  return;
}

void
reorder_seqid_charvector_from_pragma (char_vector sid, char_vector *seq_region)
{
  if ((*seq_region)->next_avail < 2) { // empty (no pragma directives) or only one sequence name
    del_char_vector (*seq_region);
    *seq_region = sid;
    sid->ref_counter++;
    return;
  }
  int i, hid, *idx, n_seq = (*seq_region)->nstrings;
  char_vector seq = *seq_region; // convenience only, hard to look at asterisks
  empfreq ef;
  hashtable ht;

  /* 1. hash table with indices of seq_region */
  ht = new_hashtable (seq->nstrings);
  for (i = 0; i < seq->nstrings; i++) insert_hashtable (ht, seq->string[i], i);

  /* 2. which seq_region name maps to each seqid */
  idx = (int*) biomcmc_malloc (sid->nstrings * sizeof (int));
  for (i = 0; i < sid->nstrings; i++) {
    hid = lookup_hashtable (ht, sid->string[i]);
    if (hid < 0) hid = n_seq++; // seqid not present in pragma, will be sorted last
    idx[i] = hid; // hid value doesnt matter, only order
  }

  /* 3. sort idx[] while keeping index i (of idx[i]) */
  ef = new_empfreq_sort_increasing (idx, sid->nstrings, 2); // 2->integer (0->char and 1->size_t)
  for (i = 0; i < sid->nstrings; i++) idx[i] = ef->i[i].idx;

  /* 4. reorder seqid preesrving preferences from pragma (seq_region) */
  char_vector_reorder_strings_from_external_order (sid, idx);

  /* 5. now that seqid is up-to-date, we don't need seq_region (but we do need its pointer!) */
  del_char_vector (*seq_region);
  *seq_region = sid;
  sid->ref_counter++;

  /* 6. free memomry */
  if (idx) free (idx);
  del_empfreq (ef);
  del_hashtable (ht);
  return;
}

void
gff3_generate_seq_vectors (gff3_t g3, hashtable hgs)
{
  int i;

  g3->seq_length = (int*) biomcmc_malloc (g3->seqname->nstrings * sizeof (int));
  g3->seq_f0_idx = (int*) biomcmc_malloc (g3->seqname->nstrings * sizeof (int));

  for (i = 0; i < g3->seqname->nstrings; i++) g3->seq_length[i] = lookup_hashtable (hgs, g3->seqname->string[i]);
  assert (g3->f0[0].seqid.id == 0);
  g3->seq_f0_idx[0] = 0; // idx of first element (on f0 vector) belonging to seqid 0
  for (i = 1; i < g3->n_f0; i++) if (g3->f0[i].seqid.id != g3->f0[i-1].seqid.id) g3->seq_f0_idx[ g3->f0[i].seqid.id ] = i;
  return;
}

void
add_fasta_to_gff3 (gff3_t g3, char_vector name, char_vector seq)
{
  int hid, i, *idx;

  /* 1. prepare to reorder fasta, with name order in gff3 (we don't care about name[], GFF3 has authoritative one) */
  idx = (int*) biomcmc_malloc (name->nstrings * sizeof (int));
  for (i = 0; i < name->nstrings; i++) {
    hid = lookup_hashtable (g3->seqname_hash, name->string[i]);
    if (hid < 0) hid = lookup_bruteforce (g3->seqname, name->string[i]);
    idx[i] = hid; // negative values if fasta sequence not present in GFF3 
  }

  /* 2. even if GFF3 has no DNA seqs, we cannot just copy char_vector b/c we need to preserve order and size (from seqname) */
  if (!g3->sequence) g3->sequence = new_char_vector (g3->seqname->nstrings);
  /* 3. Assume that the fasta file can have only info about some of the seqs */
  for (i = 0; i < name->nstrings; i++) if (idx[i] >= 0) char_vector_add_string_at_position (g3->sequence, seq->string[i], idx[i]); 
  if (idx) free (idx);
  return;
} 

int
lookup_bruteforce (char_vector haystack, const char *needle)
{
  int i, j, minlen;
  size_t n_len = strlen (needle);
  for (i = 0; i < haystack->nstrings; i++) {
    minlen = BIOMCMC_MIN(n_len, haystack->nchars[i]);
    for (j = 0; (j < minlen) && (needle[j] == haystack->string[i][j]); j++);
    if (j == minlen) return i;
  }
  return -1;
}

char *
save_fasta_from_gff3 (gff3_t g3, char *fname, bool overwrite)
{
  if (!g3->sequence) return NULL; // no fasta info 
  char *filename;
  FILE *f_exist;
  if (!fname) fname = g3->file_basename;
  size_t len = strlen (fname);
#ifdef HAVE_ZLIB
  filename = (char*) biomcmc_malloc ((len + 7) * sizeof (char));
  strcpy (filename, fname); strcat (filename, ".fa.gz");
#else
  filename = (char*) biomcmc_malloc ((len + 4) * sizeof (char));
  strcpy (filename, fname); strcat (filename, ".fa");
#endif

  if ((!overwrite) && ((f_exist = fopen(filename, "r")) != NULL)) { // file exists; don't overwrite
    fclose (f_exist);
    return filename;
  }

  save_gzfasta_from_char_vector (filename, g3->seqname, g3->sequence);
  return filename;
}

gff3_fields *
find_gff3_fields_within_position (gff3_t g3, const char *ref_genome, int location, int *n)
{
  gff3_fields *fid;
  int i, start_genome, end_genome, first, mid, last; // start and end are positions in f0[] vector of all gff3_fields_t

  *n = -1; // return value if genome does not exist (i.e. error)
  i = lookup_hashtable (g3->seqname_hash, ref_genome);
  if (i < 0) { biomcmc_warning ("no reference \"%s\" found in GFF3 file %s\n", ref_genome, g3->file_basename); return NULL; }
  start_genome = g3->seq_f0_idx[i]; //seq_f0_idx is a list of locations of first field from genome  
  if (i < g3->seqname->nstrings - 1) end_genome = g3->seq_f0_idx[i + 1] - 1;
  else                               end_genome = g3->n_f0 - 1; // n_f0 is an extra value with total number of fields

  *n = 0; // return value if no feature found (position before first or after last feature) 
  if ((g3->f0[start_genome].start > location) || (g3->f0[end_genome].end < location)) return NULL;

  fid = (gff3_fields*) biomcmc_malloc ((end_genome - start_genome + 1) * sizeof (gff3_fields));// end==start  

  first = start_genome; last = end_genome; // we might need start and end of genome again later
  while (first <= last) {
    mid = first + (last - first)/2;
    if (location < g3->f0[mid].start) last = mid - 1; // update binary search and move
    else { // location >= feature start ; now we collect all those also <= feature end
      if (location <= g3->f0[mid].end) { // start <= location <= end; we may have more than one 
        int j;
        for (j = mid; (j >= first) && (location >= g3->f0[j].start) && (location <= g3->f0[j].end); j--)
          if (g3->f0[j].type.id != GFF3_TYPE_region) fid[(*n)++] = g3->f0[j]; // add unless is "region" (i.e. whole genome)
        for (j = mid+1; (j <= last) && (location >= g3->f0[j].start) && (location <= g3->f0[j].end); j++)
          if (g3->f0[j].type.id != GFF3_TYPE_region) fid[(*n)++] = g3->f0[j];
        first = last + 1;  // getoutahere
      }
      else first = mid + 1; // update binary search for next iter (in case location is after end of f0[mid]
    } // else (i.e. if location >= feature.start)
  } // while in binary search

  if (*n) fid = (gff3_fields*) biomcmc_realloc ((gff3_fields*) fid, (*n) * sizeof (gff3_fields));
  else if (fid) { free (fid); fid = NULL; }
  return fid;
}

gff3_fields *
find_gff3_fields_within_position_all_genomes (gff3_t g3, int location, int *n)
{
  gff3_fields *fid;
  int i, start_genome, end_genome, first, mid, last; // start and end are positions in f0[] vector of all gff3_fields_t

  *n = 0; // return value if no feature found (position before first or after last feature) 
  fid = (gff3_fields*) biomcmc_malloc (g3->n_f0 * sizeof (gff3_fields));// overestimate  
  for (i = 0; i < g3->seqname->nstrings; i++) { // overall genomes/contigs/chromosomes
    start_genome = g3->seq_f0_idx[i]; 
    if (i < g3->seqname->nstrings - 1) end_genome = g3->seq_f0_idx[i + 1] - 1;
    else                               end_genome = g3->n_f0 - 1; 
    if ((g3->f0[start_genome].start > location) || (g3->f0[end_genome].end < location)) continue;

    first = start_genome; last = end_genome; 
    while (first <= last) {
      mid = first + (last - first)/2;
      if (location < g3->f0[mid].start) last = mid - 1; // update binary search and move
      else { // location >= feature start ; now we collect all those also <= feature end
        if (location <= g3->f0[mid].end) { // start <= location <= end; we may have more than one 
          int j;
          for (j = mid; (j >= first) && (location >= g3->f0[j].start) && (location <= g3->f0[j].end); j--)
            if (g3->f0[j].type.id != GFF3_TYPE_region) fid[(*n)++] = g3->f0[j]; // add unless is "region" (i.e. whole genome)
          for (j = mid+1; (j <= last) && (location >= g3->f0[j].start) && (location <= g3->f0[j].end); j++)
            if (g3->f0[j].type.id != GFF3_TYPE_region) fid[(*n)++] = g3->f0[j];
          first = last + 1;  // getoutahere
        }
        else first = mid + 1; // update binary search for next iter (in case location is after end of f0[mid]
      } // else (i.e. if location >= feature.start)
    } // while in binary search
  } // for all genomes (g3->seqname)

  if (*n) fid = (gff3_fields*) biomcmc_realloc ((gff3_fields*) fid, (*n) * sizeof (gff3_fields));
  else if (fid) { free (fid); fid = NULL; }
  return fid;
}

