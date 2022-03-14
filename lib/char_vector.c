/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file char_vector.c 
 *  \brief vector of strings (species names, leaf names, etc.) 
 */

#include "char_vector.h"


typedef struct { char *s; int idx; size_t nchars; } charvec_str;
int compare_charvecstr_decreasing (const void *a, const void *b);
int compare_charvecstr_lexicographic (const void *a, const void *b);

int
compare_charvecstr_decreasing (const void *a, const void *b)
{
  int result = (int)(((charvec_str *)b)->nchars - ((charvec_str *)a)->nchars);
  if (result) return result;
  int i;
  /* break ties by lexico order; notice order (b-a) above <bigger string first> and (a-b) below <smaller character first> */
  for (i=0; (i < (int)(((charvec_str *)b)->nchars)) && (!result); i++) result = (int)(((charvec_str *)a)->s[i] - ((charvec_str *)b)->s[i]);
  return result;
}

int
compare_charvecstr_lexicographic (const void *a, const void *b)
{
  int i, result = 0, n = (int)(((charvec_str *)b)->nchars);
  if (n < (int)(((charvec_str *)a)->nchars)) n = (int)(((charvec_str *)a)->nchars);
  for (i=0; (i < n) && (!result); i++) result = (int)(((charvec_str *)a)->s[i] - ((charvec_str *)b)->s[i]);
  if (!result) return (int)(((charvec_str *)a)->nchars - ((charvec_str *)b)->nchars); // shorter first, if otherwise identical
  return result;
}

char_vector
new_char_vector (int nstrings)
{
  char_vector vec;
  int i;

  if (nstrings < 1) biomcmc_error ("Vector of strings should have at least one string");

  vec = (char_vector) biomcmc_malloc (sizeof (struct char_vector_struct));
  vec->nstrings = nstrings;
  vec->ref_counter = 1;
  vec->next_avail = 0;

  vec->string = (char**)  biomcmc_malloc (nstrings * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_malloc (nstrings * sizeof (size_t));
  vec->alloc = NULL; // allocated only if needed (used with char_vector_append)

  for (i=0; i < nstrings; i++) {
    vec->string[i] = (char*) biomcmc_malloc (sizeof (char));
    vec->string[i][0] = '\0';
    vec->nchars[i] = 0;
  }
  return vec;
}

char_vector
new_char_vector_big (int nstrings)
{
  char_vector vec = new_char_vector (nstrings);
  vec->alloc = (size_t*) biomcmc_malloc (nstrings * sizeof (size_t));
  for (int i = 0; i < nstrings; i++) vec->alloc[i] = 0;
  return vec;
}

char_vector
new_char_vector_from_valid_strings_char_vector (char_vector vec, int *valid, int n_valid)
{
  int i;
  char_vector newvec = new_char_vector (n_valid);
  for (i=0; i < n_valid; i++) char_vector_add_string (newvec, vec->string[ valid[i] ]);
  return newvec;
} 

char_vector
new_char_vector_fixed_length (int nstrings, int nchars)
{
  char_vector vec;
  int i;

  if (nstrings < 1) biomcmc_error ("Vector of strings should have at least one string");

  vec = (char_vector) biomcmc_malloc (sizeof (struct char_vector_struct));
  vec->nstrings = nstrings;
  vec->ref_counter = 1;
  vec->next_avail = 0;

  vec->string = (char**)  biomcmc_malloc (nstrings * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_malloc (nstrings * sizeof (size_t));
  vec->alloc = NULL; // we don't worry about realloc(), careful not mix this function with 'append_big'

  for (i=0; i < nstrings; i++) {
    vec->string[i] = (char*) biomcmc_malloc ((nchars+1) * sizeof (char)); //ending '\0'
    vec->string[i][0] = '\0';
    vec->nchars[i] = (size_t) nchars;
  }
  return vec;
}

void
del_char_vector (char_vector vec)
{
  int i;
  if (!vec) return;
  if (--vec->ref_counter) return;
  if (vec->string) {
    for (i=vec->nstrings-1; i >=0; i--) if (vec->string[i]) free (vec->string[i]);
    free (vec->string);
  }
  if (vec->nchars) free (vec->nchars);
  if (vec->alloc)  free (vec->alloc);
  free (vec);
}

void
char_vector_link_string_at_position (char_vector vec, const char *string, int position)
{
  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);
  if (vec->nchars[position]) free (vec->string[position]);

  vec->nchars[position] = strlen (string); /* Actually alloc'ed memory may be larger than this (next line fix it) */
  string = (char*) biomcmc_realloc ((char*) string, (vec->nchars[position]+1) * sizeof (char));
  vec->string[position] = (char*) string;
  // next_avail may be before position; should we assume char_vector is always increasing?
  vec->next_avail = position+1;
}

void
char_vector_add_string_at_position (char_vector vec, const char *string, int position)
{
  size_t l;
  char *s = (char*) (string + strspn (string, " \t")); /* skip leading spaces */
  l = strlen (s);

  if (!l) return; /* do nothing with empty strings - like last line of some nexus alignments */
  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);

  if (l > vec->nchars[position]) {
    vec->nchars[position] = l;
    vec->string[position] = (char*) biomcmc_realloc ((char*)vec->string[position], (l+1) * sizeof (char));
  }
  strncpy (vec->string[position], s, l+1); /* l+1 will insert the ending null */
  vec->next_avail = position+1;
}

void
char_vector_add_string (char_vector vec, const char *string)
{
  char_vector_add_string_at_position (vec, string, vec->next_avail);
}

void
char_vector_append_string_at_position (char_vector vec, const char *string, int position)
{
  size_t l, this_l;
  char *s = (char*) (string + strspn (string, " \t")); /* skip leading spaces */
  l = strlen (s);

  if (!l) return; /* do nothing with empty strings - like last line of some nexus alignments */
  if (position < 0) position = 0; /* vec->next_avail works for add_string() but not here... */

  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);
  this_l = strlen (vec->string[position]); /* we assume there is an ending null at vec->string  */

  if ((l+this_l) > vec->nchars[position]) {
    vec->nchars[position] = l + this_l;
    vec->string[position] = (char*) biomcmc_realloc ((char*)vec->string[position], (l + this_l + 1) * sizeof (char));
  }
  strncpy (vec->string[position] + this_l, s, l + 1); /* l+1 will insert the ending null */
}

void
char_vector_append_string (char_vector vec, const char *string)
{
  char_vector_append_string_at_position (vec, string, vec->next_avail-1);
}

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x)) 
void
char_vector_append_string_big_at_position (char_vector vec, const char *string, int position)
{
  size_t l, this_l, new_l;
  if (!vec->alloc) char_vector_append_string_at_position (vec, string, position); // called the wrong function

  char *s = (char*) (string + strspn (string, " \t")); /* skip leading spaces */
  l = strlen (s);
  if (!l) return; /* do nothing with empty strings - like last line of some nexus alignments */
  if (position < 0) position = 0; /* vec->next_avail works for add_string() but not here... */

  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);
  this_l = vec->nchars[position];

//  printf ("DEBUG:: %4lu %4lu -> %4lu %4lu  %4lu  < - before\n", l, this_l, l + this_l + 1, vec->alloc[position], vec->nchars[position]);
  if ((l + this_l + 1) >= vec->alloc[position]) {
    new_l = l + this_l + 1;
    kroundup32(new_l); // defined above, finds next highest power of 2
    vec->string[position] = (char*) biomcmc_realloc ((char*)vec->string[position], new_l * sizeof (char));
    vec->alloc[position] = new_l;
//    printf ("DEBUG:: %4lu %4lu %4lu  -   %4lu\n", l, this_l, vec->alloc[position], vec->nchars[position]);
  }
  vec->nchars[position] += l; // have current size, not allocated size (since here they are distinct)
  strncpy (vec->string[position] + this_l, s, l + 1); /* l+1 will insert the ending null */
}
#undef kroundup32

void
char_vector_append_string_big (char_vector vec, const char *string)
{
  char_vector_append_string_big_at_position (vec, string, vec->next_avail-1);
}

void
char_vector_finalise_big (char_vector vec)
{
  for (int i = 0; i < vec->nstrings; i++) if (vec->alloc[i] > vec->nchars[i]) {
    vec->string[i] = (char*) biomcmc_realloc ((char*)vec->string[i], vec->nchars[i] * sizeof (char));
  }
  if (vec->alloc) free (vec->alloc);
  vec->alloc = NULL;
}

void
char_vector_expand_nstrings (char_vector vec, int new_size)
{
  int i;

  if (new_size < vec->nstrings) biomcmc_error ("I refuse to reduce char_vector size. This is a feature, not a bug.");

  vec->string = (char**)  biomcmc_realloc ((char**)  vec->string, new_size * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_realloc ((size_t*) vec->nchars, new_size * sizeof (size_t));

  for (i=vec->nstrings; i < new_size; i++) {
    vec->string[i] = (char*) biomcmc_malloc (sizeof (char));
    vec->string[i][0] = '\0';
    vec->nchars[i] = 0;
  }
  if (vec->alloc) {
    vec->alloc = (size_t*) biomcmc_realloc ((size_t*) vec->alloc, new_size * sizeof (size_t));
    for (i=vec->nstrings; i < new_size; i++) vec->alloc[i] = 0;
  }
  
  vec->next_avail = vec->nstrings = new_size;
}

void
char_vector_reorder_strings_from_external_order (char_vector vec, int *order)
{
  char  **ptr_c; // tmp pointer to keep track of labels (otherwise once t[2]=t[1] we loose old t[2]...)
  size_t *ptr_i; // tmp pointer to keep track of allocated memory inside char_vector_struct::
  int i;

  if (!vec->next_avail) return; /* if char_vector doesn't have any elements yet */

  ptr_c = (char**)  biomcmc_malloc (vec->nstrings * sizeof (char*));
  ptr_i = (size_t*) biomcmc_malloc (vec->nstrings * sizeof (size_t));

  /* Vector version of tmp = a; a = b; b = tmp; [geek curiosity: how would a "XOR swap" work in this case?] */
  for (i=0; i < vec->nstrings; i++) { 
    ptr_c[i] = vec->string[i]; 
    ptr_i[i] = vec->nchars[i];
  }
  for (i=0; i < vec->nstrings; i++) { 
    vec->string[i] = ptr_c[order[i]]; 
    vec->nchars[i] = ptr_i[order[i]];
  }
  if (ptr_c) free (ptr_c);
  if (ptr_i) free (ptr_i);
}

int
char_vector_remove_empty_strings (char_vector vec)
{
  int i, n_invalid = 0, n_valid=0, *valid;
  size_t length;

  if (!vec->next_avail) return 0; /* if char_vector doesn't have any elements yet */

  /* vector with indexes of valid (non-empty) strings */
  valid = (int*) biomcmc_malloc (vec->nstrings * sizeof (int));
  for (i=0; i < vec->nstrings; i++) { /* find non-empty strings and realloc them */
    length = strlen (vec->string[i]);
    if (!length) { /* empty string - this shouldn't happen and will be reported through the return value */
      n_invalid++;
      if (vec->string[i]) free (vec->string[i]); /* the condition is always true... */
    }
    else { /* valid string */
      valid[n_valid++] = i;
      if (length < vec->nchars[i]) { /* we can save some space, since nchars gives us the allocated memory */
        vec->nchars[i] = length;
        vec->string[i] = (char*) biomcmc_realloc ((char*) vec->string[i], (length+1) * sizeof (char));
      }
    }
  }

  if (!n_invalid) { if (valid) free (valid); return n_invalid; } /* all strings are valid: return zero */

  /* now we can reduce the vector of strings */
  char_vector_reduce_to_valid_strings (vec, valid, n_valid);
  
  if (valid) free (valid);
  return n_invalid;
}

int
char_vector_remove_duplicate_strings (char_vector vec)
{
  int i, j, k, *valid, n_valid = 0;
  bool equal;

  valid = (int*) biomcmc_malloc (vec->nstrings * sizeof (int));

  if (!vec->next_avail) return 0; /* if char_vector doesn't have any elements yet */

  for (i=0; i < vec->nstrings - 1; i++) if (vec->string[i]) {
    valid[n_valid++] = i; /* valid since it is leftmost, and not NULLified before */
    for (j = i + 1; j < vec->nstrings; j++) if ((vec->string[j]) && (vec->nchars[i] == vec->nchars[j])) {
      for (equal = true, k = 0; (equal == true) && (k < (int) vec->nchars[i]); k++) /* scan both strings */
        if (vec->string[i][k] != vec->string[j][k]) equal = false; /* if they are different, do nothing */
      if (equal == true) { /* if they are equal, we free the memory and NULLify superfluous element */
        if (vec->string[j]) free (vec->string[j]); /* the condition is always true... */
        vec->string[j] = NULL;
      }
    }
  }
  if (vec->string[i]) valid[n_valid++] = i; /* check if last string is unique or not (my most common mistake) */

  if (n_valid == vec->nstrings) { if (valid) free (valid); return 0; }

  i = vec->nstrings;
  /* TODO: for alignments (seq AND seqname) we want both, the original and the reduced */
  char_vector_reduce_to_valid_strings (vec, valid, n_valid);

  if (valid) free (valid);
  return (i - n_valid); /* number of chopped elements (= (old vec->nstrings) - (new vec->nstrings) ) */
}

void
char_vector_reduce_to_valid_strings (char_vector vec, int *valid, int n_valid)
{ // assumes valid[] is in increasing order s.t. valid[i] >= i always
  int i;

  for (i=0; i < n_valid; i++) {
    vec->string[i] = vec->string[ valid[i] ];
    vec->nchars[i] = vec->nchars[ valid[i] ];
  }
  vec->nstrings = vec->next_avail = n_valid;

  vec->string = (char**)  biomcmc_realloc ((char**)  vec->string, n_valid * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_realloc ((size_t*) vec->nchars, n_valid * sizeof (size_t));
}  

void
char_vector_reduce_to_trimmed_size (char_vector vec, int new_size)
{
  int i;
  if (new_size <= vec->nstrings) return;
  for (i=vec->nstrings-1; i >= new_size; i--) if (vec->string[i]) free (vec->string[i]);
  vec->nstrings = vec->next_avail = new_size;
  vec->string = (char**)  biomcmc_realloc ((char**)  vec->string, new_size * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_realloc ((size_t*) vec->nchars, new_size * sizeof (size_t));
}

void
char_vector_reorder_by_size_or_lexicographically (char_vector vec, bool lexico, int *order) 
{
  int i;
  charvec_str *cvs;

  cvs = (charvec_str*) biomcmc_malloc (vec->nstrings * sizeof (charvec_str));
  for (i = 0; i < vec->nstrings; i++) { 
    cvs[i].s = vec->string[i];
    cvs[i].idx = i;
    cvs[i].nchars = vec->nchars[i];
  }
  if (lexico) qsort (cvs, vec->nstrings, sizeof (charvec_str), compare_charvecstr_lexicographic);
  else        qsort (cvs, vec->nstrings, sizeof (charvec_str), compare_charvecstr_decreasing);

  for (i=0; i < vec->nstrings; i++) { 
    vec->string[i] = cvs[i].s; 
    vec->nchars[i] = cvs[i].nchars; 
  }
  if (order) for (i=0; i < vec->nstrings; i++) order[i] = cvs[i].idx; // if user wants to e.g. reorder tree leaves 

  if (cvs) free (cvs);
}

bool
char_vector_link_address_if_identical (char_vector *v1, char_vector *v2)
{
  int i, k;
  bool equal = true;
  if ((*v1) == (*v2)) return true; // already the same
  if ((*v1)->nstrings != (*v2)->nstrings) return false; // must be identical
  for (i = 0; (i < (*v1)->nstrings) && ((*v1)->nchars[i] == (*v2)->nchars[i]); i++);
  if (i < (*v1)->nstrings) return false; // 
  for (i = 0; equal && (i < (*v1)->nstrings); i++)
    for (k = 0; equal && (k < (int) (*v1)->nchars[i]); k++) 
      if ((*v1)->string[i][k] != (*v2)->string[i][k]) equal = false;
  if (!equal) return false; // lines below assume, then, that both char_vectors are identical
  del_char_vector (*v2);
  *v2 = *v1;
  (*v1)->ref_counter++;
  return true;
}

void // oldname: index_sptaxa_to_genetaxa
index_species_gene_char_vectors (char_vector species, char_vector gene, int *sp_idx_in_gene, int *order_external)
{ /* Don't need gene tree, just gene leaf names; alternative is to  use gene->rec instead of gene->rec->sp_id (int*) */
  int i, j, n_index = gene->nstrings, *index, *sp_order;
  
  if (!order_external) {/* Search first largest species names (so that for example "ecoli" will match only if "ecoliII" doesn't) */
    sp_order = (int*) biomcmc_malloc (species->nstrings * sizeof (int));
    for (i=0; i < species->nstrings; i++) sp_order[i] = i;
  } // assumes species char_vector is already ordered, unless order_external is present
  else sp_order = order_external;

  index = (int*) biomcmc_malloc (n_index * sizeof (int));
  for (i=0; i < n_index; i++) { 
    sp_idx_in_gene[i] = -1; /* initialize mapping */ 
    index[i] = i;  /* scan gene leaves _without_ replacement */
  }

  for (i=0; i < species->nstrings; i++) for (j=0; j < n_index; j++) /* search sp name in all unmapped gene names */
    if ((gene->nchars[index[j]] >= species->nchars[sp_order[i]]) && /* ordered species names (by nchars[]) */
        (strcasestr (gene->string[index[j]], species->string[sp_order[i]]))) { 
      /* found species name within gene name; we have a mapping */
      sp_idx_in_gene[ index[j] ] = sp_order[i];
      index[j] = index[--n_index]; // with this the whole search takes O(N ln N),
      j--; // index[j] is now a distinct element (the last)
    }

  if (n_index) {
    fprintf (stderr, "Couldn't find species for genes:\n");
    for (i=0; i < n_index; i++) fprintf (stderr, " \"%s\"\n", gene->string[index[i]]);
    biomcmc_error ("gene names should contain the name of species");
  }

  if (!order_external) if (sp_order) free (sp_order); // only if not external (o.w. will delete it)
  if (index) free (index);
}

void
update_species_count_from_gene_char_vector (char_vector species, char_vector gene, int *sp_count)
{
  int *idx_gene_to_sp, i;
  idx_gene_to_sp = (int *) biomcmc_malloc (gene->nstrings * sizeof (int));
  index_species_gene_char_vectors (species, gene, idx_gene_to_sp, NULL); /* map species names to gene names and store into idx[] */
  for (i = 0; i < gene->nstrings; i++) sp_count[ idx_gene_to_sp[i] ]++; /* update species frequencies */
  if (idx_gene_to_sp) free (idx_gene_to_sp);
  return;
}
