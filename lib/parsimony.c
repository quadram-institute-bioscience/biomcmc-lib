/*
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "parsimony.h"

int update_biparsdatmat_column_from_ones (binary_parsimony_datamatrix mrp, topology t, int *ones, bipartition bp, int *map);
void update_binary_parsimony_length (binary_parsimony pars, int new_columns_size);
void update_binary_parsimony_datamatrix_column_if_new (binary_parsimony_datamatrix mrp);
uint32_t hash_value_of_binary_parsimony_datamatrix_column (binary_parsimony_datamatrix mrp, int idx);

binary_parsimony_datamatrix
new_binary_parsimony_datamatrix (int n_sequences)
{
  int i; 
  binary_parsimony_datamatrix mrp;

  mrp = (binary_parsimony_datamatrix) biomcmc_malloc (sizeof (struct binary_parsimony_datamatrix_struct));
  mrp->ntax  = n_sequences;
  mrp->nchar = mrp->i = 0;
  mrp->ref_counter = 1;
  mrp->freq_sum = 0;
  mrp->s = (bool**) biomcmc_malloc(mrp->ntax * sizeof (bool*)); 
  for (i = 0; i < mrp->ntax; i++)  mrp->s[i] = NULL; 
  mrp->freq = NULL;
  mrp->col_hash = NULL;
  mrp->occupancy = NULL;

  return mrp;
}

binary_parsimony_datamatrix
new_binary_parsimony_datamatrix_fixed_length (int n_sequences, int n_sites)
{
  int i; 
  binary_parsimony_datamatrix mrp = new_binary_parsimony_datamatrix (n_sequences);
  mrp->nchar = n_sites;

  mrp->freq = (int*) biomcmc_malloc (mrp->nchar * sizeof (int)); 
  mrp->col_hash = (uint32_t*) biomcmc_malloc (mrp->nchar * sizeof (uint32_t)); 
  mrp->occupancy = (int*) biomcmc_malloc (mrp->nchar * sizeof (int)); 
  for (i = 0; i < mrp->ntax; i++)  mrp->s[i] = (bool*) biomcmc_malloc(mrp->nchar * sizeof (bool)); 
  for (i = 0; i < mrp->nchar; i++) mrp->freq[i] = mrp->occupancy[i] = 0;

  return mrp;
}

void 
del_binary_parsimony_datamatrix (binary_parsimony_datamatrix mrp)
{
  if (mrp) {
    int i;
    if (--mrp->ref_counter) return; /* some other place is using it, we cannot delete it yet */
    for(i=mrp->ntax-1; i>=0; i--) if (mrp->s[i]) { 
      free (mrp->s[i]); 
    }
    if (mrp->s) free (mrp->s);
    if (mrp->freq) free (mrp->freq);
    if (mrp->col_hash) free (mrp->col_hash);
    if (mrp->occupancy) free (mrp->occupancy);
    free (mrp);
  }
}

binary_parsimony
new_binary_parsimony (int n_sequences)
{
  binary_parsimony pars;
  pars = (binary_parsimony) biomcmc_malloc (sizeof (struct binary_parsimony_struct));
  pars->ref_counter = 1;
  pars->external = new_binary_parsimony_datamatrix (n_sequences);
  pars->internal = new_binary_parsimony_datamatrix (n_sequences - 1);
  pars->score = NULL; 
  return pars;
}

binary_parsimony
new_binary_parsimony_fixed_length (int n_sequences, int n_sites)
{
  binary_parsimony pars;
  pars = (binary_parsimony) biomcmc_malloc (sizeof (struct binary_parsimony_struct));
  pars->ref_counter = 1;
  pars->external = new_binary_parsimony_datamatrix_fixed_length (n_sequences, n_sites);
  pars->internal = new_binary_parsimony_datamatrix_fixed_length (n_sequences - 1, n_sites);
  if (pars->internal->freq) free (pars->internal->freq);
  if (pars->internal->col_hash) free (pars->internal->col_hash);
  if (pars->internal->occupancy) free (pars->internal->occupancy);
  pars->score = (int*) biomcmc_malloc (pars->external->nchar * sizeof (int));
  return pars;
}

void
del_binary_parsimony (binary_parsimony pars)
{
  if (pars) {
    if (--pars->ref_counter) return; /* some other place is using it, we cannot delete it yet */
    if (pars->score) free (pars->score);
    del_binary_parsimony_datamatrix (pars->internal);
    del_binary_parsimony_datamatrix (pars->external);
    free (pars);
  }
}

void
update_binary_parsimony_from_topology (binary_parsimony pars, topology t, int *map, int n_sp)
{
  int i, n1, n2, *ones = NULL;
  bipartition bp = new_bipartition (t->nleaves);
  bipartition not = new_bipartition (t->nleaves);
  binary_parsimony_datamatrix mrp = pars->external;

  if (!t->traversal_updated) update_topology_traversal (t);
  update_binary_parsimony_length (pars, 2*(t->nleaves-3));

  ones = (int*) biomcmc_malloc (t->nleaves * sizeof (int));
  for (i=0; i < t->nleaves-3; i++) { // [n-2] is root; [n-3] is leaf or redundant 
    bipartition_copy (bp, t->postorder[i]->split);
    bipartition_to_int_vector (bp, ones, bp->n_ones); // n_ones=max number of ones to check (in this case, all of them)
    n1 = update_biparsdatmat_column_from_ones (mrp, t, ones, bp, map);
    if ((n1 > 1) && (n1 < n_sp - 1)) {
      mrp->occupancy[mrp->i] = n_sp; // completeness (# species represented in bipartition)
      update_binary_parsimony_datamatrix_column_if_new (mrp);
    }
    bipartition_NOT (not, bp);
    bipartition_to_int_vector (not, ones, not->n_ones); // n_ones=max number of ones to check (in this case, all of them)
    n2 = update_biparsdatmat_column_from_ones (mrp, t, ones, not, map);
    if ((n2 > 1) && (n2 < n_sp - 1) && ((n2 + n1) != n_sp)) { 
      mrp->occupancy[mrp->i] = n_sp; // completeness (# species represented in bipartition)
      update_binary_parsimony_datamatrix_column_if_new (mrp);
    }
    if (mrp->i > mrp->nchar) biomcmc_error ("The function calling parsimony underestimated the total number of columns (tree sizes)");
  }
  pars->external->freq_sum += t->nleaves-3; 
  del_bipartition (bp);
  del_bipartition (not);
  if (ones) free (ones);
}

int
update_biparsdatmat_column_from_ones (binary_parsimony_datamatrix mrp, topology t, int *ones, bipartition bp, int *map)
{
  int j, nsp = 0;
  for (j=0; j < mrp->ntax; j++)  mrp->s[         j  ][mrp->i] = 3U; // all seqs are 'N' at first (a.k.a. {0,1}) -> absent from t in the end
  for (j=0; j < t->nleaves; j++) mrp->s[     map[j] ][mrp->i] = 1U; // species present in t start as 
  for (j=0; j < bp->n_ones; j++) mrp->s[map[ones[j]]][mrp->i] = 2U; // these species present in t are then {1} 
  for (j=0; j < mrp->ntax; j++) if (mrp->s[j][mrp->i] == 2U) nsp++;
  return nsp;
}

void
update_binary_parsimony_length (binary_parsimony pars, int new_columns_size)
{
  int i, new_size = pars->external->i + new_columns_size + 1;
  pars->external->nchar = pars->internal->nchar = new_size;
  pars->score = (int*) biomcmc_realloc ((int*) pars->score, new_size * sizeof (int));
  pars->external->freq = (int*) biomcmc_realloc ((int*) pars->external->freq, new_size * sizeof (int));
  pars->external->col_hash = (uint32_t*) biomcmc_realloc ((uint32_t*) pars->external->col_hash, new_size * sizeof (uint32_t));
  pars->external->occupancy = (int*) biomcmc_realloc ((int*) pars->external->occupancy, new_size * sizeof (int));
  // notice that internal->freq and internal->col_hash are NOT updated (should be NULL); external->ntax != internal->nchar
  for (i = 0; i < pars->external->ntax; i++) pars->external->s[i] = (bool*) biomcmc_realloc((bool*) pars->external->s[i], new_size * sizeof (bool)); 
  for (i = 0; i < pars->internal->ntax; i++) pars->internal->s[i] = (bool*) biomcmc_realloc((bool*) pars->internal->s[i], new_size * sizeof (bool));
  for (i = pars->external->i; i < new_size; i++) pars->external->freq[i] = pars->external->occupancy[i] = 0;
}

void    
update_binary_parsimony_datamatrix_column_if_new (binary_parsimony_datamatrix mrp)
{
  int i, j;
  uint32_t hashv = hash_value_of_binary_parsimony_datamatrix_column (mrp, mrp->i);

  for (i=0; i < mrp->i; i++) if ((hashv == mrp->col_hash[i]) && (mrp->occupancy[i] == mrp->occupancy[mrp->i])) { 
    // same hash value and occupancy; may be identical or may be coincidence (collision)
    for (j=0; (j < mrp->ntax) && (mrp->s[j][i] == mrp->s[j][mrp->i]); j++); // one line loop: finishes or halts when diff is found 
    if (j == mrp->ntax) { mrp->freq[i]++; return; } // premature loop end means that they're distinct; here they're the same
  }
  if (i == mrp->i) { // column loop didn't stop prematurely (i.e. column was not found and thus is unique)
    mrp->freq[mrp->i] = 1;
    mrp->col_hash[mrp->i++] = hashv;
  }
}

uint32_t 
hash_value_of_binary_parsimony_datamatrix_column (binary_parsimony_datamatrix mrp, int idx)
{ // very unoptimised, since s[i][idx] is 1,2,3 (less than 8 bits) ?
  int i;
  uint32_t hashv = biomcmc_hashint_salted ((uint32_t) mrp->s[0][idx], /*salt*/ 1);
  for (i=1; i < mrp->ntax; i++) hashv = biomcmc_hashint_mix_salted (hashv, (uint32_t) mrp->s[i][idx], /*salt*/ 2);
  return hashv;
}

int
binary_parsimony_score_of_topology (binary_parsimony pars, topology t)
{
  int i,j, pars_score = 0, incompatible = 0;
  double  incomplete = 0., complete = 0.;
  bool s1, s2, intersection;
  if (!t->traversal_updated) update_topology_traversal (t);
  for (i=0; i < pars->external->i; i++) pars->score[i] = 0;  // external->i < external->nchar since may have duplicates
#ifdef _OPENMP
#pragma omp parallel for shared(pars, t) \
  private(i,j,s1,s2,intersection) reduction (+:pars_score, incompatible, incomplete)
#endif
  for (i=0; i < pars->external->i; i++) { // pthreads would go here
    for (j=0; j < t->nleaves-2; j++) {
      /* id (0...nleaves) are leaves; (nleaves...2x nleaves-1) are internal nodes */
      if (t->postorder[j]->left->internal) s1 = pars->internal->s[t->postorder[j]->left->id - t->nleaves][i];
      else s1 = pars->external->s[t->postorder[j]->left->id][i];
      if (t->postorder[j]->right->internal) s2 = pars->internal->s[t->postorder[j]->right->id - t->nleaves][i];
      else s2 = pars->external->s[t->postorder[j]->right->id][i];
      intersection = s1 & s2; // 11, 01, 00, or 10 
      if (!intersection) { pars->score[i]++; intersection = s1|s2; } //00 only arises with 10 & 01  
      pars->internal->s[t->postorder[j]->id - t->nleaves][i] =  intersection;
    } // for j in tree node
    pars_score += (pars->score[i] * pars->external->freq[i]); // only external has freqs
    if (pars->score[i] > 1) incompatible += pars->external->freq[i];
    incomplete += (double) (pars->score[i] * pars->external->freq[i]) / (double) (pars->external->occupancy[i]); // trees w more species are less penalised
    complete += (double) (pars->score[i] * pars->external->freq[i]) / (double) (pars->external->ntax - pars->external->occupancy[i] + 1); // bigger trees more penalised
  } // for i in nchar 
  pars->costs[0] = pars_score;
  pars->costs[1] = incompatible;
  pars->costs[2] = incomplete;
  pars->costs[3] = complete;
  return pars_score;
}

void
pairwise_distances_from_binary_parsimony_datamatrix (binary_parsimony_datamatrix mrp, double **dist, int size_dist)
{
  int i,j,k, idx, dist_int;
  double dist_dbl;
  for (j = 1; j < mrp->ntax; j++) for (i = 0; i < j; i++) {
    idx = (j * (j-1)) /2 + i;
    dist_int = 0; // {10}^{01}={11}   
    dist_dbl = 0.;
    for (k = 0; k < mrp->i; k++) if (((mrp->s[i][k] ^ mrp->s[j][k]) & 3U) == 3U) { // only {10} to {01} 
      dist_int += mrp->freq[k];
      dist_dbl += ((double) (mrp->freq[k]) / biomcmc_log1p ((double) (mrp->occupancy[i])));
    }
    /* different distances/scalings: */
    dist[0][idx] = (double) (dist_int)/ (double) (mrp->freq_sum);
    if (size_dist > 0) dist[1][idx] = dist_dbl;
    if (size_dist > 1) {
      int invalid_pairs = 0;
      for (k = 0; k < mrp->i; k++) if ((mrp->s[i][k] + mrp->s[j][k]) == 6U) invalid_pairs++; 
      dist[2][idx] = (double) (dist_int)/ ((double) (mrp->freq_sum - invalid_pairs) + 1); // avoid division by zero
    }
  }
}

/* Extra ideas/todo:
 * 1. MRL (binary likelihood) for brlens; extend it to incorporate leaf uncertainty (as is flip supertrees)
 * 2. branch-wise parsimony scores 
 * 3. store columns per tree in case of jackniffing gene trees (e.g. by returning column indexes and allowing unweighted
 * parsimony)
 * 4. <done> actual score can be replaced by number of columns with perfect score (or almost perfect), as in compatible trees. 
 * 5. following above, score can be sum over 'best' columns (i.e. exclude 50% columns with worse score)
 * 6. scores if we remove some leaves 
 * 7. instead of 0/1 we can have 0/x where x is the split length (parsimony or LS?)
 * 8. each column is distance from leaf to all others (like mrd (distance) with Sankoff algo (or LS, UPGMA?). distances
 * normalised. 
 * <obs: here "weighted" parsimony is not Sankoff algo, but the pattern weights (column frequency)
 */ 
