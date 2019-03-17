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

void update_binary_mrp_matrix_column_if_new (binary_mrp_matrix mrp);

binary_mrp_matrix
new_binary_mrp_matrix (int n_sequences, int n_sites)
{
  int i; 
  binary_mrp_matrix mrp;
  
  mrp = (binary_mrp_matrix) biomcmc_malloc (sizeof (struct binary_mrp_matrix_struct));
  mrp->ntax  = n_sequences;
  mrp->nchar = n_sites;
  mrp->i= 0;  /* a.k.a. "current" in other structs; used when filling the matrix (from 0 to nchar) */
  mrp->ref_counter = 1;

  mrp->freq = (int*) biomcmc_malloc(mrp->nchar * sizeof (int)); 
  mrp->s = (bool**) biomcmc_malloc(mrp->ntax * sizeof (bool*)); 
  for (i = 0; i < mrp->ntax; i++)  mrp->s[i] = (bool*) biomcmc_malloc(mrp->nchar * sizeof (bool)); 
  for (i = 0; i < mrp->nchar; i++) mrp->freq[i] = 0;

  return mrp;
}

void 
del_binary_mrp_matrix (binary_mrp_matrix mrp)
{
  if (mrp) {
    int i;
    if (--mrp->ref_counter) return; /* some other place is using it, we cannot delete it yet */
    for(i=mrp->ntax-1; i>=0; i--) if (mrp->s[i]) { 
      free (mrp->s[i]); 
    }
    if (mrp->s) free (mrp->s);
    if (mrp->freq) free (mrp->freq);
    free (mrp);
  }
}

mrp_parsimony
new_mrp_parsimony (int n_sequences, int n_sites)
{
  int i;
  mrp_parsimony pars;
  pars = (mrp_parsimony) biomcmc_malloc (sizeof (struct mrp_parsimony_struct));
  pars->ref_counter = 1;
  pars->external = new_binary_mrp_matrix (n_sequences, n_sites);
  pars->internal = new_binary_mrp_matrix (n_sequences - 1, n_sites);
  pars->score = (int*) biomcmc_malloc (pars->external->nchar * sizeof (int));
  return pars;
}

void
del_mrp_parsimony (mrp_parsimony pars)
{
  if (pars) {
    if (--pars->ref_counter) return; /* some other place is using it, we cannot delete it yet */
    if (pars->score) free (pars->score);
    del_binary_mrp_matrix (pars->internal);
    del_binary_mrp_matrix (pars->external);
    free (pars);
  }
}

void
update_binary_mrp_matrix_from_topology (binary_mrp_matrix mrp, topology t, int *map)
{
  int i, j, k, ones[t->nleaves];
  bipartition bp = new_bipartition (t->nleaves);
  if (!t->traversal_updated) update_topology_traversal (t);

  for (i=0; i < t->nleaves-3; i++) { // [n-2] is root; [n-3] is leaf or redundant 
    for (j=0; j < mrp->ntax; j++) mrp->s[j][mrp->i] = 3U; // all seqs are 'N' at first (a.k.a. {0,1}) -> absent from t in the end
    for (j=0; j < t->nleaves; j++) mrp->s[map[j]][mrp->i] = 1U; // species present in t start as {0}
    bipartition_copy (bp, t->postorder[i]->split);
    bipartition_flip_to_smaller_set (bp); // not essencial, but helps finding same split in diff trees
    bipartition_to_int_vector (bp, ones, bp->n_ones); // n_ones=max number of ones to check (in this case, all of them)
    for (j=0; j < bp->n_ones; j++) mrp->s[map[ones[j]]][mrp->i] = 2U; // these species present in t are then {1} 
    update_binary_mrp_matrix_column_if_new (mrp);
    if (mrp->i > mrp->nchar) biomcmc_error ("The function calling parsimony underestimated the total number of columns (tree sizes)");
  }
  del_bipartition (bp);
}

void    
update_binary_mrp_matrix_column_if_new (binary_mrp_matrix mrp)
{
  int i, j;
  for (i=0; i < mrp->i; i++) {
    for (j=0; (j < mrp->ntax) && (mrp->s[i][j] == mrp->s[mrp->i][j]); j++); // one line loop: finishes or halts when diff is found 
    if (j == mrp->ntax) { mrp->freq[i]++; break; } // premature ntax loop end 
  }
  if (i ==  mrp->i) mrp->freq[mrp->i++] = 1; // column loop didn't stop prematurely (i.e. column was not found and thus is unique)
}

int
binary_mrp_parsimony_score_of_topology (mrp_parsimony pars, topology t)
{
  int i,j, pars_score = 0;
  bool s1, s2, intersection;
  if (!t->traversal_updated) update_topology_traversal (t);
  for (i=0; i < pars->external->nchar; i++) pars->score[i] = 0;
  for (i=0; i < pars->external->nchar; i++) { // pthreads would go here
    for (j=0; j < t->nleaves-2; j++) {
      /* id (0...nleaves) are leaves; (nleaves...2x nleaves-1) are internal nodes */
      if (t->postorder[j]->left->internal) s1 = pars->internal->s[t->postorder[j]->left->id - t->nleaves][i];
      else s1 = pars->external->s[t->postorder[j]->left->id][i];
      if (t->postorder[j]->right->internal) s2 = pars->internal->s[t->postorder[j]->right->id - t->nleaves][i];
      else s2 = pars->external->s[t->postorder[j]->right->id][i];
      intersection = s1 & s2; // 11, 01, 00, or 10 
      if (!intersection) { pars->score[i]++; intersection = s1|s2; } //00 only arises with 10 & 01  
      pars->external->s[t->postorder[j]->id - t->nleaves][i] =  intersection;
    } // for j in tree node
    pars_score += (pars->score[i] * pars->internal->freq[i]);
  } // for i in nchar 
  return pars_score;
}

/* Extra ideas/todo:
 * 1. MRL (binary likelihood) for brlens; extend it to incorporate leaf uncertainty (as is flip supertrees)
 * 2. branch-wise parsimony scores 
 * 3. store columns per tree in case of jackniffing gene trees (e.g. by returning column indexes and allowing unweighted
 * parsimony)
 * 4. actual score can be replaced by number of columns with perfect score (or almost perfect), as in compatible trees. 
 * 5. following above, score can be sum over 'best' columns (i.e. exclude 50% columns with worse score)
 * 6. scores if we remove some leaves 
 * 7. instead of 0/1 we can have 0/x where x is the split length (parsimony or LS?)
 * 8. each column is distance from leaf to all others (like mrd (distance) with Sankoff algo (or LS, UPGMA?). distances
 * normalised. 
 * <obs: here "weighted" parsimony is not Sankoff algo, but the pattern weights (column frequency)
 */ 
