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

binary_mrp_matrix
new_binary_mrp_matrix (int n_sequences, int n_sites)
{
  int i; 
  binary_mrp_matrix mrp;
  
  mrp = (binary_mrp_matrix) biomcmc_malloc (sizeof (struct binary_mrp_matrix_struct));
  mrp->ntax  = n_sequences;
  mrp->nchar = n_sites;
  mrp->npat  = 0; // used also when filling the matrix (from 0 to nchar)
  mrp->pattern_freq = NULL;
  mrp->ref_counter = 1;

  mrp->s = (bool**) biomcmc_malloc(mrp->ntax * sizeof (bool*)); 
  for (i = 0; i < mrp->ntax; i++) 
    mrp->s[i] = (bool*) biomcmc_malloc(mrp->nchar * sizeof (bool)); 

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
    if (mrp->pattern_freq) free (mrp->pattern_freq);
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
  int i, ones[t->nleaves];
  bipartition bp = new_bipartition (t->nleaves);
  // TODO: maybe check if bipartition exists already here?
  if (!t->traversal_updated) update_topology_traversal (t);
  for (i=0; i < t->nleaves-3; i++) { /* [n-2] is root; [n-3] is leaf or redundant */
    for (j=0; j < mrp->ntax; j++) mrp->s[j][mrp->npat] = 3U; // all seqs are 'N' at first (a.k.a. {0,1}) -> absent from t in the end
    for (j=0; j < t->nleaves; j++) mrp->s[map[j]][mrp->npat] = 1U; // species present in t start as {0}
    bipartition_copy (bp, t->postorder[i]->split);
    bipartition_flip_to_smaller_set (bp); // not essencial, but helps finding same split in diff trees
    bipartition_to_int_vector (bp, ones, bp->n_ones); // n_ones=max number of ones to check (in this case, all of them)
    for (j=0; j < bp-<n_ones; j++) mrp->s[map[ones[j]]][mrp->npat] = 2U; // these species present in t are then {1} 
    mrp->npat++;
  }

  del_bipartition (bp);
}
