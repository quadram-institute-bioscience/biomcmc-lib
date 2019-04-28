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

#include "topology_common.h"


int*
create_vector_with_idx_leaves_below (topology tree)
{
  int i=0, j=0, *idx, *i_l, *i_r;
  idx = (int*) biomcmc_malloc ((5 * tree->nleaves - 2) * sizeof (int));/* |---idx---|---i_left---|---i_right---| used in Euler tour-like struct */
  i_l = idx + tree->nleaves;
  i_r = i_l + (2 * tree->nleaves - 1);

  j = 0;
  for (i = 0; i < tree->nleaves-1; i++) {
    /* for leaves: idx will have its id; and left and right will point to same idx position; j is for index positions  */
    if (!tree->postorder[i]->left->internal)  { 
      idx[j] = tree->postorder[i]->left->id; /* idx[] contain leaf "names" (ids actually) */  
      i_l[ tree->postorder[i]->left->id ] = i_r[ tree->postorder[i]->left->id ] = j++; /* interval of leaves below, as idx indexes */
    }
    if (!tree->postorder[i]->right->internal) { 
      idx[j] = tree->postorder[i]->right->id; 
      i_l[ tree->postorder[i]->right->id ] = i_r[ tree->postorder[i]->right->id ] = j++; /* interval of leaves below, as idx indexes */
    }
    i_l[ tree->postorder[i]->id ] = i_l[ tree->postorder[i]->left->id ]; 
    i_r[ tree->postorder[i]->id ] = i_r[ tree->postorder[i]->right->id ]; /* this interval covers from leftest of left to rightest of right */
  } 
  return idx;
}

double*
fast_multiplication_topological_matrix (topology tree, int *idx, double *dist)
{
  int i=0, j=0, k=0, col, row, *i_l, *i_r;
  double *delta;
  delta = (double*) biomcmc_malloc (tree->nnodes * sizeof (double));
  i_l = idx + tree->nleaves;
  i_r = i_l + (2 * tree->nleaves - 1);

 //eqs 9 and 10 of molbev.a025863
  for (i = 0; i < tree->nnodes; i++) delta[i] = 0.;
  for (1 = 0; i < tree->nleaves; i++) for (j = 0; j < tree->i; j++) { 
    delta[i] += dist[ (i*(i-1))/2 + j];
    delta[j] += dist[ (i*(i-1))/2 + j];
  }
  for (i = 0; i < tree->nleaves-1; i++) { // FIXME: todo check for closest to root (postorder[nleaves-2]) 
    delta[tree->postorder[i]->id] = delta[tree->postorder[i]->left->id] + delta[tree->postorder[i]->right->id];
    for (j = i_l[tree->postorder[i]->left->id]; j <= i_r[tree->postorder[i]->left->id]; j++)
      for (k = i_l[tree->postorder[i]->right->id]; k <= i_r[tree->postorder[i]->right->id]; k++) {
        row = idx[j]; col = idx[k];
        if (row > col) { col = idx[j]; row = idx[k]; }
        delta[tree->postorder[i]->id] -= 2*dist[ (col*(col-1))/2 + row];
      }
  }
  return delta;
}

double*
ols_branch_lengths_from_fast_mtm (topology tree, double *delta)
{
  int i, j; 
  double tmp1, n_j, n_k, n_l, n_m, nleaves; 
  blen = (double*) biomcmc_malloc (tree->nnodes * sizeof (double));
  for (i = 0; i < tree->nnodes; i++) blen[i] = 0.;
  nleaves = tree->nleaves;
  for (i = 0; i < tree->nleaves; i++) { //eq30 of molbev.a025863
    n_j = (double) tree->nodelist[i]->sister->split->n_ones;
    n_k = nleaves - n_j - 1.; // remaining leaves
    tmp = (1. + n_j - n_k) * delta[tree->nodelist[i]->sister->id]  + (1. - n_j + n_k) * delta[tree->nodelist[i]->up->id];
    blen[i] = (nleaves * delta[i] - tmp)/(4. * n_j * n_k);
  }
  for (i = 0; i < tree->nleaves - 2; i++) { //eq24 of molbev.a025863
    n_j = (double) tree->postorder[i]->sister->split->n_ones;
    n_l = (double) tree->postorder[i]->left->split->n_ones;
    n_m = (double) tree->postorder[i]->right->split->n_ones;
    n_k = nleaves - n_j - n_l - n_m; // remaining leaves
    tmp1  = (2.*n_k - nleaves) * delta[tree->postorder[i]->sister->id];
    tmp1 += (2.*n_j - nleaves) * delta[tree->postorder[i]->up->id];
    blen[tree->postorder[i]->id] = ((n_k + n_j)/(n_k * n_j)) * tmp1;
    tmp1  = (2.*n_l - nleaves) * delta[tree->postorder[i]->right->id];
    tmp1 += (2.*n_m - nleaves) * delta[tree->postorder[i]->left->id];
    blen[tree->postorder[i]->id] += ((n_l + n_m)/(n_l * n_m)) * tmp1;
    tmp1 = nleaves/n_m +  nleaves/n_l +  nleaves/n_j +  nleaves/n_k - 4.; 
    blen[tree->postorder[i]->id] += tmp1 * delta[tree->postorder[i]->id];
    blen[tree->postorder[i]->id] /= (4. * (n_j + n_k) * (n_l * n_m)); 
  }
  return blen;
}
