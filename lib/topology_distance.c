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

#include "topology_distance.h"

/*! \brief rescales original branches and stores distance from root into #fromroot[], using up to 6 distinct rescalings */
double rescale_rooted_distances_for_patristic_distances (topology tree, double *fromroot, int mode, double tolerance);
double* fast_multiplication_topological_matrix (topology tree, int *idx, double *dist);
double* ols_branch_lengths_from_fast_mtm (topology tree, double *delta);


distance_matrix
new_distance_matrix_for_topology (int nleaves)
{
  int i;
  distance_matrix dist = new_distance_matrix (nleaves);
  dist->fromroot = (double*) biomcmc_malloc ((2 * nleaves - 1) * sizeof (double));
  /* |---idx---|---i_left---|---i_right---| used in Euler tour-like struct */
  dist->idx = (int*) biomcmc_malloc ((5 * nleaves - 2) * sizeof (int));
  dist->i_l = dist->idx + nleaves;
  dist->i_r = dist->i_l + (2 * nleaves - 1);
  for (i = 0; i < 2 * nleaves - 1; i++) dist->fromroot[i] = 0.;
  return dist;
}

void
fill_distance_matrix_from_topology (distance_matrix dist, topology tree, double *blen, bool use_upper)
{
  int i, j = 0, k, row, col;
  if (dist->size > tree->nleaves) biomcmc_error ("distance matrix is smaller than number of leaves from tree");
  if (!tree->traversal_updated) update_topology_traversal (tree);
  /* STEP 1: find distances from every node to root */
  if (!blen) for (i = 0; i < tree->nnodes; i++) dist->fromroot[i] = (double)(tree->nodelist[i]->level); /* level = nodal distance from root */
  else {
    dist->fromroot[ tree->root->id ] = 0.;
    for (i = tree->nleaves-3; i >= 0; i--)  /* internal nodes */
      dist->fromroot[ tree->postorder[i]->id ] = dist->fromroot[ tree->postorder[i]->up->id ] + blen[ tree->postorder[i]->id ];
    for (i = 0; i < tree->nleaves; i++) /* external nodes (do not belong to postorder) */
      dist->fromroot[ tree->nodelist[i]->id ] = dist->fromroot[ tree->nodelist[i]->up->id ] + blen[ tree->nodelist[i]->id ];
  }
  /* STEP 2: create tour in postorder so that we have subvectors with all leaves below it */
  j = 0;
  for (i = 0; i < tree->nleaves-1; i++) {
    /* for leaves: idx will have its id; and left and right will point to same idx position; j is for index positions  */
    if (!tree->postorder[i]->left->internal)  { 
      dist->idx[j] = tree->postorder[i]->left->id; /* idx[] contain leaf "names" (ids actually) */  
      dist->i_l[ tree->postorder[i]->left->id ] = dist->i_r[ tree->postorder[i]->left->id ] = j++; /* interval of leaves below, as idx indexes */
    }
    if (!tree->postorder[i]->right->internal) { 
      dist->idx[j] = tree->postorder[i]->right->id; 
      dist->i_l[ tree->postorder[i]->right->id ] = dist->i_r[ tree->postorder[i]->right->id ] = j++; /* interval of leaves below, as idx indexes */
    }
    dist->i_l[ tree->postorder[i]->id ] = dist->i_l[ tree->postorder[i]->left->id ]; 
    dist->i_r[ tree->postorder[i]->id ] = dist->i_r[ tree->postorder[i]->right->id ]; /* this interval covers from leftest of left to rightest of right */
  } 
  /* STEP 3: dist(A,B) = fromroot[A] + fromroot[B] - 2 * fromroot[mrca between A and B] (from STEP2 we know all A's and B's)*/
  if (use_upper) for (i = 0; i < tree->nleaves; i++) for (j = i; j < tree->nleaves; j++) dist->d[i][j] = 0.;
  else           for (i = 0; i < tree->nleaves; i++) for (j = 0; j <= i; j++)            dist->d[i][j] = 0.;

  for (i = 0; i < tree->nleaves-1; i++) 
    for (j = dist->i_l[tree->postorder[i]->left->id]; j <= dist->i_r[tree->postorder[i]->left->id]; j++)
      for (k = dist->i_l[tree->postorder[i]->right->id]; k <= dist->i_r[tree->postorder[i]->right->id]; k++) {
        row = dist->idx[j]; col = dist->idx[k];
        if (((row > col) && use_upper) || ((row < col) && !use_upper)) { col = dist->idx[j]; row = dist->idx[k]; }
        dist->d[row][col] = dist->fromroot[row] + dist->fromroot[col] - 2 * dist->fromroot[ tree->postorder[i]->id ]; 
      }
  //for (i=0;i<dist->size;i++) {for (j=0; j<dist->size;j++) printf("%12.10g ", dist->d[i][j]);printf (" DEBUG\n");}
}

void 
patristic_distances_from_topology_to_vectors (topology tree, double **dist, double *scaling, int n_dists, double tolerance)
{
  int i=0, j=0, k=0, l=0, row, col, *idx, *i_l, *i_r, onedim;
  if (!tree->traversal_updated) update_topology_traversal (tree);
  double **fromroot = NULL;

  if (tolerance < 1e-15) tolerance = 1e-15; /* any branch length shorter than this will be assumed zero (multifurcation) */
  if (n_dists > 6) { n_dists = 6; fprintf(stderr, "biomcmc WARNING: more than 6 patristic distances rescaling requested\n"); }
  fromroot = (double**) biomcmc_malloc (n_dists * sizeof (double*));
  for (i = 0; i < n_dists; i++) fromroot[i] = (double*) biomcmc_malloc ((2 * tree->nleaves - 1) * sizeof (double));

  /* STEP 1: find distances from every node to root */
  for (i = 0; i < n_dists; i++) if (dist[i]) scaling[i] = rescale_rooted_distances_for_patristic_distances (tree, fromroot[i], i, tolerance);

  /* STEP 2: create tour in postorder so that we have subvectors with all leaves below it */
  idx = create_vector_with_idx_leaves_below_for_patristic (tree);
  i_l = idx + tree->nleaves;
  i_r = i_l + (2 * tree->nleaves - 1);
  /* STEP 3: dist(A,B) = fromroot[A] + fromroot[B] - 2 * fromroot[mrca between A and B] (from STEP2 we know all A's and B's)*/
  for (i = 0; i < ((tree->nleaves - 1) * tree->nleaves)/2; i++) for (j = 0; j < n_dists; j++) dist[j][i] = 0.; 

  for (i = 0; i < tree->nleaves-1; i++) 
    for (j = i_l[tree->postorder[i]->left->id]; j <= i_r[tree->postorder[i]->left->id]; j++)
      for (k = i_l[tree->postorder[i]->right->id]; k <= i_r[tree->postorder[i]->right->id]; k++) {
        row = idx[j]; col = idx[k];
        if (row > col) { col = idx[j]; row = idx[k]; }
        onedim = (col * (col-1))/2 + row;
        for (l=0; l < n_dists; l++) if (dist[l]) dist[l][onedim] = fromroot[l][row] + fromroot[l][col] - 2 * fromroot[l][ tree->postorder[i]->id ];
      }
  if (fromroot) {
    for (i = n_dists-1; i >= 0; i--) if (fromroot[i]) free (fromroot[i]);
    free (fromroot);
  }
  if (idx) free (idx);
}

int*
create_vector_with_idx_leaves_below_for_patristic (topology tree)
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

double
rescale_rooted_distances_for_patristic_distances (topology tree, double *fromroot, int mode, double tolerance)
{
  int i;
  double *internal_d, *external_d, scale = 1.;

  internal_d = (double*) biomcmc_malloc ((tree->nleaves - 2) * sizeof (double));
  external_d = (double*) biomcmc_malloc ((tree->nleaves) * sizeof (double));
  
  fromroot[ tree->root->id ] = 0.;

  switch (mode) {
    case 0:  /*  nodal distance (number of edges > tolerance)  */
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = (tree->blength[ tree->postorder[i]->id ] > tolerance) ? 1. : 0.;
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = (tree->blength[ tree->nodelist[i]->id ] > tolerance) ? 1. : 0.; 
      break;
    case 1: /* divided by average (s.t. average blen is 1. and treelength = nnodes) */ 
      scale = 0.;
      for (i = 0; i < tree->nnodes; i++) scale += tree->blength[tree->nodelist[i]->id];
      scale /= tree->nnodes;
      if (scale < 1e-12) scale = 1e-12;
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = tree->blength[ tree->postorder[i]->id ] / scale;
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = tree->blength[ tree->nodelist[i]->id ] / scale;
      break;
    case 2:  /*  unscaled (original) distance  */
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = tree->blength[ tree->postorder[i]->id ];
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = tree->blength[ tree->nodelist[i]->id ];
      break;
    case 3: /* divided by number of nodes */ 
      scale = tree->nnodes;
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = tree->blength[ tree->postorder[i]->id ] / scale;
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = tree->blength[ tree->nodelist[i]->id ] / scale;
      break;
    case 4: /* divided by treelength (s.t. average is 1/nodes and treelength=1) */ 
      scale = 0.;
      for (i = 0; i < tree->nnodes; i++) scale += tree->blength[tree->nodelist[i]->id];
      if (scale < 1e-12) scale = 1e-12;
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = tree->blength[ tree->postorder[i]->id ] / scale;
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = tree->blength[ tree->nodelist[i]->id ] / scale;
      break;
    default: /* divided by shortest branch */ 
      scale = 1e9;
      for (i = 0; i < tree->nnodes; i++) 
        if ((scale > tree->blength[tree->nodelist[i]->id]) && (tree->blength[tree->nodelist[i]->id] > tolerance))
          scale = tree->blength[tree->nodelist[i]->id];
      if (scale > 1e8) scale = 1.;
      for (i = tree->nleaves-3; i >= 0; i--) internal_d[i] = tree->blength[ tree->postorder[i]->id ] / scale;
      for (i = 0; i < tree->nleaves; i++)    external_d[i] = tree->blength[ tree->nodelist[i]->id ] / scale;
      break;
  }
  /* main function: calculate distance from root to node in preorder */
  for (i = tree->nleaves-3; i >= 0; i--) 
    fromroot[ tree->postorder[i]->id ] = fromroot[ tree->postorder[i]->up->id ] + internal_d[i];
  for (i = 0; i < tree->nleaves; i++)  /* external nodes (do not belong to postorder) */
    fromroot[ tree->nodelist[i]->id ] = fromroot[ tree->nodelist[i]->up->id ] + external_d[i];

  if (internal_d) free (internal_d);
  if (external_d) free (external_d);
  return scale;
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
  for (i = 1; i < tree->nleaves; i++) for (j = 0; j < i; j++) { 
    delta[i] += dist[ (i*(i-1))/2 + j];
    delta[j] += dist[ (i*(i-1))/2 + j];
  }
  for (i = 0; i < tree->nleaves-1; i++) if (tree->postorder[i]->id != tree->root->left->id) { 
    delta[tree->postorder[i]->id] = delta[tree->postorder[i]->left->id] + delta[tree->postorder[i]->right->id];
    for (j = i_l[tree->postorder[i]->left->id]; j <= i_r[tree->postorder[i]->left->id]; j++)
      for (k = i_l[tree->postorder[i]->right->id]; k <= i_r[tree->postorder[i]->right->id]; k++) {
        row = idx[j]; col = idx[k];
        if (row > col) { col = idx[j]; row = idx[k]; }
        delta[tree->postorder[i]->id] -= 2*dist[ (col*(col-1))/2 + row];
      }
  }
  /* same information as right branch (and |right| <= |left| */
  delta[tree->root->left->id] = delta[tree->root->right->id]; 
  return delta;
}

double*
ols_branch_lengths_from_fast_mtm (topology tree, double *delta)
{
  int i; 
  double *blen = NULL, tmp1, n_j, n_k, n_l, n_m, nleaves; 
  blen = (double*) biomcmc_malloc (tree->nnodes * sizeof (double));
  for (i = 0; i < tree->nnodes; i++) blen[i] = 0.;
  nleaves = tree->nleaves;
  for (i = 0; i < tree->nleaves; i++) { //eq30 of molbev.a025863
    if (tree->nodelist[i]->up != tree->root) {
      n_j = (double) tree->nodelist[i]->sister->split->n_ones;
      n_k = nleaves - n_j - 1.; // remaining leaves
      tmp1 = (1. + n_j - n_k) * delta[tree->nodelist[i]->sister->id]  + (1. - n_j + n_k) * delta[tree->nodelist[i]->up->id];
      blen[i] = (nleaves * delta[i] - tmp1)/(4. * n_j * n_k);
    }
    else {
      n_j = (double) tree->nodelist[i]->sister->left->split->n_ones;
      n_k = (double) tree->nodelist[i]->sister->right->split->n_ones;
      tmp1 = (1. + n_j - n_k) * delta[tree->nodelist[i]->sister->left->id]  + (1. - n_j + n_k) * delta[tree->nodelist[i]->sister->right->id];
      blen[i] = (nleaves * delta[i] - tmp1)/(8. * n_j * n_k); // half here...
      blen[tree->nodelist[i]->sister->id] = blen[i]; // ... and half to the sister branch (skipped below)
      delta[tree->root->id] = delta[i]; // copy delta info from this leaf into root 
    }
  }
  for (i = 0; i < tree->nleaves - 3; i++) { //eq24 of molbev.a025863 (skip right daughter of root, even if internal)
    if (tree->postorder[i]->up != tree->root) {
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
      blen[tree->postorder[i]->id] /= (4. * (n_j + n_k) * (n_l + n_m)); 
    }
    else { // sister is internal node (o.w. postorder is nleaves-2)
      n_j = (double) tree->postorder[i]->sister->left->split->n_ones;
      n_k = (double) tree->postorder[i]->sister->right->split->n_ones;
      n_l = (double) tree->postorder[i]->left->split->n_ones;
      n_m = (double) tree->postorder[i]->right->split->n_ones;
      tmp1  = (2.*n_k - nleaves) * delta[tree->postorder[i]->sister->left->id];
      tmp1 += (2.*n_j - nleaves) * delta[tree->postorder[i]->sister->right->id];
      blen[tree->postorder[i]->id] = ((n_k + n_j)/(n_k * n_j)) * tmp1;
      tmp1  = (2.*n_l - nleaves) * delta[tree->postorder[i]->right->id];
      tmp1 += (2.*n_m - nleaves) * delta[tree->postorder[i]->left->id];
      blen[tree->postorder[i]->id] += ((n_l + n_m)/(n_l * n_m)) * tmp1;
      tmp1 = nleaves/n_m +  nleaves/n_l +  nleaves/n_j +  nleaves/n_k - 4.; 
      blen[tree->postorder[i]->id] += tmp1 * delta[tree->postorder[i]->id];
      blen[tree->postorder[i]->id] /= (8. * (n_j + n_k) * (n_l + n_m)); // half the true value   
      blen[tree->postorder[i]->sister->id] = blen[tree->postorder[i]->id];
    //  printf ("DEBUG::<2> %d %lf %lf %d \n", i, n_k, delta[tree->postorder[i]->up->id], tree->postorder[i]->up->id);
    }
  }
  correct_negative_branch_lengths_from_topology (tree, blen);
  return blen;
}

void
estimate_topology_branch_lengths_from_distances (topology tree, double *dist)
{
  double *blen;
  blen = new_topology_branch_lengths_from_distances (tree, dist);
  if (tree->blength) free (tree->blength);
  tree->blength = blen;
}

double*
new_topology_branch_lengths_from_distances (topology tree, double *dist)
{
  int *idx;
  double *delta, *blen;

  idx = create_vector_with_idx_leaves_below_for_patristic (tree);
  delta = fast_multiplication_topological_matrix (tree, idx, dist);
  blen = ols_branch_lengths_from_fast_mtm (tree, delta);

  if (idx) free (idx);
  if (delta) free (delta);
  return blen;
}

void
correct_negative_branch_lengths_from_topology (topology t, double *blength)
{
  int i;

  for (i=0; i < t->nleaves-1; i++) { // postorder are internal nodes only
    if (blength[t->postorder[i]->left->id] < DBL_MIN) { 
      blength[t->postorder[i]->id] -= blength[t->postorder[i]->left->id]; // left is negative number
      blength[t->postorder[i]->left->id] = 0.;
    }
    if (blength[t->postorder[i]->right->id] < DBL_MIN) { 
      blength[t->postorder[i]->id] -= blength[t->postorder[i]->right->id];
      blength[t->postorder[i]->right->id] = 0.;
    }
  }
  if (blength[t->root->id] > 0.) {
    blength[t->root->left->id]  += blength[t->root->id];
    blength[t->root->right->id] += blength[t->root->id];
    blength[t->root->id] = 0.;
  } 
}
