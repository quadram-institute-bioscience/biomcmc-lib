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

#include "upgma.h"

void
upgma_from_distance_matrix (topology tree, distance_matrix dist, bool single_linkage) 
{ /* always upper diagonal (that is, only i < j in d[i][j]) */
  int i, j, idx_i, parent = tree->nleaves, idx_j, min_row, min_col, row, col, idx_col, n_idx = tree->nleaves,
      *idx = tree->index,                            /* indexes in UPGMA */
      *idxtree = tree->index + tree->nleaves,        /* indexes in tree (since have values > nleaves) */
      *min_by_row = tree->index + 2 * tree->nleaves; /* column having min value for each row */
  double *dst_by_row, dst_row, new_dist, *gsize, *height, gs1, gs2;

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;
  
  gsize      = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* number of leaves below node */
  height     = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* distance from node to tips (ultrametric) */ 
  dst_by_row = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* min value itself for each row */

  for (i=0; i < n_idx; i++) { 
    idx[i]     = i; /* index to actual vector element for UPGMA distance matrix */
    idxtree[i] = i; /* index to actual vector element for tree nodes */
    gsize[i]   = 1.; /* group size (# elements below node in UPGMA */
    height[i]  = 0.; /* distance of node from present */
  }

  /* arbitrary initial values for vector of row-wise minimum distances and overall minimum */
  for (i=0; i < dist->size - 1; i++) dst_by_row[i] = 1.e35; 
  /* update vector with minimum distances per row */
  for (j=1; j < n_idx; j++) for (i=0; i < j; i++)
    if (dist->d[i][j] < dst_by_row[i]) { dst_by_row[i] = dist->d[i][j]; min_by_row[i] = j; }

  while (n_idx > 2) { /* draw two nodes to be connected */
    /* find vector elements with overall minimum distance */ 
    dst_row = 1.e35;
    for (i=0; i < n_idx; i++) if ((idx[i] < dist->size - 1) && (dst_by_row[idx[i]] < dst_row)) { 
      dst_row = dst_by_row[idx[i]]; min_row = i; min_col = min_by_row[idx[i]]; 
    }
    if (dst_row < 1.e-35) dst_row = 1.e-35;

    /* choose pair with smallest distance */
    i = min_row;  j = min_col; /* relative to idx[] */
    /* tree node creation */
    idx_i      = idxtree[i];
    idx_j      = idxtree[j];
    idxtree[i] = parent; /* new node we are about to create (available to sampling on next iteration) */
    idxtree[j] = idxtree[--n_idx]; /* avoid replacement (should come last since we may have i == n_idx-1 */
    create_parent_node_from_children (tree, parent, idx_i, idx_j);
    parent++; /* parent on next iteration */

    /* calculate branch lengths (UPGMA distances refer to tip, that's why we keep node UPGMA dists on  height[]) */
    gs1 = dst_row/2. - height[idx[i]]; 
    gs2 = dst_row/2. - height[idx[j]]; 
    if (gs1 < 1.e-35) gs1 = 1.e-35;
    if (gs2 < 1.e-35) gs2 = 1.e-35;
    tree->blength[idx_i] = gs1; /* UPGMA distance */ 
    tree->blength[idx_j] = gs2; /* UPGMA distance */
    height[idx[i]] = dst_row/2.;

    /* idx_j nd idx_j have indexes of actual matrix elements */
    idx_i     = idx[i]; /* unlike index of tree nodes, we don't change idx[i] which will hold dists for internal node */
    idx_j     = idx[j]; 
    idx[j]    = idx[n_idx]; /* avoid replacement */

    /* update distance matrix */
    dst_by_row[idx_i] = 1.e35; /* will need update */
    gs1 = (gsize[idx_i] + gsize[idx_j]);
    for (i=0; i < n_idx; i++) { 
      /* calculates distances to new node */
      if (single_linkage) {  /* a.k.a nearest neighbor clustering. Distance to new node is minimum between elements */
        if (idx[i] < idx_j) new_dist = dist->d[idx[i]][idx_j]; /* upper diagonal (d[row][col] --> row < col) */
        else                new_dist = dist->d[idx_j][idx[i]];
        if (idx[i] < idx_i) { row = idx[i]; col = idx_i; idx_col = min_row; }
        else                { col = idx[i]; row = idx_i; idx_col = i; }
        if ((row < col) && (new_dist < dist->d[row][col])) dist->d[row][col] = new_dist; /* skip row==col; if new > dist, then dist=dist (=MIN()) */
      }
      else { /* UPGMA: distance to new node is average between elements */
        if (idx[i] < idx_j) new_dist = gsize[idx_j] * dist->d[idx[i]][idx_j];
        else                new_dist = gsize[idx_j] * dist->d[idx_j][idx[i]];
        if (idx[i] < idx_i) { row = idx[i]; col = idx_i; idx_col = min_row; }
        else                { col = idx[i]; row = idx_i; idx_col = i; }
        if (row < col) dist->d[row][col] = (new_dist + (gsize[idx_i] * dist->d[row][col]))/gs1; /* UPGMA distance -- skip row==col */
      }

      /* check if any new minimum is found */  
      if ((dist->d[row][col] < dst_by_row[row])) { dst_by_row[row] = dist->d[row][col]; min_by_row[row] = idx_col; }

      /* rows whose minimum value was min_row need to be updated */
      if ((idx[i] < (dist->size - 1)) && ((min_by_row[idx[i]] == min_row) || (min_by_row[idx[i]] == min_col) || (min_by_row[idx[i]] >= n_idx))) {
        dst_by_row[idx[i]] = 1.e35;
        for (row = 0; row < n_idx; row++) /* "row" is just a recycled var (like i,j,k,l); no special meaning */
          if (((col=idx[row]) > idx[i]) && (dist->d[idx[i]][col] < dst_by_row[idx[i]])) { 
            dst_by_row[idx[i]] = dist->d[idx[i]][col]; min_by_row[idx[i]] = row; /* min_by_row has only row < col */
          }
      }
    } 

    gsize[idx_i] += gsize[idx_j];/* update group size of new internal node */

  } // while (n_idx > 2)
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idxtree[0], idxtree[1]);
  tree->root = tree->nodelist[parent];

  if (idx[0] < idx[1]) dst_row = dist->d[idx[0]][idx[1]];
  else                 dst_row = dist->d[idx[1]][idx[0]];
  tree->blength[idxtree[0]] = dst_row/2. - height[idx[0]]; /* UPGMA distance */ 
  tree->blength[idxtree[1]] = dst_row/2. - height[idx[1]]; /* UPGMA distance */

  update_topology_sisters   (tree);
  update_topology_traversal (tree);

  if (gsize)   free (gsize);
  if (height)  free (height);
  if (dst_by_row) free (dst_by_row);
}

void
bionj_from_distance_matrix (topology tree, distance_matrix dist) 
{ /* always use upper diagonal of distance_matrix(that is, only i < j in d[i][j]) */
  int i, j, parent = tree->nleaves, n_idx = tree->nleaves, i1, i2, b1, b2, // b1, b2 are best, b1 < b2
      *idx = tree->index,                            /* indexes in UPGMA */
      *idxtree = tree->index + tree->nleaves;        /* indexes in tree (since have values > nleaves) */
  double **delta, Q_min, Q_ij, var_1_2, diff_1_2, blen_1, blen_2, lambda;

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  delta = (double **) biomcmc_malloc (n_idx * sizeof (double*)); /* delta matrix with dists, variances and sums of distances */
  for (i=0; i < n_idx; i++) delta[i] = (double *) biomcmc_malloc (n_idx * sizeof (double));
  /* delta has distances in upper diag, var in lower diag and sum in diagonal (opposite of original BIONJ C program!) */
  for (i=0; i < n_idx; i++) for (j=i+1; j<n_idx; j++) delta[i][j] = delta[j][i] = dist->d[i][j]; // only upper diagonal of dist is used 
  for (i=0; i < n_idx; i++) delta[i][i] = 0.; 

  for (i=0; i < n_idx; i++) { 
    idx[i]     = i; /* index to actual vector element for UPGMA distance matrix */
    idxtree[i] = i; /* index to actual vector element for tree nodes */
  }

  while (n_idx > 2) { /* choose two nodes to be connected */
    /* update sums: diagonal values of delta will be the sum of distances */
    for (i=0; i < n_idx; i++) {
      delta[idx[i]][idx[i]] = 0;
      for (j=0; j < n_idx; j++) if (j!=i) { // idx(i) < idx(j) for dissimilarities
        if (idx[i] < idx[j]) delta[idx[i]][idx[i]] += delta[idx[i]][idx[j]]; 
        else                 delta[idx[i]][idx[i]] += delta[idx[j]][idx[i]];
      }
    }
    /* find pair that minimises agglomerative criterion -- matrix Q_ij */ 
    Q_min = 1.e64;
    for (i=0; i < n_idx; i++) for (j=0; j < i; j++) {
      if (idx[i] < idx[j]) { i1 = i; i2 = j; } // idx[i1] < idx[i2] always
      else                 { i1 = j; i2 = i; }
      Q_ij = (double)(n_idx - 2) * delta[idx[i1]][idx[i2]] - delta[idx[i1]][idx[i1]] - delta[idx[i2]][idx[i2]];
      if (Q_ij < Q_min - 1.e-8) { Q_min = Q_ij; b1 = i1; b2 = i2; }
    }
    //for (i=0; i < n_idx; i++) {
    //  for (j=0; j < n_idx; j++) { printf ("%9.8lf ", delta[idx[i]][idx[j]]); } printf (" <-- \n");
    //}
    //printf ("chosen: %d %d\n", b1, b2);
    diff_1_2 = (delta[idx[b1]][idx[b1]] - delta[idx[b2]][idx[b2]])/(double)(n_idx-2);
    blen_1 = 0.5 * (delta[idx[b1]][idx[b2]] + diff_1_2);
    blen_2 = 0.5 * (delta[idx[b1]][idx[b2]] - diff_1_2);
    /* calculate lambda */
    var_1_2 = delta[idx[b2]][idx[b1]];  // variance between b1 and b2
    if(var_1_2 < 1.e-12) lambda=0.5; // delta[b2][b1] is var between b1 and b2
    else {
      lambda = 0.;
      for (i=0; i< n_idx; i++) if(b1 != i && b2 != i) {
        if (idx[i] < idx[b1]) lambda += delta[idx[b1]][idx[i]]; // lambda += (var(b1,i) - var(b2,i)
        else                  lambda += delta[idx[i]][idx[b1]];
        if (idx[i] < idx[b2]) lambda -= delta[idx[b2]][idx[i]];
        else                  lambda -= delta[idx[i]][idx[b2]];
      }
      lambda = 0.5 + lambda/(2.*(double)(n_idx-2)* delta[idx[b2]][idx[b1]]);
    }
    if(lambda > 1.0) lambda = 1.0;
    if(lambda < 0.0) lambda = 0.0;
    /* update distances and variances of b1, which will be new node (b2 will be replaced) */
    for (i=0; i< n_idx; i++) if(b1 != i && b2 != i) {
      if (idx[b1] < idx[i]) {i1 = b1; i2 = i;}
      else                  {i2 = b1; i1 = i;} // idx[i1] < idx[i2] always
      /* Distance update --> i<j in delta[i][j] */
      delta[idx[i1]][idx[i2]] = lambda * (delta[idx[i1]][idx[i2]] - blen_1);
      if (idx[b2] < idx[i]) delta[idx[i1]][idx[i2]] += (1. - lambda) * (delta[idx[b2]][idx[i]] - blen_2); // distance(b2,i)
      else                  delta[idx[i1]][idx[i2]] += (1. - lambda) * (delta[idx[i]][idx[b2]] - blen_2);
      /* Variance update --> i > j in delta[i][j] */
      delta[idx[i2]][idx[i1]] =  lambda * (delta[idx[i2]][idx[i1]] - (1.-lambda) * var_1_2);
      if (idx[b2] < idx[i]) delta[idx[i2]][idx[i1]] += (1. - lambda) * (delta[idx[i]][idx[b2]]); // variance(b2,i)
      else                  delta[idx[i2]][idx[i1]] += (1. - lambda) * (delta[idx[b2]][idx[i]]); // variance(b2,i)
    }
    // do we need to make sure that b1 < b2 (not only idx[b1]<idx[b2]) otherwise we may be updating last index n_idx-1 which will  be  neglected ??
    /* tree node creation */
    create_parent_node_from_children (tree, parent, idxtree[b1], idxtree[b2]);
    tree->blength[idxtree[b1]] = blen_1;
    tree->blength[idxtree[b2]] = blen_2;
    idxtree[b1] = parent; /* new node we just created (available to sampling on next iteration) */
    idxtree[b2] = idxtree[--n_idx]; /* avoid replacement (should come last since we may have i == n_idx-1 */
    parent++; /* parent on next iteration */
    idx[b2]    = idx[n_idx]; /* avoid replacement */

  } // while (n_idx > 2)

  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idxtree[0], idxtree[1]);
  tree->root = tree->nodelist[parent];

  if (idx[0] < idx[1]) tree->blength[idxtree[0]] = tree->blength[idxtree[1]] = delta[idx[0]][idx[1]];
  else                 tree->blength[idxtree[0]] = tree->blength[idxtree[1]] = delta[idx[1]][idx[0]];

  update_topology_sisters   (tree);
  update_topology_traversal (tree);
  if (delta) {
    for (i = n_idx - 1; i >= 0; i--) if (delta[i]) free (delta[i]);
    free (delta);
  }
}
