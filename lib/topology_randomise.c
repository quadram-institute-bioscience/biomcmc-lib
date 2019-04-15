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

#include "topology_randomise.h"

/*! \brief create internal node with given children (children coalesce into parent node) */
void create_parent_node_from_children (topology tree, int parent, int lchild, int rchild);
/*! \brief recursive function that applies spr (almost always an nni) on a subtree */
bool topology_apply_shortspr_weighted_subtree (topology tree, topol_node lca, double *prob, double scale, bool update_done);

void
randomise_topology (topology tree)
{ 
  int i, lchild, rchild, parent = tree->nleaves, *idx = tree->index, n_idx = tree->nleaves;

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  for (i=0; i < n_idx; i++) idx[i] = i;

  while (n_idx > 2) { /* draw two nodes to be connected (like star-topology refinement) */
    i = biomcmc_rng_unif_int (n_idx);
    rchild = idx[i];
    idx[i] = idx[--n_idx]; /* avoid replacement */
    i = biomcmc_rng_unif_int (n_idx);
    lchild = idx[i];
    idx[i] = parent; /* new node we are about to create (available to sampling on next iteration) */
    create_parent_node_from_children (tree, parent, lchild, rchild);
    parent++; 
  }
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idx[0], idx[1]);
  tree->root = tree->nodelist[parent];
  tree->root->sister = tree->root;
  tree->root->up = NULL;

  update_topology_sisters   (tree);
  topology_apply_spr_unrooted (tree, false); // for some reason the above function does not sample uniformly from unrooted
}

void
quasi_randomise_topology (topology tree, int sample_type)
{ 
  int i, lchild, rchild, parent = tree->nleaves, swap;
  int *taxa_idx = tree->index, 
      *r_idx    = tree->index + tree->nleaves, 
      *l_idx    = tree->index + 2 * tree->nleaves,
      *idx      = tree->index + 3 * tree->nleaves; /* vectors */

  if ((!sample_type) || (!tree->quasirandom)) { sample_type = 14; tree->quasirandom = true; } /* starts a new random tree and prepare vectors */

  /* despite the weird numbering (8,2,4,1), pls do not change the order: first initialize taxa_idx, only then shuffle, etc.
   * this numbering is related to the frequency with which we could do these moves, if user calls using 
   * sample_type=1,2,3,4,5... */
  if (sample_type & 8) for (i = 0; i < tree->nleaves; i++) taxa_idx[i] = i; /* initialization */
  if (sample_type & 2) { /* change order of leaves (does not change shape) */
    for (i = tree->nleaves - 1; i > 0; i--) { /* Knuth shuffle: notice that zero is excluded from loop, and unif(0,i) is with i included */
      lchild = biomcmc_rng_unif_int (i + 1); /* unif(0,i-1) -- that is, i excluded -- would be Sattolo's algorithm */
      swap = taxa_idx[lchild]; taxa_idx[lchild] = taxa_idx[i]; taxa_idx[i] = swap;
    }
  }
  if (sample_type & 4) { /* sample a new shape (does not change leaves) */ 
    for (i = 0; i < tree->nleaves - 2; i++) {
      r_idx[i] = biomcmc_rng_unif_int (tree->nleaves - i);
      l_idx[i] = biomcmc_rng_unif_int (tree->nleaves - i - 1);
    }
  }
  if (sample_type & 1) { /* deterministic cycle over choices (may change shape) */
    for (i = 0; i < tree->nleaves - 2; i++) {
      r_idx[i]--; /* should always be between 0 and nleaves - i - 1 */
      if (r_idx[i] < 0) r_idx[i] = tree->nleaves - i - 1;
      l_idx[i]--; /* should always be between 0 and nleaves - i - 2 */
      if (l_idx[i] < 0) l_idx[i] = tree->nleaves - i - 2;
    }
  }

  /* for all sample types we must copy the indexes to the "disposable" vector */
  for (i = 0; i < tree->nleaves; i++) idx[i] = taxa_idx[i]; 

  for (i = 0; i < tree->nleaves - 2; i++) {/* draw two nodes to be connected (like star-topology refinement) */
    rchild          = idx[ r_idx[i] ];
    idx[ r_idx[i] ] = idx[tree->nleaves - i - 1]; /* avoid replacement */
    lchild          = idx[ l_idx[i] ];
    idx[ l_idx[i] ] = parent; /* new node we are about to create (available to sampling on next iteration) */
    create_parent_node_from_children (tree, parent, lchild, rchild);
    parent++; 
  }
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idx[0], idx[1]);
  tree->root = tree->nodelist[parent];
  tree->root->sister = tree->root;
  tree->root->up = NULL;

  update_topology_sisters   (tree);
  update_topology_traversal (tree);
}

void
create_parent_node_from_children (topology tree, int parent, int lchild, int rchild)
{
  tree->nodelist[parent]->left  = tree->nodelist[lchild];
  tree->nodelist[parent]->right = tree->nodelist[rchild];
  tree->nodelist[lchild]->up = tree->nodelist[parent];
  tree->nodelist[rchild]->up = tree->nodelist[parent];
  tree->nodelist[rchild]->sister = tree->nodelist[lchild];
  tree->nodelist[lchild]->sister = tree->nodelist[rchild];
}

void
topology_apply_rerooting (topology tree, bool update_done)
{
  int i, n1 = 0, n_valid = 0, n_invalid = 0,
      *valid = tree->index, *invalid = tree->index + tree->nnodes; /* index has size 4*nleaves = 2*nnodes + 2 */

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  invalid[n_invalid++] = tree->root->id; /* root and its children are the only forbidden regraft nodes */ 
  invalid[n_invalid++] = tree->root->left->id;
  invalid[n_invalid++] = tree->root->right->id;
  /* create valid[] vector of indexes by exclusion of invalid[] */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  for (i = 0; i < tree->nnodes; i++) { /* every node "i" is OR valid OR invalid */
    if ((n1 < n_invalid) && (i == invalid[n1])) n1++; /* node "i" is invalid: skip it (that's why we qsorted) */
    else valid[n_valid++] = i;                        /* node "i" is valid */
  }

  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw regraft node */

  apply_spr_at_nodes_LCAprune (tree, tree->root, tree->nodelist[n1], update_done);
  update_topology_traversal (tree);
}

void
topology_apply_shortspr (topology tree, bool update_done)
{
  bool success = false;
  double scale = 1./tree->nleaves;
  int i;
  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the sequence was destroyed */
  tree->quasirandom = false; /* IOW, quasi_randomise() will have to start from scratch */
  if (!tree->traversal_updated) update_topology_traversal (tree);

  for (i = 0; (!success) && (i < 4); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, NULL, scale, update_done);
    scale *= 2.;
  }
  if (!success) topology_apply_shortspr_weighted_subtree (tree, tree->root, NULL, 1., update_done);
  update_topology_traversal (tree);
}

void
topology_apply_shortspr_weighted (topology tree, double *prob, bool update_done)
{
  bool success = false;
  double scale = 1., *localprob = NULL;
  int i;

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the sequence was destroyed */
  tree->quasirandom = false; /* IOW, quasi_randomise() will have to start from scratch */
  if (!tree->traversal_updated) update_topology_traversal (tree);

  for (i = 0; (!success) && (i < 8); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, prob, scale, update_done);
    scale *= 2.;
  }
  if (success) { update_topology_traversal (tree); return; }
  /* something is wrong, prob values too small maybe. */
  localprob = (double*) biomcmc_malloc (tree->nnodes * sizeof(double));
  for (i = 0; i < tree->nnodes; i++) localprob[i] = prob[i] + 1.;
  scale = 1./(double) (tree->nleaves);
  for (i = 0; (!success) && (i < 8); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, localprob, scale, update_done);
    scale *= 2.;
  }
  update_topology_traversal (tree);
  if (localprob) free (localprob);
}

bool
topology_apply_shortspr_weighted_subtree (topology tree, topol_node lca, double *prob, double scale, bool update_done)
{  /* theoretically it works even if lca is root, but may leave the (unrooted version of the) tree unchanged */ 
  bool success = false;
  double prob_occurence, p_l, p_r;

  if (!lca->left->internal) return false; /* if left is a leaf then right is also a leaf (since no swaps were applied to lca yet */
  else                      success |= topology_apply_shortspr_weighted_subtree (tree, lca->left, prob, scale, update_done);

  if (lca->right->internal) success |= topology_apply_shortspr_weighted_subtree (tree, lca->right, prob, scale, update_done);
  else {  /* three leaves => only 3 possible topols => only 2 possible swaps: ((A,B),C) -> swap to (A,(B,C)) OR (B,(A,C)) */
    /* 1) the 2 possible swaps can be done by chosing one of the two leaves at subtree as regraft: A or B in ((A,B),C)
     * 2) although we called this function on children, we can still trust left subtree is larger than right */ 
    if (prob) p_l = scale * prob[lca->left->id]; 
    else      p_l = scale; /* if prob == NULL then we just use scale as common prob of swap */
    if (biomcmc_rng_unif_pos32 () < p_l) { /* branch swap will occur */
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
      return true;
    }
    return success; /* failed on this node, but may have succeeded before */
  }
  /* otherwise both children are internal */
  prob_occurence = biomcmc_rng_unif_pos32 ();
  if (prob) {
    p_l = scale * prob[lca->left->id];
    p_r = scale * prob[lca->right->id];
  }
  else p_l = p_r = scale;
  /* only left = L*(1-R) ; only right = R*(1-L) ; both = L*R ; none = (1-L)*(1-R) 
   * like in the prob line (stick break reprsentation):
   * |--------------|------*-----|-------------|--------| 
   * | only left    |    both    |  only right |  none  |        (1)
   * |---------------------*-------------------|         
   * |   left or both      |   right or both   |                 (2)
   * (note that "both" is split half between right and left) */
  if (prob_occurence < p_l + p_r - p_l * p_r) { /* at least one event happens [L(1-R)+R(1-L)+L*R] */
    topol_node np, nr;
    topol_node prunenode[2], regraftnode[3];
    if (prob_occurence < p_l - (p_l * p_r / 2.)) { /* only left subtree OR both with half prob [L(1-R)+L*R/2] -- marked (2) above */
      regraftnode[0] = lca->left; /* after SPR will be at right side */
      regraftnode[1] = lca->left->left;
      regraftnode[2] = lca->left->right;
      prunenode[0] = lca->right->left;
      prunenode[1] = lca->right->right;
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
    }
    else { /* only right subtree or both with other half prob [R(1-L) + L*R/2] -- marked (2) above  */
      regraftnode[0] = lca->right; /* after SPR will be at left side */
      regraftnode[1] = lca->right->left;
      regraftnode[2] = lca->right->right;
      prunenode[0] = lca->left->left;
      prunenode[1] = lca->left->right;
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->right->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->right->right, update_done);
    }
    if ((prob_occurence > (p_l - (p_l * p_r))) && (prob_occurence < p_l)) { /* both events happen -- marked (1) above */
      np =   prunenode[ biomcmc_rng_unif_int (2)]; 
      nr = regraftnode[ biomcmc_rng_unif_int (3)]; 
      apply_spr_at_nodes_notLCAprune (tree, np, nr, update_done); 
    }
    return true;
  }
  return success; /* failed on this node, but may have succeeded before */
}

void
topology_apply_spr_on_subtree (topology tree, topol_node lca, bool update_done)
{ 
  int i, n1, n2, *valid = tree->index, *invalid, *regraft, n_valid = 0, n_invalid = 0, n_regraft = 0;
  topol_node first_child = lca;

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  if (lca == tree->root) biomcmc_error ("root node is not eligible for SPR move (maybe root->left or root->right?)");
  if (!tree->traversal_updated) update_topology_traversal (tree);
  if (lca->split->n_ones < 3) return; /* it should not complain if subtree is a "leaf" or has only two children */
 
  if (lca->split->n_ones == 3) { /* three leaves => only 3 possible topols => only 2 possible swaps */
    /* 1) the 2 possible swaps can be done by chosing one of the two leaves at subtree as regraft: A or B in ((A,B),C)
     * 2) lca->left is internal (by design of traversal, see bipartition_is_larger() function); */
    if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
    else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
    return;
  }

  /* find number of nodes below subtree (using postorder info) */
  while (first_child->internal) first_child = first_child->left; /* walk down the subtree (left is first, in postorder) */ 

  /* all nodes below lca node (inclusive) are eligible prune nodes */
  for (i = first_child->up->mid[0]; i <= lca->mid[0]; i++) { 
    valid[n_valid++] = tree->postorder[i]->id; /* postorder[] stores only internal nodes */
    if (!tree->postorder[i]->left->internal)  valid[n_valid++] = tree->postorder[i]->left->id;
    if (!tree->postorder[i]->right->internal) valid[n_valid++] = tree->postorder[i]->right->id;
  }
#ifdef BIOMCMC_DEBUG
  if (n_valid != (2*lca->split->n_ones - 1)) biomcmc_error ("%d nodes eligible in subtree with %d leaves (SPR)",n_valid, lca->split->n_ones);
#endif
  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  /* n_valid is always smaller than nnodes so we're safe (tree->index has only 4*nleaves = 2*nnodes + 2 elements */
  invalid = tree->index + n_valid; /* invalid[] vector points to after last element of valid[] */

  /* remove prune node's possible four immediate neighbors (left, right, up and sister), lca included */
  invalid[n_invalid++] = n1; /* trivial restriction (regraft != prune) */
  if (tree->nodelist[n1]->internal) { 
    invalid[n_invalid++] = tree->nodelist[n1]->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->right->id;
  }
  if (tree->nodelist[n1] != lca) {
    invalid[n_invalid++] = tree->nodelist[n1]->up->id;
    invalid[n_invalid++] = tree->nodelist[n1]->sister->id;
  }

  regraft = invalid + n_invalid; /* regraft[] starts after last element of invalid[] */

  /* create regraft[] vector of indexes by exclusion of invalid[] from valid[] */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  qsort (  valid,   n_valid, sizeof (int), compare_int_increasing);
  n2 = 0;
  for (i = 0; i < n_valid; i++) { /* every node is OR valid OR invalid */
    if ((n2 < n_invalid) && (valid[i] == invalid[n2])) n2++; /* skip invalid */
    else regraft[n_regraft++] = valid[i];
  }

  n2 = regraft[ biomcmc_rng_unif_int (n_regraft) ]; /* regraft node (here valid[] already has idx in postorder) */

  /* this wrapper function will decide if prune node is LCA or not of regraft */
  apply_spr_at_nodes (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
  return;
}

void
topology_apply_spr (topology tree, bool update_done)
{ 
  uint32_t n_left, n_right;

  if (cant_apply_swap (tree)) return; /* checks if any subtree has more than two leaves */
  n_left  = tree->root->left->split->n_ones;
  n_right = tree->root->right->split->n_ones;

  if     (n_right < 3) topology_apply_spr_on_subtree (tree, tree->root->left, update_done);
  else if (n_left < 3) topology_apply_spr_on_subtree (tree, tree->root->right, update_done);
  else if (biomcmc_rng_unif_int (n_right + n_left) < n_left) topology_apply_spr_on_subtree (tree, tree->root->left, update_done);
  else                                                       topology_apply_spr_on_subtree (tree, tree->root->right, update_done);
}

void
topology_apply_spr_unrooted (topology tree, bool update_done)
{ /* neglects root node and applies SPR that changes eq. unrooted topology */
  int i, n1, n2, *valid = tree->index, *invalid, *regraft, n_valid = 0, n_invalid = 0, n_regraft = 0;

  if (tree->nleaves < 4) return; /* There is only one unrooted triplet */
  if (!tree->traversal_updated) update_topology_traversal (tree);
  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;
  if (tree->nleaves == 4) { /* special case, furthermore can be simplified */
    if (!tree->root->right->internal) { /* (((a,b),c),d) --> a or b regrafted to d (OR,equiv, regrafted to c) */
      valid[n_valid++] = tree->root->left->left->left->id; 
      valid[n_valid++] = tree->root->left->left->right->id;
      n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */
      apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->root->right, update_done);
      update_topology_traversal (tree);
      return;
    }
    else { /* ((a,b),(c,d)) --> a or b regrafted to c (again, since it's unrooted, it's equiv to regrafting to d) */
      valid[n_valid++] = tree->root->left->left->id; 
      valid[n_valid++] = tree->root->left->right->id;
      n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */
      apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->root->right->left, update_done);
      update_topology_traversal (tree);
      return;
    }
  }
  /* all nodes except root are eligible prune nodes */
  for (i = 0; i < tree->nnodes; i++) if (i != tree->root->id) valid[n_valid++] = i;
  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  /* tree->index has (4*nleaves) = (2*nnodes + 2) elements so we're safe since n_valid <= nnodes */
  invalid = tree->index + n_valid; /* invalid[] vector points to after last element of valid[] */

  /* remove prune node's possible four immediate neighbors (left, right, up and sister) */
  invalid[n_invalid++] = n1; /* trivial restriction (regraft != prune) */
  invalid[n_invalid++] = tree->nodelist[n1]->up->id; /* since prune cannot be root node */
  invalid[n_invalid++] = tree->nodelist[n1]->sister->id;
  if (tree->nodelist[n1]->internal) { 
    invalid[n_invalid++] = tree->nodelist[n1]->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->right->id;
  }
  if ((tree->nodelist[n1]->up == tree->root) && (tree->nodelist[n1]->sister->internal)) { /* equiv. to rerooting */ 
    invalid[n_invalid++] = tree->nodelist[n1]->sister->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->sister->right->id;
  }
  else if (tree->nodelist[n1]->up->up == tree->root) { /* equiv. to rerooting */ 
    invalid[n_invalid++] = tree->nodelist[n1]->up->sister->id;
  }

  regraft = invalid + n_invalid; /* regraft[] starts after last element of invalid[] */

  /* create regraft[] vector of indexes by exclusion of invalid[] from all nodes (including root now) */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  n2 = 0;
  for (i = 0; i < tree->nnodes; i++) { /* every node is OR valid OR invalid */
    if ((n2 < n_invalid) && (invalid[n2] == i)) n2++; /* skip invalid */
    else regraft[n_regraft++] = i;
  }

  n2 = regraft[ biomcmc_rng_unif_int (n_regraft) ]; /* regraft node (here valid[] already has idx in postorder) */

  /* this wrapper function will decide if prune node is LCA or not of regraft */
  apply_spr_at_nodes (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
  return;
}

void
topology_apply_nni (topology tree, bool update_done)
{
  int i, n1, n2, *valid = tree->index, *invalid = tree->index + tree->nnodes, n_valid = 0, n_invalid = 0, lca_idx = 0;
  /* please read comment on function topology_apply_spr() about the two vectors and the need to order them. I always
   * find this solution ugly, but at least I think is correct. The bugs didn't appear here with the tempting
   * no-replacement algorithm, but in apply_spr() so I suspect it was only with the regraft sampling (leomrtns). 
   *
   * UPDATE 2010.07.27: this function seems to be buggy (does not respect rooting) */

  /* tree->index is also used by quasi_randomise_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  invalid[n_invalid++] = tree->root->id; /* root is forbidden prune node (otherwise we have rerooting) */
  /* if a child of root is pruned, the regraft must be on subtree below (at least two nodes apart from prune) */
  if ((!tree->root->left->internal) || (!tree->root->left->left->internal && !tree->root->left->right->internal)) 
    invalid[n_invalid++] = tree->root->left->id;
  if ((!tree->root->right->internal) || (!tree->root->right->left->internal && !tree->root->right->right->internal)) 
    invalid[n_invalid++] = tree->root->right->id;

  /* create valid[] vector of indexes by exclusion of invalid[] */ 
  if (n_invalid > 1) qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  n1 = 0;
  for (i = 0; i < tree->nnodes; i++) { /* every node is OR valid OR invalid */
    if ((n1 < n_invalid) && (i == invalid[n1])) n1++; /* node "i" is invalid, skip (that's why we qsorted) */
    else valid[n_valid++] = i;                        /* node "i" is valid */
  }

  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  n_valid = 0;
  /* NNI: at most four possible regraft nodes of LCA type */
  if (tree->nodelist[n1]->internal) {
    if (tree->nodelist[n1]->left->internal) {
      valid[n_valid++] = tree->nodelist[n1]->left->left->id;
      valid[n_valid++] = tree->nodelist[n1]->left->right->id;
    }
    if (tree->nodelist[n1]->right->internal) {
      valid[n_valid++] = tree->nodelist[n1]->right->left->id;
      valid[n_valid++] = tree->nodelist[n1]->right->right->id;
    }
    lca_idx = n_valid; /* up to this index, spr is of LCA type */
  }

  if (tree->nodelist[n1]->up != tree->root) { /* root's children are only allowed to LCA-type moves */ 
    valid[n_valid++] = tree->nodelist[n1]->up->sister->id;
    if (tree->nodelist[n1]->up->up != tree->root) 
    valid[n_valid++] = tree->nodelist[n1]->up->up->id;
    if (tree->nodelist[n1]->sister->internal) {
    valid[n_valid++] = tree->nodelist[n1]->sister->left->id;
    valid[n_valid++] = tree->nodelist[n1]->sister->right->id;
    }
  }

  i = biomcmc_rng_unif_int (n_valid); /* draw regraft node */
  n2 = valid[i];

  if (i < lca_idx) apply_spr_at_nodes_LCAprune (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  else          apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
}

bool
cant_apply_swap (topology tree)
{
  if (!tree->traversal_updated) update_topology_traversal (tree);
  if ((tree->root->left->split->n_ones < 3) && (tree->root->left->split->n_ones < 3)) return true;
  return false;
}

