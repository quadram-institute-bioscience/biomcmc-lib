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

/*! \brief rescales original branches and stores distance from root into #fromroot[], using up to 6 distinct rescalings */
double rescale_rooted_distances_for_patristic_distances (topology tree, double *fromroot, int mode, double tolerance);
/*! \brief update internal bipartitions and reorder siblings by heavy child on left */
uint32_t update_subtree_bipartitions (topol_node this);
/*! \brief tree traversal with preorder and postorder node tracking */
uint32_t update_subtree_traversal (topology tree, topol_node this, int *postcount, int *undonecount);
/*! \brief Auxiliary function to topology_to_string_by_id() and toplogy_to_string_create_name(). */
void topology_subtree_to_string_by_id (char *str, const topol_node node, double *blen, bool create_name);
/*! \brief Auxiliary function to topology_to_string_by_name(). */
void topology_subtree_to_string_by_name (char *str, const topol_node node, const char **taxlabel, double *blen);
/*! \brief Flag nodes descending from this topol_node_struct as upper part undone, in pre-order. */
void undo_udone (topol_node this);
/*! \brief Flag nodes ancestral to this topol_node_struct as lower part undone. */
void undo_ddone (topol_node this);

topology
new_topology (int nleaves) 
{
  topology tree;
  int i;
  size_t sizeof_node = sizeof (struct topol_node_struct);

  tree = (topology) biomcmc_malloc (sizeof (struct topology_struct));
  tree->nleaves = nleaves;
  tree->nnodes  = 2*nleaves - 1;
  tree->n_undone = nleaves - 1;
  tree->hashID1 = tree->hashID2 = 0;
  tree->traversal_updated = false;
  tree->ref_counter = 1;
  tree->taxlabel = NULL; /* Memory allocated by topology_space:: or other place */
  tree->quasirandom = false;

  /* actual allocation */
  tree->nodelist  = (topol_node*) biomcmc_malloc (tree->nnodes * sizeof (topol_node));
  /* pointers only */
  tree->postorder = (topol_node*) biomcmc_malloc ((tree->nleaves - 1) * sizeof (topol_node));
  tree->undone    = (topol_node*) biomcmc_malloc ((tree->nleaves - 1) * sizeof (topol_node));
  tree->blength = (double*) biomcmc_malloc (tree->nnodes * sizeof(double));
  for (i = 0; i < tree->nnodes; i++) tree->blength[i] = 1.; 


  /* sandbox vector (thread-safe "global" variable) used 
   * - when drawing nodes in random spr (two vectors of size nnodes)
   * - to store accepted topology configuration in minisampler or other composite MC proposals */
  tree->index = (int*) biomcmc_malloc (4 * tree->nleaves * sizeof (int));
  for (i=0; i < (4 * tree->nleaves); i++) tree->index[i] = 0;

  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nleaves; i++) { 
    tree->nodelist[i] = (topol_node) biomcmc_malloc (sizeof_node);
    tree->nodelist[i]->u_done = tree->nodelist[i]->internal = false;
    tree->nodelist[i]->d_done = true;
    tree->nodelist[i]->left = tree->nodelist[i]->right = tree->nodelist[i]->up = NULL;
    tree->nodelist[i]->sister = tree->nodelist[i];
    tree->nodelist[i]->mid[0] = tree->nodelist[i]->mid[1] = tree->nodelist[i]->id = i;

    tree->nodelist[i]->split = new_bipartition (tree->nleaves);
    bipartition_set (tree->nodelist[i]->split, i);
  }

  for (i=tree->nleaves; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (topol_node) biomcmc_malloc (sizeof_node);
    tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = true;
    tree->nodelist[i]->left = tree->nodelist[i]->right = tree->nodelist[i]->up = NULL;
    tree->nodelist[i]->sister = tree->nodelist[i]; // root is sister of itself ;)
    tree->nodelist[i]->internal = true;
    tree->nodelist[i]->mid[0] = tree->nodelist[i]->mid[1] = tree->nodelist[i]->id = i;

    tree->nodelist[i]->split = new_bipartition (tree->nleaves);
  }
  tree->root = tree->nodelist[tree->nnodes - 1]; /* arbitrary, but usually correct */

  return tree;
}

void 
del_topology (topology tree) 
{
  if (tree) {
    if (--tree->ref_counter) return; /* free memory only after all references are deleted */
    if (tree->postorder) free (tree->postorder);
    if (tree->undone)    free (tree->undone);
    if (tree->index)     free (tree->index);
    if (tree->blength)   free (tree->blength);
    if (tree->nodelist) {
      int i;
      for (i=tree->nnodes-1; i >=0; i--) if (tree->nodelist[i]) { 
        del_bipartition (tree->nodelist[i]->split);
        free (tree->nodelist[i]);
      }
      free (tree->nodelist);
    }
    del_char_vector (tree->taxlabel);
  
    free (tree);
  }
}

void
copy_topology_from_topology (topology to_tree, topology from_tree)
{
  int i;

  if (!from_tree->traversal_updated) update_topology_traversal (from_tree);

  for (i=0; i < from_tree->nleaves; i++) {
    to_tree->nodelist[i]->up     = to_tree->nodelist[from_tree->nodelist[i]->up->id];
    to_tree->nodelist[i]->sister = to_tree->nodelist[from_tree->nodelist[i]->sister->id];
    to_tree->nodelist[i]->left   = to_tree->nodelist[i]->right = NULL; // redundant, but safe
    // to_tree->nodelist[i]->u_done = from_tree->nodelist[i]->u_done; 
    to_tree->nodelist[i]->d_done = true;
  }

  for (i = from_tree->nleaves; i < from_tree->nnodes; i++) {
    to_tree->nodelist[i]->mid[0] = from_tree->nodelist[i]->mid[0];
    to_tree->nodelist[i]->mid[1] = from_tree->nodelist[i]->mid[1];
    //to_tree->nodelist[i]->u_done = from_tree->nodelist[i]->u_done; 
    to_tree->nodelist[i]->d_done = from_tree->nodelist[i]->d_done; 

    if (from_tree->nodelist[i]->up) to_tree->nodelist[i]->up = to_tree->nodelist[from_tree->nodelist[i]->up->id];
    else {
      to_tree->nodelist[i]->up = NULL;
      to_tree->root = to_tree->nodelist[i];
    }
    to_tree->nodelist[i]->left   = to_tree->nodelist[from_tree->nodelist[i]->left->id];
    to_tree->nodelist[i]->right  = to_tree->nodelist[from_tree->nodelist[i]->right->id];
    to_tree->nodelist[i]->sister = to_tree->nodelist[from_tree->nodelist[i]->sister->id];
  } // for (nnodes)
  
  if (from_tree->blength) { 
    if (!to_tree->blength) to_tree->blength = (double*) biomcmc_malloc (from_tree->nnodes * sizeof (double));
    for (i = 0; i < from_tree->nnodes; i++) to_tree->blength[i] = from_tree->blength[i];
  }
  update_topology_traversal (to_tree);
  if (from_tree->taxlabel) { // in case from_tree is a dummy/temp but to_tree is important
    del_char_vector (to_tree->taxlabel);
    to_tree->taxlabel = from_tree->taxlabel;
    to_tree->taxlabel->ref_counter++;
  }
}

void
debug_topol (topology tree)
{
  int i;
  //for (i=0; i < tree->nleaves; i++)
  //  printf ("%2d.%2d %2d] ", i, tree->nodelist[i]->up->id, tree->nodelist[i]->sister->id);
  //printf ("\n");
  for (i=tree->nleaves; i < tree->nnodes; i++) {
    if (tree->nodelist[i]->up) printf ("%2d.%2d *", i, tree->nodelist[i]->up->id);
    else                       printf ("%2d.root*", i);
    printf ("%2d %2d %2d| ", tree->nodelist[i]->left->id, tree->nodelist[i]->right->id, tree->nodelist[i]->sister->id);
  }
  printf ("   (index . Up * Left Right Sister |\n");
  for (i=0; i < tree->nleaves - 1; i++) 
    printf ("[%2d. %1d] ", tree->postorder[i]->id, tree->postorder[i]->d_done); 
  printf ("   (postorder\n");
  for (i=0; i < tree->n_undone; i++) 
    printf ("[%2d. %1d] ", tree->undone[i]->id, tree->undone[i]->d_done); 
  if (tree->undo_prune) printf ("   (undone {%3d %3d}\n\n", tree->undo_prune->id, tree->undo_regraft->id);
  else printf ("   (undone {}\n\n");
}

void 
update_topology_sisters (topology tree)
{
  int i;

  for (i=0; i < tree->nnodes; i++) {
    if (tree->nodelist[i]->up) {
      if (tree->nodelist[i]->up->left == tree->nodelist[i]) tree->nodelist[i]->sister = tree->nodelist[i]->up->right;
      else                                                  tree->nodelist[i]->sister = tree->nodelist[i]->up->left;
    }
    else tree->nodelist[i]->sister = tree->nodelist[i]; /* root node */
  }
}

void
update_topology_traversal (topology tree)
{
  int i = 0, j = 0;
 
  tree->hashID1 = update_subtree_bipartitions (tree->root);
  tree->hashID2 = update_subtree_traversal (tree, tree->root, &i, &j); 
  tree->n_undone = j;

  /* preorder scan == postorder[nleaves-3 -> 0] since postorder[nleaves-2] == root */
  tree->root->level = 0;
  for (i = tree->nleaves-3; i >= 0; i--) tree->postorder[i]->level = tree->postorder[i]->up->level + 1; /*internal nodes */
  for (i = 0; i < tree->nleaves; i++)    tree->nodelist[i]->level  = tree->nodelist[i]->up->level  + 1;

  tree->traversal_updated = true;
}

uint32_t
update_subtree_bipartitions (topol_node this)
{
  uint32_t hash1 = 0x55555UL, hash2 = 0xefc6dUL; /* each tree has a (ideally) unique hash value */
  if (this->left->internal)  hash1 = update_subtree_bipartitions (this->left);
  else                       hash1 = this->left->id;
  if (this->right->internal) hash2 = update_subtree_bipartitions (this->right);
  else                       hash2 = this->right->id;

  bipartition_OR (this->split, this->left->split, this->right->split, false);
  if (bipartition_is_larger (this->right->split, this->left->split)) {
    // heavy child (more leaves - or leaves with larger ids in case of a tie) at left 
    topol_node tmp = this->left; 
    this->left = this->right; 
    this->right = tmp;
    hash1 = biomcmc_hashint_salted (hash1, /*salt*/ 4);
  }
  else hash2 = biomcmc_hashint_salted (hash2, /*salt*/ 4); // only one of them is hashed (the lighter subtree, or the right leaf for cherries)
  if (hash1 > hash2) return hash1 - hash2 + 1; /* avoid numbers out of range; the +1 is arbitrary here, just to avoid zero below */
  else               return hash2 - hash1 + 1;
}

uint32_t
update_subtree_traversal (topology tree, topol_node this, int *postcount, int *undonecount)
{
  uint32_t hash1 = 0x9d2dUL, hash2 = 0x8cab4UL; /* each tree has a (ideally) unique hash value */
  if (this->left->internal)  hash1 = update_subtree_traversal (tree, this->left, postcount, undonecount);
  else                       hash1 = this->left->id;
  if (this->right->internal) hash2 = update_subtree_traversal (tree, this->right, postcount, undonecount);
  else                       hash2 = this->right->id;

  this->mid[0] = *postcount;
  tree->postorder[(*postcount)++] = this;
  if (!this->d_done) { this->mid[1] = *undonecount; tree->undone[(*undonecount)++] = this; }
  return biomcmc_hashint_mix_salted (hash1, hash2 + (uint32_t) (*postcount), /*salt*/ 1);
}

bool
topology_is_equal (topology t1, topology t2)
{ // this is a simple (lowlevel) function, doesn't check if taxlabels are equivalent when both are present
  int i; 
  if (t1->taxlabel && t2->taxlabel && (t1->taxlabel != t2->taxlabel)) return false;
  if (!t1->traversal_updated) update_topology_traversal (t1);
  if (!t2->traversal_updated) update_topology_traversal (t2);

  if (t1->hashID1 != t2->hashID1) return false;
  if (t1->hashID2 != t2->hashID2) return false;
  for (i=0; i < t1->nleaves-1; i++) // only if by chance both hash values are the same (double collision) 
    if (!bipartition_is_equal (t1->postorder[i]->split, t2->postorder[i]->split)) return false;
  return true;
}

bool
topology_is_equal_unrooted (topology t1, topology t2, bool use_root_later) 
{ /* check if structures are same once root is removed; do not check induced tree (i.e. excluding unique leaves) */
  int i, n = t1->nleaves - 3;
  bipartition *b1, *b2;
  if (t1->taxlabel && t2->taxlabel && (t1->taxlabel != t2->taxlabel)) return false;
  if (t1->nleaves != t2->nleaves) return false;
  if (!t1->traversal_updated) update_topology_traversal (t1);
  if (!t2->traversal_updated) update_topology_traversal (t2);

  b1 = (bipartition*) biomcmc_malloc (2 * n * sizeof (bipartition));
  b2 = b1 + n; // first half from t1 and second half from t2
  for (i=0; i < n; i++) { // excludes root node and one before to avoid redundancy
    b1[i] = t1->postorder[i]->split;
    b2[i] = t2->postorder[i]->split;
    bipartition_flip_to_smaller_set (b1[i]); /* WARNING: this changes original splits, leading to problems if root is then used */
    bipartition_flip_to_smaller_set (b2[i]);
  }
  qsort (b1, n, sizeof (bipartition), compare_bipartitions_increasing);
  qsort (b2, n, sizeof (bipartition), compare_bipartitions_increasing);
  for (i=0; (i < n) && (bipartition_is_equal (b1[i], b2[i])); i++);

  if (use_root_later) for (i=0; i < n; i++) {
    bipartition_OR (t1->postorder[i]->split, t1->postorder[i]->left->split, t1->postorder[i]->right->split, false);
    bipartition_OR (t2->postorder[i]->split, t2->postorder[i]->left->split, t2->postorder[i]->right->split, false);
  }
  if (b1) free (b1);
  if (i == n) return true;
  return false;
}

void
reorder_topology_leaves (topology tree)
{
  int i, *order;
  topol_node *pivot;
  double *blen;

  if (!tree->taxlabel) return; 
  order = (int*) biomcmc_malloc (tree->taxlabel->nstrings * sizeof (int));
  char_vector_reorder_by_size_or_lexicographically (tree->taxlabel, false, order); // false/true -> by size/lexico

  pivot = (topol_node*) biomcmc_malloc (tree->nleaves * sizeof (topol_node));
  blen = (double*) biomcmc_malloc (tree->nnodes *sizeof (double*));
  for (i=0; i < tree->nleaves; i++) {
    pivot[i] = tree->nodelist[i];
    blen[i] = tree->blength[i];
  }
  for (i=0; i < tree->nleaves; i++) {
    tree->nodelist[i] = pivot[ order[i] ];
    tree->nodelist[i]->id = i;
    bipartition_zero (tree->nodelist[i]->split);
    bipartition_set  (tree->nodelist[i]->split, i);
    tree->blength[i] = blen[ order[i] ];
  }
  if (pivot) free (pivot);
  if (order) free (order);
  if (blen) free (blen);
  update_topology_traversal (tree);
}

bool 
node1_is_child_of_node2 (topol_node node1, topol_node node2)
{
  topol_node this = node1;
  while (this->up && (this != node2)) this = this->up;
  if (this == node2) return true;
  return false;
}

char *
topology_to_string_by_id (const topology tree, double *blen) 
{
  char *str;
  /* allocate space for str (overestimate size) */
  int size = tree->nnodes * 3 + tree->nleaves * 8;
  if (blen) size += 16 * tree->nnodes;

  str = (char *) biomcmc_malloc (sizeof (char) * size);
  memset (str, '\0', sizeof (char) * size);
  topology_subtree_to_string_by_id (str, tree->root, blen, false);
  str = strcat (str, ";");
  return str;
}

char *
topology_to_string_create_name (const topology tree, double *blen) 
{
  char *str;
  /* allocate space for str (overestimate size) */
  int size = tree->nnodes * 3 + tree->nleaves * 12;
  if (blen) size += 16 * tree->nnodes;

  str = (char *) biomcmc_malloc (sizeof (char) * size);
  memset (str, '\0', sizeof (char) * size);
  topology_subtree_to_string_by_id (str, tree->root, blen, true);
  str = strcat (str, ";");
  return str;
}

void
topology_subtree_to_string_by_id (char *str, const topol_node node, double *blen, bool create_name)
{
  char s1[64];
  if (node->internal) { /* internal node */
    str = strcat (str, "(");
    topology_subtree_to_string_by_id (str, node->left, blen, create_name);
    str = strcat (str, ",");
    topology_subtree_to_string_by_id (str, node->right, blen, create_name);
    if (blen) { sprintf (s1, "):%.12g", blen[node->id]); str = strcat (str, s1); }
    else str = strcat (str,  ")");
  } else {
    if (create_name) { /* taxa names will be s1, s2 etc. */
      if (blen) sprintf (s1, "T%d:%.12g", node->id+1, blen[node->id]);
      else      sprintf (s1, "T%d", node->id+1);
    } else {
      if (blen) sprintf (s1, "%d:%.12g", node->id+1, blen[node->id]);
      else      sprintf (s1, "%d", node->id+1);
    }
    str = strcat (str, s1); 
  } // else (not internal)
}


char *
topology_to_string_by_name (const topology tree, double *blen)
{
  char *str;
  int size = 1 + tree->nnodes * 3, i;

  if (!tree->taxlabel) return topology_to_string_by_id (tree, blen);

  /* allocate space for str (overestimate size) */
  for (i=0; i < tree->nleaves; i++) size += strlen (tree->taxlabel->string[i]) + 1;
  if (blen) size += 16 * tree->nnodes;

  str = (char *) biomcmc_malloc (sizeof (char) * size);
  memset (str, '\0', sizeof (char) * size);
  topology_subtree_to_string_by_name (str, tree->root, (const char **) tree->taxlabel->string, blen);
  str = strcat (str, ";");
  return str;
}

void
topology_subtree_to_string_by_name (char *str, const topol_node node, const char **taxlabel, double *blen)
{
  char s1[64];
  if (node->internal) { /* internal node */
    str = strcat (str, "(");
    topology_subtree_to_string_by_name (str, node->left, taxlabel, blen);
    str = strcat (str, ",");
    topology_subtree_to_string_by_name (str, node->right, taxlabel, blen);
    if (blen) { sprintf (s1, "):%.12g", blen[node->id]); str = strcat (str, s1); }
    else str = strcat (str,  ")");
  }
  else {
    str = strcat (str, taxlabel[node->id]);
    if (blen) { sprintf (s1, ":%.12g", blen[node->id]); str = strcat (str, s1); }
  }
}

void
graphviz_file_topology (FILE * fout, char *label, const topology tree) 
{
  int i;
  fprintf (fout, "graph G {\n");
  fprintf (fout, "  graph [ size=\"7,9\" page=\"8.5,11\" center=\"\" ]\n");
  fprintf (fout, "  node  [ fontsize = \"8\" width=.08, hight=.08 ]\n");
  fprintf (fout, "  edge  [ fontsize = \"6\" len=1.5 ]\n");
  for (i = 0; i < tree->nnodes; i++)
   {
    if (!tree->nodelist[i]->internal)
      fprintf (fout,
               "  %d\t[ label = \"%d\" width=.16, hight=.16 ];\n",
               tree->nodelist[i]->id, tree->nodelist[i]->id);
    if (tree->nodelist[i]->left)
      fprintf (fout, "  %d -- %d\t[ label = \"%f\" ];\n",
               tree->nodelist[i]->id, tree->nodelist[i]->left->id,
               tree->blength[tree->nodelist[i]->left->id]);
    if (tree->nodelist[i]->right)
      fprintf (fout, "  %d -- %d\t[ label = \"%f\" ];\n",
               tree->nodelist[i]->id,
               tree->nodelist[i]->right->id,
               tree->blength[tree->nodelist[i]->right->id]);
   }
  if (label)
    fprintf (fout, "  label =\"%s\";\n", label);
  fprintf (fout, "\n}\n");
  fflush (fout);
}

void
apply_spr_at_nodes (topology tree, topol_node prune, topol_node regraft, bool update_done)
{
  if (node1_is_child_of_node2 (regraft, prune)) apply_spr_at_nodes_LCAprune (tree, prune, regraft, update_done);
  else                                       apply_spr_at_nodes_notLCAprune (tree, prune, regraft, update_done);
}

void
apply_spr_at_nodes_LCAprune (topology tree, topol_node prune, topol_node regraft, bool update_done)
{
  /* <b> prune is lca</b>:
   * Algorithm equivalent to rerooting, regraft node climbs up till it finds prune.
   *
   * \verbatim
   *
   *                        /prune.up        prune.up\
   * regraft_________ prune/            ==>           \prune___________C
   *          |   |        \                          /        |   |  
   *          A   B         \C                regraft/         A   B
   *
   * \endverbatim
   */
  topol_node r           = regraft, 
             rup         = regraft->up, 
             tmp         = rup->up, 
             newchild    = regraft->up, 
             prunesister = prune->sister;
  bool regraft_is_left = false;

  if (!node1_is_child_of_node2 (regraft, prune)) printf ("gotcha! This is a BUG, not your fault\n"); //DEBUG
  r->up = prune;
  rup->up = prune;
  if (rup->left == r) { /* regraft is left child */
    regraft_is_left = true;
    rup->left = tmp;   
    rup->right->sister = tmp; tmp->sister = rup->right;
  }
  else {
    rup->right = tmp; 
    rup->left->sister = tmp; tmp->sister = rup->left;
  }
  rup->sister = r; r->sister = rup;

  r = rup;
  rup = tmp;
  tmp = tmp->up;
  while (rup != prune) {
    if (rup->left == r) { /* this won't work for last iteration, see fix [1,2] below */ 
      rup->left = tmp; 
      rup->right->sister = tmp; tmp->sister = rup->right;
    }
    else {
      rup->right = tmp;
      rup->left->sister = tmp; tmp->sister = rup->left;
    }

    rup->up = r;
    r = rup;
    rup = tmp;
    tmp = tmp->up;
  }
  /* rup == prune */
  if (prune->left == r) tmp = prune->right; /* sister */
  else tmp = prune->left;
  tmp->up = r;
  if (r->left == prune) { /* [1] fix commented above */
    r->left = tmp;
    r->right->sister = tmp; tmp->sister = r->right;
  }
  else {
    r->right = tmp;
    r->left->sister = tmp; tmp->sister = r->left;
  }
  prune->sister = prunesister; /* [2] recover original sister */

  /* It could be arbitrary, but this way original positions are recovered when SPR is undone (I think) */
  if (regraft_is_left) { prune->left  = regraft; prune->right = newchild; }
  else                 { prune->right = regraft; prune->left  = newchild; }

  /* how to undo this move */
  tree->undo_prune = prune;
  tree->undo_regraft = tmp;
  tree->undo_lca = true;

  /* update u_done and d_done information */
  if (update_done) {
    if (prune->left->internal)  undo_udone (prune->left); 
    if (prune->right->internal) undo_udone (prune->right);
    if (!(tmp->internal)) tmp = tmp->up;
    undo_ddone (tmp);
  }

  /* calling function should update traversal (since it may update d_done before doing that) */
  tree->traversal_updated = false;
}

void
apply_spr_at_nodes_notLCAprune (topology tree, topol_node prune, topol_node regraft, bool update_done)
{
  /*! <b>prune is not lca</b>:
   *  Detach the prune subtree and reinsert it just above the regraft node (regraft node may be root).
   *
   * \verbatim
   *  Prune: 
   *
   *  p.left\              /p.up.up                       p.left\                |p.up.up                  
   *         \prune___p.up/                       ==>            \p_______prune  |                         
   *         /            \                                      /               |                         
   * p.right/              \p.up.left || p.up.right      p.right/                |p.up.left || p.up.right  
   *
   * \endverbatim
   */
  topol_node psister, p = prune;

  /* update u_done and d_done information */
  if (update_done) {
    undo_ddone (prune->up);

    /* must be done by hand (no undo_ddone()) since we must find the LCA */
    if (regraft->up) {
      prune = regraft->up;
      while ((prune->up) && (prune->d_done)) { prune->d_done = false; prune = prune->up; }
    }
    else { /* regraft is root node */
      prune = regraft; 
      prune->d_done = false;
    }

    /* prune == LCA between original prune and regraft nodes */
    if (prune->left->internal)  undo_udone (prune->left);
    if (prune->right->internal) undo_udone (prune->right);
    prune = p;
  }

  psister = prune->sister;
  prune = prune->up;
  if (prune != tree->root) {
    psister->sister = prune->sister; 
    prune->sister->sister = psister;
    psister->up = prune->up;
    if (prune->up->left == prune) prune->up->left = psister;
    else prune->up->right = psister;
  }
  else {
    tree->root = psister->sister = psister; /* root's sister is itself */
    prune->sister->sister = psister;
    psister->up = prune->up;
  }

  /*! \verbatim
   *  Regraft: 
   *
   *  p.left\                |r.up        p.left\               /prune.up (=r.up.up)
   *         \p_______prune  |      ==>          \p_______prune/ 
   *         /               |                   /             \
   * p.right/                |r          p.right/               \r
   *
   * \endverbatim
   */
  if (prune->left == psister) prune->left = regraft;
  else prune->right = regraft;

  p->sister = regraft;
  regraft->sister = p;

  prune->up = regraft->up;
  
  if (regraft != tree->root) {
    if (regraft->up->left == regraft) {
      regraft->up->left = prune;
      regraft->up->right->sister = prune;
      prune->sister = regraft->up->right;
    }
    else {
      regraft->up->right = prune;
      regraft->up->left->sister = prune;
      prune->sister = regraft->up->left;
    }
  }
  else { /* prune is the new root node */
    tree->root = prune->sister = prune;  
  }
  
  regraft->up = prune;

  /* how to undo this move */
  tree->undo_prune = p;
  tree->undo_regraft = psister;
  tree->undo_lca = false;

  /* calling function should update traversal (since it may update d_done before doing that) */
  tree->traversal_updated = false;
}

void
topology_undo_random_move (topology tree, bool update_done)
{
  if (tree->undo_lca) apply_spr_at_nodes_LCAprune (tree, tree->undo_prune, tree->undo_regraft, update_done);
  else             apply_spr_at_nodes_notLCAprune (tree, tree->undo_prune, tree->undo_regraft, update_done);
}

void
undo_udone (topol_node this)
{
  this->u_done = false;
  if (this->left->internal) undo_udone (this->left);
  if (this->right->internal) undo_udone (this->right);
}

void
undo_ddone (topol_node this)
{
  while ((this) && (this->d_done)) { this->d_done = false; this = this->up; }
}

void
clear_topology_flags (topology tree)
{
  int i;
  for (i = tree->nleaves; i < tree->nnodes; i++) tree->nodelist[i]->d_done = tree->nodelist[i]->u_done = true;
}

void
raise_topology_flags (topology tree)
{
  int i;
  for (i = tree->nleaves; i < tree->nnodes; i++) tree->nodelist[i]->d_done = tree->nodelist[i]->u_done = false;
}

void
topology_reset_random_move (topology tree)
{
  topology_undo_random_move (tree, false);
  clear_topology_flags (tree);
}

int
copy_topology_to_intvector_by_postorder (topology tree, int *ivec)
{
  int j, k = 0;
  if (!tree->traversal_updated) update_topology_traversal (tree);

  for (j = 0; j < tree->nleaves; j++) /* leaves first */
    ivec[k++] = tree->nodelist[j]->up->mid[0] + tree->nleaves;
  for (j = 0; j < tree->nleaves - 2; j++) /* position here is the same as mid[0] + nleaves */
    ivec[k++] = tree->postorder[j]->up->mid[0] + tree->nleaves;

  return k;
}

int
copy_intvector_to_topology_by_postorder  (topology tree, int *ivec)
{
  int j;

  for (j = 0; j < tree->nnodes; j++) 
    tree->nodelist[j]->up = tree->nodelist[j]->left = tree->nodelist[j]->right = NULL;
  
  for (j = 0; j < tree->nnodes - 1; j++) {
    tree->nodelist[j]->up = tree->nodelist[ ivec[j] ];
    if (tree->nodelist[ ivec[j] ]->left == NULL) tree->nodelist[ ivec[j] ]->left  = tree->nodelist[j];
    else                                         tree->nodelist[ ivec[j] ]->right = tree->nodelist[j];
  }

  tree->root = tree->nodelist[j]; /* now we have  j = tree->nnodes - 1 */

  update_topology_sisters (tree);
  update_topology_traversal (tree);

  return tree->nnodes - 1;
}

void
copy_topology_to_intvector_by_id (topology tree, int *ivec)
{
  int j;

  for (j = 0; j < tree->nnodes; j++) {
    if (tree->nodelist[j]->up) ivec[j] = tree->nodelist[j]->up->id;
    else if (tree->nodelist[j] == tree->root) ivec[j] = -1; /* root node */
    else biomcmc_error ("orphan node is not root, cannot copy it to int vector");
  }
}

void
copy_intvector_to_topology_by_id (topology tree, int *ivec)
{
  int j; 

  for (j = tree->nleaves; j < tree->nnodes; j++) 
    tree->nodelist[j]->up = tree->nodelist[j]->left = tree->nodelist[j]->right = NULL;
  
  for (j = 0; j < tree->nnodes; j++) {
    if (ivec[j] >= 0) {
      tree->nodelist[j]->up = tree->nodelist[ ivec[j] ];
      if (tree->nodelist[j]->up->left == NULL) tree->nodelist[j]->up->left  = tree->nodelist[j];
      else                                     tree->nodelist[j]->up->right = tree->nodelist[j];
    }
    else tree->root = tree->nodelist[j];
  }
//  for (j = 0; j < tree->nnodes; j++) printf ("%4d ", ivec[j]);  printf (" || DEBUG\n"); fflush(stdout);

  update_topology_sisters (tree);
  update_topology_traversal (tree);
}

