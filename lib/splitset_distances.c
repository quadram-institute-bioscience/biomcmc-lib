/* 
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

#include "splitset_distances.h"

splitset new_splitset (int nleaves, int sp_nleaves); 
void prepare_genetree_sptree_split (topology gene, topology species, splitset split);
int rf_hdist_topology_lowlevel (splitset split, bool exit_at_rf);
int dSPR_topology_lowlevel (splitset split);

/*! \brief vector with pointers to bipartitions that are identical on both trees */
void split_create_agreement_list (splitset split);
void split_remove_agree_edges (splitset split, bipartition *b, int *nb);
void split_remove_duplicates (bipartition *b, int *nb);
void split_compress_agreement (splitset split);
void split_create_disagreement_list (splitset split);
void split_disagreement_assign_match (splitset split);
void split_find_small_disagreement (splitset split);
void split_remove_small_disagreement (splitset split);
void split_minimize_subtrees (splitset split);
void split_remove_redundant_bit (splitset split, int id);
void split_replace_bit (splitset split, int to, int from);
void split_new_size (splitset split, int size, bool update_bipartitions);
void split_swap_position (bipartition *b, int i1, int i2);

splitset 
new_splitset (int nleaves, int sp_nleaves)
{
  splitset split;
  int i;

  split = (splitset) biomcmc_malloc (sizeof (struct splitset_struct));
  split->size = split->spsize = nleaves - 1; /* actually can be nleaves-3, but we use disagree[] as temporary for gene/species comparison */
  split->n_agree = split->n_disagree = 0;
  split->match = false;
  split->spr = split->spr_extra = split->rf = split->hdist = split->hdist_reduced = 0;

  /* species tree before removal of missing species has at most 2*sp_nleaves-1 bipartitions, and after their removal _but_
   * considering subtrees spanned by gene tree (multree) has at most 2*gene_nleaves-1 bipartitions. We don't need
   * the leaves (trivial bipartitions) -- therefore the 2 x nleaves is overkill --, but sp0[] needs sp_nleaves. */
  if (nleaves > sp_nleaves) split->spsize = 2 * sp_nleaves + nleaves; 
  else                      split->spsize = 3 * sp_nleaves; /* +sp leaves for the first elems */

  split->s_split  = (bipartition*) biomcmc_malloc (split->spsize * sizeof (bipartition));
  split->g_split  = (bipartition*) biomcmc_malloc (split->size * sizeof (bipartition));
  split->agree    = (bipartition*) biomcmc_malloc (split->size * sizeof (bipartition));
  split->disagree = (bipartition*) biomcmc_malloc (split->size * split->size * sizeof (bipartition));

  split->s_split[0]  = new_bipartition (nleaves);
  split->g_split[0]  = new_bipartition (nleaves);
  split->agree[0]    = new_bipartition (nleaves); 
  split->disagree[0] = new_bipartition (nleaves); 
  for (i = 1; i < split->size; i++) {
    split->g_split[i] = new_bipartition_from_bipsize (split->g_split[0]->n); /* use same bipsize */
    split->agree[i]   = new_bipartition_from_bipsize (split->agree[0]->n);
  }
  for (i = 1; i < split->size * split->size; i++) split->disagree[i] = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (i = 1; i < split->spsize; i++) split->s_split[i] = new_bipartition_from_bipsize (split->s_split[0]->n);

  split->n_g = split->n_s = 0; /* number of elements actively in use */
  split->prune = new_bipartition_from_bipsize (split->disagree[0]->n);
  split->h = new_hungarian (split->size, false); // true/false -> double/integer

  split->sp0 = split->s_split; /* E.g.  sp0 => [0|1|2|3|4|5|6|7] while s_split => [4|5|6|7] */
  split->s_split = split->s_split + sp_nleaves; /* pointer to nleaves-th element */

  return split;
}

void
del_splitset (splitset split)
{
  int i;
  if (!split) return;

  del_bipartition (split->prune);
  if (split->disagree) {
    for (i = split->size * split->size - 1; i >= 0; i--) del_bipartition (split->disagree[i]);
    free (split->disagree);
  }
  if (split->agree) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->agree[i]);
    free (split->agree);
  }
  if (split->g_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->g_split[i]);
    free (split->g_split);
  }
  if (split->sp0) {
    for (i = split->spsize - 1; i >= 0; i--) del_bipartition (split->sp0[i]);
    free (split->sp0);
  }
  del_hungarian (split->h);
  free (split);
}

splitset
new_splitset_genespecies (topology gene, topology species, reconciliation rec)
{
  int i, j; 
  splitset split = new_splitset (gene->nleaves, species->nleaves);

  for (i=0; i < species->nleaves; i++) { /* initialize bipartitions at leaves with gene leave information */
    bipartition_zero (split->sp0[i]);
    for (j = 0; j < gene->nleaves; j++) if (rec->sp_id[j] == i) bipartition_set (split->sp0[i], j);
  }
  /* temporarily store (total, not reduced) number of bipartitions on both trees (to available to max_distance calcs etc.)*/
  split->n_g = gene->nleaves - 3;
  split->n_s = rec->sp_size - 3;
  for (i=0; i < species->nleaves; i++) if (rec->sp_count[i] > 1) split->n_s++; // species not a leaf, but a cherry actually
  /* we could exclude species absent from gene tree (more memory efficient) [1], but it's easier to leave them since indexes
   * are preserved (when calling postorder[]->left->id for instance), otherwise we should need a index vector to map tree
   * IDs to positions here.  [1] means swap_bipartition, s_split -= removed_leaves, realloc (sp0, smaller_size)  
   * Everything works fine here since bitstrings always refer to gene labels (that is, absent species are all zero or
   * repeated bitstrings, when copying sptree structure) */
  return split;
}

void
prepare_genetree_sptree_split (topology gene, topology species, splitset split)
{
  int i;
  if (!gene->traversal_updated)    update_topology_traversal (gene);
  if (!species->traversal_updated) update_topology_traversal (species);

  bipsize_resize (split->g_split[0]->n, split->g_split[0]->n->original_size); 
  bipsize_resize (split->s_split[0]->n, split->s_split[0]->n->original_size); 

  /* use sp0[] that has info about the species leaves (coded for gene ids as a multifurcation) */ 
  for (i=0; i < species->nleaves - 1; i++) bipartition_zero (split->s_split[i]); /* some might be skipped below */
  for (i=0; i < species->nleaves - 1; i++) {
    bipartition_OR (split->sp0[species->postorder[i]->id], /* postorder[i]->id is always larger than nleaves */
                    split->sp0[species->postorder[i]->left->id], /* sp0 has elems not accessible by s_split (w/ leaves info) */
                    split->sp0[species->postorder[i]->right->id], false);
    /* TODO this is where tripartition is updated */
  }
  split->n_s = species->nleaves - 1;
  /* some spleaves may be zero (absent on gene) leading to internal nodes with < 2 elems, thus we need to "minimize" the species subtrees 
   * (not to mention some other skipped node in postorder) */ 
  for (i = 0; i < split->n_s; i++) { /* here we use s_split[], that points to internal bipartitions only */
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; split_swap_position (split->s_split, i, split->n_s); i--; }
  }

  split_remove_duplicates (split->s_split, &(split->n_s));
  for (i=0; i < species->nleaves; i++) if (split->sp0[i]->n_ones > 1) { /* add multifurcations representing species on gene tree */
    bipartition_copy (split->s_split[split->n_s], split->sp0[i]); /* idea from RF distance (arxiv.1210.2665) */
    bipartition_flip_to_smaller_set (split->s_split[split->n_s++]); /* this species has many many copies, that is most of leaves are same species... */
  }
  split_remove_duplicates (split->s_split, &(split->n_s));

  /* gene tree bipartitions are simpler */
  split->n_g = gene->nleaves - 3;
  for (i=0; i < split->n_g; i++) {
    bipartition_copy (split->g_split[i], gene->postorder[i]->split);
    /* TODO this is where tripartition is updated */
    bipartition_flip_to_smaller_set (split->g_split[i]);
  }
}

/* empirical observation: split->spr below tends to overestimate (very rarely it overestimates), while 
 * spr + spr_extra may overestimate, have a higher variability but has the impressive property of recognizing even a
 * lot of SPRs. A good compromise seems to weight the contribution of spr_extra. BTW this variable counts the number
 * of "swapped" prune edges -- that is, a prune bipartition representing subtrees that are in opposite sides of the
 * edge, that usually mean 2 SPRs like (A,B)-+-(C,D)  <-> (A,C)-+-(B,D) */  
/* split->spr += (split->spr_extra/2); // I'll comment it out and leave it to the calling function */

int
dSPR_gene_species (topology gene, topology species, splitset split)
{
  // first calculate Hdist on original (not reduced) trees, then prepare again (to use reduced trees)
  prepare_genetree_sptree_split (gene, species, split);
  rf_hdist_topology_lowlevel (split, false); // hdist, rf
  if (!split->rf) return 0;
  prepare_genetree_sptree_split (gene, species, split);
  return dSPR_topology_lowlevel (split); // hdist_reduced, spr, spr_extra
}

int
dSPR_gene_species_rf (topology gene, topology species, splitset split)
{
  prepare_genetree_sptree_split (gene, species, split);
  return rf_hdist_topology_lowlevel (split, true); // true -> exit as soon as RF is calculated
}

int
dSPR_gene_species_hdist (topology gene, topology species, splitset split)
{
  prepare_genetree_sptree_split (gene, species, split);
  return rf_hdist_topology_lowlevel (split, false);
}

int
rf_hdist_topology_lowlevel (splitset split, bool exit_at_rf)
{
  split->hdist_reduced = split->hdist = split->rf = split->spr = split->spr_extra = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits); 
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits); 
  split_create_agreement_list (split);  // vector of identical bipartitions
  // Importantly, here we do NOT call split_compress_agreement()
  split->rf = split->n_g + split->n_s;
  if (exit_at_rf) return split->rf;
  if (!split->rf) return 0; // all edges were in agreement
  split->match = true; // only calculate hdist_reduced if match == True (first time) 
  split_create_disagreement_list (split); // vector of smallest disagreements
  split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist_reduced after first time) */
  split->hdist = split->hdist_reduced;
  return split->hdist; /* do not calculate SPR, exit now */
}

int
dSPR_topology_lowlevel (splitset split)
{
  int i = 0, mismatch = -1;
  split->match = true; 
  split->hdist_reduced = split->spr = split->spr_extra = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits); 
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits); 
  //for (i = 0; i < split->n_g; i++) { bipartition_print_to_stdout (split->g_split[i]); } printf ("G   start::DEBUG::\n");
  //for (i = 0; i < split->n_s; i++) { bipartition_print_to_stdout (split->s_split[i]); } printf ("S\n");

  i++; /* to trick -Werror, since we don't use it unless for debug */
  while (mismatch) {
    split_create_agreement_list (split);  // vector of identical bipartitions
    split_compress_agreement (split);     // iterative replacement of cherry by new leaf

    //for (i = 0; i < split->n_g; i++)     { bipartition_print_to_stdout (split->g_split[i]); } printf ("G   ::DEBUG::\n");
    //for (i = 0; i < split->n_s; i++)     { bipartition_print_to_stdout (split->s_split[i]); } printf ("S\n");
    //for (i = 0; i < split->n_agree; i++) { bipartition_print_to_stdout (split->agree[i]);   } printf ("A\n");

    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
    if (!mismatch) return split->spr;
    
    split_create_disagreement_list (split); // vector of smallest disagreements
    split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist after first time) */
    
    split_remove_duplicates (split->disagree, &(split->n_disagree)); // some elements are equal; this function also qsorts 
    split_find_small_disagreement (split);  // could also be one leaf only 
    
    //for (i = 0; i < split->n_disagree; i++) { bipartition_print_to_stdout (split->disagree[i]); printf ("\n"); }
    //printf ("{%d} prune: ", split->n_disagree); bipartition_print_to_stdout (split->prune); printf ("\n");

    split->spr++;
    split_remove_small_disagreement (split);

    split_minimize_subtrees (split);
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
  }
  return split->spr;
}

void
split_create_agreement_list (splitset split)
{
  int s, g;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    if (bipartition_is_equal (split->g_split[g], split->s_split[s])) {
      bipartition_copy (split->agree[split->n_agree++], split->g_split[g]); 
      split->n_g--; split_swap_position (split->g_split, g, split->n_g); /* if we don't swap them, we lose ref to "old" value on g_split[] */
      split->n_s--; split_swap_position (split->s_split, s, split->n_s); 
      g--; s = split->n_s; /* pretend loop finished, examine again with new values */
    }
  split_remove_agree_edges (split, split->g_split, &(split->n_g));
  split_remove_agree_edges (split, split->s_split, &(split->n_s));
}

void
split_remove_agree_edges (splitset split, bipartition *b, int *nb)
{
  int i, a;
  for (i = 0; i < (*nb); i++) for (a = 0; a < split->n_agree; a++) 
    if (bipartition_is_equal (b[i], split->agree[a])) {
      (*nb)--; 
      split_swap_position (b, i, (*nb));
      i--;
      a = split->n_agree; /* loop again over new value */
    }
}

void
split_remove_duplicates (bipartition *b, int *nb)
{
  int i, j;
  bipartition pivot;
  if ((*nb) < 2) return; /* only if we have a vector with > 1 element */
  qsort (b, (*nb), sizeof (bipartition), compare_bipartitions_increasing);

  for (i = (*nb) - 1; i >= 1; i--)
    if (bipartition_is_equal (b[i], b[i-1])) {
      pivot = b[i]; /* do not lose a pointer to this element */
      for (j = i; j < (*nb)-1; j++) b[j] = b[j+1];
      b[j] = pivot; /* j = (*nb) - 1, which will become obsolete through next line -->  (*nb)-- */
      (*nb)--;
    }
}

void
split_compress_agreement (splitset split)
{
  int i, j, pair[2];

  for (i = 0; i < split->n_agree; i++) if (split->agree[i]->n_ones == 2) { /* cherry in common, can be represented by just one leaf */
    bipartition_to_int_vector (split->agree[i], pair, 2);
    split_remove_redundant_bit (split, pair[1]);
    split_new_size (split,split->agree[0]->n->bits - 1, false); /* false = do not recalculate every bipartition's last elem */
    bipartition_resize_vector (split->agree, split->n_agree);
    for (j = 0; j < split->n_agree; j++) { /* minimize subtree size and remove single leaves */
      bipartition_flip_to_smaller_set (split->agree[j]); /* agree only */
      if (split->agree[j]->n_ones < 2) split_swap_position (split->agree, j--, --split->n_agree);
    }
    i = -1; /* redo all iterations, with new info (agree[] will be smaller) */
  }
  bipartition_resize_vector (split->g_split, split->n_g);
  bipartition_resize_vector (split->s_split, split->n_s);
}

void
split_create_disagreement_list (splitset split)
{
  int g, s;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++) { 
    bipartition_XOR (split->disagree[g * split->n_s + s], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[g * split->n_s + s]);
  }
  split->n_disagree = split->n_g * split->n_s;
}

void
split_disagreement_assign_match (splitset split)
{ /* also calculates split->hdist_reduced */ 
  int g, s, max_n, sum = 0;
  
  if (split->n_g > split->n_s) max_n = split->n_g;
  else                         max_n = split->n_s;
  if (max_n < 2) return;
  
  hungarian_reset (split->h);
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)  
    hungarian_update_cost (split->h, g, s, &(split->disagree[g * split->n_s + s]->n_ones));
  hungarian_solve (split->h, max_n);
  /* now split->h->col_mate will have the pairs */
  /* if we do the matching below it becomes much faster, but we may miss the best prune subtrees in a few cases (do not compromise the algo) */
  split->n_disagree = 0;
  for (g = 0; g < max_n; g++) if ((g < split->n_g) && ( split->h->col_mate[g] < split->n_s)) { /* some matchings might be to dummy edges */
    bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[split->h->col_mate[g]], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    sum += split->disagree[split->n_disagree-1]->n_ones;
  }
  if (split->match) { split->hdist_reduced = split->h->final_cost+split->h->initial_cost; split->match = false; }
}

void
split_find_small_disagreement (splitset split)
{
  bipartition dis;
  int a, d;

  bipartition_copy (split->prune, split->disagree[0]); /* smallest, in case we don't find a better one in loop below */ 
  if (split->prune->n_ones < 2) return;

  dis = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (d = 0; d < split->n_disagree; d++) for (a = 0; a < split->n_agree; a++) {
    if ((split->disagree[d]->n_ones == split->agree[a]->n_ones) || 
        (split->disagree[d]->n_ones == (split->agree[a]->n->bits - split->agree[a]->n_ones))) {
      bipartition_XOR (dis, split->disagree[d], split->agree[a], true);
      if      (!dis->n_ones)               { bipartition_copy (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      else if (dis->n_ones == dis->n->bits) { bipartition_NOT (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      //      bipartition_XORNOT (dis, split->disagree[d], split->agree[a], true); /* equiv to above XOR + NOT */
      //      if (!dis->n_ones) { bipartition_NOT (split->prune, split->disagree[d]); del_bipartition (dis); return; }
    }
  }
  /* check if prune nodes are all on same side of a tree or if they are actually two SPRs (one from each tree) */
  for (d = 0; d < split->n_g; d++) {
    if (!bipartition_contains_bits (split->g_split[d], split->prune)) {
      bipartition_NOT (dis, split->g_split[d]);
      if (!bipartition_contains_bits (dis, split->prune)) { split->spr_extra++; d = split->n_g; }
    }
  } 
  del_bipartition (dis); 
}

void
split_remove_small_disagreement (splitset split)
{
  int *index, i, j = split->prune->n_ones - 1, k = 0, size = split->agree[0]->n->bits; 

  index = (int*) biomcmc_malloc (split->prune->n_ones * sizeof (int));
  bipartition_to_int_vector (split->prune, index, split->prune->n_ones);

  for (i = size - 1; i >= (size - split->prune->n_ones); i--) {
    if (index[k] >= (size - split->prune->n_ones)) i = -1;
    else {
      if (i == index[j]) j--;
      else split_replace_bit (split, index[k++], i); 
    }
  }

  split_new_size (split,size - split->prune->n_ones, true); 
  if (index) free (index);
}

void
split_minimize_subtrees (splitset split)
{
  int i;

  for (i = 0; i < split->n_s; i++) {      
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; split_swap_position (split->s_split, i, split->n_s); i--; }
  }
  for (i = 0; i < split->n_g; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    if (split->g_split[i]->n_ones < 2) { split->n_g--; split_swap_position (split->g_split, i, split->n_g); i--; }
  }
  for (i = 0; i < split->n_agree; i++) {
    bipartition_flip_to_smaller_set (split->agree[i]);
    if (split->agree[i]->n_ones < 2) { split->n_agree--; split_swap_position (split->agree, i, split->n_agree); i--; }
  }
}

void
split_remove_redundant_bit (splitset split, int id)
{
  int  last = split->agree[0]->n->bits-1;
  if (id < last) split_replace_bit (split, id, last);
}

void
split_replace_bit (splitset split, int to, int from)
{
  if (from <= to) return;
  /*not needed for disagree[] */
  bipartition_replace_bit_in_vector (split->agree,   split->n_agree, to, from, true);
  bipartition_replace_bit_in_vector (split->g_split, split->n_g,     to, from, true);
  bipartition_replace_bit_in_vector (split->s_split, split->n_s,     to, from, true);
}

void
split_new_size (splitset split, int size, bool update_bipartitions)
{
  bipsize_resize (split->g_split[0]->n, size);
  bipsize_resize (split->s_split[0]->n, size);
  bipsize_resize (split->agree[0]->n, size);
  bipsize_resize (split->disagree[0]->n, size);
  if (update_bipartitions) {
    bipartition_resize_vector (split->g_split, split->n_g);
    bipartition_resize_vector (split->s_split, split->n_s);
    bipartition_resize_vector (split->agree, split->n_agree);
  }
}

void
split_swap_position (bipartition *b, int i1, int i2)
{
  bipartition pivot = b[i1];
  b[i1] = b[i2];
  b[i2] = pivot;
}

