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

#include "reconciliation.h"

topol_node mrca_between_nodes (speciestree sptre, int i, int j);
void gene_tree_reconcile_unrooted (topology gene, topology species);
void prepare_for_loss_calculation (topology gene, topology species);

reconciliation
new_reconciliation (int gene_nleaves, int sp_nleaves)
{
  int i, nnodes = 2 * gene_nleaves - 1, sizeofnode = sizeof (topol_node);

  reconciliation r = (reconciliation) biomcmc_malloc (sizeof (struct reconciliation_struct));
  r->ndups = -1;
  r->nloss = -1;
  r->ndcos = -1;
  r->sp_size = 0; /* number of species represented */
  r->size_diff = 0; /* 2 X (gene_nleaves - sp_size) */ 

  r->map_d  = (topol_node*) biomcmc_malloc (nnodes * sizeofnode); /* sptree node "below" edge */ 
  r->map_u  = (topol_node*) biomcmc_malloc (nnodes * sizeofnode); /* sptree node "above" edge */
  r->ndup_d = (int*)        biomcmc_malloc (nnodes * sizeof (int)); /* partial number of dups below */
  r->ndup_u = (int*)        biomcmc_malloc (nnodes * sizeof (int)); /* partial number of dups above */
  r->nlos_d = (int*)        biomcmc_malloc (nnodes * sizeof (int)); /* partial number of losses below */
  r->nlos_u = (int*)        biomcmc_malloc (nnodes * sizeof (int)); /* partial number of losses above */
  r->sp_id  = (int*)        biomcmc_malloc (gene_nleaves * sizeof (int));
  r->sp_count = (int*) biomcmc_malloc (sp_nleaves * sizeof (int)); /* frequency of each species (set only once) */ 
  r->map_d[0] = r->map_u[0] = NULL;
  for (i = 0; i < gene_nleaves; i++) r->nlos_d[i] = r->ndup_d[i] = 0; /* number of duplosses on terminal branches (below) */

  return r;
}

reconciliation
new_reconciliation_from_reconciliation (int gene_nleaves, int sp_nleaves, reconciliation from)
{
  int i;
  reconciliation r;

  r = new_reconciliation (gene_nleaves, sp_nleaves);
  r->ndups = from->ndups;
  r->nloss = from->nloss;
  r->ndcos = from->ndcos;
  r->sp_size = from->sp_size; 
  r->size_diff = from->size_diff;
  for (i = 0; i < gene_nleaves; i++) {
     r->nlos_d[i] = from->nlos_d[i];
     r->ndup_d[i] = from->ndup_d[i];
     r->sp_id[i] = from->sp_id[i];
  }
  for (i = 0; i < sp_nleaves; i++) r->sp_count[i] = from->sp_count[i];

  return r;
}

void
del_reconciliation (reconciliation r)
{
  if (!r) return;
  if (r->map_d)  free (r->map_d);
  if (r->map_u)  free (r->map_u);
  if (r->sp_id)  free (r->sp_id);
  if (r->ndup_d) free (r->ndup_d);
  if (r->ndup_u) free (r->ndup_u);
  if (r->nlos_d) free (r->nlos_d);
  if (r->nlos_u) free (r->nlos_u);
  if (r->sp_count) free (r->sp_count);
  free (r);
}

void
reconciliation_index_sptaxa_to_genetaxa (char_vector species, char_vector gene, int *sp_idx_in_gene, empfreq ef_external)
{ /* Don't need gene tree, just gene leaf names; alternative is to  use gene->rec instead of gene->rec->sp_id (int*) */
  int i, j, n_index = gene->nstrings, *index;
  empfreq ef;

  /* Search first largest species names (so that for example "ecoli" will match only if "ecoliII" doesn't) */
  if (ef_external) ef = ef_external; /* empfreq initialized outside (by calling function, e.g) */
  else ef = new_empfreq_sort_decreasing (species->nchars, species->nstrings, 1); /* local ("1" means size_t) */

  index = (int*) biomcmc_malloc (n_index * sizeof (int));
  for (i=0; i < n_index; i++) { 
    sp_idx_in_gene[i] = -1; /* initialize mapping */ 
    index[i] = i;  /* scan gene leaves _without_ replacement */
  }

  for (i=0; i < species->nstrings; i++) for (j=0; j < n_index; j++) /* search sp name in all unmapped gene names */
    if ((gene->nchars[index[j]] >= species->nchars[ ef->i[i].idx ]) && /* ordered species names (by nchars[]) */
        (strcasestr (gene->string[index[j]], species->string[ ef->i[i].idx ]))) { 
      /* found species name within gene name; we have a mapping */
      sp_idx_in_gene[ index[j] ] = ef->i[i].idx;
      index[j] = index[--n_index]; // with this the whole search takes O(N ln N),
      j--; // index[j] is now a distinct element (the last)
    }

  if (n_index) {
    fprintf (stderr, "Couldn't find species for genes:\n");
    for (i=0; i < n_index; i++) fprintf (stderr, " \"%s\"\n", gene->string[index[i]]);
    biomcmc_error ("gene names should contain the name of species");
  }

  if (!ef_external) del_empfreq (ef); 
  if (index) free (index);
}

void
initialize_reconciliation_sp_count (reconciliation rec, int n_sp, int n_idx)
{
  int i; /* rec->sp_id[i] is the species index for gene i */
  for (i = 0; i < n_sp; i++)  rec->sp_count[i] = 0; /* representativity of each species in gene family */
  for (i = 0; i < n_idx; i++) rec->sp_count[ rec->sp_id[i] ]++; /* update species frequencies */

  rec->sp_size = 0;
  for (i = 0; i < n_sp; i++) if (rec->sp_count[i]) rec->sp_size++;
  rec->size_diff = 2 * (n_idx - rec->sp_size); /* term to calculate deepcoals */ 
}

void
initialize_reconciliation_from_new_species_tree (genetree gtre, speciestree sptre)
{
  int i, n_mrca = (sptre->t->nnodes * (sptre->t->nnodes-1))/2;
  if (sptre == gtre->sptre) return; // already up-to-date species tree
  del_speciestree (gtre->sptre);
  gtre->speciestree = sptre; gtre->speciestree->ref_counter++;
  for (i=0; i < gtre->t->nleaves; i++) gtre->rec->map_d[i] = sptre->t->nodelist[ gtre->rec->sp_id[i] ];
  
  if (!sptre->t->traversal_updated) {
    for (i=0; i < n_mrca; i++) sptre->mrca[i] = NULL;
    update_topology_traversal (sptre->t);
  }
}

topol_node
mrca_between_nodes (speciestree sptre, int i, int j)
{
  topol_node p;
  bipartition lowsplit;
  int i1, j1, index;

  if (i == j) return sptre->t->nodelist[i];
  if (j > i) { int tmp = i; i = j; j = tmp; } /* i should always be larger than j */
  index = (i * (i-1))/2 + j;
  if (sptre->mrca[index]) return sptre->mrca[index];


  /* Minimize sequencial search by choosing node closer to root (topol_node::level = distance from root) */
  if (sptre->t->nodelist[i]->level > sptre->t->nodelist[j]->level) 
   { i1 = i; j1 = j; } /* node[i1] will be further from root */
  else 
   { i1 = j; j1 = i; } /* node[i1] will be further from root */
 
  p = sptre->t->nodelist[j1]; /* here we start to actually search for the lca */

  /* if node[i1] - which is further from root - is internal we must compare fully both bipartitions - O(n) */
  if ((sptre->t->nodelist[i1]->internal) && (lowsplit = sptre->t->nodelist[i1]->split))
    while ((p) && (!bipartition_contains_bits (p->split, lowsplit))) p = p->up; /* climb up the tree */
  else /* if node[i1] is a leaf the comparison is faster - O(1) - for very large trees */
    while ((p) && (!bipartition_is_bit_set (p->split, i1))) p = p->up; /* climb up the tree */

  if (!p) biomcmc_error ("Couldn't find the MRCA. Possible bug related to root node."); /* this shouldn't happen(R) */

  return sptre->t->mrca[index] = p;
}

void
reconciliation_gene_tree_reconcile (genetree gtre, speciestree sptre)
{
  int i, g_id;
  topol_node map_lchild, map_rchild;

  initialize_reconciliation_from_species_tree (gtre, stre);

  prepare_for_loss_calculation (gtre->t, sptre->t);

  for (i=0; i < gtre->t->nleaves-1; i++) {
    g_id = gtre->t->postorder[i]->id; /* node ID on gene tree */
    map_lchild = gtre->t->rec->map_d[ gtre->t->postorder[i]->left->id ]; /* gene->map[] are nodes on species tree */
    map_rchild = gtre->t->rec->map_d[ gtre->t->postorder[i]->right->id ];
    gtre->rec->map_d[g_id] = mrca_between_nodes (sptre, map_lchild->id, map_rchild->id);

    /* cummulative number of duplications below node, following e.g. Bioinformatics.2001.821 */
    gtre->rec->ndup_d[g_id] = gtre->rec->ndup_d[gtre->t->postorder[i]->left->id] + gtre->rec->ndup_d[gtre->t->postorder[i]->right->id]; 
    /* cummul. number of losses below node, following SIAM.2000.729 (not over edges like ACMTransComputBiolBioinfo.2010.14) */
    gtre->rec->nlos_d[g_id] = gtre->rec->nlos_d[gtre->t->postorder[i]->left->id] + gtre->rec->nlos_d[gtre->t->postorder[i]->right->id];

    if ((gtre->rec->map_d[g_id] == map_lchild) || (gtre->rec->map_d[g_id] == map_rchild)) {
      gtre->rec->ndup_d[g_id]++;
      if (map_lchild != map_rchild) { /* if all three are the same there are no losses */
        if (map_lchild == gtre->rec->map_d[g_id]) /* number of intermediate nodes in sptree + 1 */
          gtre->rec->nlos_d[g_id] += (map_rchild->mid[4] - gtre->rec->map_d[g_id]->mid[4]);
        else
          gtre->rec->nlos_d[g_id] += (map_lchild->mid[4] - gtre->rec->map_d[g_id]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's children ("-2" since difference in levels=1 means NO interm) */
      gtre->rec->nlos_d[g_id] += (map_lchild->mid[4] + map_rchild->mid[4] - 2 * gtre->rec->map_d[g_id]->mid[4] - 2);
  }

  gene_tree_reconcile_unrooted (gtre, sptre);
}

void
gene_tree_reconcile_unrooted (genetree gtre, speciestree sptre)
{ /* total n_dups rooted at node should be ndup_u + ndup_d + Indicator{ mrca(map_d, map_u) } */
  int i, g_id, r_left = gtre->t->root->left->id, r_right = gtre->t->root->right->id, thisloss, thisdups, thiscoal, 
      min_coal = 0xffffff, min_dups = 0xffffff, min_loss = 0xffffff; /* large number */
  topol_node map_root, map_up, map_sister;

  gtre->rec->map_u[r_left]   = gtre->rec->map_d[r_right]; /* neglect root node; r_left and r_right have same info */
  gtre->rec->map_u[r_right]  = gtre->rec->map_d[r_left];
  gtre->rec->ndup_u[r_left]  = gtre->rec->ndup_d[r_right];
  gtre->rec->ndup_u[r_right] = gtre->rec->ndup_d[r_left];
  gtre->rec->nlos_u[r_left]  = gtre->rec->nlos_d[r_right];
  gtre->rec->nlos_u[r_right] = gtre->rec->nlos_d[r_left];

  /* as "rooted" version, but replacing left and right by up and sister (obs: postorder[gene->nleaves-2] => root) */
  for (i = gtre->t->nleaves-3; i >= 0; i--) if ((gtre->t->postorder[i]->id != r_left) && (gtre->t->postorder[i]->id != r_right)) {
    g_id = gtre->t->tpostorder[i]->id; /* node ID on gene tree */
    map_up     = gtre->rec->map_u[ gtre->t->postorder[i]->up->id ]; /* gene->map[] are nodes on species tree */
    map_sister = gtre->rec->map_d[ gtre->t->postorder[i]->sister->id ];
    gtre->rec->map_u[g_id] = mrca_between_nodes (sptre, map_up->id, map_sister->id);

    gtre->rec->ndup_u[g_id] = gtre->rec->ndup_u[gtre->t->postorder[i]->up->id] + gtre->rec->ndup_d[gtre->t->postorder[i]->sister->id]; 
    gtre->rec->nlos_u[g_id] = gtre->rec->nlos_u[gtre->t->postorder[i]->up->id] + gtre->rec->nlos_d[gtre->t->postorder[i]->sister->id]; 
    
    if ((gtre->rec->map_u[g_id] == map_up) || (gtre->rec->map_u[g_id] == map_sister)) {
      gtre->rec->ndup_u[g_id]++;
      if (map_up != map_sister) { /* if all three are the same there are no losses */
        if (map_sister == gtre->rec->map_u[g_id]) /* then map_up is distinct from mrca(up,sister) */ 
          gtre->rec->nlos_u[g_id] += (map_up->mid[4] - gtre->rec->map_u[g_id]->mid[4]);
        else
          gtre->rec->nlos_u[g_id] += (map_sister->mid[4] - gtre->rec->map_u[g_id]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's "children" ("-2" since difference in levels=1 means NO interm) */
      gtre->rec->nlos_u[g_id] += (map_sister->mid[4] + map_up->mid[4] - 2 * gtre->rec->map_u[g_id]->mid[4] - 2);
  }

  /* we went in preorder over internal nodes; at last, we go over leaves */
  for (i = 0; i < gtre->t->nleaves; i++) if ((gtre->t->nodelist[i]->id != r_left) && (gtre->t->nodelist[i]->id != r_right)) {
    map_up     = gtre->rec->map_u[ gtre->t->nodelist[i]->up->id ]; /* gene->map[] are nodes on species tree */
    map_sister = gtre->rec->map_d[ gtre->t->nodelist[i]->sister->id ];
    gtre->rec->map_u[i] = mrca_between_nodes (sptre, map_up->id, map_sister->id);

    gtre->rec->ndup_u[i] = gtre->rec->ndup_u[gtre->t->nodelist[i]->up->id] + gtre->rec->ndup_d[gtre->t->nodelist[i]->sister->id]; 
    gtre->rec->nlos_u[i] = gtre->rec->nlos_u[gtre->t->nodelist[i]->up->id] + gtre->rec->nlos_d[gtre->t->nodelist[i]->sister->id];

    if ((gtre->rec->map_u[i] == map_up) || (gtre->rec->map_u[i] == map_sister)) {
      gtre->rec->ndup_u[i]++;

      if (map_up != map_sister) { /* if all three are the same there are no losses */
        if (map_sister == gtre->rec->map_u[i]) /* then map_up is distinct from mrca(up,sister) */ 
          gtre->rec->nlos_u[i] += (map_up->mid[4] - gtre->rec->map_u[i]->mid[4]);
        else
          gtre->rec->nlos_u[i] += (map_sister->mid[4] - gtre->rec->map_u[i]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's "children" ("-2" since difference in levels=1 means NO interm) */
      gtre->rec->nlos_u[i] += (map_sister->mid[4] + map_up->mid[4] - 2 * gtre->rec->map_u[i]->mid[4] - 2);
  }

  /* create virtual roots at every edge to calc the total number of dups (above + below) */
  for (i = 0; i < gtre->t->nnodes; i++) if ((i != r_right) && (i != gtre->t->root->id)) { /* only r_left (or r_right) is needed */
    map_root = mrca_between_nodes (sptre, gtre->rec->map_u[i]->id, gtre->rec->map_d[i]->id);
    /* duplications */ 
    thisdups = gtre->rec->ndup_u[i] + gtre->rec->ndup_d[i]; /* g_id now is the number of dups at virtual root */
    if ((map_root == gtre->rec->map_u[i]) || (map_root == gtre->rec->map_d[i])) thisdups++;
    /* losses */
    thisloss = gtre->rec->nlos_u[i] + gtre->rec->nlos_d[i]; /* number of losses at virtual root */
    if      ((map_root == gtre->rec->map_u[i]) && (map_root != gtre->rec->map_d[i]))
      thisloss += (gtre->rec->map_d[i]->mid[4] - map_root->mid[4]); // "d(a(g),g) + 1" in SIAM.2000.729 
    else if ((map_root != gtre->rec->map_u[i]) && (map_root == gtre->rec->map_d[i]))
      thisloss += (gtre->rec->map_u[i]->mid[4] - map_root->mid[4]); // "d(a(g),g) + 1" in SIAM.2000.729  
    else if ((map_root != gtre->rec->map_u[i]) && (map_root != gtre->rec->map_d[i]))
      thisloss += (gtre->rec->map_u[i]->mid[4] + gene->rec->map_d[i]->mid[4] - 2 * map_root->mid[4] - 2);
    // else "loss(g) = 0" following SIAM.2000.729
    /* deep coalescences = loss - 2 x dups + 2 x |leaf difference between gene and species trees| */
    thiscoal = thisloss - 2 * thisdups + gtre->rec->size_diff;

    if (thisdups < min_dups) min_dups = thisdups;
    if (thisloss < min_loss) min_loss = thisloss; 
    if (thiscoal < min_coal) min_coal = thiscoal;
  }
  gtre->rec->ndups = min_dups;
  gtre->rec->nloss = min_loss;
  gtre->rec->ndcos = min_coal;
}

void
prepare_for_loss_calculation (topology gene, topology species)
{ 
  int i, c_l, c_r;

  for (i = 0; i < species->nleaves; i++) species->nodelist[i]->mid[2] = gene->rec->sp_count[i];
  for (i = 0; i < species->nleaves - 1; i++) { 
    /* mid[2] is the "effective" subtree cardinality: corrects for duplicates and absent species */
    c_l = species->postorder[i]->left->mid[2]; 
    c_r = species->postorder[i]->right->mid[2];
    species->postorder[i]->mid[2] = c_l + c_r; 

    /* mid[3] indicates if node is active or not (0 = pruned; 1 = normal; 0xffff = dummy node (only one child) */
    if ((!c_l) && (!c_r)) species->postorder[i]->mid[3] = 0; /* both children absent from gene: exclude this node */ 
    else if ((c_l) && (c_r)) species->postorder[i]->mid[3] = 1; /* both children present: this node is eligible */
    /* only one active child: this node will simply duplicate values from single valid child */
    //else species->postorder[i]->mid[3] = (c_l? 2: 3); /* id=2 if active child is left, id=3 if active is right node */
    else species->postorder[i]->mid[3] = 0xffff; 
  }

  /* in preorder: mid[4] has level (distance from root) taking into account only active species */
  if (species->root->mid[3] == 1) species->root->mid[4] = 0;
  else                            species->root->mid[4] = -1;

  for (i = species->nleaves-3; i >= 0; i--) {
    if (species->postorder[i]->mid[3] == 1)
      species->postorder[i]->mid[4] = species->postorder[i]->up->mid[4] + 1;
    else /* works if mid[3] is larger than one; if mid[3] is zero this node won't be mapped into anyway */
      species->postorder[i]->mid[4] = species->postorder[i]->up->mid[4];
  }
  for (i = 0; i < species->nleaves; i++) if (species->nodelist[i]->mid[2]) /* only leaves with one or more homologs */
    species->nodelist[i]->mid[4] = species->nodelist[i]->up->mid[4] + 1;
}

 /* OLD NOTE about DEEPCOAL (now we use nloss - 2xndups):
  * following ThanNahleh.PLoSComputBiol.2009.preprint, for "regular nodes" extra lineages above a node are equal to
  * |subtree rooted at node| - # coalescences below subtree (including root node of subtree) - 1 
  * No details are given for leaves, but on MolPhylEvol.1997.349 they mention the need for creating artificial nodes
  * for duplicated species and removal of subtrees with species unrepresented in gene family. For deep coal, we look
  * only at the species tree [1], where the mapping represents the coalescences. 
  * [1] the exception are maybe the species with more than one copy, but I didn't look at that. */ 
