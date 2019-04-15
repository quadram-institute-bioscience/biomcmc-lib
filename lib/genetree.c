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

#include "genetree.h"
#include "reconciliation.h"
#include "splitset_distances.h"
#define NDISTS 6

/* not used yet, probably more useful in downstream programs */
const char *distance_names[][2] = {  // short names and long names
   {(char*) "dup", (char*) "duplication"},
   {(char*) "los", (char*) "loss"},
   {(char*) "dco", (char*) "deep coalescence"},
   {(char*) "rfd", (char*) "RF"},
   {(char*) "hdi", (char*) "Hungarian"},
   {(char*) "spr", (char*) "approx. SPR"}
};

genetree
new_genetree_speciestree_pair (topology gene, topology species)
{
  speciestree sptre = new_speciestree (species, NULL);
  genetree gtre = new_genetree (gene, sptre);
  return gtre;
}

genetree
new_genetree (topology gene, speciestree sptre)
{
  int i;
  genetree gtre;

  gtre = (genetree) biomcmc_malloc (sizeof (struct genetree_struct));
  gtre->ref_counter = 1;
  gtre->t = gene;
  gtre->t->ref_counter++;
  gtre->sptre = NULL; 

  gtre->distance = (int*) biomcmc_malloc (NDISTS * sizeof (int));
  gtre->minmax   = (int*) biomcmc_malloc (2 * NDISTS * sizeof (int));
  for (i=0; i < NDISTS; i++) { gtre->minmax[i] = INT32_MAX;  gtre->minmax[i + NDISTS] = INT32_MIN; }
 
  gtre->rec = new_reconciliation (gtre->t->nleaves, sptre->t->nleaves);
  reconciliation_index_sptaxa_to_genetaxa (sptre->t->taxlabel, gene->taxlabel, gtre->rec->sp_id, sptre->spnames_order);
  initialize_reconciliation_sp_count (gtre->rec, sptre->t->taxlabel->nstrings, gene->nleaves);
  initialize_reconciliation_from_new_species_tree (gtre, sptre); // points to current sptree and updates node pointers
  gtre->split = new_splitset_genespecies (gene, sptre->t, gtre->rec);
  return gtre;
}

void
del_genetree (genetree gtre)
{
  if (!gtre) return;
  if (--gtre->ref_counter) return;
  del_topology (gtre->t);
  del_speciestree (gtre->sptre);
  del_reconciliation (gtre->rec);
  if (gtre->distance) free (gtre->distance);
  if (gtre->minmax) free (gtre->minmax);
  free (gtre);
  return;
}

speciestree
new_speciestree (topology species, int *order_of_species_names)
{
  speciestree sptre;
  int i, n_mrca = (species->nnodes * (species->nnodes-1))/2; // mrca has all _nodes_ not only leaves
  sptre = (speciestree) biomcmc_malloc (sizeof (struct speciestree_struct));
  sptre->ref_counter = 1;
  sptre->t = species;
  sptre->t->ref_counter++;
  sptre->mrca = (topol_node*) biomcmc_malloc (n_mrca * sizeof (topol_node));
  for (i=0; i < n_mrca; i++) sptre->mrca[i] = NULL;
  if (order_of_species_names) sptre->spnames_order = order_of_species_names; // CAUTION: no pointer checking or ref_counter 
  else { // default is to reorder char_vector, unless user wants specific order (and thus must provide order, from largest to smallest)
    sptre->spnames_order = NULL;
    reorder_topology_leaves (sptre->t);
  }
  return sptre;
}

void
del_speciestree (speciestree sptre)
{
  if (!sptre) return;
  if (--sptre->ref_counter) return;
  if (sptre->mrca) free (sptre->mrca);
  if (sptre->spnames_order) free (sptre->spnames_order);
  del_topology (sptre->t);
  free (sptre);
  return;
}

void
genetree_speciestree_distances (genetree gtre, speciestree sptre)
{
  int k;
  reconciliation_gene_tree_reconcile (gtre, sptre);
  dSPR_gene_species (gtre->t, sptre->t, gtre->split);
  gtre->distance[0] = gtre->rec->ndups;
  gtre->distance[1] = gtre->rec->nloss;
  gtre->distance[2] = gtre->rec->ndcos;
  gtre->distance[3] = gtre->split->rf;
  gtre->distance[4] = gtre->split->hdist;
  gtre->distance[5] = gtre->split->spr + gtre->split->spr_extra;
  for (k=0; k < NDISTS; k++) if (gtre->minmax[k] > gtre->distance[k]) gtre->minmax[k] = gtre->distance[k]; // min 
  for (k=0; k < NDISTS; k++) if (gtre->minmax[k+NDISTS] < gtre->distance[k]) gtre->minmax[k+NDISTS] = gtre->distance[k]; // max
}

void  // these should call lower level local functions
genetree_reconcile_speciestree (genetree gtre, speciestree sptre)
{ // updates gtre->sptre since uses pointers to sptre nodes
  reconciliation_gene_tree_reconcile (gtre, sptre);
  return;
}

void
genetree_dSPR_speciestree (genetree gtre, speciestree sptre, int level)
{ // does not update gtre->sptre
  if (level == 2) dSPR_gene_species (gtre->t, sptre->t, gtre->split);
  else if (level == 1) dSPR_gene_species_hdist (gtre->t, sptre->t, gtre->split);
  else dSPR_gene_species_rf (gtre->t, sptre->t, gtre->split);
}
