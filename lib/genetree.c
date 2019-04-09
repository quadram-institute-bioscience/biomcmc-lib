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

genetree
new_genetree (topology gene, speciestree sptre)
{
  genetree gtre;

  gtre = (genetree) biomcmc_malloc (sizeof (struct genetree_struct));
  gtre->ref_counter = 1;
  gtre->t = gene;
  gtre->t->ref_counter++;
  gtre->sptre = NULL; 
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
  free (gtre);
  return;
}

speciestree
new_speciestree (topology species, empfreq order_of_species_names)
{
  speciestree sptre; 
  int i, n_mrca = (species->nnodes * (species->nnodes-1))/2; // mrca has all _nodes_ not only leaves
  sptre = (speciestree) biomcmc_malloc (sizeof (struct speciestree_struct));
  sptre->ref_counter = 1;
  sptre->t = species;
  sptre->t->ref_counter++;
  sptre->mrca = (topol_node*) biomcmc_malloc (n_mrca * sizeof (topol_node));
  for (i=0; i < n_mrca; i++) sptre->mrca[i] = NULL;
  if (order_of_species_names) {
    sptre->spnames_order = order_of_species_names;
    sptre->spnames_order->ref_counter++;
  }
  else sptre->spnames_order = new_empfreq_sort_decreasing (species->taxlabel->nchars, species->taxlabel->nstrings, 1);
  return sptre;
}

void
del_speciestree (speciestree sptre)
{
  int i;
  if (!sptre) return;
  if (--sptre->ref_counter) return;
  if (sptre->mrca) free (sptre->mrca);
  del_empfreq (sptre->spnames_order);
  del_topology (sptre->t);
  free (sptre);
  return;
}

void  // these should call lower level local functions
genetree_reconcile_speciestree (genetree gtre, speciestree sptre)
{
  reconciliation_gene_tree_reconcile (gtre, sptre);
  return;
}
