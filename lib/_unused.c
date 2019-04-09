

void
init_tree_recon_from_species_topology (topology gene, topology species)
{  /* safe, since uses index_sptaxa_to_genetaxa() which does not assume sorted spnames */
  if (!gene->rec) gene->rec = new_reconciliation (gene->nleaves, species->nleaves);
  if (!species->mrca) new_mrca_for_topology (species);
  index_sptaxa_to_genetaxa (species->taxlabel, gene->taxlabel, gene->rec->sp_id, NULL);
  initialize_reconciliation_sp_count (gene->rec, species->taxlabel->nstrings, gene->nleaves);
  gene_tree_reconcile (gene, species); 
}

void
init_tree_recon_from_species_names (topology gene, char_vector sptaxlabel)
{  /* safe, since uses index_sptaxa_to_genetaxa() which does not assume sorted spnames */
  /* used by dSPR (topology_splitset) only (for now) */ 
  if (!gene->rec) gene->rec = new_reconciliation (gene->nleaves, sptaxlabel->nstrings);
  index_sptaxa_to_genetaxa (sptaxlabel, gene->taxlabel, gene->rec->sp_id, NULL);
  initialize_reconciliation_sp_count (gene->rec, sptaxlabel->nstrings, gene->nleaves);
}

/*! \brief find occurences of ordered species->string[] inside gene->string[] filling indexes in rec->sp_id[] and updating
 * rec->sp_count[]. 
 *
 *  This function can only be used whe wen know that the species names are ordered from longer to shorter. This happens
 *  only when there are no topologies (where the leaves are the indexes of taxlabel) associated with the species'
 *  taxlabels. So they cannot be used within a topology_space_struct but can be used if our primary data is only these
 *  taxlabels (like in the guenomu MCMC algrithm). Calls initialize_reconciliation_sp_count() automatically. */
void 
index_sptaxa_to_reconciliation (char_vector species, char_vector gene, reconciliation rec)
{/* ASSUMES that largest species names already sorted (so that e.g. "ecoli" will match only if "ecoliII" didn't) */
  int i, j, n_index = gene->nstrings, *index;

  index = (int*) biomcmc_malloc (n_index * sizeof (int));
  for (i=0; i < n_index; i++) { 
    rec->sp_id[i] = -1; /* initialize mapping */ 
    index[i] = i;  /* scan gene leaves _without_ replacement */
  }

  for (i=0; i < species->nstrings; i++) for (j=0; j < n_index; j++) /* search sp name in all unmapped gene names */
    if ((gene->nchars[index[j]] >= species->nchars[i]) && /* ordered species names (by nchars[]) */
        (strcasestr (gene->string[index[j]], species->string[i]))) { 
      /* found species name within gene name; we have a mapping */
      rec->sp_id[ index[j] ] = i;
      index[j] = index[--n_index]; // with this the whole search takes O(N ln N),
      j--; // index[j] is now a distinct element (the last)
    }

  if (n_index) {
    fprintf (stderr, "Couldn't find species for genes:\n");
    for (i=0; i < n_index; i++) fprintf (stderr, " \"%s\"\n", gene->string[index[i]]);
    biomcmc_error ("gene names should contain the name of ordered species");
  }
  if (index) free (index);

  initialize_reconciliation_sp_count (rec, species->nstrings, gene->nstrings);
}

// splitset_distances.c


void
prepare_split_from_topologies (topology t1, topology t2, splitset split, int recycle_t1)
{
  int i;
  if ((!recycle_t1) && (!t1->traversal_updated)) update_topology_traversal (t1); /* recycle_t1 > 0 --> we just used t1 in prev iteration */
  if (!t2->traversal_updated) update_topology_traversal (t2);
  /* the vector elements share a single bitstring size, that is modified by e.g. dSPR calculation */
  bipsize_resize (split->g_split[0]->n, split->g_split[0]->n->original_size); 
  bipsize_resize (split->s_split[0]->n, split->s_split[0]->n->original_size); 

  /* heavy child (more leaves or leaf with larger ID in case of tie) at left, after traversal */
  for (i=0; i < t1->nleaves-3; i++) {
    /* Q: why i from 0 to nleaves-4 ?       A: [nleaves-2] is root; [nleaves-3] is OR
     * left child of root  - in which case right child is leaf and left bipartition will be singleton; OR
     * right child of root - which will have same info as (previously visited) left child and therefore redundant */
    bipartition_copy (split->s_split[i], t2->postorder[i]->split);
    /* TODO this is where tripartition is updated BUT must go to nleaves-2 (both children of root are used)*/
    bipartition_flip_to_smaller_set (split->s_split[i]);
  }
  if (!recycle_t1) for (i=0; i < t1->nleaves-3; i++) {  /* Q: why i from 0 to nleaves-4 ? --> same as above */ 
    bipartition_copy (split->g_split[i], t1->postorder[i]->split);
    /* TODO this is where tripartition is updated BUT must go to nleaves-2 (both children of root are used)*/
    bipartition_flip_to_smaller_set (split->g_split[i]);
  }

  split->n_g = split->n_s = i;
  //for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("::DEBUG::  :: ::S\n");
  if (!recycle_t1) qsort (split->g_split, split->n_g, sizeof (bipartition), compare_splitset_bipartition_increasing);
  qsort (split->s_split, split->n_s, sizeof (bipartition), compare_splitset_bipartition_increasing);
}

void /* this function recreates the subtree spanned by given species on gene tree */
split_add_gene_subtree (splitset split, int taxa)
{
  int j, ndis = 0, size = split->sp0[taxa]->n_ones, last_elem = split->spsize - split->size;
  /* special case 1: this cherry is on leaf of sptree (accessible through sp0 but not trough s_split) */
  bipartition_copy (split->s_split[split->n_s++], split->sp0[taxa]); /* RF distance (arxiv.1210.2665) would stop here */
  if (size < 4) return;
  /* temporarily use agree[] and disagree[] */
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits); 
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits); 
  /* 1) create subtree spanned from gene tree elements from taxa */
  for (j = 0; j < split->n_g; j++) {
    bipartition_AND (split->agree[0], split->sp0[taxa], split->g_split[j], true);
    if ( (split->agree[0]->n_ones > 1) && (split->agree[0]->n_ones < (size-1)) ) { /* internal node of this subtree */
      bipartition_ANDNOT (split->agree[1], split->sp0[taxa], split->g_split[j], true);
      if (bipartition_is_larger (split->agree[0], split->agree[1])) bipartition_copy (split->disagree[ndis++], split->agree[1]);
      else                                                          bipartition_copy (split->disagree[ndis++], split->agree[0]);
    }
  }
  /* 2) maintain only distinct bipartitions, being careful not to apply unmasked operations (the mask is sp0[taxa]) */
  split_remove_duplicates (split->disagree, &(ndis));
  /* 3) copy distinct bipartitions to s_split[]; if more than last_elem, something's wrong (I miscalculated vector size?)*/
  for (j = 0; (j < ndis) && (j < last_elem); j++) bipartition_copy (split->s_split[split->n_s++], split->disagree[j]); 
}

int
dSPR_topology (topology t1, topology t2, splitset split)
{
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  rf_hdist_topology_lowlevel (split, false); // calculate Hdist on full trees (without pruning common subtrees)
  prepare_split_from_topologies (t1, t2, split, false); // prepare bipartitions again
  return dSPR_topology_lowlevel (split);
}

int
dSPR_topology_rf (topology t1, topology t2, splitset split)
{  // will calculate only RF, not Hdist or Hdist_reduced
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  return rf_hdist_topology_lowlevel (split, true); // true -> exit as soon as RF is calculated
}

int
dSPR_topology_hdist (topology t1, topology t2, splitset split)
{
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  return rf_hdist_topology_lowlevel (split, false);
}


