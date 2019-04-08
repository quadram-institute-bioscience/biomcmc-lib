

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
