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

// OBS: topology_space object should always have <distinct> set, and <freq> scaled to one 

#include "topology_space.h"

/*! \brief Allocates memory for topology_space_struct (set of trees present in nexus file).  */
topology_space new_topology_space (void);

/*! \brief Reads tree from file and store in treespace. 
 *
 * Comment to myself about the external hashtable:
 * Given a hashtable with names, we want to number tree taxlabels according to these hash values. 
 * After checking if all names from taxlabel have a corresponding hash key, this function will create a vector with 
 * the position, in hash, of elements of taxlabel. This is what I call a mapping (!!). For instance, if we have the 
 * hashtable and taxlabels, the mapping is given as in the following:
 * \verbatim
 *     hash["A"] = 0       taxlabel[0] = "C"             order[0] = 2 
 *     hash["B"] = 1       taxlabel[1] = "B"   lead to   order[1] = 1
 *     hash["C"] = 2  and  taxlabel[2] = "D"   mapping   order[2] = 3
 *     hash["D"] = 3       taxlabel[3] = "E"             order[3] = 4
 *     hash["E"] = 4       taxlabel[4] = "A"             order[4] = 0
 * \endverbatim
 * Using this ordering, all trees belonging to a topology_space_struct will be relabeled by the mapping. 
 * Notice than even in the same topology_space_struct distinct trees may have distinct leaf label orders, 
 * despite the taxlabel[] vector with leaf names is shared.
 * This is necessary since despite newick_tree_struct has freedom about the order of nodes (including leaves), in 
 * topology_struct the order is defined. */
void add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, hashtable external_hash, int **order, double tree_weight, bool use_root_location);

/*! \brief Reads translation table (one line) of the form "number = taxa name" in tree file. */
void translate_taxa_topology_space (topology_space tsp, char *string, hashtable external_hash);

/*! \brief string with original file name, with extension stripped -- be caredul not to overwrite it on program */
void store_filename_in_topology_space (topology_space tre, char *filename);

// TODO: create from_newick_space_to_topol_space

bool
is_file_nexus_tree_file (char *seqfilename) 
{
  FILE *seqfile;
  char *line = NULL, *line_read = NULL;
  size_t linelength = 0;
  int is_nexus = 0, i;

  seqfile = biomcmc_fopen (seqfilename, "r");

  /* Search for evidence of nexus file; if obligatory syntax "#NEXUS" is found, see if it is an alignment */
  for (i=0; (is_nexus < 3) && (i < 256) && (biomcmc_getline (&line_read, &linelength, seqfile) != -1);) {
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_string (line)) { 
      if ((!is_nexus) && strcasestr (line, "#NEXUS")) is_nexus++; 
      else if ((is_nexus == 1) && strcasestr (line, "BEGIN") && strcasestr (line, "TREES")) is_nexus++; 
      else if ((is_nexus > 1) && (strcasestr (line, "TRANSLATE") || strcasestr (line, "TREE"))) is_nexus++; 
    }
    if (!is_nexus) i++; /* give up if "#NEXUS" is not found in first 256 lines; otherwise keep looking */
  }

  fclose (seqfile);
  if (line_read) free (line_read);

  if (is_nexus == 3) return true;
  return false;
}

/*! \brief Auxiliary function for the python module */
void
add_string_with_size_to_topology_space (topology_space *tsp, char *long_string, size_t string_size, bool use_root_location)
{
  char *local_string;
  int i, index;
  newick_tree tree;
  topology topol;
  /* read string into a nexus tree */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (string_size + 1));
  strncpy (local_string, long_string, string_size + 1); /* adds '\0' only when long_string is smaller!! */
  local_string[string_size] = '\0'; /* not null-terminated by default */
  tree  = new_newick_tree_from_string (local_string);
  if (local_string) free (local_string);

  /* create taxlabels (shared across trees) or prepare nexus tree leaves to be reordered */
  if (!(*tsp)) {
    (*tsp) = new_topology_space();
    /* tsp->taxlabel will point to names of first tree. This info will be also available at the hashtable */ 
    (*tsp)->taxlabel = new_char_vector (tree->nleaves);
    (*tsp)->taxlabel_hash = new_hashtable (tree->nleaves);
    for (i=0; i< tree->nleaves; i++) { 
      char_vector_link_string_at_position ((*tsp)->taxlabel, tree->leaflist[i]->taxlabel, i); 
      insert_hashtable ((*tsp)->taxlabel_hash, (*tsp)->taxlabel->string[i], i);
      tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
      tree->leaflist[i]->id = i;
    }
  }
  else {
    if ((*tsp)->taxlabel->nstrings != tree->nleaves) {
      del_newick_tree (tree); del_topology_space (*tsp); 
      biomcmc_error ( "All trees from nexus file must have same number of leaves\n");
    }
    for (i=0; i< tree->nleaves; i++) { /* nexus tree leaflist IDs must follow hashtable */
      index = lookup_hashtable ((*tsp)->taxlabel_hash, tree->leaflist[i]->taxlabel);
      if (index < 0) {
        del_newick_tree (tree); del_topology_space (*tsp);
        biomcmc_error ( "Leaf names are not the same across all trees in nexus file\n");
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }
  /* now topology is ready to receive information from nexus tree */
  topol = new_topology (tree->nleaves);
  copy_topology_from_newick_tree (topol, tree, false); // false=don't copy taxlabels from newick_tree 
  topol->taxlabel = (*tsp)->taxlabel; /* taxlabel is shared among all topologies */
  (*tsp)->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */
  del_newick_tree (tree);
  /* comparisons below _assume_ that trees share a char_vector so lines above important */
  add_topology_to_topology_space_if_distinct (topol, (*tsp), 1., use_root_location); // 1. is this tree's weight
  return;
}

topology_space
read_topology_space_from_file (char *seqfilename, hashtable external_taxhash, bool use_root_location)
{
  return read_topology_space_from_file_with_burnin_thin (seqfilename, external_taxhash, 0, 1, use_root_location);
}

topology_space
read_topology_space_from_file_with_burnin_thin (char *seqfilename, hashtable external_taxhash, int burnin, int thin, bool use_root_location)
{
  topology_space treespace=NULL;
  FILE *seqfile;
  char *line=NULL, *line_read=NULL, *needle_tip=NULL;
  bool option_begin_trees    = false,
       option_translate_perm = false,
       option_translate_temp = false,
       option_include_tree = false;
  size_t linelength = 0;
  double freq_sum = 0., this_tree_weight;  /* posterior frequency per tree (if present) */
  int i=0, j=0, iteration = 1, *order_external = NULL; // leaves will follow external_taxhash if exists (malloc'ed by add_tree_to_topology_space())
  if (burnin < 0) burnin = 0;
  if (thin < 1) thin = 1;

  seqfile = biomcmc_fopen (seqfilename, "r");
  
  /* the variable *line_read should point always to the same value (no line++ or alike) */
  biomcmc_getline (&line_read, &linelength, seqfile);
  line = remove_nexus_comments (&line_read, &linelength, seqfile);
  while (!nonempty_string (line)) {
    /* skip (possibly not following NEXUS format) initial comments and blank lines */
    if (biomcmc_getline (&line_read, &linelength, seqfile) < 0) 
      biomcmc_error ("Premature end of NEXUS tree file %s\n", seqfilename);
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
  }
  if (!strcasestr (line, "NEXUS")) 
    biomcmc_error ( "%s is not ot a Nexus tree file (first line should be \"#NEXUS\")\n", seqfilename);

  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    /* read frequency ('posterior distribution) information, in mrbayes' .trprobs format -> before remove_comments */
    needle_tip = line_read;
    this_tree_weight = 1.;
    if ((iteration > burnin) && !(iteration%thin)) option_include_tree = true;
    else option_include_tree = false;

    if ( option_include_tree && (needle_tip = strcasestr (needle_tip, "TREE")) && 
         (needle_tip = strrchr (needle_tip, '=')) &&
         (sscanf (needle_tip, "= [ &W %lf ]", &this_tree_weight) != 1) ) this_tree_weight = 1.; 

    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (nonempty_string (line)) { /* don't do anything untill 'BEGIN TREES' block */
      if ((!option_begin_trees) && (strcasestr (line, "BEGIN TREES"))) {
        option_begin_trees = true;
        treespace = new_topology_space ();
      } 

      else if (!option_translate_temp) {/* check if we need to translate taxon names; in any case see if we have trees to read */
        if (strcasestr (line, "TRANSLATE")) {
          option_translate_perm = true;
          option_translate_temp = true;
        }
        else if (strcasestr (line, "TREE") && (needle_tip = strcasestr (line, "="))) {
          needle_tip++; /* remove "=" from string */
          iteration++;
          if (option_include_tree)
            add_tree_to_topology_space (treespace, needle_tip, option_translate_perm, external_taxhash, &order_external, this_tree_weight, use_root_location);
        }
      }
    
      if (option_translate_temp) {/* we are reading translation table token <-> taxlabel */  
        translate_taxa_topology_space (treespace, line, external_taxhash);
        if (strchr (line, ';')) option_translate_temp = false;
      }
    } // if (line)
  } //while (biomcmc_getline)

  if (external_taxhash) {
    treespace->taxlabel_hash = external_taxhash;
    external_taxhash->ref_counter++; /* since we are sharing the hashfunction */
    // reorder taxlabels to conform to hashtable 
    char_vector_reorder_strings_from_external_order (treespace->taxlabel, order_external);
  }

  fclose (seqfile);
  if (order_external) free (order_external);
  if (line_read) free (line_read);

  /* char_vector_remove_empty_strings() also updates string lengths into taxlabel->nchars so we call it anyway */
  if (char_vector_remove_empty_strings (treespace->taxlabel)) /* this problem should have been detected before... */
    biomcmc_error ("empty taxon names in nexus tree file (reading problem or wrong/duplicate numbers in translate)");
  /* each branch value is sum over all trees with same topol; must be scaled to mean value */
  for (i=0; i < treespace->ndistinct; i++) for (j=0; j < treespace->distinct[i]->nnodes; j++)
    treespace->distinct[i]->blength[j] /= treespace->freq[i]; /* weighted avge = sum{W.x} / sum{W} */
  /* Now we store frequency values (over all trees) in treespace->freq[], for weighted or unweighted files */
  for (i=0; i < treespace->ndistinct; i++) freq_sum += treespace->freq[i];
  for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] /= freq_sum; /* normalize to one */

  store_filename_in_topology_space (treespace, seqfilename);
  return treespace;
}

void
merge_topology_spaces (topology_space ts1, topology_space ts2, double weight_ts1, bool use_root_location)
{ /* ts1->tree is not correct anymore, should not be used after calling this function */
  int i, j, *idx, n_idx = ts1->ndistinct; 
  double total_freq = 0.;

  // TODO: check if taxlabel_hash is the same; branch lengths; mark where convergence could go.
  if (weight_ts1 <= 0.) weight_ts1 = 1.; // usually ts1->ntrees/ts2->ntrees

  idx = (int*) biomcmc_malloc (ts1->ndistinct * sizeof (int)); 
  for (i=0; i < ts1->ndistinct; i++) {
    ts1->freq[i] *= weight_ts1; // weight (number of trees) of ts1 relative to ts2
    idx[i] = i; /* index of trees from ts1 not compared to ts2 yet */
  }

  for (j=0; j < ts2->ndistinct; j++) {
    int found_id = -1;
    /* comparison includes root location (faster than unrooted since uses hash) */ 
    for (i=0; (i < n_idx) && (found_id < 0); i++) if (topology_is_equal (ts2->distinct[j], ts1->distinct[idx[i]] )) found_id = i;
    if ((!use_root_location) && (found_id < 0)) { /* if they look distinct (different root), then do slower unrooted calculation */ 
      for (i=0; (i < n_idx) && (found_id < 0); i++) if (topology_is_equal_unrooted (ts2->distinct[j], ts1->distinct[idx[i]], true)) found_id = i;
    }
    if (found_id >= 0) {
      ts1->freq[ idx[found_id] ] += ts2->freq[j];
      idx [found_id] = idx[--n_idx]; // assuming each tree from ts1 can be found at most once in ts2 
    } // if tree not found 
    else { // tree ts2->distinct[j] is unique to ts2 
      int new_id = ts1->ndistinct++;
      ts2->distinct[j]->id = new_id;
      ts1->freq =     (double*)   biomcmc_realloc ((double*)   ts1->freq,     sizeof (double) * (ts1->ndistinct));
      ts1->distinct = (topology*) biomcmc_realloc ((topology*) ts1->distinct, sizeof (topology) * (ts1->ndistinct));
      ts1->freq[new_id] = ts2->freq[j];
      ts1->distinct[new_id] = ts2->distinct[j];
      ts2->distinct[j] = NULL;
      for (i=0; i < ts1->distinct[new_id]->nleaves; i++) {
        /* the leaf bipartitions never change, so can be shared among all topologies */
        ts1->distinct[new_id]->nodelist[i]->split = ts1->distinct[0]->nodelist[i]->split;
        ts1->distinct[0]->nodelist[i]->split->ref_counter++;
      }
    }
  } // for j in ts2->ndistinct

  for (i=0; i < ts1->ndistinct; i++) total_freq += ts1->freq[i];
  for (i=0; i < ts1->ndistinct; i++) ts1->freq[i] /= total_freq;
  if (idx) free (idx);
}

/*
void 
sort_topology_space_by_frequency(topology_space tsp, double *external_freqs) 
{ 
  double part_sum = 0., freq, *local_freqs = tsp->freq, *pivot_d;
  topology pivot_t;
  char *stree;
  empfreq_double efd;
  if (external_freqs) local_freqs = external_freqs; 
  efd = new_empfreq_double_sort_decreasing (local_freqs, tsp->ndistinct);
  // FIXME: stopped here: must change and freq[] if external is NULL. tree[i] is ponter to distinct, wont change
} */

void
save_topology_space_to_trprobs_file (topology_space tsp, char *filename, double credible)
{ // works even if tsp->freq is unscaled (i.e. total counts and not frequencies summing to one)
  int i, idx;
  FILE *stream;
  double part_sum = 0., freq, *scaled_freqs;
  char *stree;
  empfreq_double efd;

  if (credible > 1.) credible = 1.;

  stream = biomcmc_fopen (filename, "w");
  fprintf (stream, "#NEXUS\n[While frequency 'p' is unscaled, 'P' and 'W' are scaled by credible=%.4lf]\n", credible);
  fprintf (stream, "\n\nBegin trees;\n Translate\n");
  fprintf (stream, "\t1  %s", tsp->taxlabel->string[0]);
  for (i=1; i < tsp->distinct[0]->nleaves; i++) fprintf (stream, ",\n\t%d  %s", i+1, tsp->taxlabel->string[i]);
  fprintf (stream, "\n;\n");

  scaled_freqs = (double*) biomcmc_malloc (sizeof (double) * tsp->ndistinct);
  for (i = 0; i < tsp->ndistinct; i++) part_sum += tsp->freq[i];
  for (i = 0; i < tsp->ndistinct; i++) scaled_freqs[i] = tsp->freq[i]/part_sum;

  efd = new_empfreq_double_sort_decreasing (scaled_freqs, tsp->ndistinct);

  part_sum = 0.;
  for (i = 0; (i < tsp->ndistinct) && (part_sum < 1.); i++) {
    idx = efd->d[i].idx; /* element ordered by frequency */
    freq = scaled_freqs[idx] / credible; /* rescaling s.t. new frequencies sum to one (and not to "credible") */
    part_sum += freq;
    stree = topology_to_string_by_id (tsp->distinct[idx], false); /* from topol to newick */
    fprintf (stream, "tree tree_%d \t[p= %.5lf, P= %.5lf] = [&W %.8lf] %s;\n", i, tsp->freq[idx], part_sum, freq, stree);
    free (stree);
  }
  fprintf (stream, "\nEnd;\n");

  fclose (stream);
  del_empfreq_double (efd);
  if (scaled_freqs) free (scaled_freqs);
}

int
estimate_treesize_from_file (char *seqfilename)
{
  FILE *seqfile;
  char *line=NULL, *line_read=NULL, *needle_tip=NULL;
  size_t linelength = 0;
  int this_size, size = 0, ntrees = 0; 

  seqfile = biomcmc_fopen (seqfilename, "r");
  /* the variable *line should point always to the same value (no line++ or alike) */
  while ((biomcmc_getline (&line_read, &linelength, seqfile) != -1) && (ntrees < 10)) {
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (strcasestr (line, "TREE") && (needle_tip = strcasestr (line, "="))) {
      needle_tip++; /* remove "=" from string */
      this_size  = number_of_leaves_in_newick (&needle_tip); 
      if (this_size) { size += this_size; ntrees++; }
    }
  }
  fclose (seqfile);
  if (line_read) free (line_read);

  if (!ntrees) return -1;
  return size/ntrees;
}

topology_space
new_topology_space (void) 
{
  topology_space tsp;
  
  tsp = (topology_space) biomcmc_malloc (sizeof (struct topology_space_struct));

  tsp->ndistinct = tsp->ntrees = 0;
  tsp->is_rooted = true;  /* if true then tree comparison is straightforward; o.w. we must reroot to same leaf (to compare branches) */
  tsp->tree = tsp->distinct = NULL;
  tsp->taxlabel = NULL;
  tsp->taxlabel_hash = NULL;
  tsp->filename = NULL;
  tsp->distinct = NULL;
  tsp->freq = NULL; /* will be created online (as we read more trees) or when reading "[&W 0.01]" posterior probability data */
  /* tsp->distinct vector is increased by add_tree_to_topology_space() while tsp->tree are pointers to dsp->distinct
   * tsp->taxlabel vector is setup by translate_taxa_topology_space() or add_tree_to_topology_space(), whichever is used
   * tsp->taxlabel_hash is setup by one of the above in absence of global external_hash -- will be replaced by a pointer to an external hastable (from
   * alignment or another tree file, for instance) */
  return tsp;
}

void
del_topology_space (topology_space tsp) 
{
  if (tsp) {
    int i;
    if (tsp->distinct) { 
      for (i=tsp->ndistinct-1; i>=0; i--) del_topology (tsp->distinct[i]);
      free (tsp->distinct);
    }
    if (tsp->tree)     free (tsp->tree);
    if (tsp->freq)     free (tsp->freq);
    if (tsp->filename) free (tsp->filename);
    del_hashtable (tsp->taxlabel_hash);
    del_char_vector (tsp->taxlabel);
    free (tsp);
  }
}

void
add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, hashtable external_hash, int **order, double tree_weight, bool use_root_location)
{
  int i, index, *original_order;
  char *local_string;
  newick_tree tree;
  topology topol;

  /* use local copy of string to avoid problems with biomcmc_getline() */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (strlen (string) + 1));
  strcpy (local_string, string);
  local_string[strlen(string)] = '\0'; /* not null-terminated by default */
  tree  = new_newick_tree_from_string (local_string); // also copies string, but relies on '\0' which may not be present 
  if (local_string) free (local_string);

  if ((tsp->ntrees == 0) && (!translate)) { /* CASE 1: first tree read and no TRANSLATE command in nexus file */
    /* tsp->taxlabel will point to names of first tree. This info will be also available at the hashtable */ 
    tsp->taxlabel = new_char_vector (tree->nleaves); 

    for (i=0; i< tree->nleaves; i++) { 
      // leaf names of first newick_tree will be shared by topol_space
      char_vector_link_string_at_position (tsp->taxlabel, tree->leaflist[i]->taxlabel, i); 
      tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
      tree->leaflist[i]->id = i;
    }

    if (external_hash) {
      /* leaves should be renumbered according to external taxlabel_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int)); /* two vecs, second is temporary */
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          del_newick_tree (tree); del_topology_space (tsp);
          biomcmc_error ( "tree label %s not found in sequence data\n", tsp->taxlabel->string[i]); 
        }
      }
    }
    else { /* no global (external) hash, must create local one */
      tsp->taxlabel_hash = new_hashtable (tree->nleaves);
      for (i=0; i< tree->nleaves; i++) insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
    }
  }

  else if ((tsp->ntrees == 0) && (translate)) { /* CASE 2: first tree read and TRANSLATE command in nexus file */
    for (i=0; i< tree->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index < 0 || index >= tree->nleaves) {
        del_newick_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE \n", index);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); // we only care about the node ID 
      tree->leaflist[i]->id = index;
    }

    if (external_hash) { //same as before, remembering that taxlabel_hash doesn't exist 
      /* leaves should be renumbered according to external_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int));
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          del_newick_tree (tree); del_topology_space (tsp);
          biomcmc_error ( "tree label %s not found in external hash table with mapped names (from alignment, generally)\n", tsp->taxlabel[i]); 
        }
      }
    }
  }

  else if ((tsp->ntrees > 0) && (!translate)) { /* CASE 3: not the first tree read and no TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_newick_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between trees of same file\n");
    }
    /* use hashtable to check if names are consistent and point all leaves to taxlabel vector */
    for (i=0; i< tree->nleaves; i++) {
      if (external_hash) index = lookup_hashtable (external_hash, tree->leaflist[i]->taxlabel);
      else          index = lookup_hashtable (tsp->taxlabel_hash, tree->leaflist[i]->taxlabel);
      if (index < 0) {
        del_newick_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf names disagree between trees of same file\n");
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  else {             /*  CASE 4: not the first tree read and TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_newick_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between tree and TRANSLATE command\n");
    }
    for (i=0; i< tree->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index < 0 || index >= tree->nleaves) {
        int nleaveslocal = tree->nleaves;
        del_newick_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE -- NTAX = %d\n", index + 1, nleaveslocal);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  if (external_hash) {
    original_order = (*order) + tsp->taxlabel->nstrings;
    for (i=0; i < tree->nleaves; i++) original_order[i] = tree->leaflist[i]->id;
    for (i=0; i < tree->nleaves; i++) tree->leaflist[i]->id = (*order)[ original_order[i] ];
  }

  topol = new_topology (tree->nleaves);
  copy_topology_from_newick_tree (topol, tree, false);
  topol->taxlabel = tsp->taxlabel; /* taxlabel is shared among all topologies */
  tsp->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */
  del_newick_tree (tree);
  add_topology_to_topology_space_if_distinct (topol, tsp, tree_weight, use_root_location);
  return;
}

void
add_topology_to_topology_space_if_distinct (topology topol, topology_space tsp, double tree_weight, bool use_root_location)
{
  int i, found_id = -1;
  tsp->tree = (topology*) biomcmc_realloc ((topology*) tsp->tree, sizeof (topology) * (tsp->ntrees+1));

  /* distinct trees might have same unrooted info, but rooted calculation is faster */
  for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal (topol, tsp->distinct[i])) found_id = i;
  if ((!use_root_location) && (found_id < 0)) /* if they look distinct (different root), then do slower unrooted calculation */ 
    for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal_unrooted (topol, tsp->distinct[i], true)) found_id = i;
  // better solution is to have vector of splits for each tree, since slowest part is to copy/order bipartitions 

  if (found_id >= 0) { 
    tsp->tree[tsp->ntrees] = tsp->distinct[found_id];
    tsp->freq[found_id] += tree_weight;
    for (i = 0; i < topol->nnodes; i++) { 
      tsp->distinct[found_id]->blength[i] += (tree_weight * topol->blength[i]); /* weighted mean length */
    }
    del_topology (topol);
  }
  else {
    topol->id = tsp->ndistinct++;
    for (i = 0; i < topol->nnodes; i++) topol->blength[i] *= tree_weight; // weighted average 
    tsp->distinct = (topology*) biomcmc_realloc ((topology*) tsp->distinct, sizeof (topology) * (tsp->ndistinct));
    tsp->freq =     (double*)   biomcmc_realloc ((double*)   tsp->freq,     sizeof (double) * (tsp->ndistinct));
    tsp->distinct[topol->id] = topol;
    tsp->freq[topol->id] = tree_weight;
    tsp->tree[tsp->ntrees] = tsp->distinct[topol->id];
    if (topol->id > 0) for (i=0; i < tsp->distinct[topol->id]->nleaves; i++) {/* leaf bipartitions never change -> shared among all topologies */
      del_bipartition (tsp->distinct[topol->id]->nodelist[i]->split);
      tsp->distinct[topol->id]->nodelist[i]->split = tsp->distinct[0]->nodelist[i]->split;
      tsp->distinct[0]->nodelist[i]->split->ref_counter++;
    }
    if      (!(tsp->ndistinct%10000)) fprintf (stderr, "+");
    else if (!(tsp->ndistinct%1000))  fprintf (stderr, "."); /* the "else" is to avoid printing both */
    fflush (stdout); 
  }
  tsp->ntrees++;
}

/* TODO: create unroot_topol_space() */
void
translate_taxa_topology_space (topology_space tsp, char *string, hashtable external_hash) 
{
  int i, index;
  char *c, *s, *last, *comma_position, seqname[MAX_NAME_LENGTH]="";
  bool good_scan;

  /* the file may have the first token<->taxon_name at the same line as the "TRANSLATE" command. */
  if  ( (s = strcasestr (string, "TRANSLATE")) ) s += strlen ("TRANSLATE");
  else s = (string);
  last = s + strlen (s) +1;

  if (!tsp->taxlabel) tsp->taxlabel = new_char_vector (1); /* we don't know beforehand how many taxa */

  /* one or more "token1 taxon_id1, token2 taxon_id2," */
  while ((comma_position = strchr (s, ',')) && (s<last)) {
    if (strchr (s, '\"')) good_scan = (sscanf (s, " %d \"%[^\"]\",", &index, seqname) == 2); /* double quotes in name */
    else good_scan = (sscanf (s, " %d %[^,]", &index, seqname) == 2); /* leaf name is everything before a comma */ 
    if (!good_scan) biomcmc_error ("could not scan leaf info in TRANSLATE command");
    
    index--; /* in nexus we have index = 1...NTAX */
    if ((!strlen (seqname)) || (index<0)) biomcmc_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    char_vector_add_string_at_position (tsp->taxlabel, seqname, index);
    
    s = comma_position+1;
  }

  /* maybe only one "token taxon_id" (the last line, e.g.) */
  if (strchr (s, '\"')) good_scan = (sscanf (s, " %d \"%[^\"]\" ", &index, seqname) == 2); /* double quotes in name */
  else good_scan = (sscanf (s, " %d %[^,;] ", &index, seqname) == 2); /* leaf name is everything before a comma */
  if (good_scan) {
    good_scan = false; /* recicling the variable (use it to see if we found the end of translation) */
    while ( (c = strrchr (seqname, ';')) ) { 
      *c = '\0';  /* remove possible semicolon at end (even a freaky case of several ones!) */
      good_scan = true; /* A semicolon (even within seqname) means that translation table is completely read */
    }
    
    index--; /* in nexus we have index = 1...NTAX */
    if ((!strlen (seqname)) || (index<0)) biomcmc_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    char_vector_add_string_at_position (tsp->taxlabel, seqname, index);
  }
  
  /* when we finished reading the leaf names (semicolon separated (line below) or not (good_scan) by space) */ 
  if ((strchr (string, ';') || good_scan) && !external_hash) { /* if external_has is present, do not create local one */
    tsp->taxlabel_hash = new_hashtable (tsp->taxlabel->nstrings);
    for (i=0; i<tsp->taxlabel->nstrings; i++) insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
  }
}

void
store_filename_in_topology_space (topology_space tre, char *filename)
{
  char *first_char = NULL, *last_char = NULL;
  size_t size;

  last_char = strrchr (filename, '.'); /* point to last occurence of "." */
  if (last_char) size = last_char - filename + 1; /* plus one for terminating EOL */
  else           size = strlen (filename) + 1;
  first_char = strrchr (filename, '/'); /* point to last occurence of "/" -- the "beginning" of file name removes all directories */
  if (first_char) size -= (first_char - filename); 
  else            first_char = filename;

  tre->filename = biomcmc_malloc (size * sizeof (char));
  memcpy (tre->filename, first_char, size);
  memcpy (tre->filename + size - 1, "\0", 1);
}

