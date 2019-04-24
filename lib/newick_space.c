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

#include "newick_space.h"

newick_space
new_newick_space ()
{
  newick_space nwk;
  nwk = (newick_space) biomcmc_malloc (sizeof (struct newick_space_struct));
  nwk->ntrees = 0;
  nwk->t = NULL;
  nwk->ref_counter = 1;
  return nwk;
}

void
del_newick_space (newick_space nwk)
{
  if (!nwk) return;
  if (--nwk->ref_counter) return;
  int i;
  for (i = nwk->ntrees - 1; i >= 0; i--) del_topology (nwk->t[i]);
  free (nwk);
}

newick_space
new_newick_space_from_file (char *filename)
{
  newick_space nwk = new_newick_space ();
  update_newick_space_from_file (nwk, filename);
  return nwk;
}

topology
new_single_topology_from_newick_file (char *filename)
{
  FILE *file;
  char *line_read=NULL, *str_tree_start=NULL, *str_tree_end=NULL;
  size_t linelength = 0, str_tree_size = 0;
  newick_tree tree;
  topology topol;

  file = biomcmc_fopen (filename, "r");
  biomcmc_getline (&line_read, &linelength, file);
  str_tree_start = strchr (line_read, '(');
  str_tree_end   = strchr (str_tree_start, ';');
  if (str_tree_end == NULL) str_tree_size = strlen (str_tree_start); 
  else str_tree_size = str_tree_end - str_tree_start;
  line_read[str_tree_size] = '\0'; /* not null-terminated by default */
  tree  = new_newick_tree_from_string (line_read); // this might increase size, but not decrease (I hope)
  if (line_read) free (line_read);
  topol = new_topology (tree->nleaves);
  copy_topology_from_newick_tree (topol, tree);
  del_newick_tree (tree);
  return topol;
}

void
update_newick_space_from_file (newick_space nwk, char *filename)
{
  FILE *file;
  char *line=NULL, *line_read=NULL, *str_tree_start=NULL, *str_tree_end=NULL;
  size_t linelength = 0, str_tree_size = 0;
 // TODO: newick may span several lines
  if (!nwk) biomcmc_error ("The calling function should allocate memory for newick_space_struct"); 
  file = biomcmc_fopen (filename, "r");
  /* line_read can't be changed (no line_read++ etc.) otherwise we lose first position and can't free() it */
  while (biomcmc_getline (&line_read, &linelength, file) != -1) {
    line = remove_nexus_comments (&line_read, &linelength, file); // will remove comments inside tree
    str_tree_start = strchr (line, '(');
    while (str_tree_start) { /* possible to store several newick files in one line */
      str_tree_end = strchr (str_tree_start, ';');
      if (str_tree_end == NULL) str_tree_size = strlen (str_tree_start); /* function strchrnul() does this, but may not be portable? */
      else str_tree_size = str_tree_end - str_tree_start;
      update_newick_space_from_string (nwk, str_tree_start, str_tree_size);
      if (str_tree_end) str_tree_start = strchr (str_tree_end, '(');
      else str_tree_start = NULL;
    }
  }
  fclose (file);
  if (line_read) free (line_read);
}

void
update_newick_space_from_string (newick_space nwk, char *tree_string, size_t string_size)
{
  char *local_string;
  newick_tree tree;
  topology topol;

  /* read string into a newick_tree; better to work on copy since we could have several trees per line */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (string_size + 1));
  strncpy (local_string, tree_string, string_size + 1); /* adds '\0' only when long_string is smaller!! */
  local_string[string_size] = '\0'; /* not null-terminated by default */
  tree  = new_newick_tree_from_string (local_string);
  if (local_string) free (local_string);
  topol = new_topology (tree->nleaves);
  copy_topology_from_newick_tree (topol, tree);
  del_newick_tree (tree);

  nwk->t = (topology*) biomcmc_realloc ((topology*) nwk->t, sizeof (topology) * (nwk->ntrees + 1));
  nwk->t[nwk->ntrees++] = topol;
  return;
}

void
update_newick_space_from_topology (newick_space nwk, topology topol)
{
  nwk->t = (topology*) biomcmc_realloc ((topology*) nwk->t, sizeof (topology) * (nwk->ntrees + 1));
  nwk->t[nwk->ntrees++] = topol;
  topol->ref_counter++;
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
      this_size  = number_of_leaves_in_newick (&needle_tip); /* false = do not attempt to change tree string */
      if (this_size) { size += this_size; ntrees++; }
    }
  }
  fclose (seqfile);
  if (line_read) free (line_read);

  if (!ntrees) return -1;
  return size/ntrees;
}

