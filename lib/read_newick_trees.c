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

#include "read_newick_trees.h"

#define DEFAULTBLENGTH 1. /*!< \brief Default branch length. */

typedef struct newick_node_struct* newick_node;
typedef struct newick_tree_struct* newick_tree; /*! \brief newick trees have minimal information, unlike topology_struct */

struct newick_node_struct
{
  newick_node up, right, left; /*! \brief Parent and children nodes. */
  int id;               /*! \brief Initial pre-order numbering of node. */
  double branch_length; /*! \brief Branch length from node to node->up. */
  char *taxlabel;       /*! \brief Leaf sequence name */
};

struct newick_tree_struct
{
  newick_node *nodelist; /*! \brief Vector with pointers to every internal node. */
  newick_node *leaflist; /*! \brief Vector with pointers to tree leaves. */
  newick_node root;      /*! \brief Pointer to root node. */
  bool has_branches;     /*! \brief Boolean saying if tree has branch lengths or not. (topology alsways has, even if one-zero) */
  int nnodes, nleaves;   /*! \brief Number of nodes (including leaves), and number of leaves */
};

/*! \brief Allocates memory for newick_tree_struct. */
newick_tree new_newick_tree (int nleaves);
/*! \brief Frees memory used by tree. */
void del_newick_tree (newick_tree T);
/*! \brief Copy information from newick_tree struct to topology_struct  */
void copy_topology_from_newick_tree (topology tree, newick_tree nwk_tree);
/*! \brief Creates newick_tree structure. */ 
newick_tree new_newick_tree_from_string (char *external_string);
/*! \brief Recursive function that creates a node based on parenthetic structure. */
newick_node subtree_newick_tree (newick_tree tree, char *lsptr, char *rsptr, int *node_id, newick_node up);
/*! \brief Reads leaf name (or number, if translation table is present). */
char* read_taxlabel ( const char *name_start, const char *name_end);
/*! \brief Preorder initialization of leaves. */
void create_leaflist_newick_tree (newick_tree tree, newick_node this, int *id);
/*! \brief Preorder initialization of _internal_ nodes; 'id' should be >= nleaves. */
void create_node_id_newick_tree (newick_node this, int *id);
/*! \brief Searches for (last reported) branch length on string or return default value.  */
double read_branch_length (char *right_string_ptr);
/*! \brief Returns position of innermost comma (divides string into two subtrees). */
int find_branch_split_newick (char *left_string_ptr, char *right_string_ptr);
/*! \brief Counts the number of leaves and resolves (one) trifurcation of tree string. */
int number_of_leaves_in_newick (char **string);
/*! \brief Tries to resolve multifurcations on string, by replacing (,,) for (,(,)) ; heap_depth to avoid infinite recursion */
size_t remove_multifurcations_newick (char **string, size_t i_left, size_t i_right, int heap_depth);
size_t create_new_bifurcation_newick (char **string, size_t i_left, size_t comma_location);

/** global functions **/

newick_tree
new_newick_tree (int nleaves) 
{
  newick_tree tree;
  int i;
  size_t sizeof_node = sizeof (struct newick_node_struct);

  tree = (newick_tree) biomcmc_malloc (sizeof (struct newick_tree_struct));
  tree->nleaves = nleaves;
  tree->nnodes  = 2 * nleaves - 1;
  tree->has_branches = false;
  
  tree->nodelist = (newick_node*) biomcmc_malloc (tree->nnodes * sizeof (newick_node));
  tree->leaflist = (newick_node*) biomcmc_malloc (tree->nleaves * sizeof (newick_node));
  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (newick_node) biomcmc_malloc (sizeof_node);
    tree->nodelist[i]->up = tree->nodelist[i]->right = tree->nodelist[i]->left = NULL;
    tree->nodelist[i]->taxlabel = NULL;
  }
  return tree;
}

void 
del_newick_tree (newick_tree T) 
{
  if (!T) return;
  int i;
  for (i=T->nnodes-1; i >=0; i--) if (T->nodelist[i]) free (T->nodelist[i]);
  if (T->nodelist) free (T->nodelist);
  if (T->leaflist) free (T->leaflist);
  free (T);
}

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

/** local functions (declared here) **/

void
copy_topology_from_newick_tree (topology tree, newick_tree nwk_tree)
{
  int i, node_id;

  tree->taxlabel = new_char_vector (nwk_tree->nleaves);
  for (i=0; i< nwk_tree->nleaves; i++) { 
    char_vector_link_string_at_position (tree->taxlabel, nwk_tree->leaflist[i]->taxlabel, i); 
    nwk_tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
    nwk_tree->leaflist[i]->id = i;
  }

  for (i = 0; i < tree->nnodes; i++) {
    node_id = nwk_tree->nodelist[i]->id;
  
    tree->nodelist[node_id]->mid[0] = tree->nodelist[node_id]->mid[1] = tree->nodelist[node_id]->id = node_id;
    if (tree->blength) tree->blength[node_id] = nwk_tree->nodelist[i]->branch_length; /* should be malloc'ed beforehand */
    if (nwk_tree->nodelist[i]->up)    tree->nodelist[node_id]->up    = tree->nodelist[nwk_tree->nodelist[i]->up->id];
    else                              tree->nodelist[node_id]->up    = NULL;
    if (nwk_tree->nodelist[i]->left)  tree->nodelist[node_id]->left  = tree->nodelist[nwk_tree->nodelist[i]->left->id];
    else                              tree->nodelist[node_id]->left  = NULL; 
    if (nwk_tree->nodelist[i]->right) tree->nodelist[node_id]->right = tree->nodelist[nwk_tree->nodelist[i]->right->id];    
    else                              tree->nodelist[node_id]->right = NULL;
  } // for (nnodes)
  tree->root = tree->nodelist[nwk_tree->root->id];
  for (i = 0; i < tree->nleaves; i++) { tree->nodelist[i]->u_done = false; tree->nodelist[i]->d_done = true; }
  for (i = tree->nleaves; i < tree->nnodes; i++) tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = false;
  update_topology_sisters (tree);
  update_topology_traversal (tree);
}


newick_tree
new_newick_tree_from_string (char *external_string) 
{
  char *lptr, *rptr, *string;
  int id, nleaves;
  size_t string_size = strlen (external_string);
  newick_tree T; 
  
  string = (char*) biomcmc_malloc (sizeof (char) * (string_size + 1));
  strncpy (string, external_string, string_size + 1); /* adds '\0' only when long_string is smaller!! */
  string[string_size] = '\0'; /* not null-terminated by default */
  string = remove_space_from_string (string);
  nleaves = number_of_leaves_in_newick (&string);
  T = new_newick_tree (nleaves);
  if (strchr (string, ':')) T->has_branches = true;
  
  /* begin & end of string */
  lptr = string;  
  rptr = string + strlen (string) - 1;
  
  id = 0; /* This function does the actual creation of the tree */
  T->root = subtree_newick_tree (T, lptr, rptr, &id, NULL);
  id = 0; /* vector of pointers to the tree leaves */
  create_leaflist_newick_tree (T, T->root, &id); 
  create_node_id_newick_tree (T->root, &id); 

  return T;
}

newick_node
subtree_newick_tree (newick_tree tree, char *lsptr, char *rsptr, int *node_id, newick_node up) 
{
  newick_node thisnode;
  
  thisnode = tree->nodelist[*node_id];
  thisnode->up = up;
  thisnode->id = -1; 
  thisnode->branch_length = tree->has_branches ? read_branch_length (rsptr) : DEFAULTBLENGTH;
  thisnode->left = NULL;
  thisnode->right = NULL;
  thisnode->taxlabel = NULL;

  (*node_id)++;

  if (*lsptr == '(') { /* internal node */
    char *newend = rsptr;
    int comma_pos = find_branch_split_newick (lsptr, rsptr);

    thisnode->left = subtree_newick_tree (tree, lsptr+1, lsptr+comma_pos-1, node_id, thisnode);
    while ((newend != lsptr) && (newend != NULL) && (*newend != ')')) newend--;
    if (newend == lsptr) newend = rsptr;
    thisnode->right = subtree_newick_tree (tree, lsptr+comma_pos+1, newend-1, node_id, thisnode);
  }

  else thisnode->taxlabel = read_taxlabel (lsptr, rsptr); /* leaf */

  return thisnode;
}

int
find_branch_split_newick (char *left_string_ptr, char *right_string_ptr) 
{
  int nLevel = 0;
  int comma_position = 0;
  char *treeptr = left_string_ptr;
  
  while ((*treeptr != ',') || (nLevel != 1)) { /* stop when *treeptr == ',' AND nLevel == 1 */
    if (treeptr == right_string_ptr) biomcmc_error ("unbalanced tree: couldn't find innermost comma for subtree");
    if (*treeptr == '(') nLevel++;
    if (*treeptr == ')') nLevel--;
    treeptr++;
    comma_position++;
  }
  return comma_position;
}

char *
read_taxlabel ( const char *name_start, const char *name_end) 
{
  size_t seqsize;
  size_t i;
  char *this_end, *this_start, *label=NULL;
  this_start = (char*) name_start;
  while (isspace(*this_start)) this_start++;
  this_end = this_start;
  while ((this_end <= name_end) && (*this_end != ',') && (*this_end != ')') && (*this_end != ':')) this_end++;
  while (isspace(*this_end)) this_end--; // go backwards only after reaching end (since spaces in middle are fine)
  seqsize = this_end - this_start;
  i = sizeof (char)*(seqsize+1);
  label = (char*) biomcmc_malloc (i);
  label[0] = '\0';
  strncat (label, this_start, seqsize);
  return label;
}

void
create_leaflist_newick_tree (newick_tree tree, newick_node this, int *id) 
{
  if (this->taxlabel != NULL) {
    this->id = (*id);
    tree->leaflist[(*id)++] = this;
  }
  else {
    if (this->left)  create_leaflist_newick_tree (tree, this->left, id);
    if (this->right) create_leaflist_newick_tree (tree, this->right, id);
  }
}

void
create_node_id_newick_tree (newick_node this, int *id) 
{
  this->id = (*id)++;
  if ( (this->left)  && (this->left->id < 0)  ) create_node_id_newick_tree (this->left, id);
  if ( (this->right) && (this->right->id < 0) ) create_node_id_newick_tree (this->right, id);
}

double
read_branch_length (char *right_string_ptr) 
{
   char *backwards = right_string_ptr;
   double branch = DEFAULTBLENGTH;
  
  if ((*backwards == ')') || (*backwards == ',')) return DEFAULTBLENGTH;
  while (*backwards != ':') backwards--;
  sscanf (backwards, ": %lf", &branch);
  if (branch < 0.) return 0.;
  else return branch;
}

void prtdbg (char *string, int start, int end) { // not declared above since just debugging tool
  int i; 
  printf(" <<");
  for (i=start; i<=end;i++) printf ("%c",string[i]); 
  printf(">> ");
}

int 
number_of_leaves_in_newick (char **string)
{
  int nopen = 0, nclose = 0, ncommas = 0, i;
  int has_branches = 0; /* could be a bool, but I want to debug number of branches */
  int len = strlen (*string);
  char current;

  if (*(*string + len - 1) == ';') *(*string + len - 1) = '\0';
  //printf ("\nDEBUG1:: <<%s>>\n", *string); 
  remove_multifurcations_newick (string, 0, (size_t) len, 0); // will resolve all (<2048) polytomies
  //printf ("\nDEBUG2:: <<%s>>\n", *string); 
  len = strlen (*string); // updated
  for (i = 0; i < len; i++) {
    current = (*string)[i];
    if (current == ',' && (nopen - nclose) == 1) ncommas++;
    else if (current == '(') nopen++;
    else if (current == ')') nclose++; 
    else if (current == ':') has_branches++;
  }
  if (nopen != nclose || ncommas > 2 || ncommas < 1) biomcmc_error ("%d %d %d Invalid tree structure n_leaves_newick(): <<%s>>", nopen, nclose, ncommas, *string);
  return nopen + 1;
}

size_t 
remove_multifurcations_newick (char **string, size_t i_left, size_t i_right, int heap_depth)
{ /* returns how many chars it shifted */
  char current;
  size_t i, shifted_chars = 0, *nsplit;
  int nopen = 0, nclose = 0, ncommas = 0;
  if ((*string)[i_left] != '(') return 0;
  if (heap_depth > 2048) { fprintf (stderr, "biomcmc WARNING: Too many multifurcations, I give up!\n"); return 0; }

  nsplit = (size_t*) biomcmc_malloc (2 * sizeof (size_t));
  for (i = i_left; i <= i_right; i++) {
    current = (*string)[i];
    if ((current == ',') && ((nopen - nclose) == 1) && (ncommas < 2)) nsplit[ncommas++] = i; // store two first commas 
    else if (current == '(') nopen++;
    else if (current == ')') nclose++;
  }
  //printf ("\n\tDEBUG::start  this_left_right: %d %d __ %d :: %d %d %d\n", i_left, i_right, heap_depth, nopen, nclose, ncommas);  //DEBUG
  if ((!nopen) || (!ncommas)) {free (nsplit); return 0; } // reached a leaf, although the != '(' above should have catched it...
  if (nopen != nclose || ncommas < 1) biomcmc_error ("<%d %d | %d> Invalid tree remove_multifurcations() : %s\n", nopen, nclose, ncommas, *string);
  if (ncommas > 1) {
    //prtdbg (*string, i_left, nsplit[1]); printf ("\t::DEBUG::%d %d|\t", i_left, nsplit[1]); prtdbg (*string, nsplit[1], i_right);
    shifted_chars = create_new_bifurcation_newick (string, i_left, nsplit[1]);
    //printf ("\nDEBUG::After resolving:: "); prtdbg (*string, i_left, i_right + shifted_chars); printf ("\n");
    i = remove_multifurcations_newick (string, i_left, i_right + shifted_chars, heap_depth+1);
    shifted_chars += i;
  }
  else {// Updating left side will change string size, so when we update right size we must add these new bifurcations 
    size_t new_right = i_right;
    while ((new_right != i_left) && ((*string)[new_right] != ')')) new_right--;
    if (new_right == i_left) new_right = i_right;
    //printf ("::DEBUG::regular:: new_left and right: %d  %d \tcomma: %d   (below are subtrees)\n", i_left, new_right, nsplit[0]);
    //prtdbg (*string, i_left+1, nsplit[0]); printf("\t\t"); prtdbg (*string, nsplit[0], new_right-1); printf ("\n");
    shifted_chars = remove_multifurcations_newick (string, i_left + 1, nsplit[0] - 1, heap_depth); // nsplits[0] has usual comma even in binary trees
    i = remove_multifurcations_newick (string, nsplit[0] + 1 + shifted_chars, new_right - 1 + shifted_chars, heap_depth);
    shifted_chars += i;
  }
  free (nsplit);
  return shifted_chars;
}

size_t
create_new_bifurcation_newick (char **string, size_t i_left, size_t comma_location)
{ 
  char *pivot = NULL, *tstring = NULL, *next = NULL, *fixed[2]={"(", "):0.0"};
  size_t len = strlen (*string), flen0, flen1; /* how many extra chars we add with this function */
  flen0 = strlen(fixed[0]);
  flen1 = strlen(fixed[1]);

  //*string = (char *) biomcmc_realloc ((char*) (*string), sizeof (char) * (len + flen0 + flen1));
  tstring = (char *) biomcmc_malloc (sizeof (char) * (len + flen0 + flen1 + 1));
  memset (tstring, '\0',  sizeof (char) * (len + flen0 + flen1 + 1)); 
  next = tstring;
  memcpy (next, *string, i_left); next += i_left;
  memcpy (next, fixed[0], flen0); next += flen0;
  memcpy (next, *string + i_left, comma_location - i_left); next += comma_location - i_left; 
  memcpy (next, fixed[1], flen1); next += flen1;
  memcpy (next, *string + comma_location, len - comma_location);

  //memcpy (*string, tstring, len + flen0 + flen1); free (tstring);
  pivot = *string;  
  *string = tstring; 
  free (pivot);
  return flen0 + flen1;
}
