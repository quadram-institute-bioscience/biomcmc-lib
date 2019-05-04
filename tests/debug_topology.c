#include <biomcmc.h> 
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif

char filename[2048] = TEST_FILE_DIR; // now we can memcpy() file names _after_ prefix_size
size_t prefix_size = strlen(TEST_FILE_DIR); // all modifications to filename[] come after prefix_size

void debug_simple_tree_functions (void);
void check_distance_estimation (topology tree);

int main(void)
{
  printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour\n");

  debug_simple_tree_functions ();

  return TEST_SKIPPED; 
}

void
debug_simple_tree_functions (void)
{
  newick_space nwk;
  speciestree sptre;
  genetree gtre;
  int i,j;
  char *s; 
  bool is_it_true = false;

  printf ("\n<newick trees with first containing multifurcations>\n");
  memcpy(filename + prefix_size, "multifurcation.nwk", 18);
  nwk = new_newick_space_from_file (filename);
  for (i=0; i<nwk->ntrees; i++) reorder_topology_leaves (nwk->t[i]);
  char_vector_link_address_if_identical (&(nwk->t[0]->taxlabel), &(nwk->t[1]->taxlabel)); 
  printf ("RF distance between multifurcating and binary trees: %d\t",(int)(topology_is_equal_unrooted (nwk->t[0], nwk->t[1], false)));
  s = topology_to_string_by_name (nwk->t[0], nwk->t[0]->blength); 
  printf("for file 'multifurcation.nwk'\n Resolved tree: %s\n\n", s); free (s);
  check_distance_estimation (nwk->t[1]);
  del_newick_space (nwk);

  printf ("\n<orthologous newick file>\n");
  memcpy(filename + prefix_size, "ortho.nwk", 10);
  nwk = new_newick_space_from_file (filename);
  for (i=0; i<nwk->ntrees; i++) reorder_topology_leaves (nwk->t[i]);
  printf("check for char_vector equality (if leaves are identical)::\n");
  for (i=0; i<nwk->ntrees; i++) for (j=0; j<i; j++) {
    is_it_true = char_vector_link_address_if_identical (&(nwk->t[i]->taxlabel), &(nwk->t[j]->taxlabel)); 
    printf ("leaves from trees %i and %i are %sidentical; ", i,j, (is_it_true?"":"not "));
  }
  //printf("\nDEBUG::trees::\n");char *s; for (i=0; i<nwk_spc->ntrees; i++){ s = topology_to_string_by_name (nwk_spc->t[i], NULL); printf ("tree %d\t %s\n", i, s); free (s); }
  printf("\nAnd now unchecked topol equality matrix::\n");
  for (i=0; i<nwk->ntrees; i++) {
    for (j=0; j<nwk->ntrees; j++) printf ("%d\t",(int)(topology_is_equal_unrooted (nwk->t[i], nwk->t[j], false)));
    printf ("\n");
  }

  printf ("\n<genetree/speciestree from this ortho trees>\n");
  sptre = new_speciestree (nwk->t[0], NULL); 
  for (i=1; i< nwk->ntrees; i++) {
    gtre = new_genetree (nwk->t[i], sptre);
    genetree_reconcile_speciestree (gtre, sptre); 
    genetree_dSPR_speciestree (gtre, sptre, 2); 
    printf ("%d -> %d %d %d %d %d %d\n", i, gtre->rec->ndups, gtre->rec->nloss, gtre->rec->ndcos, gtre->split->rf, gtre->split->hdist, gtre->split->spr+gtre->split->spr_extra);
    del_speciestree (sptre);
  }

  del_newick_space (nwk);
  return;
}

void
check_distance_estimation (topology tree)
{
  topology this = new_topology (tree->nleaves);
  double **dist, *scale;
  int i, n_pairs = (tree->nleaves * (tree->nleaves -1))/2;
  int j;
  char *s;

  copy_topology_from_topology (this, tree);
  dist = (double**) biomcmc_malloc (6 * sizeof (double*));
  for (i=0;i<6;i++) dist[i] = (double*) biomcmc_malloc (n_pairs * sizeof (double));
  scale = (double*) biomcmc_malloc (6 * sizeof (double));
  patristic_distances_from_topology_to_vectors (tree, dist, scale, 6, 1e-9);

  for (i=0; i < 6; i++) {
    estimate_topology_branch_lengths_from_distances (this, dist[i]);
    for (j=0;j<tree->nnodes;j++) this->blength[j] *= scale[i];
    s = topology_to_string_by_name (this, this->blength); 
    printf("%s\n", s); free (s);
  }

  if (dist) {
    for (i=5; i >=0; i--) if (dist[i]) free (dist[i]);
    free (dist);
  }
  if (scale) free (scale);
  del_topology (this);
}

