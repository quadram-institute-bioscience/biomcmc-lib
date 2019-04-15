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

int main(void)
{
  topology t;
  newick_space nwk_spc;
  speciestree sptre;
  genetree gtre;
  int i,j;
  bool is_it_true = false;
  
  memcpy(filename + prefix_size, "iqtree.nwk", 11); // 11 is sizeof file (or larger)
  t = new_single_topology_from_newick_file (filename);
  printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour\n");
  memcpy(filename + prefix_size, "ortho.nwk", 10);
  nwk_spc = new_newick_space_from_file (filename);
  for (i=0; i<nwk_spc->ntrees; i++) reorder_topology_leaves (nwk_spc->t[i]);
  printf("::iqtree has nleaves =  %d\n Check for char_vector equality (if leaves are identical)::\n", t->nleaves);
  for (i=0; i<nwk_spc->ntrees; i++) for (j=0; j<i; j++) {
    is_it_true = char_vector_link_address_if_identical (&(nwk_spc->t[i]->taxlabel), &(nwk_spc->t[j]->taxlabel)); 
    printf ("leaves from trees %i and %i are %s identical; ", i,j, (is_it_true?"":"not"));
  }
  //printf("\nDEBUG::trees::\n");char *s; for (i=0; i<nwk_spc->ntrees; i++){ s = topology_to_string_by_name (nwk_spc->t[i], NULL); printf ("tree %d\t %s\n", i, s); free (s); }

  printf("\nAnd now unchecked topol equality matrix::\n");
  for (i=0; i<nwk_spc->ntrees; i++) {
    for (j=0; j<nwk_spc->ntrees; j++) printf ("%d\t",(int)(topology_is_equal_unrooted (nwk_spc->t[i], nwk_spc->t[j], false)));
    printf ("\n");
  }

  sptre = new_speciestree (nwk_spc->t[0], NULL); 
  for (i=1; i< nwk_spc->ntrees; i++) {
    gtre = new_genetree (nwk_spc->t[i], sptre);
    genetree_reconcile_speciestree (gtre, sptre); 
    genetree_dSPR_speciestree (gtre, sptre, 2); 
    printf ("%d -> %d %d %d %d %d %d\n", i, gtre->rec->ndups, gtre->rec->nloss, gtre->rec->ndcos, gtre->split->rf, gtre->split->hdist, gtre->split->spr+gtre->split->spr_extra);
    del_speciestree (sptre);
  }

  return TEST_SKIPPED; 
}
