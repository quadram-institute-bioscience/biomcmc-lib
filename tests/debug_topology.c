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
  int i,j;
  memcpy(filename + prefix_size, "iqtree.nwk", 11); // 11 is sizeof file (or larger)
  t = new_single_topology_from_newick_file (filename);
  printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour");
  printf("::iqtree has nleaves =  %d\n And now unchecked topol equality matrix::\n", t->nleaves);
  memcpy(filename + prefix_size, "ortho.nwk", 10);
  nwk_spc = new_newick_space_from_file (filename);
  for (i=0; i<nwk_spc->ntrees; i++) reorder_topology_leaves (nwk_spc->t[i]);
  for (i=0; i<nwk_spc->ntrees; i++) {
    for (j=0; j<nwk_spc->ntrees; j++) printf ("%d\t",(int)(topology_is_equal_unrooted (nwk_spc->t[i], nwk_spc->t[j], false)));
    printf ("\n");
  }

  for (i=0; i< nwk_spc->ntrees; i++) {
    sptre = new_speciestree (nwk_spc->t[i], NULL);
    printf ("tree %5d longest name: %s\n", i, sptre->t->taxlabel->string[ sptre->spnames_order->i[0].idx ]);
    del_speciestree (sptre);
  }


  return TEST_SKIPPED; 
}
