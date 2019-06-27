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

void debug_simple_minhash_functions (void);

int main(void)
{

 printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour\n");
 debug_simple_minhash_functions();
  return TEST_SKIPPED;
}

void
debug_simple_minhash_functions (void)
{
  int i,j,k;
  cm_sketch *cm;
  double dist[8];
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();
  memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
  aln = read_alignment_from_file (filename);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  cm = (cm_sketch*) biomcmc_malloc (aln->ntax * sizeof (cm_sketch));
  for (i=0; i < aln->ntax; i++) cm[i] = new_fixedhash_sketch_from_dna (aln->character->string[i], aln->character->nchars[i], 1024);
  time1 = clock (); printf ("  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) {
    compare_cm_sketches (cm[i], cm[j], dist);
    printf ("\n%40s %40s ", aln->taxlabel->string[j], aln->taxlabel->string[i]); 
    for (k=0;k<8;k++) printf ("%12.8lf ", dist[k]);
  }
  time1 = clock (); printf ("\n  time to compare sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i= aln->ntax -1; i >- 0; i--) del_cm_sketch (cm[i]); 
  if (cm) free (cm);
  del_alignment (aln);
}
