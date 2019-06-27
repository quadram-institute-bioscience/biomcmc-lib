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

START_TEST(fixedhash_sketch_small_function)
{
  int i,j;
  char dna[] = "AAGGCCTTAGTCTGTGTCACACGTGTGTGTGTGTACACACACACACACACACCCCTCTCTCTCTCTCTC";
  cm_sketch cm = new_fixedhash_sketch_from_dna (dna, strlen(dna), 16);
  for (i = 0; i < cm->size; i++) {
    for (j = 0; j < 8; j++) printf ("%12.8lf ", (double)(cm->freq[j][i])/(double)(cm->count));
    printf ("\n");
  }
  del_cm_sketch(cm);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

START_TEST(fixedhash_sketch_alignment_function)
{
  int i,j;
  cm_sketch *cm;
  double dist[8];
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();
  memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
  aln = read_alignment_from_file (filename);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  cm = (cm_sketch*) biomcmc_malloc (aln->ntax * sizeof (cm_sketch));
  for (i=0; i < aln->ntax; i++) cm[i] = new_fixedhash_sketch_from_dna (aln->character->string[i], aln->character->nchars[i], 64);
  time1 = clock (); printf ("  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) compare_cm_sketches (cm[i], cm[j], dist);

  time1 = clock (); printf ("  time to compare sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i= aln->ntax -1; i >- 0; i--) del_cm_sketch (cm[i]); 
  if (cm) free (cm);
  del_alignment (aln);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

Suite * minhash_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("MinHash");

  tc_case = tcase_create("fixed hash skecth");
  tcase_add_test(tc_case, fixedhash_sketch_small_function);
  tcase_add_test(tc_case, fixedhash_sketch_alignment_function);
  suite_add_tcase(s, tc_case);

  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (minhash_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}
