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

START_TEST(file_match_function)
{
  int i,j;
  clock_t time0, time1;
  char *txt[] = {"CCTACAAAGATTAAA\0", "AAAACTATAC\0", "AAAACTATACAAAA\0"};
  suffix_tree st;
  st_matches match;
  alignment aln;

  time0 = clock ();
  memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
  aln = read_alignment_from_file (filename);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  st = new_suffix_tree (aln->character->string[0], aln->character->nchars[0], false); // true=copy text

  for (i = 0; i < 3; i++) {
    match = new_st_matches_from_pattern (txt[i], st);
    printf ("[%s]\nmatch %d\tpartial: %d\tn_matches = %d \t length = %d\n", txt[i], i, match->is_partial, match->n_idx, match->length);
    for (j = 0; j < match->n_idx; j++) printf ("[%7d]\t %.20s\n", match->idx[j], aln->character->string[0] + match->idx[j]);
    del_st_matches (match);
  }
  time1 = clock (); printf ("  time to find matches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_suffix_tree (st);
//  del_alignment (aln); // FIXME
  if (false) ck_abort_msg ("dummy");
}
END_TEST

Suite * suffix_tree_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("SuffixTree");

  tc_case = tcase_create("suffix tree-based string search");
  tcase_add_test(tc_case, file_match_function);
  suite_add_tcase(s, tc_case);

  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (suffix_tree_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}
