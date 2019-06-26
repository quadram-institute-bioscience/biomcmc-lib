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


START_TEST(fixedhash_sketch_function)
{
  char dna[] = "AAGGCCTTAGTCTGTGTCACACGTGTGTGTGTGTACACACACACACACACACCCCTCTCTCTCTCTCTC";
  fixedhash_sketch_from_dna (dna);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

Suite * minhash_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("MinHash");

  tc_case = tcase_create("fixed hash skecth");
  tcase_add_test(tc_case, fixedhash_sketch_function);
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
