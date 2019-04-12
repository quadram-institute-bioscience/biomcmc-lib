#include <biomcmc.h> 
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif

START_TEST(test_should_work)
{
  ck_assert_int_eq(5, 5);
}
END_TEST

START_TEST(test_should_not_work)
{
  if (0)
    ck_abort_msg ("This never fails hahaha");
}
END_TEST

START_TEST(test_should_not_work2)
{
  ck_assert_msg (1 > 0, "one is larger than zero, duh");

}
END_TEST

Suite * money_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("Money");

  tc_case = tcase_create("Core");
  tcase_add_test(tc_case, test_should_work);
  tcase_add_test(tc_case, test_should_not_work);
  suite_add_tcase(s, tc_case);
  tc_case = tcase_create("Case2");
  tcase_add_test(tc_case, test_should_not_work2);
  suite_add_tcase(s, tc_case);

  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  fprintf (stderr, "TEST_FILE_DIR:\n%s\n", TEST_FILE_DIR); 
  sr = srunner_create (money_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
