#include <biomcmc.h> 
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
  xz_file_t *f;
  char *line_read = NULL;
  size_t linelength = 0;

  printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour\n");
  if (argc != 2) { fprintf (stderr, "I need one argument\n"); return TEST_SKIPPED; }
  f = biomcmc_xz_open (argv[1], "r", 1024);
  while (biomcmc_getline_xz (&line_read, &linelength, f) != -1) printf ("<<%s", line_read);
  biomcmc_xz_close (f);
  return TEST_SKIPPED; 
}
