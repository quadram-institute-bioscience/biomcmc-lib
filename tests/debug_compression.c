#include <biomcmc.h> 
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
 file_compress_t fc; 
  char *line_read = NULL;
  size_t linelength = 0;

  printf("DEBUG:: This program does not perform true tests, and serves to debug functions and to show expected behaviour\n");
  if (argc != 2) { fprintf (stderr, "I need one argument\n"); return TEST_SKIPPED; }
  fc = biomcmc_open_compress (argv[1], "r");
  while (biomcmc_getline_compress (&line_read, &linelength, fc) != -1) printf ("<<%s", line_read);
  biomcmc_close_compress (fc);
  return TEST_SKIPPED; 
}
