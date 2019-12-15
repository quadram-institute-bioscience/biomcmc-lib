#include <biomcmc.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99


int main(int argc, char **argv)
{
  clock_t time0, time1;
  int i, n_samples;

  if (argc == 1) return TEST_SKIPPED;
  time0 = clock ();

  biomcmc_random_number_init (0ULL);
  sscanf (argv[1], " %d ", &n_samples);
  for (i=0; i < n_samples; i++) printf ("%u\n", biomcmc_rng_get_32());
  time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  biomcmc_random_number_finalize();
  return TEST_SKIPPED;
}

