#include <biomcmc.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99


int main(int argc, char **argv)
{
  clock_t time0, time1;
  int i, n_samples;
  uint64_t x;

  uint8_t algo;

  if (argc == 1) return TEST_SKIPPED;
  time0 = clock ();

  biomcmc_random_number_init (0ULL);
  sscanf (argv[1], " %d ", &n_samples);
  sscanf (argv[2], " %hhd ", &algo);
//  printf ("type: d\ncount: %d\nnumbit: 32\n", n_samples);
  biomcmc_rng_set_algorithm (algo);

//  for (i=0; i < n_samples; i++) printf ("%u\n", (uint32_t)(biomcmc_rng_get () & 0xffffffffU));
  for (i=0; i < n_samples; i++) {
    x = biomcmc_rng_get ();
    fwrite (&x, sizeof (uint64_t), 1, stdout);
  } 
  time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
//  biomcmc_fprintf_fortune (stderr);
//  biomcmc_fprintf_bofh (stderr);

  biomcmc_random_number_finalize();
  return TEST_SKIPPED;
}

