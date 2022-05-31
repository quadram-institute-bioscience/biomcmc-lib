#include <biomcmc.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99


int main(int argc, char **argv)
{
  clock_t time0, time1;
  uint32_t i, x, salt, salt_final;
  uint32_t *counter;
  uint32_t n_zeroes = 0, n_collisions = 0, large_collision = 0;

  if (argc == 1) return TEST_SKIPPED;

  sscanf (argv[1], " %u ", &salt);
  if (argc == 2) {  // DIEHARDER test
    for (i = 0; i < UINT_MAX; i++) {
      x = biomcmc_hashint_salted (i, salt);
      fwrite (&x, sizeof (uint64_t), 1, stdout);
    }
    return TEST_SKIPPED;  
  }
  /* else; no dieharder */
  sscanf (argv[2], " %u ", &salt_final); // run from salt to salt_final
  if (salt_final < salt) { x = salt_final; salt_final = salt; salt = x; }

  FILE *fp = fopen ("hash_stats.txt", "a");
  if (!fp) biomcmc_error ("could not open stats file");
  counter = (uint32_t*) biomcmc_malloc (UINT_MAX * sizeof (uint32_t));

  for (; salt <= salt_final; salt++) {
    for (i = 0; i < UINT_MAX; i++) counter[i] = 0;
    time0 = clock ();
    for (i = 0; i < UINT_MAX; i++) {
      x = biomcmc_hashint_salted (i, salt);
      counter[x]++;
    }
    time1 = clock (); 
    fprintf (fp, "salt = %5u timing: %.8f secs\t", salt, (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    n_zeroes = 0; n_collisions = 0; large_collision = 0;
    for (i = 0; i < UINT_MAX; i++) {
      if (!counter[i]) n_zeroes++;
      if (counter[i] > 1) n_collisions++;
      if (counter[i] > large_collision) large_collision = counter[i];
    }

    fprintf (fp, "n_zeroes = %6u n_collisions = %6u largest_collision = %6u\n", n_zeroes, n_collisions, large_collision);
    fflush (fp);
  }
  fclose (fp);
  if (counter) free (counter);
  return TEST_SKIPPED;
}

