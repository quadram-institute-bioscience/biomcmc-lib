#include <biomcmc.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99


int main(int argc, char **argv)
{
  int64_t time0[2];
  int i, j, k, n;
  char *s;
  gff3_t g3;
  gff3_fields *hit;
  if (argc == 1) return TEST_SKIPPED;

  biomcmc_get_time (time0); 

  g3 = read_gff3_from_file(argv[1]);
  for (i = 0; i < g3->n_f0; i++) printf ("id=%5d seqid=%46s id=%4d [%5d -%5d] attr=%16s type=%8s\n", i, g3->f0[i].seqid.str, g3->f0[i].seqid.id,
                                          g3->f0[i].start, g3->f0[i].end, g3->f0[i].attr_id.str, g3->f0[i].type.str);

  fprintf (stderr, "timing: %.8f secs\n",biomcmc_update_elapsed_time (time0)); // updates time0 with internal call to biomcmc_get_time

  for (i = 0; i < g3->seqname->nstrings; i++) 
    for (j = 0; j < (int) g3->sequence->nchars[i]; j += 100) {
      hit = find_gff3_fields_within_position (g3, g3->seqname->string[i], j, &n);
      for (k = 0; k < n; k++) printf ("genome:%4d j:%5d start:%5d end:%5d attr:%16s\n",i, j, hit[k].start, hit[k].end, hit[k].attr_id.str);
      if (hit) free (hit); 
    }
  fprintf (stderr, "timing: %.8f secs\n",biomcmc_update_elapsed_time (time0)); // updates time0 with internal call to biomcmc_get_time

  printf ("number of seqnames (contigs/genomes/chromosomes) = %d\n", g3->seqname->nstrings);
  for (i = 0; i < g3->seqname->nstrings; i++) printf ("FASTA: %s\t%ld\n", g3->seqname->string[i], (g3->sequence==NULL? 0:g3->sequence->nchars[i]));
  s = save_fasta_from_gff3 (g3, NULL, true); // true means to overwrite file even if exists
  if (s) free (s);

  fprintf (stderr, "timing: %.8f secs\n",biomcmc_update_elapsed_time (time0)); // updates time0 with internal call to biomcmc_get_time

  return TEST_SKIPPED;
}

