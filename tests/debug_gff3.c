#include <biomcmc.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99


int main(int argc, char **argv)
{
  clock_t time0, time1;
  int i;
  gff3_t g3;
  if (argc == 1) return TEST_SKIPPED;
  time0 = clock ();

  g3 = read_gff3_from_file(argv[1]);
  for (i = 0; i < g3->n_cds; i++) printf ("%5d %32s %32s : %5d\n", i, g3->cds[i]->seqid.str, g3->cds[i]->attr_id.str, g3->cds[i]->start);

  if (g3->seqname) { 
    printf ("n = %d\n", g3->seqname->nstrings);
    for (i = 0; i < g3->seqname->nstrings; i++) printf ("FASTA: %s\n", g3->seqname->string[i]);
  }
//  else printf ("no contig names\n");
  time1 = clock ();  fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  return TEST_SKIPPED;
}

