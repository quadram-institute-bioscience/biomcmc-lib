void
accumulate_kmers_from_dna (char *dna, int dna_length, void (*reduce)(kmerhash, void *), void *reduce_params)
{ 
  int i;
  kmerhash kmer = new_kmerhash();
  uint64_t hash_f = 0UL, hash_r = 0UL;

//  for (i = 0; i < 15; i++) fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r); 
  // must call two functions: once at beginning, using all bits, and another only newest kmers every step of for()
  for (i = 15; i < dna_length; i++) {
    fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r);
    (*reduce)(kmer, reduce_params); // external function decides what to do with new kmers
  }
  del_kmerhash (kmer);
}

static void
fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr)
{
  *hf = *hf << 4 | dna_in_4_bits[dnachar][0]; // forward
  *hr = *hf >> 4 | ((uint64_t)(dna_in_4_bits[dnachar][1]) << 60UL); // reverse
}
