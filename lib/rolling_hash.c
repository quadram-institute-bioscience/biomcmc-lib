/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 */

#include "rolling_hash.h"
#include "random_number_lists.h" // predetermined rnadomised DNA values and semi-ordered "salt"

uint32_t**
new_dna_salted_hash_encoding (uint32_t salt) 
{
  uint8_t i=255, j, tbl1, tbl2; 
  uint32_t **shash = (uint32_t**) biomcmc_malloc (2 * sizeof (uint32_t*)); // opposite order as dna_in_2_bits[]
  for (i = 0; i < 2; ++i) shash[i] = (uint32_t*) biomcmc_malloc (256 * sizeof (uint32_t));
  /* notice do{}while() instead of for() since i is always < 256 */
  i = 255; do { shash[0][i] = shash[1][i] = 4;} while (i--); // anything else is fifth state

  shash[0]['A'] = shash[0]['a'] = 0; shash[1]['A'] = shash[1]['a'] = 3;  /*  A  <-> T  */
  shash[0]['C'] = shash[0]['c'] = 1; shash[1]['C'] = shash[1]['c'] = 2;  /*  C  <-> G  */
  shash[0]['G'] = shash[0]['g'] = 2; shash[1]['G'] = shash[1]['g'] = 1;  /*  G  <-> C  */
  shash[0]['T'] = shash[0]['t'] = 3; shash[1]['T'] = shash[1]['t'] = 0;  /*  T  <-> A  */
  shash[0]['U'] = shash[0]['u'] = 3; shash[1]['U'] = shash[1]['u'] = 0;  /*  U  <-> A  */

  tbl1 = (salt & ((rnd_salt_h16_list_size >> 2) - 1)) << 2;  // 4 * ( salt % (list_size/4) )
  tbl2 = (int) (salt / (rnd_salt_h16_list_size >> 2)) % prime_salt_list_size;
  /** now we transform the indexes for their equiv. hash values; all have same salt */
  for (j = 0; j < 2; ++j) {
    i = 255; do { shash[j][i] = rnd_salt_h16_list[shash[j][i] + tbl1] + prime_salt_list[tbl2]; } while (i--);
  }
  return shash;
}

void
del_dna_salted_hash_encoding (uint32_t** shash) 
{
  if (!shash) return;
  if (shash[1]) free (shash[1]);
  if (shash[0]) free (shash[0]);
  free (shash);
}

#define RoL(val, numbits) ((val) << (numbits)) | ((val) >> (32 - (numbits)))
#define RoR(val, numbits) ((val) >> (numbits)) | ((val) << (32 - (numbits)))
void
roll_hash_add (uint32_t *h, const char dna_base, const uint8_t rol, const uint32_t* shashcode)
{
  *h = RoL(*h, rol);
  *h ^= shashcode[(uint8_t) dna_base]; // XOR h and ki
}

void // kmer_size can't be 0 or 32 
roll_hash_replace_f (uint32_t *h, const char old_base, const char new_base, const uint8_t kmer_size, const uint8_t rol, const uint32_t* shashcode)
{
  uint8_t remain = (rol * (kmer_size-1)) & 31U; // since x % y = x & (y-1) // this can be calculated outside
  *h ^= RoL(shashcode[(uint8_t) old_base], remain); // remove "leftmost" base 
  *h  = RoL(*h, rol);
  *h ^= shashcode[(uint8_t) new_base];
}

void // kmer_size can't be 0 or 32 
roll_hash_replace_r (uint32_t *h, const char old_base, const char new_base, const uint8_t kmer_size, const uint8_t rol, const uint32_t* shashcode)
{
  uint8_t remain = (rol * (kmer_size-1)) & 31U; // since x % y = x & (y-1) // this can be calculated outside 
  *h ^= shashcode[(uint8_t) old_base];
  *h  = RoR(*h, rol);
  *h ^= RoL(shashcode[(uint8_t) new_base], remain);  
}
#undef RoL
#undef RoR

rolling_hash
new_rolling_hash (uint8_t kmer_size, uint32_t salt)
{
  rolling_hash rh = (rolling_hash) biomcmc_malloc (sizeof (struct rolling_hash_struct));
  rh->salted_hashcode = new_dna_salted_hash_encoding (salt); 
  rh->kmer[0] = rh->kmer[1] = 0;
  rh->kmer_size = kmer_size;
  if (kmer_size > 30) rh->kmer_size = 30;
  if (kmer_size < 4) rh->kmer_size = 4;
  rh->rol = (salt >> 3) & 15U + 1; /* rolling window size between 1 and 16 */
// STOPHERE:: create rh->remain 
  rh->canonical= 0; /* default direction is forward (zero), */
  kmer->dna = NULL; /* pointer to DNA sequence */
  kmer->n_dna = 0;  /* DNA sequence length */
  kmer->i = 0;      /* current position on DNA sequence */
  kmer->ref_counter = 1;  /* just in case this struct is shared */
  return kmer;
}

del_rolling_hash (rolling_hash rh)
{
  if (!rh) return;
  if (--rh->ref_counter) return;
  del_dna_salted_hash_encoding (rh->salted_hashcode);
  free (rh);
}

void
link_rolling_hash_to_dna_sequence (rolling_hash rh, char *dna, size_t dna_length)
{
  int j;
  rh->dna = dna;
  rh->n_dna = dna_length; 
  /* generate first k-mer */
  rh->kmer[0] = rh->salted_hashcode[0][(uint8_t) dna[0]];
  rh->kmer[1] = rh->salted_hashcode[1][(uint8_t) dna[rh->kmer_size-1]];
  for (j = 1; j < rh->kmer_size; ++j) {
    roll_hash_add (&(rh->kmer[0]), dna[j],                 rh->rol, rh->salted_hashcode[0]);
    roll_hash_add (&(rh->kmer[1]), dna[rh->kmer_size-1-j], rh->rol, rh->salted_hashcode[1]);
  }
  rh->i = rh->kmer_size - 1; // FIXME special number; now the first iteration already has the k-mers
}

// TODO: created linked list if kmer_size > 30, and possibly with spacer region between them

bool
rolling_hash_iterator (rolling_hash rh)
{

  if (rh->i == rh->n_dna) return false;

  if (rh->i >= rh->kmer_size) {
    roll_hash_replace_f (&(rh->kmer[0]), rh->dna[rh->i - rh->kmer_size], rh->dna[j], rh->kmer_size, rh->rol, shcode2[0]);
    roll_hash_replace_r (&h1R, my_dna_seq[j - kmer_size], my_dna_seq[j], kmer_size, rol, shcode1[1]);
  }

  rh->i = kmer_size; // (or size+1?)

// STOPHERE
  // assume for/rev are full, solve initial case later
  if (kmer->p->dense == 2) { // AT vs GC comparison
    while ((kmer->i < kmer->n_dna) && (dna_in_1_bits[(int)(kmer->dna[kmer->i])][0] > 1)) kmer->i++;
    if (kmer->i == kmer->n_dna) return false;
    dnachar = (int) kmer->dna[kmer->i];
    //ABCD->BCDE: forward [C D] [A B] --> [D E] [B C] , reverse: [b a] [d c] -->  [c b] [e d]
    if (kmer->p->n2) { // only if we need 2 uint64_t elements forward[1] and reverse[0] are used
      kmer->forward[1] = kmer->forward[1] << 1 | kmer->forward[0] >> 63UL;  // forward[1] at left of forward[0]
      kmer->reverse[0] = kmer->reverse[0] >> 1 | kmer->reverse[1] << 63UL;  // reverse[0] is at left
    }
    kmer->forward[0] = kmer->forward[0] << 1 | dna_in_1_bits[dnachar][0]; // forward[0] receives new value
    kmer->reverse[1] = kmer->reverse[1] >> 1 | ((uint64_t)(dna_in_1_bits[dnachar][1]) << 63UL); // reverse[1] receives new value
  }
  else if (kmer->p->dense == 1) {
    while ((kmer->i < kmer->n_dna) && (dna_in_2_bits[(int)(kmer->dna[kmer->i])][0] > 3)) kmer->i++;
    if (kmer->i == kmer->n_dna) return false;
    dnachar = (int) kmer->dna[kmer->i];
    if (kmer->p->n2) { 
      kmer->forward[1] = kmer->forward[1] << 2 | kmer->forward[0] >> 62UL;  // forward[1] at left of forward[0]
      kmer->reverse[0] = kmer->reverse[0] >> 2 | kmer->reverse[1] << 62UL;  // reverse[0] is at left
    }
    kmer->forward[0] = kmer->forward[0] << 2 | dna_in_2_bits[dnachar][0]; // forward[0] receives new value
    kmer->reverse[1] = kmer->reverse[1] >> 2 | ((uint64_t)(dna_in_2_bits[dnachar][1]) << 62UL); // reverse[1] receives new value
  }
  else { // 4 bits (dense==false) can handle non-orthodox DNA characters
    dnachar = (int) kmer->dna[kmer->i];
    if (kmer->p->n2) { 
      kmer->forward[1] = kmer->forward[1] << 4 | kmer->forward[0] >> 60UL; 
      kmer->reverse[0] = kmer->reverse[0] >> 4 | kmer->reverse[1] << 60UL;
    }
    kmer->forward[0] = kmer->forward[0] << 4 | dna_in_4_bits[dnachar][0];
    kmer->reverse[1] = kmer->reverse[1] >> 4 | ((uint64_t)(dna_in_4_bits[dnachar][1]) << 60UL);
  } // else if (dense)
  kmer->i++;
  if (kmer->i < kmer->p->size[0]) kmerhash_iterator (kmer); //do not update hashes, just fill forward[] and reverse[] 

  for (i=0; i < kmer->p->n1; i++) { // fit into one uint64_t 
    if (kmer->i >= kmer->p->size[i]) { 
      hf = kmer->forward[0] & kmer->p->mask1[i]; hr = kmer->reverse[1] >> kmer->p->shift1[i];
      if (hr < hf) hf = hr;
      kmer->kmer[i] = hf;
      kmer->hash[i] = kmer->p->hashfunction (&hf, kmer->p->nbytes[i], kmer->p->seed[i]); // 313 is just a seed
    }
  }

  for (i=0; i < kmer->p->n2; i++) { // need to compare two uint64_t; no kmer->kmer[] possible 
    j = kmer->p->n1 + i;
    if (kmer->i >= kmer->p->size[j]) {
      //ABCDE : forward[0][1]=[DE][BC] reverse[0][1]=[cb][ed] ; 1. compare [DE] and [ed]
      if (kmer->forward[0] < kmer->reverse[1])
        kmer->hash[j] = kmer->p->hashfunction (kmer->forward, kmer->p->nbytes[j], kmer->p->seed[j]);
      else if ((kmer->forward[0] == kmer->reverse[1]) && // 2. compare [_C] and [c_] (where '_' means masked out)
               ((kmer->forward[1] & kmer->p->mask2[i]) < (kmer->reverse[0] >> kmer->p->shift2[i])) )
        kmer->hash[j] = kmer->p->hashfunction (kmer->forward, kmer->p->nbytes[j], kmer->p->seed[j]);
      else { // 3. skip b from reverse[0]
        rev8b = (uint8_t*) kmer->reverse; // so that rev8b[3] is 3 bytes after start of reverse
        kmer->hash[j] = kmer->p->hashfunction (rev8b + kmer->p->shift2[i], kmer->p->nbytes[j], kmer->p->seed[j]);
      }
    } // if dna[i] > min size for this operation
  }
  return true;
}

