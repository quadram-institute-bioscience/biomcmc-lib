/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 */

#include "rolling_hash.h"

static uint8_t  rand_hash_list_size = 128; /*! \brief hardcoded table size (*must* be power of 4) */ 
static uint32_t rand_hash_list[] = { // 16 bits 
  0xe02f, 0x2076, 0xea7b, 0x8547, 0x1c49, 0x211b, 0x7af3, 0x5460, 0x3e49, 0xc657, 0xa0e7, 0x169c, 0x26c2, 0x04e9, 0xcaa4, 0x88d0,
  0x5ce8, 0xa00c, 0x5c21, 0xcdf2, 0x024a, 0xbac6, 0xc8ac, 0x0a76, 0x973d, 0x5fd7, 0x79aa, 0x99cf, 0xbe46, 0x2e28, 0x4ff0, 0x4c33,
  0xbee5, 0x7ef3, 0xd911, 0x7b59, 0xa574, 0xe5d7, 0xecf4, 0xcada, 0x79ac, 0xea92, 0xbcb8, 0x19b3, 0x0998, 0xab8f, 0xfdac, 0x399b,
  0x3dea, 0x30c4, 0x00a9, 0xbe33, 0x32b8, 0x46df, 0x2931, 0x99e5, 0x5dc6, 0xe750, 0xc8cf, 0x2a15, 0xffdb, 0x2f1f, 0x99a1, 0x8c84, 
  0xcc11, 0xb0d5, 0x2123, 0x981f, 0x50c2, 0x7afd, 0x6300, 0xa28e, 0x9c25, 0xbfc7, 0x3f26, 0x32b5, 0xc3f1, 0x95f9, 0x8acd, 0x800e, 
  0x8c64, 0xd492, 0xeba4, 0xda95, 0x8e6c, 0x56f3, 0x74d8, 0x75b8, 0x9b7a, 0xf584, 0x4a39, 0x86ca, 0x879a, 0x3bea, 0x92f6, 0xdd71, 
  0xd93c, 0x8a10, 0x0166, 0xb109, 0x269e, 0xe172, 0xc726, 0x4bcf, 0xc317, 0xa53c, 0x6a31, 0x616c, 0x9acf, 0x4ba0, 0x5519, 0xc128, 
  0x2f8b, 0xc0d5, 0x2721, 0x1335, 0xdd9e, 0x1593, 0x37d1, 0xd22f, 0x82e5, 0x612f, 0x2a5c, 0xc027, 0x3b54, 0xb025, 0x7d9b, 0xf47e};
  
static uint8_t  prime_salt_list_size = 256; /*! \brief hardcoded table size */
static uint32_t prime_salt_list[] = { // 128 x 16 bits  +  128 x 24 bits
  0x0,    0x6d,   0x7877, 0xde7d, 0xa937, 0x14cb, 0xdea7, 0xfe8f, 0xb99,  0x8ad7, 0xb7ef, 0x5387, 0xf475, 0x94d,  0x191b, 0xe37d,
  0x6511, 0xd7cb, 0xc307, 0xdf2b, 0xaf6d, 0x521f, 0xeb61, 0x59db, 0xc2e3, 0x8fab, 0x575f, 0x9e0b, 0xd321, 0xb07d, 0x399b, 0xb435, 
  0x5897, 0xa111, 0xd459, 0x4549, 0xca75, 0x3b7b, 0x20c5, 0x900d, 0x658f, 0x321d, 0x16cf, 0xe34d, 0x2acd, 0x2cf,  0x3bb7, 0xdab, 
  0xf7eb, 0xa535, 0x6983, 0xee09, 0x3b9f, 0xf8ab, 0xa321, 0xe84f, 0x3fd3, 0xcb03, 0x263b, 0x377f, 0x56d5, 0xff8b, 0x1c19, 0x3eb9, 
  0xfef3, 0xb633, 0x607,  0xae05, 0x9715, 0x95b7, 0xcae5, 0x97eb, 0x4f67, 0x21d,  0xc133, 0x6109, 0xf697, 0x7b01, 0x72a3, 0xdf13, 
  0xd12f, 0x7e2f, 0x3cdd, 0xdd5,  0xd7b3, 0xa93,  0x9793, 0x2329, 0x58ff, 0x1e25, 0xa6c1, 0x759d, 0xa931, 0xc7f9, 0x29d5, 0x4bbf, 
  0xbc77, 0x892d, 0x72ef, 0x277f, 0x7669, 0x908b, 0x4d6b, 0x293f, 0x10f3, 0x7d63, 0x7243, 0x5047, 0x1ea1, 0xf47,  0x8c27, 0x838f, 
  0x7d3f, 0xb417, 0x169f, 0xfcfb, 0xd9b5, 0x1d9f, 0xe4b,  0x9185, 0xc775, 0xd939, 0x380b, 0xc767, 0x4c55, 0x8ce3, 0x70ab, 0x606b,
  0x3df1c9, 0x424859, 0x8fea9,  0x1ec5e3, 0x2c6c5,  0x8fce75, 0x5fa37d, 0x87223,  0x96ea4d, 0x5d106d, 0x13321d, 0x200e39, 0x3a3975, 0x60d5d,  0x184901, 0x550273, 
  0x32b1a1, 0x83ffe7, 0x4cc11d, 0x71741b, 0x3fe293, 0x1a7b9b, 0x463bd7, 0x68ec85, 0x37ab87, 0x15c2bd, 0x319103, 0x392567, 0x285fe3, 0x748a9b, 0x56ae9,  0x1a835b, 
  0x808169, 0x92fa0b, 0x1005c5, 0x7b656b, 0x464c39, 0x3623ff, 0x87c8c9, 0x1fd6bd, 0x4d6df7, 0x785dd,  0x91d051, 0x2b0643, 0x3eb669, 0x65a3e3, 0x441f3b, 0x381a6d, 
  0x293893, 0x895a09, 0x7d8d5d, 0x4a0345, 0x613c9f, 0x696e09, 0x1d6cbf, 0x1f0f87, 0x949603, 0x5991cb, 0x581ec3, 0x170033, 0x7801df, 0x195d49, 0x783d4b, 0x74ca4f, 
  0x2f7f05, 0xe884d,  0xe68c9,  0x5163c7, 0x2d4dbf, 0x246ddb, 0x2a749,  0xe4e05,  0x1e42dd, 0x21c391, 0x9166e1, 0x280efd, 0x6a05d3, 0x5a4adb, 0x631d03, 0x48e8f7, 
  0x11ad07, 0x9114a3, 0xf3f41,  0x3cf47,  0x3cab41, 0x10b0d7, 0x62e40b, 0x469949, 0x31f9bd, 0x21a567, 0x5fe02b, 0x1bd897, 0x15805,  0x37b547, 0x57dfdf, 0x8fd5c9, 
  0x24119f, 0x855979, 0x339853, 0x792a21, 0x6f7cd,  0x17a703, 0x805c63, 0xaad01,  0x8fa569, 0x4d1dc1, 0x2a809f, 0x61d611, 0x83578f, 0x446e79, 0x7f7f4d, 0x20ea5d, 
  0x47c545, 0x3e07b9, 0x2d85f3, 0x1ffa3f, 0x267e59, 0x4e212b, 0x16bd43, 0x68acf7, 0x92b48d, 0x80aca7, 0x65f191, 0x3164ad, 0xd1f6b,  0x106951, 0x30a847, 0x493d7b};

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

  tbl1 = (salt & ((rand_hash_list_size >> 2) - 1)) << 2;  // 4 * ( salt % (list_size/4) )
  tbl2 = (int) (salt / (rand_hash_list_size >> 2)) % prime_salt_list_size;
  /** now we transform the indexes for their equiv. hash values; all have same salt */
  for (j = 0; j < 2; ++j) {
    i = 255; do { shash[j][i] = rand_hash_list[shash[j][i] + tbl1] + prime_salt_list[tbl2]; } while (i--);
  }
  return shash;
}

void
del_dna_salted_hash_encoding (uint32_t** shash) {
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

