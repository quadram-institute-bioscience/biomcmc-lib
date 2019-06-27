/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "hashtable.h"

/* Aux functions for hashtable of strings */
uint32_t hash(hashtable ht, const char* str);
inline uint32_t hash1 (hashtable ht);
inline uint32_t hash2 (hashtable ht);
/* Aux functions for hashtable of bipartitions */
inline uint32_t biphash1 (bip_hashtable ht);
inline uint32_t biphash2 (bip_hashtable ht);


hashtable 
new_hashtable (int size)
{
  int i; 
  double x;
  hashtable ht;
  
  ht = (hashtable) biomcmc_malloc (sizeof (struct hashtable_struct));

  /* Setting hashtable to be size of power of 2. This is required so that proper open addressing can occur, i.e. 
   * eventually every entry in the table is considered.
   */
  x = ceil (log (size) / log (2));
  size = pow (2, x+1);

  ht->size = size;
  
  /* allocate memory for table */
  ht->table = (hashtable_item*) biomcmc_malloc(ht->size * sizeof (hashtable_item)); 

  /* initialize to NULLS */
  for(i = 0; i < ht->size; i++) ht->table[i] = 0; 
  ht->P = 2147483647; /* initialize P (large prime)*/
  ht->probelength = 0;
  
  /*initialize hash1 and hash2 variables*/
  srand (time(0)); /* the GSL library would be overkill... */
  ht->a1 = rand() * (ht->P - 1) + 1;
  ht->a2 = rand() * (ht->P - 1) + 1;
  ht->b1 = rand() * ht->P;
  ht->b2 = rand() * ht->P;

  ht->ref_counter = 1; /* at least one place (the calling structure/function) is using this hashtable */

  return ht;
}

void 
del_hashtable (hashtable ht)
{
  if (ht) {
    int i;
    if (--ht->ref_counter) return; /* some other place is using this hashtable, we cannot delete it yet */
    for(i=ht->size-1; i>=0; i--) if (ht->table[i]) { 
      if(ht->table[i]->key) free (ht->table[i]->key); 
      free (ht->table[i]); 
    }
    if (ht->table) free (ht->table);
    free (ht);
  }
}

/* This is called before hash1() and hash2(). It sets the  value of h to be used by hash1() and hash2() and saves 
 * excessive calls to the hash function. */
uint32_t 
hash (hashtable ht, const char* key) 
{
  uint32_t g;
  
  ht->h = 0;
  while (*key) {
    ht->h = (ht->h << 4) + *key++;
    g = ht->h & 0xF0000000L;
    if (g) ht->h ^= g >> 24;
    ht->h &= ~g;
  }
  return (ht->h % ht->P); /*P_ is a large prime*/
}

uint32_t 
hash1 (hashtable ht) 
{
  return ((((ht->a1 * ht->h) + ht->b1) % ht->P) % ht->size) ;
}

uint32_t 
hash2 (hashtable ht) 
{
  return ((((ht->a2 * ht->h + ht->b2) % ht->P) % (ht->size - 3)) | 1);
}

void 
insert_hashtable (hashtable ht, char* key, int value) 
{
  uint32_t i;
  int h1, h2;
  
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);

  ht->probelength = 0;
  
  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!strcmp (ht->table[i]->key, key)) /* key was already inserted */
      return; 
  }
  /* alloc space for new key */
  ht->table[i] = biomcmc_malloc (sizeof (struct hashtable_item_struct));
  ht->table[i]->key = (char*) biomcmc_malloc((strlen (key)+1) * sizeof (char));
  strcpy(ht->table[i]->key, key);
  ht->table[i]->value = value;
  return;
}

int 
lookup_hashtable (hashtable ht, char* key) 
{
  uint32_t i;
  int h1, h2;
  
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);
  
  ht->probelength = 0;

  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!(ht->table[i])) return -1;
    else if ( (ht->table[i]) && (!strcmp (ht->table[i]->key, key)) ) 
      return ht->table[i]->value;
  }
  return -2;
}

/* * Functions to work with hashtable of bipartitions * */

bip_hashtable 
new_bip_hashtable (int size)
{
  int i; 
  double x;
  bip_hashtable ht;
  
  ht = (bip_hashtable) biomcmc_malloc (sizeof (struct bip_hashtable_struct));

  /* Setting hashtable to be size of power of 2. This is required so that proper open addressing can occur, i.e. 
   * eventually every entry in the table is considered.  */
  x = ceil (log (size) / log (2));
  size = pow (2, x+1);
  ht->size = size;
  
  /* allocate memory for table */
  ht->table = (bip_hashitem*) biomcmc_malloc(ht->size * sizeof (bip_hashitem)); 
  /* initialize to NULLS */
  for(i = 0; i < ht->size; i++) ht->table[i] = NULL; 
  ht->P = 2147483647; /* initialize P (large prime)*/
  ht->probelength = 0;
  ht->maxfreq = 1; // will store count of most frequent bipartition, to normalize frequency 
  
  /*initialize hash1 and hash2 variables*/
  srand (time(0)); /* the GSL library would be overkill... */
  ht->a1 = rand() * (ht->P - 1) + 1;
  ht->a2 = rand() * (ht->P - 1) + 1;
  ht->b1 = rand() * ht->P;
  ht->b2 = rand() * ht->P;

  ht->ref_counter = 1; /* at least one place (the calling structure/function) is using this hashtable */

  return ht;
}

void 
del_bip_hashtable (bip_hashtable ht)
{
  if (ht) {
    int i;
    if (--ht->ref_counter) return; /* some other place is using this hashtable, we cannot delete it yet */
    for(i=ht->size-1; i>=0; i--) if (ht->table[i]) { 
      if(ht->table[i]->key) del_bipartition (ht->table[i]->key); 
      free (ht->table[i]); 
    }
    if (ht->table) free (ht->table);
    free (ht);
  }
}

uint32_t biphash1 (bip_hashtable ht) { return ((((ht->a1 * ht->h) + ht->b1) % ht->P) % ht->size) ; }

uint32_t biphash2 (bip_hashtable ht) { return ((((ht->a2 * ht->h + ht->b2) % ht->P) % (ht->size - 3)) | 1); }

void 
bip_hashtable_insert (bip_hashtable ht, bipartition key) 
{
  uint32_t i;
  int h1, h2;
  
  ht->h = bipartition_hash (key) % ht->P; /*P_ is a large prime*/
  h1 = biphash1 (ht);
  h2 = biphash2 (ht);

  ht->probelength = 0;
  
  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (bipartition_is_equal (ht->table[i]->key, key)) {  // already observed
      if (ht->maxfreq < ++ht->table[i]->count) ht->maxfreq = ht->table[i]->count; // notice the 'count++' doing main work 
      return; 
    }
  }
  /* alloc space for new key */
  ht->table[i] = biomcmc_malloc (sizeof (struct hashtable_item_struct));
  ht->table[i]->key = new_bipartition_copy_from (key); 
  ht->table[i]->count = 1;
  return;
}

double
bip_hashtable_get_frequency (bip_hashtable ht, bipartition key) 
{
  uint32_t i;
  int h1, h2;
  
  ht->h = bipartition_hash (key) % ht->P; /*P_ is a large prime*/
  h1 = biphash1 (ht);
  h2 = biphash2 (ht);
  
  ht->probelength = 0;

  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!(ht->table[i])) return 0.;
    else if ( (ht->table[i]) && (bipartition_is_equal (ht->table[i]->key, key)) ) return (double)(ht->table[i]->count)/(double)(ht->maxfreq);
  }
  return -1.;
}

/* * * collection of hash functions from internet 
 * Sources:
 * http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
 * http://burtleburtle.net/bob/hash/integer.html 
 * http://www.cse.yorku.ca/~oz/hash.html
 * https://lemire.me/blog/2018/08/15/fast-strongly-universal-64-bit-hashing-everywhere/
 * http://www.concentric.net/~Ttwang/tech/inthash.htm * * */

uint32_t 
biomcmc_hashint_1 (uint32_t a)
{
  a -= (a<<6); a ^= (a>>17); a -= (a<<9); a ^= (a<<4); a -= (a<<3); a ^= (a<<10); a ^= (a>>15);
  return a;
}

uint32_t 
biomcmc_hashint_2 (uint32_t a)
{
  a += ~(a<<15); a ^= (a>>10); a += (a<<3); a ^= (a>>6); a += ~(a<<11); a ^= (a>>16);
  return a;
}

uint32_t 
biomcmc_hashint_3 (uint32_t a)
{
  a = (a+0x479ab41d) + (a<<8); a = (a^0xe4aa10ce) ^ (a>>5); a = (a+0x9942f0a6) - (a<<14); 
  a = (a^0x5aedd67d) ^ (a>>3); a = (a+0x17bea992) + (a<<7);
  return a;
}

uint32_t 
biomcmc_hashint_4 (uint32_t a)
{
  a = (a^0xdeadbeef) + (a<<4); a = a ^ (a>>10); a = a + (a<<7); a = a ^ (a>>13);
  return a;
}

uint32_t 
biomcmc_hashint_5 (uint32_t a)
{ /* must use use at least the 17 lowest bits */
  a = a ^ (a>>4); a = (a^0xdeadbeef) + (a<<5); a = a ^ (a>>11);
  return a;
}

/* hashCodes that differ only by constant multiples at each bit position have a bounded number of collisions 
 * (approximately 8 at default load factor). */
uint32_t 
biomcmc_hashint_6 (uint32_t a)
{
  a ^= (a >> 20) ^ (a >> 12); a = a ^ (a >> 7) ^ (a >> 4);
  return a;
}

uint32_t 
biomcmc_hashint_7 (uint32_t key)
{
  /* (~a + (a << K)) means ((a << K) - a - 1) and (a * 2057) means (a + (a << 3) + (a << 11)) */ 
  key = ~key + (key << 15); key = key ^ (key >> 12); key = key + (key << 2);
  key = key ^ (key >> 4);   key = key * 2057;        key = key ^ (key >> 16);
  return key;
}

uint32_t 
biomcmc_hashint_8 (uint32_t  a)
{
  a = (a+0x7ed55d16) + (a<<12); a = (a^0xc761c23c) ^ (a>>19); a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);  a = (a+0xfd7046c5) + (a<<3);  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

uint32_t 
biomcmc_hashint_9 (uint32_t a)
{
  a = (a ^ 61) ^ (a >> 16); a = a + (a << 3); a = a ^ (a >> 4); a = a * 0x27d4eb2d; a = a ^ (a >> 15);
  return a;
}

int ulx_size = 12;
uint64_t ulx[] = {
  0x65d200ce55b19ad8UL, 0x4f2162926e40c299UL, 0x162dd799029970f8UL, // 0...2
  0x68b665e6872bd1f4UL, 0xb6cfcf9d79b51db2UL, 0x7a2b92ae912898c2UL, // 3...5
  0xff51afd7ed558ccdUL, 0xc4ceb9fe1a85ec53UL, 0x87c37b91114253d5UL, 0x4cf5ad432745937fUL, // 6...9 (murmurhash) 
  0x52dce729UL, 0x38495ab5UL}; // 10...11 

uint64_t 
biomcmc_hashint64_1 (uint64_t key) /* 64 bits */
{
  /* ((key + (key << 3)) + (key << 8)) = key * 265 and ((key + (key << 2)) + (key << 4)) = key * 21 */
  key = (~key) + (key << 21); key = key ^ (key >> 24); key = (key + (key << 3)) + (key << 8);
  key = key ^ (key >> 14); key = (key + (key << 2)) + (key << 4); key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

uint64_t 
biomcmc_hashint64_2 (uint64_t x) /* 64 bits, splits into 2 32bits blocks */
{
  uint32_t low = x;
  uint32_t high = x >> 32UL;
  return ((ulx[0] * low + ulx[1] * high + ulx[2]) >> 32) | ((ulx[3] * low + ulx[4] * high + ulx[5]) & 0xFFFFFFFF00000000UL);
}

void
biomcmc_hashint64_to_vector (uint64_t x, uint32_t *out) /* 64 bits, splits into two 32bits blocks */
{
  uint32_t low = x;
  uint32_t high = x >> 32UL;
  out[0] = (ulx[0] * low + ulx[1] * high + ulx[2]) >> 32;
  out[1] = (ulx[3] * low + ulx[4] * high + ulx[5]) >> 32;
  out[2] =  (low ^ 61) ^ (low >> 16);  out[2] += (out[2] << 3); out[2] ^= (out[2] >> 4); out[2] *= 0x27d4eb2d; out[2] ^= (out[2] >> 15);
  out[3] = (high ^ 61) ^ (high >> 16); out[3] += (out[3] << 3); out[3] ^= (out[3] >> 4); out[3] *= 0x27d4eb2d; out[3] ^= (out[3] >> 15);
}

uint64_t 
biomcmc_hashint64_seed (uint64_t h, int seed) /* based on murmur (mixer algo); original uses seed = 6 */
{
  int i = seed % (ulx_size - 1); 
  h ^= h >> 33; h *= ulx[i];
  h ^= h >> 33; h *= ulx[i+1];
  h ^= h >> 33;
  return h;
}

uint32_t 
biomcmc_hashint_64to32_seed (uint64_t x, int seed) /* 64 bits */
{ // published algo uses seed=0 or seed=3
  uint32_t low = x;
  uint32_t high = x >> 32UL;
  int i = seed % (ulx_size - 2); 
  return (uint32_t) ((ulx[i] * low + ulx[i+1] * high + ulx[i+2]) >> 32);
}

uint32_t 
biomcmc_hashint_64to32 (uint64_t key) /* input is 64 bits but output is 32 bits (good for hash, not for RNG...) */
{
  key = (~key) + (key << 18); key = key ^ (key >> 31); key = key * 21;
  key = key ^ (key >> 11);    key = key + (key << 6);  key = key ^ (key >> 22);
  return (uint32_t) key;
}

uint32_t 
biomcmc_hashint_mix (uint32_t a, uint32_t b, uint32_t c) 
{
  a=a-b;  a=a-c;  a=a^(c >> 13);  b=b-c;  b=b-a;  b=b^(a << 8);   c=c-a;  c=c-b;  c=c^(b >> 13);
  a=a-b;  a=a-c;  a=a^(c >> 12);  b=b-c;  b=b-a;  b=b^(a << 16);  c=c-a;  c=c-b;  c=c^(b >> 5);
  a=a-b;  a=a-c;  a=a^(c >> 3);   b=b-c;  b=b-a;  b=b^(a << 10);  c=c-a;  c=c-b;  c=c^(b >> 15);
  return c;
}

uint32_t
biomcmc_hashstring_1 (unsigned char *str) /* CRC algortihm */
{ /* ASCII chars use 5 or 6 bits out of the 8 available, thus this function worries more about the lowest 5bits */
  uint32_t c, highorder, hash = 0;
  while ((c = *str++)) {
    highorder = hash & 0xf8000000;  /*extract high-order 5 bits from h */
    /*shift left by 5 bits, move the highorder 5 bits to the low-order end and XOR into h */
    hash = ((hash << 5) ^ (highorder >> 27)) ^ c; 
  }
  return hash;
}

uint32_t 
biomcmc_hashstring_2 (unsigned char *str) /* PJW algorithm */
{ /* ASCII chars use 5 or 6 bits out of the 8 available, thus this function worries more about the lowest 5bits */
  uint32_t c, g, hash = 0;
  while ((c = *str++)) {
    // The top 4 bits of h are all zero
    hash = (hash << 4) + c; // shift h 4 bits left, add c
    g = hash & 0xf0000000;  // get the top 4 bits of h
    if (g != 0) hash = (hash ^ (g >> 24)) ^ g; // if the top 4 bits aren't zero, move them to the low end of h 
  }
  return hash;
}

uint32_t 
biomcmc_hashstring_3 (unsigned char *str) /*dbj2 algorithm */
{
  uint32_t c, hash = 5381;
  while ((c = *str++)) hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  return hash;
}

uint32_t 
biomcmc_hashstring_4 (unsigned char *str) /* dbj2 algorithm with XOR */
{
  uint32_t c, hash = 5381;
  while ((c = *str++)) hash = ((hash << 5) + hash) ^ c; /* (hash * 33) XOR c */
  return hash;
}

uint32_t 
biomcmc_hashstring_5 (unsigned char *str) /* sdbm algorithm */
{ /* hash(i) = hash(i - 1) * 65599 + str[i] */
  uint32_t c, hash = 0;
  while ((c = *str++)) hash = c + (hash << 6) + (hash << 16) - hash;
  return hash;
}

uint32_t 
bipartition_hash (bipartition bip) 
{ // assumes bipartition is flipped to smaller set
  int i; 
  uint32_t h = biomcmc_hashint_9 ((uint32_t) bip->n_ones); // inspired by FNV hash
  for (i=0; i < bip->n->ints; i++) h = biomcmc_hashint_5 (h) ^ biomcmc_hashint_64to32 ( bip->bs[i]);
  return h;
}

/*** MurmurHash3 from https://github.com/PeterScott/murmur3/  written originally by Austin Appleby and CC0 ***/

void biomcmc_murmurhash3 ( const void * key, const int len, const uint32_t seed, void * out )
{  /* out[] is 128 bits (4 x uint_32, for instance) */
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 16;
  int i;
  uint64_t h1 = seed, h2 = seed;
  const uint64_t * blocks = (const uint64_t *)(data);

  for(i = 0; i < nblocks; i++) {
    uint64_t k1 = blocks[i*2+0];
    uint64_t k2 = blocks[i*2+1];
    k1 *= ulx[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx[9]; h1 ^= k1;
    h1 = (h1 << 27) | (h1 >> 37); h1 += h2; h1 = h1 * 5 + ulx[10];
    k2 *= ulx[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx[8]; h2 ^= k2;
    h2 = (h2 << 31) | (h2 >> 33); h2 += h1; h2 = h2 * 5 + ulx[11];
  }

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  uint64_t k1 = 0, k2 = 0;

  switch(len & 15) {
    case 15: k2 ^= (uint64_t)(tail[14]) << 48; break;
    case 14: k2 ^= (uint64_t)(tail[13]) << 40; break;
    case 13: k2 ^= (uint64_t)(tail[12]) << 32; break;
    case 12: k2 ^= (uint64_t)(tail[11]) << 24; break;
    case 11: k2 ^= (uint64_t)(tail[10]) << 16; break;
    case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;  break;
    case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;  break;
             k2 *= ulx[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx[8]; h2 ^= k2; break;

    case  8: k1 ^= (uint64_t)(tail[ 7]) << 56; break;
    case  7: k1 ^= (uint64_t)(tail[ 6]) << 48; break;
    case  6: k1 ^= (uint64_t)(tail[ 5]) << 40; break;
    case  5: k1 ^= (uint64_t)(tail[ 4]) << 32; break;
    case  4: k1 ^= (uint64_t)(tail[ 3]) << 24; break;
    case  3: k1 ^= (uint64_t)(tail[ 2]) << 16; break;
    case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;  break;
    default: k1 ^= (uint64_t)(tail[ 0]) << 0; // equiv to "case 1"  
             k1 *= ulx[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx[9]; h1 ^= k1; break;
  };

  h1 ^= len; h2 ^= len;
  h1 += h2; h2 += h1;

  h1 ^= h1 >> 33; h1 *= ulx[6]; h1 ^= h1 >> 33; h1 *= ulx[7]; h1 ^= h1 >> 33;
  h2 ^= h2 >> 33; h2 *= ulx[6]; h2 ^= h2 >> 33; h2 *= ulx[7]; h2 ^= h2 >> 33;

  h1 += h2; h2 += h1;
  ((uint64_t*)out)[0] = h1;
  ((uint64_t*)out)[1] = h2;
}
 // see also https://github.com/torvalds/linux/blob/master/lib/xxhash.c
