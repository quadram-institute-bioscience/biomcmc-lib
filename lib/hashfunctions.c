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

/** collection of hash functions from internet. whenever possible I use a 'salt' number which will indicate the hash
 * function to be used (specially usefull for one-liners). This seems more readable than the hash_1()... hash_10()
 * version implemented in my older libraries (may a bit slower than hardcoding each function?)
 * Sources:
 * http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
 * http://burtleburtle.net/bob/hash/integer.html 
 * http://www.cse.yorku.ca/~oz/hash.html
 * https://lemire.me/blog/2018/08/15/fast-strongly-universal-64-bit-hashing-everywhere/
 * https://github.com/Cyan4973/smhasher
 * https://github.com/PeterScott/murmur3/ 
 * http://www.concentric.net/~Ttwang/tech/inthash.htm 
 **/

#include "hashfunctions.h"
#include "constant_random_lists.h" 

/** \brief large deterministic list of random numbers, that can be used to "salt" hashes etc. index[] can be all zeroes
 * (since idea of this list is to explore all vector combinations, hashing etc should be done outside of it) */
uint32_t biomcmc_get_salt_set_from_spice_table (uint32_t index[], uint32_t *salt, uint32_t salt_length);

#define RoL(val, numbits) (((val) << (numbits)) | ((val) >> (32 - (numbits))))
#define RoR(val, numbits) (((val) >> (numbits)) | ((val) << (32 - (numbits))))
uint32_t 
biomcmc_hashint_salted (uint32_t a, unsigned int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint32_t h = a;
  switch(salt & 15) { // 4 last bits
    case 10: // murmur3 avalanche (from this version)
      a *= 0xcc9e2d51u; a = (a << 15) | (a >> 17); a *= 0x1b873593u; 
      a = 0x499606EDu ^ a; a = (a << 13) | (a >> (19)); a = (a * 5) + 0xe6546b64u; break; // 499606ED is large prime added by leomrtns (should be a seed)
    case 9: // murmur-like  avalanche
      a += 0x5851f4; a ^= (a >> 16); a *= 0x85ebca6b; a ^= (a >> 13); a *= 0xc2b2ae35; a ^= (a >> 16); break;
    case 8: 
      a += 0xe6543b; a -= (a<<6); a ^= (a>>17); a -= (a<<9); a ^= (a<<4); a -= (a<<3); a ^= (a<<10); a ^= (a>>15); break;
    case 7:
      a += ~(a<<15); a ^= (a>>10); a += (a<<3); a ^= (a>>6); a += ~(a<<11); a ^= (a>>16); break;
    case 6:// half-avalanche: Every input bit affects itself and all higher output bits, plus a few lower output bits
      a = (a+0x479ab41d) + (a<<8); a = (a^0xe4aa10ce) ^ (a>>5); a = (a+0x9942f0a6) - (a<<14); 
      a = (a^0x5aedd67d) ^ (a>>3); a = (a+0x17bea992) + (a<<7); break;
    case 5:
      a = (a^0xdeadbeef) + (a<<4); a = a ^ (a>>10); a = a + (a<<7); a = a ^ (a>>13); break;
    case 4: // must use use at least the 17 lowest bits
      a = a ^ (a>>4); a = (a^0xdeadbeef) + (a<<5); a = a ^ (a>>11); break;
    case 3: /* Java hashCodes that differ by multiples at each bit have bounded number of collisions (~8) <- but bad with small numbers */
      a = a * 0x27d9ab + 0xdca2; a ^= (a >> 20) ^ (a >> 12); a = a ^ (a >> 7) ^ (a >> 4); break; /* must be a large value */
    case 2: // (~a + (a << K)) means ((a << K) - a - 1) and (a * 2057) means (a + (a << 3) + (a << 11)) 
      a = ~a + (a << 15); a ^= (a >> 12); a+= (a << 2); a ^= (a >> 4); a *= 2057; a ^= (a >> 16); break;
    case 1:  // full avalanche
      a = (a+0x7ed55d16) + (a<<12); a = (a^0xc761c23c) ^ (a>>19); a = (a+0x165667b1) + (a<<5);
      a = (a+0xd3a2646c) ^ (a<<9);  a = (a+0xfd7046c5) + (a<<3);  a = (a^0xb55a4f09) ^ (a>>16); break;
    default:
      a = (a ^ 61) ^ (a >> 16); a = a + (a << 3); a = a ^ (a >> 4); a = a * 0x27d4eb2d; a = a ^ (a >> 15); break;
  };
  return a;
}

/* originally for char (1 byte) but can be used for uint32_t (4 bytes) or uint64_t (8 bytes) */
uint32_t 
biomcmc_hashbyte_salted (const void *str, size_t size, unsigned int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint8_t *c = (uint8_t*) str;
  uint32_t g, hash = 0;
  switch(salt & 7) { // last 3 bits
    case 6: /* FNV1 hash */
      hash = 0x811C9DC5;
      for (; (*c) && (size > 0); c++, size--) { 
        hash += (hash<<1) + (hash<<4) + (hash<<7) + (hash<<8) + (hash<<24); hash ^= *c;
      }; break;

    case 5: // Jenkins's one_at_a_time hash
      for (; (*c) && (size > 0); c++, size--) { hash += *c; hash += hash << 10; hash ^= hash >> 6; }
      hash += hash << 3; hash ^= hash >> 11; hash += hash << 15; break;

    case 4: /* CRC algortihm: ASCII chars use 5-6 bits out of 8, this function worries more about the lowest 5 bits */
      /* extract high-order 5 bits (g), shift left by 5 bits, move g to low-order end and XOR into h */
      for (; (*c) && (size > 0); c++, size--) {
        g = hash & 0xf8000000; hash = ((hash << 5) ^ (g >> 27)) ^ *c; 
      }; break;

    case 3:/* PJW algorithm: like CR, worries more about lowest 5 */
      // top 4 bits are zero; shift 4 bits left, add c, get top 4 bits, if top 4 bits aren't zero, move them to low end 
      for (; (*c) && (size > 0); c++, size--) { 
        hash = (hash << 4) + *c; g = hash & 0xf0000000; 
        if (g != 0) hash = (hash ^ (g >> 24)) ^ g; 
      }; break;

    case 2: /* dbj2 algorithm */
      hash = 5381;
      for (; (*c) && (size > 0); c++, size--) hash = ((hash << 5) + hash) + *c; /* hash * 33 + c */
      break;

    case 1: /* dbj2 algorithm with XOR */
      hash = 5381;
      for (; (*c) && (size > 0); c++, size--) hash = ((hash << 5) + hash) ^ *c; /* (hash * 33) XOR c */
      break;

    default: /* sdbm algorithm */
      for (; (*c) && (size > 0); c++, size--) hash = *c + (hash << 6) + (hash << 16) - hash;
      break;
  };
  return hash;
}

uint64_t
biomcmc_hashint64_salted (uint64_t k, unsigned int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint32_t low = k, high = k >> 32UL;
  uint64_t b = k >> 32UL; // only 32bit of each key and b used 

  if (salt > 7) k = (k << 27) | (k >> 37);  // rotate key
  switch (salt & 7) { // last 3 bits
    case 7: // https://github.com/drobilla/zix/blob/main/src/digest.c
    k ^= k >> 23u; k *= 0x2127599BF4325C37ull; k ^= k >> 47u; break;
    case 6: // https://github.com/vigna/MRG32k3a/blob/master/MRG32k3a.c
      k |= 1ULL; k = (k ^ (k >> 30)) * 0xbf58476d1ce4e5b9; k = (k ^ (k >> 27)) * 0x94d049bb133111eb; k = (k >> 1) ^ (k >> 32); break;
    case 5: // Wang Yi https://github.com/Cyan4973/smhasher/blob/master/wyhash.h
      k = biomcmc_hashint64_mix_salted (biomcmc_hashint64_mix_salted (k * 0x60bee2bee120fc15ull, 0xa3b195354a39b70dull, 0), 0x1b03738712fad5c9ull, 0);
      break;
    case 4: // xxhash avalanche
      k ^= k >> 33; k *= ulx_h64[13]; k ^= k >> 29; k *= ulx_h64[14]; k ^= k >> 32; break;
    case 3:  /* ((key + (key << 3)) + (key << 8)) = key * 265 and ((key + (key << 2)) + (key << 4)) = key * 21 */
      k = (~k) + (k << 21); k = k ^ (k >> 24); k = (k + (k << 3)) + (k << 8);
      k = k ^ (k >> 14); k = (k + (k << 2)) + (k << 4); k = k ^ (k >> 28);
      k = k + (k << 31); break;
    case 2: // Lemire's blog post about concatenating two 32 bits
      k = ((ulx_h64[0] * low + ulx_h64[1] * high + ulx_h64[2]) >> 32) | ((ulx_h64[3] * low + ulx_h64[4] * high + ulx_h64[5]) & 0xFFFFFFFF00000000UL);
      break;
    case 1: // two 32bits
      k = (k+0x479ab41d)+(k<<8); k = (k^0xe4aa10ce)^(k>>5);  k = (k+0x9942f0a6)-(k<<14); k = (k^0x5aedd67d)^(k>>3); k = (k+0x17bea992)+(k<<7);
      b = (b+0x7ed55d16)+(b<<12); b = (b^0xc761c23c)^(b>>19); b = (b+0x165667b1)+(b<<5); // k=half-avalanche, b=full-avalanche
      b = (b+0xd3a2646c)^(b<<9);  b = (b+0xfd7046c5)+(b<<3);  b = (b^0xb55a4f09)^(b>>16);
      k =  (k << 32) | b; break;
    default: // mixer algo of murmur
      k ^= k >> 33; k *= ulx_h64[6]; k ^= k >> 33; k *= ulx_h64[7]; k ^= k >> 33; break;
  };
  return k;
}

uint32_t
biomcmc_hashint_mix_salted (uint32_t a, uint32_t b, unsigned int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint64_t la = a^0x7b16763u, lb = b^0xe4f5a905u;
  uint32_t c = (salt^0xdeadbeef) + (salt<<4); 

  switch (salt & 3) { // last 2 bits
    case 3: // wyhash32 from Wang Yi (2 mixes)
      la = la * lb; la ^= (la>>32); la *= 0x4a9e6939u; a = (uint32_t)(la ^ (la>>32)); break;
    case 2:
      c = c ^ (c>>10); c = c + (c<<7); c = c ^ (c>>13);
      a=a-b; a=a-c; a=a^(c >> 13); b=b-c; b=b-a; b=b^(a << 8);  c=c-a; c=c-b; c=c^(b >> 13);
      a=a-b; a=a-c; a=a^(c >> 12); b=b-c; b=b-a; b=b^(a << 16); c=c-a; c=c-b; c=c^(b >> 5);
      a=a-b; a=a-c; a=a^(c >> 3);  b=b-c; b=b-a; b=b^(a << 10); c=c-a; c=c-b; c=c^(b >> 15); a = c; break;
    case 1: // modified FNV hash with 64 bits expansion
      a = (uint32_t) (ulx_h64[0] * (uint64_t) a + ulx_h64[1] * ((uint64_t)(b) << 3) + ulx_h64[2]);
     a ^= b; a *= 16777619; break;
    default: // FNV hash
      a = 2166136261U ^ a; a *= 16777619; a ^= b; a *= 16777619;
  };
  return a;
}

uint64_t
biomcmc_hashint64_mix_salted (uint64_t a, uint64_t b, unsigned int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint64_t ha = a >> 32, hb = b >> 32, la = a & 0xFFFFFFFF00000000UL, lb = b & 0xFFFFFFFF00000000UL;
  uint64_t rh = ha * hb, rm_0 = ha * lb, rm_1 = hb * la, rl =  la * lb, t = rl + (rm_0<<32), c = t < rl;
  switch (salt & 1) { // last bit
    case 1: // mumhash high and low
      a = t + (rm_1 << 32); b = rh + (rm_0 >> 32) + (rm_1 >> 32); a += b; break;
    default: // from  Wang Yi wyhash, very similar to mumhash
      b ^= 0x60bee2bee120fc15ull;
      a = t + (rm_1 << 32); c += a<t; b = rh + (rm_0 >> 32) + (rm_1 >> 32) + c; a ^= b; break;
  };
  return a;
}

uint32_t 
biomcmc_hashint_64to32_seed (uint64_t x, unsigned int seed) /* 64 bits */
{ // published algo uses seed=0 or seed=3
  uint32_t low = x;
  uint32_t high = x >> 32UL;
  int i = seed % (ulx_h64_size - 2); 
  return (uint32_t) ((ulx_h64[i] * low + ulx_h64[i+1] * high + ulx_h64[i+2]) >> 32);
}

uint32_t 
biomcmc_hashint_64to32 (uint64_t key) /* input is 64 bits but output is 32 bits (good for hash, not for RNG...) */
{
  key = (~key) + (key << 18); key = key ^ (key >> 31); key = key * 21;
  key = key ^ (key >> 11);    key = key + (key << 6);  key = key ^ (key >> 22);
  return (uint32_t) key;
}

uint32_t 
bipartition_hash (bipartition bip) 
{ // assumes bipartition is flipped to smaller set
  int i; 
  uint32_t a = bip->n_ones; // inspired by FNV hash
  a = ~a + (a << 15); a ^= (a >> 12); a+= (a << 2); a ^= (a >> 4); a *= 2057; a ^= (a >> 16);
  for (i=0; i < bip->n->ints; i++) a = (a ^ biomcmc_hashint_64to32 (bip->bs[i])) * 16777619;
  return a;
}

/*** functions using constant_random_lists.h values ***/

uint32_t
biomcmc_get_salt_set_from_spice_table (uint32_t seeds[], uint32_t *salt, uint32_t salt_length)
{
  uint16_t id;
  uint32_t index, si = 0;
  union { uint32_t i; double f; } un_fi;  // same memory space shared 

  if (!salt_length) return 0;

  index = seeds[0];
  // below, if y is power of 2, then x & (y-1) == x % y --> (index[] & size-1) or (index%size)
  id = index & (rnd_salt_h64_list_size - 1); index /= rnd_salt_h64_list_size;  // size = 256 (total_0 8 bits) 
  salt[si++] = (uint32_t) rnd_salt_h64_list[id]; 

  if (salt_length > 1) {
    id = index & (rnd_salt_h64_list_size - 1); index /= rnd_salt_h64_list_size; // size = 256 (total_0 16 bits)
    salt[si++] = (rnd_salt_h64_list[id] >> 32); 
  }
  if (salt_length > 2) {
    id = index & (rnd_salt_h16_list_size - 1); index /= rnd_salt_h16_list_size; // size = 256 (total_0 24 bits) 
    salt[si++] += rnd_salt_h16_list[id]; 
  }
  if (salt_length > 3) {
    index = seeds[1];
    id = index & (prime_salt_list_size - 1); index /= prime_salt_list_size;  // size = 512 (total_1 9 bits)
    salt[si++] = prime_salt_list[id]; 
  }
  if (salt_length > 4) {
    id = index % marsaglia_constants_size;  index /= marsaglia_constants_size; // size = 81 (total_1 15 bits)
    salt[si++] = (marsaglia_constants[id] << 16) - 1;
  }
  if (salt_length > 5) {
    id = index % marsaglia_constants_size;  index /= marsaglia_constants_size; // size = 81 (total_1 21 bits) 
    salt[si++] = (marsaglia_constants[id] << 15) - 1;
  }
  if (salt_length > 6) {
    id = index % ulx_h64_size; index /= ulx_h64_size;      // size = 185 (total_1 29 bits) 
    salt[si++] = (uint32_t) ulx_h64[id];
  }
  if (salt_length > 7) {
    index = seeds[2];
    id = index % ulx_h64_size; index /= ulx_h64_size; // size = 185 (total_2 7 bits)
    salt[si++] = ulx_h64[id] >> 16; // cannot be 32 since some value (ulx_h64[10]) has only lower 32 bits
  }
  if (salt_length > 8) {
    id = index % lgamma_algmcs_size; index /= lgamma_algmcs_size; // size = 15 (total_2 11 bits)
    un_fi.f = lgamma_algmcs[id];
    salt[si++] = un_fi.i;
  }
  if (salt_length > 9) {
    id = index % lgamma_coeffs_size; index /= lgamma_coeffs_size; // size = 40 (total_2 16 bits)
    un_fi.f = lgamma_coeffs[id];
    salt[si++] = un_fi.i;
  }
  if (salt_length > 10) {
    id = index % stirl_sferr_halves_size; index /= stirl_sferr_halves_size; // size = 31 (total_2 21 bits)
    un_fi.f = stirl_sferr_halves[id] + 2.7; // first value is zero
    salt[si++] = un_fi.i;
  }
  /* combinations using consecutive element */
  if (salt_length > 11) { 
    id = index % (ulx_h64_size - 1);  index /= (ulx_h64_size-1); // size = 184 (total_2 29 bits)
    salt[si++] = (uint32_t) ((ulx_h64[id] >> 9) + (ulx_h64[id+1] >> 7));
  }
  if (salt_length > 12) {
    index = seeds[3];
    id = index % (rnd_salt_h64_list_size - 1); index /= (rnd_salt_h64_list_size-1); // size = 255 (total_3 8 bits)
    salt[si++] = (uint32_t)((rnd_salt_h64_list[id] >> 28) + (rnd_salt_h64_list[id+1] >> 25));
  }
  if (salt_length > 13) {
    id = index % (rnd_salt_h16_list_size - 1); index /= (rnd_salt_h16_list_size-1); // size = 255 (total_3 16 bits)
    salt[si++] = rnd_salt_h16_list[id] + ((uint32_t)(rnd_salt_h16_list[id+1])  << 4);
  }
  if (salt_length > 14) {
    id = index % (prime_salt_list_size - 1); index /= (prime_salt_list_size-1); // size = 511 (total_3 25 bits)
    salt[si++] = (prime_salt_list[id] << 2) + (prime_salt_list[id+1]  >> 2);
  }
  if (salt_length > 15) {
    id = index % (prime_salt_list_size - 1); index /= (prime_salt_list_size-1); // size = 511 (total_3 25 bits)
    salt[si++] = (prime_salt_list[id] << 2) + (prime_salt_list[id+1]  >> 2);
  }
  return si; 
}

void
biomcmc_invert_bits32_by_address (uint32_t *n)
{
  *n = ((*n >>  1) & 0x55555555) | ((*n <<  1) & 0xaaaaaaaa);
  *n = ((*n >>  2) & 0x33333333) | ((*n <<  2) & 0xcccccccc);
  *n = ((*n >>  4) & 0x0f0f0f0f) | ((*n <<  4) & 0xf0f0f0f0);
  *n = ((*n >>  8) & 0x00ff00ff) | ((*n <<  8) & 0xff00ff00);
  *n = ((*n >> 16) & 0x0000ffff) | ((*n << 16) & 0xffff0000);
}

void
biomcmc_salt_vector32_from_spice_table (uint32_t *a, uint32_t n_a, uint32_t seed[])
{
  uint32_t i,j;
  for (i=0; i < n_a;) {
    i += biomcmc_get_salt_set_from_spice_table (seed, a + i, n_a - i);
    for (j=0;j<4;j++) seed[j] = biomcmc_hashint_salted (seed[j], j); // given same seed, salts are the same
  }
}

void
biomcmc_salt_vector64_from_spice_table (uint64_t *a, uint32_t n_a, uint32_t seed[])
{
  uint32_t i, seed2[4], *x = (uint32_t*)biomcmc_malloc (sizeof(uint32_t) * n_a);
  uint64_t tmp;
  for (i=0; i < n_a; i++) x[i] = (uint32_t) a[i]; 
  // right-most bits 
  biomcmc_salt_vector32_from_spice_table (x, n_a, seed);
  for (i=0; i < n_a; i++) { tmp = a[i]; a[i] = x[i]; x[i] = tmp >> 32; }
  // left-most bits are in backwards order from table, and with inverted bits
  for (i=0;i<4;i++) seed2[i] = RoL(seed[i], i+3); // given same seed, salts are the same
  biomcmc_salt_vector32_from_spice_table (x, n_a, seed2);
  for (i=0; i < n_a; i++) { biomcmc_invert_bits32_by_address (&x[n_a-i-1]); a[i] |= ((uint64_t)(x[n_a-i-1]) << 32); }
  if (x) free (x);
}
#undef RoL
#undef RoR


/*** MurmurHash3 from https://github.com/PeterScott/murmur3/  written originally by Austin Appleby and CC0 ***/

uint64_t
biomcmc_murmurhash3_64bits (const void *key, const size_t len, const uint32_t seed)
{
  return biomcmc_murmurhash3_128bits (key, len, seed, NULL);
}

uint64_t
biomcmc_murmurhash3_128bits (const void *key, const size_t len, const uint32_t seed, void *out)
{  /* out[] is 128 bits (4 x uint_32, for instance), assumes 64bits machine  */
  const uint8_t * data = (const uint8_t*) key;
  const int nblocks = len / 16;
  int i;
  uint64_t h1 = seed, h2 = seed;
  const uint64_t * blocks = (const uint64_t *) (data);

  for(i = 0; i < nblocks; i++) {
    uint64_t k1 = blocks[i*2+0];
    uint64_t k2 = blocks[i*2+1];
    k1 *= ulx_h64[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx_h64[9]; h1 ^= k1;
    h1 = (h1 << 27) | (h1 >> 37); h1 += h2; h1 = h1 * 5 + ulx_h64[10];
    k2 *= ulx_h64[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx_h64[8]; h2 ^= k2;
    h2 = (h2 << 31) | (h2 >> 33); h2 += h1; h2 = h2 * 5 + ulx_h64[11];
  }

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  uint64_t k1 = 0, k2 = 0;

  switch(len & 15) {
    case 15: k2 ^= (uint64_t)(tail[14]) << 48; attribute_FALLTHROUGH // macro defined in lowlevel.h since clang doesn't understand this __atribute
    case 14: k2 ^= (uint64_t)(tail[13]) << 40; attribute_FALLTHROUGH
    case 13: k2 ^= (uint64_t)(tail[12]) << 32; attribute_FALLTHROUGH
    case 12: k2 ^= (uint64_t)(tail[11]) << 24; attribute_FALLTHROUGH
    case 11: k2 ^= (uint64_t)(tail[10]) << 16; attribute_FALLTHROUGH
    case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;  attribute_FALLTHROUGH
    case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;  
             k2 *= ulx_h64[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx_h64[8]; h2 ^= k2; attribute_FALLTHROUGH
    case  8: k1 ^= (uint64_t)(tail[ 7]) << 56; attribute_FALLTHROUGH
    case  7: k1 ^= (uint64_t)(tail[ 6]) << 48; attribute_FALLTHROUGH
    case  6: k1 ^= (uint64_t)(tail[ 5]) << 40; attribute_FALLTHROUGH
    case  5: k1 ^= (uint64_t)(tail[ 4]) << 32; attribute_FALLTHROUGH
    case  4: k1 ^= (uint64_t)(tail[ 3]) << 24; attribute_FALLTHROUGH
    case  3: k1 ^= (uint64_t)(tail[ 2]) << 16; attribute_FALLTHROUGH
    case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;  attribute_FALLTHROUGH 
    case  1: k1 ^= (uint64_t)(tail[ 0]) << 0; // equiv to "case 1"  
             k1 *= ulx_h64[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx_h64[9]; h1 ^= k1; break;
  };

  h1 ^= len; h2 ^= len;
  h1 += h2; h2 += h1;

  h1 ^= h1 >> 33; h1 *= ulx_h64[6]; h1 ^= h1 >> 33; h1 *= ulx_h64[7]; h1 ^= h1 >> 33;
  h2 ^= h2 >> 33; h2 *= ulx_h64[6]; h2 ^= h2 >> 33; h2 *= ulx_h64[7]; h2 ^= h2 >> 33;

  h1 += h2; h2 += h1;
  if (out) {
    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
  }
  return biomcmc_hashint64_mix_salted (h1, h2, /*salt*/ 1);
}

// https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
uint32_t 
biomcmc_murmurhash3_32bits (const void *data, const size_t nbytes, const uint32_t seed)
{ /* assumes 32 bits machine, and returns 32 bits */
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const int nblocks = nbytes / 4;
  const uint32_t *blocks = (const uint32_t *) (data);
  const uint8_t *tail = (const uint8_t *) (data + (nblocks * 4));
  uint32_t k, h = seed;
  int i;

  for (i = 0; i < nblocks; i++) {
    k = blocks[i]; k *= c1; k = (k << 15) | (k >> (32 - 15)); k *= c2;
    h ^= k; h = (h << 13) | (h >> (32 - 13)); h = (h * 5) + 0xe6546b64;
  }

  k = 0;
  switch (nbytes & 3) {
    case 3: k ^= tail[2] << 16; attribute_FALLTHROUGH
    case 2: k ^= tail[1] << 8;  attribute_FALLTHROUGH
    case 1: k ^= tail[0]; k *= c1; k = (k << 15) | (k >> (32 - 15)); k *= c2; h ^= k; break;
  };
  h ^= nbytes; h ^= h >> 16; h *= 0x85ebca6b; h ^= h >> 13; h *= 0xc2b2ae35; h ^= h >> 16;
  return h;
}

// https://github.com/torvalds/linux/blob/master/lib/xxhash.c and https://github.com/Cyan4973/smhasher/blob/master/xxhash.c
static uint64_t XXH_read64 (const void* memPtr);
static uint32_t XXH_read32 (const void* memPtr);
static uint64_t XXH_swap64 (uint64_t x);
static uint32_t XXH_swap32 (uint32_t x);
static int XXH_isLittleEndian(void);
static uint64_t xxh_get_unaligned_le64 (const void* ptr);
static uint64_t xxh_get_unaligned_le32 (const void* ptr);
static uint64_t xxh64_round (uint64_t acc, const uint64_t input);
static uint64_t xxh64_merge_round (uint64_t acc, uint64_t val);

/*!XXH_FORCE_MEMORY_ACCESS :
 * By default, access to unaligned memory is controlled by `memcpy()`, which is safe and portable.
 * Unfortunately, on some target/compiler combinations, the generated assembly is sub-optimal.
 * The below switch allow to select different access method for improved performance.
 * Method 0 (default) : use `memcpy()`. Safe and portable.
 * Method 1 : `__packed` statement. It depends on compiler extension (ie, not portable).
 *            This method is safe if your compiler supports it, and *generally* as fast or faster than `memcpy`.
 * Method 2 : direct access. This method doesn't depend on compiler but violate C standard.
 *            It can generate buggy code on targets which do not support unaligned memory accesses.
 *            But in some circumstances, it's the only known way to get the most performance (ie GCC + ARMv6)
 * See http://stackoverflow.com/a/32095106/646947 for details.
 * Prefer these methods in priority order (0 > 1 > 2)
 */
#ifndef XXH_FORCE_MEMORY_ACCESS   /* can be defined externally, on command line for example */
#  if defined(__GNUC__) && ( defined(__ARM_ARCH_6__) || defined(__ARM_ARCH_6J__) \
                        || defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) \
                        || defined(__ARM_ARCH_6ZK__) || defined(__ARM_ARCH_6T2__) )
#    define XXH_FORCE_MEMORY_ACCESS 2
#  elif (defined(__INTEL_COMPILER) && !defined(_WIN32)) || \
  (defined(__GNUC__) && ( defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) \
                    || defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) \
                    || defined(__ARM_ARCH_7S__) ))
#    define XXH_FORCE_MEMORY_ACCESS 1
#  endif
#endif

#if (defined(XXH_FORCE_MEMORY_ACCESS) && (XXH_FORCE_MEMORY_ACCESS==2))
/* Force direct memory access. Only works on CPU which support unaligned memory access in hardware */
static uint64_t XXH_read64 (const void* memPtr) { return *(const uint64_t*) memPtr; }
static uint32_t XXH_read32 (const void* memPtr) { return *(const uint32_t*) memPtr; }
#elif (defined(XXH_FORCE_MEMORY_ACCESS) && (XXH_FORCE_MEMORY_ACCESS==1))
/* __pack instructions are safer but compiler specific, potentially problematic for some compilers (currently only defined for gcc and icc) */
typedef union { uint32_t u32; uint64_t u64; } __attribute__((packed)) unalign64;
typedef union { uint32_t u32; } __attribute__((packed)) unalign;
static uint64_t XXH_read64 (const void* ptr) { return ((const unalign64*)ptr)->u64; }
static uint32_t XXH_read32 (const void* ptr) { return ((const unalign*)ptr)->u32; }
#else
/* portable and safe solution. Generally efficient. see : http://stackoverflow.com/a/32095106/646947   */
static uint64_t XXH_read64(const void* memPtr) { uint64_t val;  memcpy(&val, memPtr, sizeof(val)); return val; }
static uint32_t XXH_read32(const void* memPtr) { uint32_t val;  memcpy(&val, memPtr, sizeof(val)); return val; }
#endif /* XXH_FORCE_DIRECT_MEMORY_ACCESS */
 
static uint64_t XXH_swap64 (uint64_t x)
{
  return  
    ((x << 56) & 0xff00000000000000ULL) | ((x << 40) & 0x00ff000000000000ULL) | ((x << 24) & 0x0000ff0000000000ULL) | ((x << 8)  & 0x000000ff00000000ULL) |
    ((x >> 8)  & 0x00000000ff000000ULL) | ((x >> 24) & 0x0000000000ff0000ULL) | ((x >> 40) & 0x000000000000ff00ULL) | ((x >> 56) & 0x00000000000000ffULL);
}

static uint32_t XXH_swap32 (uint32_t x)
{
  return  
    ((x << 24) & 0xff000000 ) | ((x <<  8) & 0x00ff0000 ) | ((x >>  8) & 0x0000ff00 ) | ((x >> 24) & 0x000000ff );
}

#ifndef XXH_CPU_LITTLE_ENDIAN
static int XXH_isLittleEndian(void)
{
    const union { uint32_t u; uint8_t c[4]; } one = { 1 };   /* don't use static : performance detrimental  */
    return one.c[0];
}
#   define XXH_CPU_LITTLE_ENDIAN   XXH_isLittleEndian()
#endif

static uint64_t xxh_get_unaligned_le64 (const void* ptr)
{
  return XXH_CPU_LITTLE_ENDIAN ? XXH_read64(ptr) : XXH_swap64(XXH_read64(ptr));
}

static uint64_t xxh_get_unaligned_le32 (const void* ptr)
{
  return XXH_CPU_LITTLE_ENDIAN ? XXH_read32(ptr) : XXH_swap32(XXH_read32(ptr));
}

static uint64_t xxh64_round (uint64_t acc, const uint64_t input)
{
	acc += input * ulx_h64[13];
  acc = (acc << 31) | (acc >> 33);
	acc *= ulx_h64[12];
	return acc;
}

static uint64_t xxh64_merge_round (uint64_t acc, uint64_t val)
{
	val = xxh64_round (0, val);
	acc ^= val;
	acc = acc * ulx_h64[12] + ulx_h64[15];
	return acc;
}

uint64_t 
biomcmc_xxh64 (const void *input, const size_t len, const uint32_t seed)
{
	const uint8_t *p = (const uint8_t *) input;
	const uint8_t *const b_end = p + len;
	uint64_t h64 = (uint64_t) seed;

  if (len > 31) {
    const uint8_t *const limit = b_end - 32;
    uint64_t v1 = h64 + ulx_h64[12]+ ulx_h64[13];
    uint64_t v2 = h64 + ulx_h64[13];
    uint64_t v3 = h64 + 0;
    uint64_t v4 = h64 - ulx_h64[12];

    do {
      v1 = xxh64_round(v1, xxh_get_unaligned_le64(p)); p += 8;
      v2 = xxh64_round(v2, xxh_get_unaligned_le64(p)); p += 8;
      v3 = xxh64_round(v3, xxh_get_unaligned_le64(p)); p += 8;
      v4 = xxh64_round(v4, xxh_get_unaligned_le64(p)); p += 8;
    } while (p <= limit);

    h64 = ((v1<<1)|(v1>>63)) + ((v2<<7)|(v2>>57)) + ((v3<<12)|(v3>>52)) + ((v4<<18)|(v4>>46)); // rotations
    h64 = xxh64_merge_round (h64, v1);
    h64 = xxh64_merge_round (h64, v2);
    h64 = xxh64_merge_round (h64, v3);
    h64 = xxh64_merge_round (h64, v4);
  } else {
    h64  += ulx_h64[16];
  }

  h64 += (uint64_t) len;

  while (p + 8 <= b_end) {
    const uint64_t k1 = xxh64_round (0, xxh_get_unaligned_le64(p));
    h64 ^= k1;
    h64 = ((h64 << 27) | (h64 >> 37)) * ulx_h64[12] + ulx_h64[15];
    p += 8;
  }
  if (p + 4 <= b_end) {
    h64 ^= (uint64_t) (xxh_get_unaligned_le32(p)) * ulx_h64[12];
    h64 = ((h64 << 23)|(h64 >> 41)) * ulx_h64[13] + ulx_h64[14];
    p += 4;
  }
  while (p < b_end) {
    h64 ^= (*p) * ulx_h64[16];
    h64 = ((h64<<11)|(h64>>53)) * ulx_h64[12];
    p++;
  }
  h64 ^= h64 >> 33; h64 *= ulx_h64[13]; h64 ^= h64 >> 29; h64 *= ulx_h64[14]; h64 ^= h64 >> 32; // avalanche
  return h64;
}

/*** google Highway Hash ***/

static void ghh_ZipperMergeAndAdd (const uint64_t v1, const uint64_t v0, uint64_t* add1, uint64_t* add0);
static void ghh_Update (const uint64_t lanes[4], HighwayHashState* state);
static uint64_t ghh_Read64 (const uint8_t* src);
static void ghh_Rotate32By (uint64_t count, uint64_t lanes[4]);
static void ghh_Permute (const uint64_t v[4], uint64_t* permuted);
void ghh_PermuteAndUpdate (HighwayHashState* state);
static void ghh_ModularReduction (uint64_t a3_unmasked, uint64_t a2, uint64_t a1, uint64_t a0, uint64_t* m1, uint64_t* m0);
static void ghh_ProcessAll (const uint8_t* data, size_t size, const uint64_t key[4], HighwayHashState* state); 
 
/* Internal implementation          */
void HighwayHashReset (const uint64_t key[4], HighwayHashState* state) {
  state->mul0[0] = 0xdbe6d5d5fe4cce2full;
  state->mul0[1] = 0xa4093822299f31d0ull;
  state->mul0[2] = 0x13198a2e03707344ull;
  state->mul0[3] = 0x243f6a8885a308d3ull;
  state->mul1[0] = 0x3bd39e10cb0ef593ull;
  state->mul1[1] = 0xc0acf169b5f18a8cull;
  state->mul1[2] = 0xbe5466cf34e90c6cull;
  state->mul1[3] = 0x452821e638d01377ull;
  state->v0[0] = state->mul0[0] ^ key[0];
  state->v0[1] = state->mul0[1] ^ key[1];
  state->v0[2] = state->mul0[2] ^ key[2];
  state->v0[3] = state->mul0[3] ^ key[3];
  state->v1[0] = state->mul1[0] ^ ((key[0] >> 32) | (key[0] << 32));
  state->v1[1] = state->mul1[1] ^ ((key[1] >> 32) | (key[1] << 32));
  state->v1[2] = state->mul1[2] ^ ((key[2] >> 32) | (key[2] << 32));
  state->v1[3] = state->mul1[3] ^ ((key[3] >> 32) | (key[3] << 32));
}

static void ghh_ZipperMergeAndAdd (const uint64_t v1, const uint64_t v0, uint64_t* add1, uint64_t* add0) {
  *add0 += (((v0 & 0xff000000ull)   | (v1 & 0xff00000000ull)) >> 24)    | (((v0 & 0xff0000000000ull)          | (v1 & 0xff000000000000ull)) >> 16) |
           (  v0 & 0xff0000ull)     | ((v0 & 0xff00ull) << 32)          | ((v1 & 0xff00000000000000ull) >> 8) | (v0 << 56);
  *add1 += (((v1 & 0xff000000ull)   | (v0 & 0xff00000000ull)) >> 24)    | (v1 & 0xff0000ull)                  | ((v1 & 0xff0000000000ull) >> 16) |
           ((v1 & 0xff00ull) << 24) | ((v0 & 0xff000000000000ull) >> 8) | ((v1 & 0xffull) << 48)              | (v0 & 0xff00000000000000ull);
}

static void ghh_Update (const uint64_t lanes[4], HighwayHashState* state) {
  int i;
  for (i = 0; i < 4; ++i) {
    state->v1[i] += state->mul0[i] + lanes[i];
    state->mul0[i] ^= (state->v1[i] & 0xffffffff) * (state->v0[i] >> 32);
    state->v0[i] += state->mul1[i];
    state->mul1[i] ^= (state->v0[i] & 0xffffffff) * (state->v1[i] >> 32);
  }
  ghh_ZipperMergeAndAdd(state->v1[1], state->v1[0], &state->v0[1], &state->v0[0]);
  ghh_ZipperMergeAndAdd(state->v1[3], state->v1[2], &state->v0[3], &state->v0[2]);
  ghh_ZipperMergeAndAdd(state->v0[1], state->v0[0], &state->v1[1], &state->v1[0]);
  ghh_ZipperMergeAndAdd(state->v0[3], state->v0[2], &state->v1[3], &state->v1[2]);
}

static uint64_t ghh_Read64 (const uint8_t* src) {
  return (uint64_t)src[0] | ((uint64_t)src[1] << 8) |
      ((uint64_t)src[2] << 16) | ((uint64_t)src[3] << 24) |
      ((uint64_t)src[4] << 32) | ((uint64_t)src[5] << 40) |
      ((uint64_t)src[6] << 48) | ((uint64_t)src[7] << 56);
}

void HighwayHashUpdatePacket (const uint8_t* packet, HighwayHashState* state) {
  uint64_t lanes[4];
  lanes[0] = ghh_Read64(packet + 0);
  lanes[1] = ghh_Read64(packet + 8);
  lanes[2] = ghh_Read64(packet + 16);
  lanes[3] = ghh_Read64(packet + 24);
  ghh_Update (lanes, state);
}

static void ghh_Rotate32By (uint64_t count, uint64_t lanes[4]) {
  int i;
  for (i = 0; i < 4; ++i) {
    uint32_t half0 = lanes[i] & 0xffffffff;
    uint32_t half1 = (lanes[i] >> 32);
    lanes[i] = (half0 << count) | (half0 >> (32 - count));
    lanes[i] |= (uint64_t)((half1 << count) | (half1 >> (32 - count))) << 32;
  }
}

void HighwayHashUpdateRemainder (const uint8_t* bytes, const size_t size_mod32, HighwayHashState* state) {
  int i;
  const size_t size_mod4 = size_mod32 & 3;
  const uint8_t* remainder = bytes + (size_mod32 & ~3);
  uint8_t packet[32] = {0};
  for (i = 0; i < 4; ++i) state->v0[i] += ((uint64_t)size_mod32 << 32) + size_mod32;
  ghh_Rotate32By (size_mod32, state->v1);
  for (i = 0; i < remainder - bytes; i++) packet[i] = bytes[i];
  if (size_mod32 & 16) {
    for (i = 0; i < 4; i++) packet[28 + i] = remainder[i + size_mod4 - 4];
  } else {
    if (size_mod4) {
      packet[16 + 0] = remainder[0];
      packet[16 + 1] = remainder[size_mod4 >> 1];
      packet[16 + 2] = remainder[size_mod4 - 1];
    }
  }
  HighwayHashUpdatePacket (packet, state);
}

static void ghh_Permute (const uint64_t v[4], uint64_t* permuted) {
  permuted[0] = (v[2] >> 32) | (v[2] << 32);
  permuted[1] = (v[3] >> 32) | (v[3] << 32);
  permuted[2] = (v[0] >> 32) | (v[0] << 32);
  permuted[3] = (v[1] >> 32) | (v[1] << 32);
}

void ghh_PermuteAndUpdate (HighwayHashState* state) {
  uint64_t permuted[4];
  ghh_Permute (state->v0, permuted);
  ghh_Update (permuted, state);
}

static void ghh_ModularReduction (uint64_t a3_unmasked, uint64_t a2, uint64_t a1, uint64_t a0, uint64_t* m1, uint64_t* m0) {
  uint64_t a3 = a3_unmasked & 0x3FFFFFFFFFFFFFFFull;
  *m1 = a1 ^ ((a3 << 1) | (a2 >> 63)) ^ ((a3 << 2) | (a2 >> 62));
  *m0 = a0 ^ (a2 << 1) ^ (a2 << 2);
}

uint64_t HighwayHashFinalize64 (HighwayHashState* state) {
  int i;
  for (i = 0; i < 4; i++) ghh_PermuteAndUpdate(state);
  return state->v0[0] + state->v1[0] + state->mul0[0] + state->mul1[0];
}

void HighwayHashFinalize128 (HighwayHashState* state, uint64_t hash[2]) {
  int i;
  for (i = 0; i < 6; i++) ghh_PermuteAndUpdate (state);
  hash[0] = state->v0[0] + state->mul0[0] + state->v1[2] + state->mul1[2];
  hash[1] = state->v0[1] + state->mul0[1] + state->v1[3] + state->mul1[3];
}

void HighwayHashFinalize256 (HighwayHashState* state, uint64_t hash[4]) {
  int i;
  /* We anticipate that 256-bit hashing will be mostly used with long messages because storing and using the 256-bit hash (in contrast to 128-bit)
     carries a larger additional constant cost by itself. Doing extra rounds here hardly increases the per-byte cost of long messages. */
  for (i = 0; i < 10; i++)  ghh_PermuteAndUpdate(state);
  ghh_ModularReduction (state->v1[1] + state->mul1[1], state->v1[0] + state->mul1[0], 
                        state->v0[1] + state->mul0[1], state->v0[0] + state->mul0[0], &hash[1], &hash[0]);
  ghh_ModularReduction (state->v1[3] + state->mul1[3], state->v1[2] + state->mul1[2], 
                        state->v0[3] + state->mul0[3], state->v0[2] + state->mul0[2], &hash[3], &hash[2]);
}

/* Non-cat API: single call on full data                                      */
static void ghh_ProcessAll (const uint8_t* data, size_t size, const uint64_t key[4], HighwayHashState* state) {
  size_t i;
  HighwayHashReset (key, state);
  for (i = 0; i + 32 <= size; i += 32) HighwayHashUpdatePacket (data + i, state);
  if ((size & 31) != 0) HighwayHashUpdateRemainder (data + i, size & 31, state);
}

uint64_t HighwayHash64 (const uint8_t* data, size_t size, const uint64_t key[4]) {
  HighwayHashState state;
  ghh_ProcessAll (data, size, key, &state);
  return HighwayHashFinalize64 (&state);
}

void HighwayHash128 (const uint8_t* data, size_t size, const uint64_t key[4], uint64_t hash[2]) {
  HighwayHashState state;
  ghh_ProcessAll (data, size, key, &state);
  HighwayHashFinalize128 (&state, hash);
}

void HighwayHash256 (const uint8_t* data, size_t size, const uint64_t key[4], uint64_t hash[4]) {
  HighwayHashState state;
  ghh_ProcessAll (data, size, key, &state);
  HighwayHashFinalize256 (&state, hash);
}

/* Cat API: allows appending with multiple calls                              */
void HighwayHashCatStart (const uint64_t key[4], HighwayHashCat* state) {
  HighwayHashReset (key, &state->state);
  state->num = 0;
}

void HighwayHashCatAppend (const uint8_t* bytes, size_t num, HighwayHashCat* state) {
  size_t i;
  if (state->num != 0) {
    size_t num_add = num > (32u - state->num) ? (32u - state->num) : num;
    for (i = 0; i < num_add; i++) state->packet[state->num + i] = bytes[i];
    state->num += num_add;
    num -= num_add;
    bytes += num_add;
    if (state->num == 32) { HighwayHashUpdatePacket (state->packet, &state->state); state->num = 0; }
  }
  while (num >= 32) { HighwayHashUpdatePacket(bytes, &state->state); num -= 32; bytes += 32; }
  for (i = 0; i < num; i++) { state->packet[state->num] = bytes[i]; state->num++; }
}

uint64_t HighwayHashCatFinish64 (const HighwayHashCat* state) {
  HighwayHashState copy = state->state;
  if (state->num) HighwayHashUpdateRemainder (state->packet, state->num, &copy);
  return HighwayHashFinalize64 (&copy);
}

void HighwayHashCatFinish128 (const HighwayHashCat* state, uint64_t hash[2]) {
  HighwayHashState copy = state->state;
  if (state->num) HighwayHashUpdateRemainder (state->packet, state->num, &copy);
  HighwayHashFinalize128 (&copy, hash);
}

void HighwayHashCatFinish256 (const HighwayHashCat* state, uint64_t hash[4]) {
  HighwayHashState copy = state->state;
  if (state->num) HighwayHashUpdateRemainder (state->packet, state->num, &copy);
  HighwayHashFinalize256 (&copy, hash);
}

/*  google highwayhash Usage examples:
void Example64() {
  uint64_t key[4] = {1, 2, 3, 4};
  const char* text = "Hello world!";
  size_t size = strlen(text);
  uint64_t hash = HighwayHash64((const uint8_t*)text, size, key);
  printf("%016"PRIx64"\n", hash);
}

void Example64Cat() {
  uint64_t key[4] = {1, 2, 3, 4};
  HighwayHashCat state;
  uint64_t hash;
  HighwayHashCatStart(key, &state);
  HighwayHashCatAppend((const uint8_t*)"Hello", 5, &state);
  HighwayHashCatAppend((const uint8_t*)" world!", 7, &state);
  hash = HighwayHashCatFinish64(state);
  printf("%016"PRIx64"\n", hash);
}
*/
