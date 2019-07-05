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

static int ulx_size = 12;
static uint64_t ulx[] = {
  0x65d200ce55b19ad8UL, 0x4f2162926e40c299UL, 0x162dd799029970f8UL, // 0...2
  0x68b665e6872bd1f4UL, 0xb6cfcf9d79b51db2UL, 0x7a2b92ae912898c2UL, // 3...5
  0xff51afd7ed558ccdUL, 0xc4ceb9fe1a85ec53UL, 0x87c37b91114253d5UL, 0x4cf5ad432745937fUL, // 6...9 (murmurhash) 
  0x52dce729UL, 0x38495ab5UL, // 10...11 
  11400714785074694791ULL, 14029467366897019727ULL, 1609587929392839161ULL, // 12..14 used by xxhash
  9650029242287828579ULL, 2870177450012600261ULL}; // 15..16 used by xxhash

uint32_t 
biomcmc_hashint_salted (uint32_t a, int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  switch(salt & 15) { // 4 last bits
    case 8: 
      a -= (a<<6); a ^= (a>>17); a -= (a<<9); a ^= (a<<4); a -= (a<<3); a ^= (a<<10); a ^= (a>>15); break;
    case 7:
      a += ~(a<<15); a ^= (a>>10); a += (a<<3); a ^= (a>>6); a += ~(a<<11); a ^= (a>>16); break;
    case 6:// half-avalanche: Every input bit affects itself and all higher output bits, plus a few lower output bits
      a = (a+0x479ab41d) + (a<<8); a = (a^0xe4aa10ce) ^ (a>>5); a = (a+0x9942f0a6) - (a<<14); 
      a = (a^0x5aedd67d) ^ (a>>3); a = (a+0x17bea992) + (a<<7); break;
    case 5:
      a = (a^0xdeadbeef) + (a<<4); a = a ^ (a>>10); a = a + (a<<7); a = a ^ (a>>13); break;
    case 4: // must use use at least the 17 lowest bits
      a = a ^ (a>>4); a = (a^0xdeadbeef) + (a<<5); a = a ^ (a>>11); break;
    case 3: /* hashCodes that differ only by constant multiples at each bit have  bounded number of collisions (~8 at default load factor) */
      a ^= (a >> 20) ^ (a >> 12); a = a ^ (a >> 7) ^ (a >> 4); break;
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
biomcmc_hashbyte_salted (const void *str, size_t size, int salt)
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
biomcmc_hashint64_salted (uint64_t k, int salt)
{ // salt != seed, since main usage is to determine which hash function is used
  uint32_t low = k, high = k >> 32UL;
  uint64_t b = k >> 32UL; // only 32bit of each key and b used 

  if (salt > 7) k = (k << 27) | (k >> 37);  // rotate key
  switch (salt & 7) { // last 3 bits
    case 5: // Wang Yi https://github.com/Cyan4973/smhasher/blob/master/wyhash.h
      k = biomcmc_hashint64_mix_salted (biomcmc_hashint64_mix_salted (k * 0x60bee2bee120fc15ull, 0xa3b195354a39b70dull, 0), 0x1b03738712fad5c9ull, 0);
      break;
    case 4: // xxhash avalanche
      k ^= k >> 33; k *= ulx[13]; k ^= k >> 29; k *= ulx[14]; k ^= k >> 32; break;
    case 3:  /* ((key + (key << 3)) + (key << 8)) = key * 265 and ((key + (key << 2)) + (key << 4)) = key * 21 */
      k = (~k) + (k << 21); k = k ^ (k >> 24); k = (k + (k << 3)) + (k << 8);
      k = k ^ (k >> 14); k = (k + (k << 2)) + (k << 4); k = k ^ (k >> 28);
      k = k + (k << 31); break;
    case 2: // Lemire's blog post about concatenating two 32 bits
      k = ((ulx[0] * low + ulx[1] * high + ulx[2]) >> 32) | ((ulx[3] * low + ulx[4] * high + ulx[5]) & 0xFFFFFFFF00000000UL);
      break;
    case 1: // two 32bits
      k = (k+0x479ab41d)+(k<<8); k = (k^0xe4aa10ce)^(k>>5);  k = (k+0x9942f0a6)-(k<<14); k = (k^0x5aedd67d)^(k>>3); k = (k+0x17bea992)+(k<<7);
      b = (b+0x7ed55d16)+(b<<12); b = (b^0xc761c23c)^(b>>19); b = (b+0x165667b1)+(b<<5); // k=half-avalanche, b=full-avalanche
      b = (b+0xd3a2646c)^(b<<9);  b = (b+0xfd7046c5)+(b<<3);  b = (b^0xb55a4f09)^(b>>16);
      k =  (k << 32) | b; break;
    default: // mixer algo of murmur
      k ^= k >> 33; k *= ulx[6]; k ^= k >> 33; k *= ulx[7]; k ^= k >> 33; break;
  };
  return k;
}

uint32_t
biomcmc_hashint_mix_salted (uint32_t a, uint32_t b, int salt)
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
      a = (uint32_t) (ulx[0] * (uint64_t) a + ulx[1] * ((uint64_t)(b) << 3) + ulx[2]);
      a ^= b; a *= 16777619; break;
    default: // FNV hash
      a = 2166136261U ^ a; a *= 16777619; a ^= b; a *= 16777619;
  };
  return a;
}

uint64_t
biomcmc_hashint64_mix_salted (uint64_t a, uint64_t b, int salt)
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
bipartition_hash (bipartition bip) 
{ // assumes bipartition is flipped to smaller set
  int i; 
  uint32_t a = bip->n_ones; // inspired by FNV hash
  a = ~a + (a << 15); a ^= (a >> 12); a+= (a << 2); a ^= (a >> 4); a *= 2057; a ^= (a >> 16);
  for (i=0; i < bip->n->ints; i++) a = (a ^ biomcmc_hashint_64to32 (bip->bs[i])) * 16777619;
  return a;
}

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
    k1 *= ulx[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx[9]; h1 ^= k1;
    h1 = (h1 << 27) | (h1 >> 37); h1 += h2; h1 = h1 * 5 + ulx[10];
    k2 *= ulx[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx[8]; h2 ^= k2;
    h2 = (h2 << 31) | (h2 >> 33); h2 += h1; h2 = h2 * 5 + ulx[11];
  }

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  uint64_t k1 = 0, k2 = 0;

  switch(len & 15) {
    case 15: k2 ^= (uint64_t)(tail[14]) << 48;  __attribute__ ((fallthrough));
    case 14: k2 ^= (uint64_t)(tail[13]) << 40;  __attribute__ ((fallthrough));
    case 13: k2 ^= (uint64_t)(tail[12]) << 32;  __attribute__ ((fallthrough));
    case 12: k2 ^= (uint64_t)(tail[11]) << 24;  __attribute__ ((fallthrough));
    case 11: k2 ^= (uint64_t)(tail[10]) << 16;  __attribute__ ((fallthrough));
    case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;   __attribute__ ((fallthrough));
    case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;  
             k2 *= ulx[9]; k2 = (k2 << 33) | (k2 >> 31); k2 *= ulx[8]; h2 ^= k2;
             __attribute__ ((fallthrough));
    case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;  __attribute__ ((fallthrough));
    case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;  __attribute__ ((fallthrough));
    case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;  __attribute__ ((fallthrough));
    case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;  __attribute__ ((fallthrough));
    case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;  __attribute__ ((fallthrough));
    case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;  __attribute__ ((fallthrough));
    case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;   __attribute__ ((fallthrough));
    case  1: k1 ^= (uint64_t)(tail[ 0]) << 0; // equiv to "case 1"  
             k1 *= ulx[8]; k1 = (k1 << 31) | (k1 >> 33); k1 *= ulx[9]; h1 ^= k1; break;
  };

  h1 ^= len; h2 ^= len;
  h1 += h2; h2 += h1;

  h1 ^= h1 >> 33; h1 *= ulx[6]; h1 ^= h1 >> 33; h1 *= ulx[7]; h1 ^= h1 >> 33;
  h2 ^= h2 >> 33; h2 *= ulx[6]; h2 ^= h2 >> 33; h2 *= ulx[7]; h2 ^= h2 >> 33;

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
    case 3: k ^= tail[2] << 16; __attribute__ ((fallthrough));
    case 2: k ^= tail[1] << 8;  __attribute__ ((fallthrough));
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
	acc += input * ulx[13];
  acc = (acc << 31) | (acc >> 33);
	acc *= ulx[12];
	return acc;
}

static uint64_t xxh64_merge_round (uint64_t acc, uint64_t val)
{
	val = xxh64_round (0, val);
	acc ^= val;
	acc = acc * ulx[12] + ulx[15];
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
    uint64_t v1 = h64 + ulx[12]+ ulx[13];
    uint64_t v2 = h64 + ulx[13];
    uint64_t v3 = h64 + 0;
    uint64_t v4 = h64 - ulx[12];

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
    h64  += ulx[16];
  }

  h64 += (uint64_t) len;

  while (p + 8 <= b_end) {
    const uint64_t k1 = xxh64_round (0, xxh_get_unaligned_le64(p));
    h64 ^= k1;
    h64 = ((h64 << 27) | (h64 >> 37)) * ulx[12] + ulx[15];
    p += 8;
  }
  if (p + 4 <= b_end) {
    h64 ^= (uint64_t) (xxh_get_unaligned_le32(p)) * ulx[12];
    h64 = ((h64 << 23)|(h64 >> 41)) * ulx[13] + ulx[14];
    p += 4;
  }
  while (p < b_end) {
    h64 ^= (*p) * ulx[16];
    h64 = ((h64<<11)|(h64>>53)) * ulx[12];
    p++;
  }
  h64 ^= h64 >> 33; h64 *= ulx[13]; h64 ^= h64 >> 29; h64 *= ulx[14]; h64 ^= h64 >> 32; // avalanche
  return h64;
}

