/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULLAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/* leo@uvigo: In most PRNGs which need a struct (vector and counter) the initial state of the vector is assumed 
 * "in equilibrium" (in a good initial state). The idea from sprng parallel implementation is to replace the "poor"
 * vector initialization (implemented in GSL, for instance) by an explicit sampling from another PRNG, since 
 * distinct vectors (of size 2^14+1 in case of GFSR4 e.g.) will lead to independent streams. Besides the weak point
 * of shift register generators is the possible correlation of vector elements given the single seed. 
 *
 * My idea was to increase the randomness by combining several PRNGs, in the initialization and the higher level
 * functions (in random_number.c). The algorithms are approximately in order of complexity, so that functions at the
 * bottom should not call functions at the top (at the risk of overflowing the stack ;). */

#include "random_number_gen.h"
#include "constant_random_lists.h" 

void rng_set_marsaglia_constants (uint32_t *m, uint32_t s);
#define RoL64(val, numbits) (((val) << (numbits)) | ((val) >> (64 - (numbits))))
#define RoL(val, numbits) (((val) << (numbits)) | ((val) >> (32 - (numbits))))

inline uint64_t
rng_get_taus (rng_taus_struct *r)
{ /* r->x[30] (6 vectors of size 4 or 5) */
  int i;
  uint64_t combined = 0ULL, 
           *A  = r->x, /* different (non-overlapping) elements of same vector */
           *q  = r->x + (2 * r->n), /* for this generator only r->n doesn't mean current iteration   */
           *rr = r->x + (3 * r->n), /* it means the number of elements (4 or 5, depending on stream) */
           *C  = r->x + (4 * r->n), 
           *s  = r->x + (5 * r->n);
  /* obs: "r" from the article is called "rr" here since I already called "r" the RNG... */

  for (i=0; i < r->n; i++) { /* order in TAUSWORTHE define: (A,q,r,C,s) where C=2^64 - 2^{64-k} and r=k-s */
    A[i] = ((((A[i] & C[i]) << s[i]) & 0xffffffffffffffffULL) ^ 
                   ((((A[i] << q[i]) & 0xffffffffffffffffULL)^A[i]) >> rr[i]));
    combined ^= A[i];
  }
  return combined;
}

void
rng_set_taus (rng_taus_struct *r, uint64_t seed, int stream)
{ /* r->x[30] (6 vectors of size 4 or 5) */
  int i;
  uint64_t *A, *k, *q; 
 
  rng_set_stream_taus (r, stream);
 
  A = r->x; /* different (non-overlapping) elements of same vector */
  k = r->x +      r->n;
  q = r->x + (2 * r->n);

  if (!seed) seed = 0x2f72b5f978acb838ULL; 

  seed = rng_randomize_array_64bits (r->x, r->n, seed, true); /* initialize vector using fixed (but quite unique to seed) table */
  seed = rng_randomize_array_64bits (r->x, r->n, seed, false); /* increase randomness */ 
  rng_twist_array_64bits (r->x, r->n, seed, 3);

  for (i=0; i < r->n; i++) /* Initial values should be larger or equal to 2^(64-k) */
    if (A[i] < ( 1ULL << (64 - k[i]))) A[i] += (1ULL << (64 - k[i]));

  for (i=0; i < r->n; i++) /* may not be necessary for our stream parameters */
    A[i] = ((((A[i] <<  q[i]) ^ A[i]) >> k[i]) ^A[i]);

  /* warm up */
  for (i=0; i < 10; i++) rng_get_taus (r);
}

void
rng_set_stream_taus (rng_taus_struct *r, int stream_number)
{
  int i;
  uint64_t *k, *q, *rr, *C, *s;
  
  stream_number %= 150; /* we have 150 distinct streams, 44 with five components and 106 with four */
  if (stream_number < 44) r->n = 5;
  else r->n = 4;

  k  = r->x +       r->n;/* different (non-overlapping) elements of same vector */
  q  = r->x + (2 * r->n);
  rr = r->x + (3 * r->n);
  C  = r->x + (4 * r->n);
  s  = r->x + (5 * r->n);

  /* number of elements of each stream: 20, 24, 4, 2, 92, 8 */
  if (stream_number < 20) { /* table 7 */
    for (i=0; i < r->n; i++) { q[i] = qTable76[0][i]; k[i] = kTable76[0][i]; s[i] = sTable76[stream_number][i]; }
  }
  else if (stream_number < 44) { /* table 6 */
    for (i=0; i < r->n; i++) { q[i] = qTable76[1][i]; k[i] = kTable76[1][i]; s[i] = sTable76[stream_number][i]; }
  }
  else if (stream_number < 48) { /* table 5 (1 to 4) */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[0][i]; k[i] = kTable543[0][i]; s[i] = sTable543[sn][i]; }
  }
  else if (stream_number < 50) { /* table 5 (5 to 6) */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[1][i]; k[i] = kTable543[1][i]; s[i] = sTable543[sn][i]; }
  }
  else if (stream_number < 142) { /* table 4 */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[2][i]; k[i] = kTable543[2][i]; s[i] = sTable543[sn][i]; }
  }
  else { /* (stream_number < 150) table 3 */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[3][i]; k[i] = kTable543[3][i]; s[i] = sTable543[sn][i]; }
  }

  for (i=0; i < r->n; i++) { /* auxiliary variables, common to all streams ("rr" is original "r" in article) */
    rr[i] = k[i] - s[i];
    C[i] = Cmask[k[i] - 36]; /* we start at k=36 and k < 64 always */
  }
}

inline uint64_t
rng_get_xorshift (rng_xorshift_struct *r)
{ /* x[65]  but last one (r->x[64]) is not part of the stream (is aux variable) */
  int i;
  uint64_t t, v;

  i = r->n = (r->n + 1) & 63;
  t = r->x[i];
  v = r->x[(i+11) & 63];     /* Index is (i-53) mod 64 */
  t ^= t << 33;  t ^= t >> 26; /* (I + L^a)(I + R^b) */
  v ^= v << 27;  v ^= v >> 29; /* (I + L^c)(I + R^d) */
  r->x[i] = (v ^= t) & 0xffffffffffffffffULL; /* Update circular array */
  /* 0x61c8864680b583ebULL = odd approximation to 2**64*(3-sqrt(5))/2. */
  r->x[64] += 0x61c8864680b583ebULL;/* Update Weyl generator */

  return (v + (r->x[64] ^ (r->x[64] >> 27))) & 0xffffffffffffffffULL;
}

void
rng_set_xorshift (rng_xorshift_struct *r, uint64_t seed)
{ /* x[65]  but last one (r->x[64]) is not part of the stream (is aux variable) */
  int i, j;
  uint64_t t;

  if (!seed) seed = 0x1db9b83a20cc6503ULL; 

  seed = rng_randomize_array_64bits (r->x, 64, seed, true); // initialise array with table of random numbers
  rng_twist_array_64bits (r->x, 64, seed, 4);

  /* Avoid correlations for close seeds; Recurrence has period 2**64-1 */
  for (r->x[64] = seed, i = 0; i < 64; i++) { /* Initialise circular array */
    rng_get_brent_64bits (&seed);
    r->x[i] = (r->x[i] ^ seed) + (r->x[64] += 0x61c8864680b583ebULL);
  }

  for (i = 63, j = 0; j < 256; j++) { /* Discard first 256 results */
    i = (i+1) & 63;
    t = r->x[i];
    seed = r->x[(i+11) & 63];               /* Index is (i-53) mod 64 */
    t ^= t << 33; t ^= t >> 26;             /* (I + L^a)(I + R^b) */
    seed ^= seed << 27; seed ^= seed >> 29; /* (I + L^c)(I + R^d) */
    r->x[i] = (seed ^ t) & 0xffffffffffffffffULL; /* Update circular array */
  }
  r->n = i;
  return;
}

inline uint64_t
rng_get_mt19937 (rng_mt19937_struct *r)
{
  static const uint64_t mag01[2]={ 0ULL, 0xB5026F5AA96619E9ULL}; /* this is magic vector, don't change */
  uint64_t x;

  if (r->n >= 312) { /* generate all 312 words at once */
    int i;
    for (i = 0; i < 156; i++) {
      x = (r->x[i] & 0xFFFFFFFF80000000ULL)| (r->x[i+1] & 0x7FFFFFFFULL);
      r->x[i] = r->x[i+156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    for (; i < 311; i++) {
      x = (r->x[i] & 0xFFFFFFFF80000000ULL) | (r->x[i+1] & 0x7FFFFFFFULL);
      r->x[i] = r->x[i-156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    x = (r->x[311] & 0xFFFFFFFF80000000ULL) | (r->x[0] & 0x7FFFFFFFULL);
    r->x[311] = r->x[155] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    r->n = 0; 
  }
  x = r->x[r->n++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x;
}

void
rng_set_mt19937 (rng_mt19937_struct *r, uint64_t seed)
{
  if (!seed) seed = 0x33cba2d924f83a89ULL;
  r->n = 313;
  seed = rng_randomize_array_64bits (r->x, 312, seed, true); 
  rng_twist_array_64bits (r->x, 312, seed, 8);
}

inline uint32_t
rng_get_mt19937ar (rng_mt19937ar_struct *r)
{
  static const uint32_t mag01[2]={ 0ULL, 0x9908b0dfUL}; /* this is magic vector, don't change */
  uint32_t i, y;
  if (r->n >= 624) { /* generate N words at one time */
    for (i = 0; i < 227; i++) {
      y = (r->x[i] & 0x80000000UL) | (r->x[i+1] & 0x7fffffffUL);
      r->x[i] = r->x[i+397] ^ (y >> 1) ^ mag01[y & 1UL];
    }
    for (; i < 623; i++) {
      y = (r->x[i] & 0x80000000UL) | (r->x[i+1] & 0x7fffffffUL);
      r->x[i] = r->x[i-277] ^ (y >> 1) ^ mag01[y & 1UL];
    }
    y = (r->x[623] & 0x80000000UL) | (r->x[0] & 0x7fffffffUL);
    r->x[623] = r->x[396] ^ (y >> 1) ^ mag01[y & 1UL];
    r->n = 0;
  }
  y = r->x[r->n++];
  /* Tempering */
  y ^= (y >> 11); y ^= (y << 7) & 0x9d2c5680UL; y ^= (y << 15) & 0xefc60000UL; y ^= (y >> 18);
  return y;
}

void
rng_set_mt19937ar (rng_mt19937ar_struct *r, uint64_t seed)
{
  if (!seed) seed = 0x22b417f12d5f1072ULL;
  r->n = 625;
  seed = rng_randomize_array_32bits (r->x, 624, seed, true); 
  rng_twist_array_32bits (r->x, 624, seed, 8);
}

inline uint32_t
rng_get_gfsr4 (rng_gfsr4_struct *r)
{ /* r->x[16384] */
  r->n = ((r->n)+1) & 16383;
  return r->x[r->n] =
  (r->x[(r->n + 15913) & 16383] ^ r->x[(r->n + 14798) & 16383] ^ r->x[(r->n + 9396)  & 16383] ^ r->x[(r->n + 6695)  & 16383]) & 0xffffffffUL;
}

void
rng_set_gfsr4 (rng_gfsr4_struct *r, uint64_t seed)
{ /* r->x[16384] */
  if (!seed) seed = 0x06e346963311b7e3ULL;
  seed = rng_randomize_array_32bits (r->x, 16384, seed, true); 
  rng_twist_array_32bits (r->x, 16384, seed, 4);
  r->n = 0;
}

/* * * below, generators whose initial states depend only on simple functions (single-variable PRNGs etc.)  */

uint32_t
rng_get_diaconis (rng_diaconis_struct *r)
{ /*r->x[128] */
  /* two implicid bits version. Period lenght: (2**R + 1) * 2**31. */
  uint32_t b0, sr, ss, br, bs;
  /* fib(n) = fib(n-R) * fib(n-S); with all fib() odd. */
  r->n--;
  br = r->x[(r->n + 127) & 127];
  bs = r->x[(r->n + 30)  & 127];
  sr = br & 1; br ^= sr;
  ss = bs & 1; bs ^= ss;
  b0 = 4 * br * bs;
  if (sr) br *= 3;
  if (ss) bs *= 3;
  b0 += br + bs + sr + ss;
  r->x[r->n & 127] = b0;
  return b0 + (b0 >> 16); /* low bit improvement */
}

uint32_t
rng_get_diaconis_onebit (rng_diaconis_struct *r)
{ /*r->x[128] */
/* one implicid bit version. Period lenght: (2**R + 1) * 2**30. */
  uint32_t b0, br, bs;
  /* fib(n) = fib(n-R) * fib(n-S); with all fib() odd. */
  r->n--;
  br = r->x[(r->n + 127) & 127];
  bs = r->x[(r->n + 30)  & 127];
  b0 = br + bs + 2*br*bs;
  r->x[r->n & 127] = b0;
  return b0 + (b0 >> 16); /* low bit improvement */
}

void
rng_set_diaconis (rng_diaconis_struct *r, uint64_t seed)
{ /*r->x[128] */
  r->n = 0;
  int i;
  if (!seed) seed = 0x1c9cc5643af25686ULL;
  seed = rng_randomize_array_32bits (r->x, 127, seed, true);  /* initialise vector */
  seed = rng_randomize_array_32bits (r->x, 127, seed, false); /* increase randomness */
  for (i=0; i < 128; i++) r->x[i] |= 1; // initial state must be odd 
}

/* Makoto Matsumoto & Y. Kurita, Twisted GFSR Generators II, ACM Trans. Model. Comput. Simul., 4 (1994) 254-266 */
uint32_t
rng_get_tt800 (rng_tt800_struct *r)
{ /* r->x[25] */
  static const uint32_t mag01[2]={ 0x0, 0x8ebfd028}; /* this is magic vector, don't change */
  uint32_t y;

  if (r->n >= 25) { /* generate N words at one time */
    int i;
    for (i = 0; i < 18; i++) r->x[i] = r->x[i+7]  ^ (r->x[i] >> 1) ^ mag01[r->x[i] % 2];
    for (; i < 25; i++)      r->x[i] = r->x[i-18] ^ (r->x[i] >> 1) ^ mag01[r->x[i] % 2];
    r->n = 0;
  }
  y = r->x[r->n];
  y ^= (y << 7)  & 0x2b5b2500UL; /* s and b, magic vectors */
  y ^= (y << 15) & 0xdb8b0000UL; /* t and c, magic vectors */
  r->n++;
  return (y ^ (y >> 16)) & 0xffffffffUL;
}

void
rng_set_tt800 (rng_tt800_struct *r, uint64_t seed)
{ /* r->x[25] */
  if (!seed) seed = 0x273a3292263c330eULL;
  seed = rng_randomize_array_32bits (r->x, 25, seed, true); 
  rng_twist_array_32bits (r->x, 25, seed, 10);
  r->n = 26;
}

uint32_t
rng_get_lfib4 (rng_lfib4_struct *r)
{ /* r->x[256] */  /* set = vector filling */
  r->n = (r->n + 1) & 255;
  return r->x[r->n] = r->x[r->n & 255] + r->x[(r->n + 58) & 255] + r->x[(r->n + 119) & 255] + r->x[(r->n + 178) & 255];
}

void
rng_set_lfib4 (rng_lfib4_struct *r, uint64_t seed)
{
  if (!seed) seed = 0x395894461ab4c493ULL;
  rng_randomize_array_32bits (r->x, 256, seed, true); /* increase randomness by concatenation */
  rng_twist_array_32bits (r->x, 256, seed, 7);
  r->n = 0;
}

uint32_t
rng_get_swb (rng_swb_struct *r)
{ /* r->x[258] (x[256] + two aux variables) */
  r->x[256] = r->x[(r->n + 15) & 255]; /* x in original algorithm */
  r->x[(r->n + 237) & 255] = r->x[256] - (r->x[257] = r->x[(r->n+1) & 255] + (r->x[256] < r->x[257]));
  r->n = (r->n + 1) & 255;
  return r->x[r->n];
}

void
rng_set_swb (rng_swb_struct *r, uint64_t seed)
{
  r->n = 0;
  r->x[256] = r->x[257] = 0UL;
  if (!seed) seed = 0x123733ca1b72b747ULL;
  rng_randomize_array_32bits (r->x, 256, seed, true); 
  rng_twist_array_32bits (r->x, 256, seed, 6);
}

/* Panneton, L'Ecuyer, and Matsumoto WELLRNG1024a */
uint32_t
rng_get_well1024 (rng_well1024_struct *r)
{ /* r->x[32] */ /* set = vector filling */
  uint32_t z0, z1, z2;
  z0 = r->x[(r->n + 31) & 0x0000001fU];
  z1 = r->x[r->n] ^ (r->x[(r->n + 3) & 0x0000001fU] ^ (r->x[(r->n + 3) & 0x0000001fU] >> 8));
  z2 = ((r->x[(r->n + 24) & 0x0000001fU] ^ (r->x[(r->n + 24) & 0x0000001fU] << 19)) ^
        (r->x[(r->n + 10) & 0x0000001fU] ^ (r->x[(r->n + 10) & 0x0000001fU] << 14)));
  r->x[r->n] = z1 ^ z2;
  r->x[(r->n + 31) & 0x0000001fU] = (z0 ^ (z0 << 11)) ^ (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
  r->n = (r->n + 31) & 0x0000001fU;
  return r->x[r->n];
}

void
rng_set_well1024 (rng_well1024_struct *r, uint64_t seed)
{
  r->n = 0;
  if (!seed) seed = 0x165da3840ea4bba9ULL;
  seed = rng_randomize_array_32bits (r->x, 32, seed, true); 
  rng_twist_array_32bits (r->x, 32, seed, 5);
}

/* * * Simple generators (single value or simple vectors of size 2 or 4) */
/* better not to call mix() on these functions to avoid loops, since mix call these */ 

uint64_t
rng_get_gamerand64 (uint64_t *game)
{ /* game[2] */
  game[0] = (game[0] << 32) + (game[0] >> 32); game[0] += game[1]; game[1] += game[0];
  return game[0];
}

/* http://prng.di.unimi.it/xoroshiro128plusplus.c; commented out is http://xoroshiro.di.unimi.it/xoroshiro128plus.c */
uint64_t 
rng_get_xoroshiro128 (uint64_t *s) 
{
  uint64_t s1 = s[1], result;
  // result = s[0] + s1; // 128+ 
  // result = RoL64(s[0] * 5, 7) * 9; // 128*
	result = RoL64(s[0] + s1, 17) + s[0]; // 128++
  s1 ^= s[0];
  // s[0] = RoL64(s[0], 24) ^ s1 ^ (s1 << 16); s[1] = RoL64(s1, 37);  // 128+ V 2018
  // s[0] = RoL64(s[0], 55) ^ s1 ^ (s1 << 14); s[1] = RoL64(s1, 36);  // 128+ V 2016
  // s[0] = RoL64(s[0], 24) ^ s1 ^ (s1 << 16); s[1] = RoL64(s1, 37);  // 128* 
     s[0] = RoL64(s[0], 49) ^ s1 ^ (s1 << 21); s[1] = RoL64(s1, 28);  // 128++
  return result;
}

uint64_t 
rng_get_xoroshiro128p (uint64_t *s) 
{
  uint64_t s1 = s[1], result;
  result = s[0] + s1; // 128+ 
  // result = RoL64(s[0] * 5, 7) * 9; // 128*
  s1 ^= s[0];
  s[0] = RoL64(s[0], 24) ^ s1 ^ (s1 << 16); s[1] = RoL64(s1, 37);  // 128+ V 2018
  // s[0] = RoL64(s[0], 24) ^ s1 ^ (s1 << 16); s[1] = RoL64(s1, 37);  // 128* 
  return result;
}

uint64_t 
rng_get_xoroshiro128s (uint64_t *s) 
{
  uint64_t s1 = s[1], result;
  result = RoL64(s[0] * 5, 7) * 9; // 128*
  s1 ^= s[0];
  s[0] = RoL64(s[0], 24) ^ s1 ^ (s1 << 16); s[1] = RoL64(s1, 37);  // 128* 
  return result;
}

/* Equivalent to 2^64 calls (e.g. generate 2^64 non-overlapping subsequences for parallel computations) */
void
rng_jump_64_xoroshiro128 (uint64_t *s) 
{
  uint64_t s0 = 0, s1 = 0;
  int b;

  for (b = 0; b < 64; b++) { if (0xdf900294d8f554a5ULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; } rng_get_xoroshiro128 (s); }
  for (b = 0; b < 64; b++) { if (0x170865df4b3201fcULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; } rng_get_xoroshiro128 (s); }

  s[0] = s0; s[1] = s1;
}

/* Equivalent to 2^96 calls (e.g. generate 2^32 non-overlapping subsequences for parallel computations) */
void
rng_jump_96_xoroshiro128 (uint64_t *s) 
{
  uint64_t s0 = 0, s1 = 0;
  int b;
  for (b = 0; b < 64; b++) { if (0xd2a98b26625eee7bULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; } rng_get_xoroshiro128 (s); }
  for (b = 0; b < 64; b++) { if (0xdddf9b1090aa7ac1ULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; } rng_get_xoroshiro128 (s); }
  s[0] = s0; s[1] = s1;
}

/* http://prng.di.unimi.it/xoshiro256plusplus.c */
uint64_t 
rng_get_xoroshiro256 (uint64_t *s) // 4 x 64bits 
{
  uint64_t result, t = s[1] << 17;
  result = RoL64(s[0] + s[3], 23) + s[0]; // 256++
  // result = RoL64(s[1] * 5, 7) * 9; // 256*
	s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
  s[2] ^= t; s[3] = RoL64(s[3], 45);
	return result;
}

/* Equivalent to 2^128 calls (e.g. generate 2^128 non-overlapping subsequences for parallel computations) */
void
rng_jump_128_xoroshiro256 (uint64_t *s) 
{	
  uint64_t s0 = 0,	s1 = 0, s2 = 0, s3 = 0;
  int b;
  for (b = 0; b < 64; b++) { if (0x180ec6d33cfd0abaULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0xd5a61266f0c9392cULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0xa9582618e03fc9aaULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0x39abdc4529b1661cULL & (1ULL << b)) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3;
}

/* Equivalent to 2^192 calls (e.g. generate 2^64 non-overlapping subsequences for parallel computations) */
void
rng_jump_192_xoroshiro256 (uint64_t *s) 
{
  uint64_t s0 = 0,	s1 = 0, s2 = 0, s3 = 0;
  int b;
  for (b = 0; b < 64; b++) { if (0x76e15d3efefdcbbfULL & 1ULL << b) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0xc5004e441c522fb3ULL & 1ULL << b) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0x77710069854ee241ULL & 1ULL << b) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  for (b = 0; b < 64; b++) { if (0x39109bb02acbe635ULL & 1ULL << b) { s0 ^= s[0]; s1 ^= s[1]; s2 ^= s[2]; s3 ^= s[3]; } rng_get_xoroshiro256 (s); }
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3;
}

uint32_t
rng_get_gamerand (uint32_t *game)
{ /* game[2] */
  game[0] = (game[0] << 16) + (game[0] >> 16); game[0] += game[1]; game[1] += game[0];
  return game[0];
}

void
rng_set_gamerand (uint32_t *game, uint64_t seed)
{ /* game[2] */
  int i;

  if (!seed) seed = 0x2a6256952c056553ULL; 

  game[0] = seed; 
  game[1] = biomcmc_hashint64_salted (seed, /*salt*/ 0);
  for (i = 0; i < 4; i++) rng_get_brent (game + 1);
}

/* Marsaglia's Super-Duper (two-components multiply-with-carry) algortihm */
uint32_t
rng_get_marsaglia (uint32_t *m)
{ /* m[4] */
  return (((m[0] = m[2] * (m[0] & 0xffffUL) + (m[0] >> 16)) << 16) + 
          ((m[1] = m[3] * (m[1] & 0xffffUL) + (m[1] >> 16 )) & 0xffffUL)) & 0xffffffffUL;
}

void
rng_set_marsaglia (uint32_t *m, uint64_t seed)
{ /* m[4] */
  if (!seed) seed = 0x2f3e89e73907c3f8ULL;
  m[0] = (uint32_t) (seed); 
  m[1] = 1UL + (uint32_t) (biomcmc_hashint64_salted (seed + 1, /*salt*/ 0) >> 4);
  
  rng_set_marsaglia_constants (m, seed); // seed can work as 'stream', since all prime pairs are explored in order
}

void
rng_set_marsaglia_constants (uint32_t *m, uint32_t s)
{
  uint16_t idx1, idx2;
  /* choose two *distinct* marsaglia_constants[], from a pool of 81. We have therefore 80 x 81 possible streams */
  idx1 = s % marsaglia_constants_size; // 0...80
  idx2 = (s/marsaglia_constants_size) % (marsaglia_constants_size - 1); // 0...79
  if (idx1 == idx2) idx2 = marsaglia_constants_size - 1; 
  m[2] = marsaglia_constants[idx1];
  m[3] = marsaglia_constants[idx2];
  if (m[0] == m[1]) m[1] *= 69069UL; // this function is called even if rng_set is not, so check is here
}

uint64_t
rng_get_std61 (uint64_t *x)
{
  (*x) = ((*x) >> 31) + (((*x) << 30) & 0x1fffffffffffffffULL) - ((*x) >> 42) - (((*x) << 19) & 0x1fffffffffffffffULL);
  if ((int64_t) (*x) < 0) (*x) += 0x1fffffffffffffffULL;
  return (*x);
}

/* The Minimal Portable Random Number Generator (32 bits)
   a = 7^5 = 16807   m = 2^31 - 1 = 2147483647 = 0x7fffffff ;   x[n+1] = a * x[n] (mod m) */
uint32_t
rng_get_std31 (uint32_t *x)
{
  uint32_t zh, zl, z;
  z = (*x) << 1; 
  zl = z & 0xffffUL; 
  zh = z >> 16; 
  zl *= 48271; zh *= 48271; 
  zh += zl >> 16; 
  zl = (zl & 0xffffUL) + (zh << 16);
  zh = (zh >> 16) << 1UL; 
  zl += zh;
  if (zh > zl) zl += 2UL;
  (*x) = zl >> 1;
  return (*x);
}
/* http://prng.di.unimi.it/xoroshiro64starstar.c */
uint32_t 
rng_get_xoroshiro64 (uint32_t *x)  // s[2]
{
	uint32_t result, s1 = x[1];
	result = RoL(x[0] * 0x9E3779BB, 5) * 5; // 64**
  //result = x[0] * 0x9E3779BB; // 64* lowest six bits have low linear complexity
	s1 ^= x[0];
	x[0] = RoL(x[0], 26) ^ s1 ^ (s1 << 9); x[1] = RoL(s1, 13); 
	return result;
}

uint32_t
rng_get_shr (uint32_t *x)
{
  if (!(*x)) (*x) = 0x12585408ULL;
  (*x) ^= ((*x) << 17); (*x) ^= ((*x) >> 13); (*x) ^= ((*x) << 5);
  return (*x);
}

uint32_t
rng_get_brent (uint32_t *x)
{
  if (!(*x)) (*x) = 0x0eba937eULL;
  (*x) ^= (*x) << 10; (*x) ^= (*x) >> 15; (*x) ^= (*x) << 4;  (*x) ^= (*x) >> 13;
  return *x;
}

uint64_t
rng_get_brent_64bits (uint64_t *x)
{
  if (!(*x)) (*x) = 0x20a4f71433a9481fULL; 
  (*x) ^= (*x) << 10; (*x) ^= (*x) >> 15; (*x) ^= (*x) << 4;  (*x) ^= (*x) >> 13;
  return *x;
}

// https://github.com/FastFilter/xor_singleheader/blob/master/include/xorfilter.h
uint64_t 
rng_get_splitmix64 (uint64_t *x) 
{
  uint64_t z = (*x += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

uint32_t
rng_get_cong (uint32_t *x)
{
  uint32_t b1, b2; /* we use only the leading half bits (the last ones are too regular */
  b1 = ((*x) = (69069UL * (*x)) + 1234567) >> 16;
  b2 = ((*x) = (69069UL * (*x)) + 1234567) & 0xffff0000UL;
  return (b1 | b2); /* notice that here we don't update x with (b1|b2) */
}

uint32_t
rng_get_cong_many (uint32_t *x)
{ /* quick-and-dirty (well, not so quick)  from GSL */
  uint32_t j, t = 0UL, bit =  0x80000000UL;
  for (j = 0; j < 32; j++) {
    rng_get_cong (x);
    if ((*x) & 0x80000000UL) t |= bit;
    bit >>= 1;
  }
  return t; /* we don't return x, despite it's being being updated */
}

uint64_t
rng_twist_array_32bits (uint32_t *a, uint32_t n_a, uint64_t seed, uint64_t stream)
{
  uint32_t i, im, s32, mars[4];
  uint64_t x64, sx[2];
  if (!seed) seed = 0x085764f60bc8797eULL;

  /* this shuffling runs every time (notice that xoroshiro is 64bits) */
  sx[0] = (uint64_t) a[0] | 1ULL; sx[1] = seed;
  for (i = 0; i < n_a - 1; i++) {
    x64 = rng_get_xoroshiro128 (sx);
    a[i]   ^= (uint32_t) (x64); 
    a[i+1] ^= (uint32_t) (x64 >> 32);
  }
  seed = sx[0];

  /* initialize Marsaglia's Super-Duper generator */
  mars[0] = s32 = (uint32_t) (sx[1]);
  mars[1] = (uint32_t) (sx[1] >> 32);
  rng_set_marsaglia_constants (mars, (uint32_t) (seed));

  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 1
    im = (i + n_a - 1) % n_a; /* im = i - 1 but warping around (so that zero minus one is last element) */ 
    a[i] ^= ((a[im] ^ (a[im] >> 30)) * 0x071a9b16UL) ^ rng_get_marsaglia (mars);
  }
  mars[0] = RoL(seed, 6);
  stream >>= 1;
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 2
    im = (i + n_a - 2) % n_a; /* im = i - 2 but with warp-around */ 
    a[i] ^= ((a[im] ^ (a[im] >> 28)) * 0x09f68db7UL) ^ rng_get_std31 (mars); // std31 uses just mars[0]
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 3 = 4
    a[i] ^= rng_get_marsaglia (mars) ^ rng_get_shr (&s32); 
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 1; i < n_a; i++) { // bit 4 = 8
    a[i] ^= (1812433253UL * (a[i-1] ^ (a[i-1] >> 30)) + i); // used by mt19937
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 5 = 16
    a[i] ^= rng_get_cong_many (mars) ^ rng_get_shr (&s32); // cong_many uses just mars[0] 
  }
  return seed;
}

uint64_t
rng_twist_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed, uint64_t stream)
{ /* modified from MT19937 (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) */
  uint32_t i, im;
  uint64_t s1, sx[2];

  if (!seed) seed = 0x1b422e75022494afULL;
  s1 = biomcmc_hashint64_salted (seed, 4); 

  /* this shuffling runs every time */
  sx[0] = a[0]; sx[1] = seed;
  for (i = 0; i < n_a; i++) a[i] ^= rng_get_xoroshiro128 (sx);
  seed = sx[0];

  /* the ones below only if bitmask is present */
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 1
    im = (i + n_a - 1) % n_a; /* im = i - 1 but warping around (so that zero minus one is last element) */ 
    a[i] ^= ((a[im] ^ (a[im] >> 62)) * 0x72b5f90702b838ULL) + rng_get_std61 (&seed);
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 2
    im = (i + n_a - 2) % n_a; /* im = i - 2 but with warp-around */ 
    a[i] ^= ((a[im] ^ (a[im] >> 58)) * 0x548ba82e1b6ce1ULL) ^ rng_get_std61 (&seed);
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 0; i < n_a; i++) { // bit 3
    a[i] ^= rng_get_std61 (&seed) ^ rng_get_brent_64bits (&s1);
  }
  stream >>= 1;
  if (stream & 1UL) for (i = 1; i < n_a; i++) { // bit 4
    a[i] ^= (6364136223846793005ULL * (a[i-1] ^ (a[i-1] >> 62)) + i); // used by mt19937
  }
  
  return seed; /* so that this function can be called several times */
}

uint64_t
rng_randomize_array_32bits (uint32_t *a, uint32_t n_a, uint64_t seed, bool first_time)
{
  uint32_t i, t[2], m[4];

  if (!seed) seed = (69069 * n_a) + 69069;
  
  if (first_time) { 
    /* spice table needs 4 32bit ints */
    m[0] = (uint32_t) (seed); m[1] = (uint32_t)(seed >> 32);
    seed = biomcmc_hashint64_salted (seed, 3);
    m[2] = (uint32_t) (seed); m[3] = (uint32_t)(seed >> 32);
    biomcmc_salt_vector32_from_spice_table (a, n_a, m);
  }
  t[0] = a[0] + 1ull; t[1] = a[1] + 1ull;
  for (i = 0; i < n_a; i++) a[i] ^= (rng_get_shr (t) ^ rng_get_std31 (t+1));

  return ((uint64_t)(t[0]) << 32) | t[1];
}

uint64_t
rng_randomize_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed, bool first_time)
{
  uint64_t t[3]; 
  uint32_t i, m[4];

  if (!seed) seed = 69069 * n_a;
  
  if (first_time) { 
    /* spice table needs 4 32bit ints (even the 64bits one) */
    m[0] = (uint32_t) (seed); m[1] = (uint32_t)(seed >> 32);
    seed = biomcmc_hashint64_salted (seed, 3);
    m[2] = (uint32_t) (seed); m[3] = (uint32_t)(seed >> 32);
    biomcmc_salt_vector64_from_spice_table (a, n_a, m);
  }
  t[0] = a[0] + 1ull; t[1] = a[1] + 1ull; t[2] = a[2] + 1ull;
  for (i = 0; i < n_a; i++) a[i] ^= rng_get_brent_64bits (t) ^ rng_get_xoroshiro128 (t+1);
  return t[0];
}

#undef RoL64
#undef RoL

