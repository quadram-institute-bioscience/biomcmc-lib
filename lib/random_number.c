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

#include "random_number.h"

/*! \brief pointer to pseudo-random number generator (should point to real stream, even when there are several) */
biomcmc_rng biomcmc_random_number; /* Actual definition here with external declaration on random_number.h (to 
                                      guarantee a single copy) */


void
biomcmc_random_number_init (uint64_t seed)
{
  if (biomcmc_random_number) { biomcmc_random_number->ref_counter++; return;} // could return error, but let's asume calling function is overzealous 

  /* seed != 0 only when debugging; defaults to larger stream (=0) seeded by current date */
  if (seed) biomcmc_random_number = new_biomcmc_rng (seed, 0);
  else      biomcmc_random_number = new_biomcmc_rng (biomcmc_rng_get_initial_seed (), 0);
}

void
biomcmc_random_number_finalize (void)
{
  if (!biomcmc_random_number) return;
  del_biomcmc_rng (biomcmc_random_number);
  biomcmc_random_number = NULL;
}

void
del_biomcmc_rng (biomcmc_rng r)
{
  if(--r->ref_counter) return;
  if (r) free (r);
}

uint64_t
biomcmc_rng_get_initial_seed (void)
{
  int64_t timeseed[2];
  uint64_t low, high, top, pid1, pid2;

  biomcmc_get_time (timeseed);
  pid1 = (uint64_t)  getpid (); /* may be a small number; return type is pid_t (=int) */
  pid2 = (uint64_t) getppid (); /* parent process ID */
  if (!pid1) pid1 = 1ULL;
  if (!pid2) pid2 = 1ULL;

  /* get first 16 bits since usec may be < 1e6 (~20 bits) in systems without posix timers */
  low  = (((uint64_t)(timeseed[1]) * pid1) & 0xffffULL);
  /* use the lowest 32 bits of time in seconds as highest 48 bits (with lowest 16 masked) */
  high = (((uint64_t)(timeseed[0]) * pid2) << 16) & 0xffffffff0000ULL;
  /* highest 16 bits */
  top  = (uint64_t) (biomcmc_hashint_salted ((int)(timeseed[0] + timeseed[1]), /*salt*/ 6)) << 48;
//  fprintf (stderr, "seed %ju\n", (uintmax_t) low|high|top);

  return (uint64_t) (low | high | top);
}

/*
uint64_t *seed_rng (void)// check also arc4random (from BSD)
{
  int fp = open("/dev/urandom", O_RDONLY);
  uint64_t seed[4]; unsigned pos = 0; // BUGGY: should malloc if we really want to return address
  if (fp>=0) while (pos < sizeof(seed)) { int amt = read(fp, (char *) &seed + pos, sizeof(seed) - pos); if (amt > 0) pos += amt; else break; }
  close(fp); return seed;
}
*/

biomcmc_rng
new_biomcmc_rng (uint64_t seed, int stream_number)
{
  int i;
  uint64_t useed = seed;
  biomcmc_rng r = (biomcmc_rng) biomcmc_malloc (sizeof (struct biomcmc_rng_struct));
  
  rng_set_taus (&(r->taus), useed, stream_number);
  rng_get_brent_64bits (&useed); /* fast (one step) PRNG to change seed */
  rng_set_mt19937 (&(r->mt), useed); /* receive modified seed (but even same seed should work) */

  for (i = 0; i < 32; i++) { rng_get_taus (&(r->taus)); rng_get_mt19937 (&(r->mt)); }

  /* rnorm 64 and 32 bits (since each normal draw produces two values), 32 bits and 16 bits temp values */
  r->have_rnorm32 = r->have_rnorm64 = r->have_bit32 = false;
  r->algorithm = 0;
  return r;
}

biomcmc_rng
new_biomcmc_rng_with_parallel_seeds (uint64_t seed, int stream_number)
{
  int i;
  biomcmc_rng r; /* pointer (which will be visible globally through biomcmc_random_number) */
  rng_xorshift_struct xor; /* local variable */
  uint64_t this_seed, useed = seed;
  /* in parallel environments this function will initalize the local PRNG after receiving the seed */

  rng_set_xorshift (&xor, useed); /* initialize XOR-shift (our seed vector) and skip first elements */
  for (i=0; i < 32; i++)            rng_get_xorshift (&xor); /* "tempering" */
  for (i=0; i < stream_number; i++) rng_get_xorshift (&xor); /* avoid seeding two streams with same seed */
  
  this_seed = rng_get_xorshift (&xor); /* each stream will receive a distinct value from xorshift */

  r = new_biomcmc_rng (this_seed, stream_number);
  /* global variable ::random_number points to actual stream; some calling function can change it later */
  biomcmc_random_number = r; 

  return r;
}

inline double
biomcmc_rng_snorm32 (void)
{
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  if (biomcmc_random_number->have_rnorm32) {
    biomcmc_random_number->have_rnorm32 = false;
    return biomcmc_random_number->rnorm32;
  }
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    u = -1. + 2. * (((double) biomcmc_rng_get_32 ()) / 4294967295.0); /* 2^32 - 1 => 1 included */
    v = -1. + 2. * (((double) biomcmc_rng_get_32 ()) / 4294967295.0); /* 2^32 - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  biomcmc_random_number->rnorm32 = u * s;
  biomcmc_random_number->have_rnorm32 = true;
  return v * s;
}

inline double
biomcmc_rng_snorm (void)
{
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  if (biomcmc_random_number->have_rnorm64) {
    biomcmc_random_number->have_rnorm64 = false;
    return biomcmc_random_number->rnorm64;
  } 
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    /*    u and v are between -1 and +1 => -1 + 2 * U       */
    u = -1. + 2. * (biomcmc_rng_get_52 () / 4503599627370495.0); /* (2^52) - 1 => 1 included */
    v = -1. + 2. * (biomcmc_rng_get_52 () / 4503599627370495.0); /* (2^52) - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  biomcmc_random_number->rnorm64 = u * s;
  biomcmc_random_number->have_rnorm64 = true;
  return v * s;
}

inline double
biomcmc_rng_unif32 (void)
{
  return (biomcmc_rng_get_32 () / 4294967295.0); /* 2^32 - 1 => 1 included */
}

inline double
biomcmc_rng_unif (void)
{
  return (biomcmc_rng_get_52 () / 4503599627370495.0); /* 2^52 - 1 => 1 included */
}

inline double
biomcmc_rng_unif_pos32 (void)
{
  double x;
  do { x = (biomcmc_rng_get_32 () / 4294967295.0); } while (x < 2 * DBL_MIN);
  return x;
}

inline double
biomcmc_rng_unif_pos (void)
{
  double x;
  do { x = (biomcmc_rng_get_52 () / 4503599627370495.0); } while (x < 2 * DBL_MIN);
  return x;
}

inline uint32_t
biomcmc_rng_unif_int (uint32_t n)
{
  uint32_t scale;
  uint32_t k;
  if (!n) biomcmc_error ("n must be larger than zero in uniform random number generator [32bits]");
  scale = 0xffffffffU/n;

  do { k = biomcmc_rng_get_32 ()/scale; } while (k >= n);
  return k;
}

inline uint64_t
biomcmc_rng_unif_int64 (uint64_t n)
{
  uint64_t scale;
  uint64_t k;
  if (!n) biomcmc_error ("n must be larger than zero in uniform random number generator [64bits]");
  scale = 0xffffffffffffffffULL/n;

  do { k = biomcmc_rng_get ()/scale; } while (k >= n);
  return k;
}

void
biomcmc_rng_set_next_algorithm (void)
{
   biomcmc_random_number->algorithm = (biomcmc_random_number->algorithm+1)%10;
}

void
biomcmc_rng_set_algorithm (uint8_t algo)
{
   biomcmc_random_number->algorithm = algo%10; 
}

uint8_t
biomcmc_rng_get_algorithm (void)
{
   return biomcmc_random_number->algorithm; 
}

inline uint64_t
biomcmc_rng_get (void)
{
  switch (biomcmc_random_number->algorithm) {
    case 0:
      return rng_get_mt19937 (&(biomcmc_random_number->mt)); // best dieharder results
    case 1:
      return rng_get_taus (&(biomcmc_random_number->taus));
    case 2:
      return (rng_get_taus (&(biomcmc_random_number->taus)) ^ rng_get_mt19937 (&(biomcmc_random_number->mt)));
    case 3: 
      return rng_get_xoroshiro128p (&(biomcmc_random_number->mt.x[0])); // 2 vars 
    case 4: 
      return rng_get_xoroshiro128s (&(biomcmc_random_number->mt.x[4])); // 2 vars 
    case 5:
      return rng_get_xoroshiro128 (&(biomcmc_random_number->mt.x[8])); // 2 vars 
    case 6:
      return rng_get_brent_64bits (&(biomcmc_random_number->mt.x[12])); // 1 var 
    case 7:
      return rng_get_splitmix64 (&(biomcmc_random_number->mt.x[16])); // 1 var 
    case 8:
      return rng_get_xoroshiro256 (&(biomcmc_random_number->mt.x[20])); // 4 vars 
    default:
      return rng_get_std61 (&(biomcmc_random_number->mt.x[0])) ^ rng_get_gamerand64 (&(biomcmc_random_number->mt.x[1]));
  }
}

inline double
biomcmc_rng_get_52 (void)
{
  /* In Matsumoto's MT19937 code they use 53 bits (total double precision) but I think that the integer fraction of a
   * double is only 52 - so the integer-to-double conversion should use only the first 52 bits */
  return (double) (biomcmc_rng_get () >> 12);
}

inline uint32_t
biomcmc_rng_get_32 (void)
{ 
  /* Any subsequence has the same equidistribution properties of the most significant bits according to L'Ecuyer 
   * (for the Tausworthe algorithm) and Matsumoto (for the Mersenne twister); thus we can have 32 bit sequences by 
   * the independent sets of the 32 most and 32 least significant  bits. Each call to biomcmc_rng_get () thus 
   * generates two 32 bits pseudo-random numbers */
  if (biomcmc_random_number->have_bit32) {
    biomcmc_random_number->have_bit32 = false;
    return (uint32_t) ((biomcmc_random_number->bit32 >> 32) & 0xffffffffU);
  }
  biomcmc_random_number->bit32 = biomcmc_rng_get ();
  biomcmc_random_number->have_bit32 = true;
  return (uint32_t) (biomcmc_random_number->bit32 & 0xffffffffU);
}

/* * * * extra functions: gettimeofday() with maximum precision * * * */

void
biomcmc_get_time (int64_t *time)
{
#if _POSIX_TIMERS > 0
  struct timespec now;
  clock_gettime (CLOCK_REALTIME, &now);
  time[1] = now.tv_nsec; // always less than 1billion thus 32bits is enough
#else
  struct timeval now;
  gettimeofday (&now, NULL);
  time[1] = now.tv_usec; // always less than 1million thus 32bits is enough
#endif
  time[0] = now.tv_sec;
  return;
}

#ifdef _POSIX_TIMERS
#define TIMEWARP 1.e9
#else
#define TIMEWARP 1.e6
#endif

double
biomcmc_elapsed_time (int64_t *now, int64_t *past)
{
  if (now[1] < past[1]) return (((double)(past[1] - now[1]) / (double)(TIMEWARP)) - 1. + (double)(now[0] - past[0]));
  else                  return (((double)(now[1] - past[1]) / (double)(TIMEWARP))      + (double)(now[0] - past[0]));
//  return (double)(now[0] + now[1]/TIMEWARP) - (double)(past[0] + past[1]/TIMEWARP);
}

double
biomcmc_update_elapsed_time (int64_t *past)
{
  int64_t now[2];
  double seconds;
  biomcmc_get_time (now);
  if (now[1] < past[1]) seconds = (((double)(past[1] - now[1]) / (double)(TIMEWARP)) - 1. + (double)(now[0] - past[0]));
  else                  seconds = (((double)(now[1] - past[1]) / (double)(TIMEWARP))      + (double)(now[0] - past[0]));
  past[0] = now[0]; past[1] = now[1];
  return seconds;
}
