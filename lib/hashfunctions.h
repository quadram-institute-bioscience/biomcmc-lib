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

/*! \file hashfunctions.h
 *  \brief Collections of hash functions for 32 and 64 bits, including one-liners, murmurhash, and xxhash 
 */

#ifndef _biomcmc_hashfunctions_h_ 
#define _biomcmc_hashfunctions_h_ 

#include "bipartition.h" // includes lowlevel.h

uint32_t biomcmc_hashint_salted (uint32_t a, unsigned int salt);
uint32_t biomcmc_hashbyte_salted (const void *str, size_t size, unsigned int salt);
uint64_t biomcmc_hashint64_salted (uint64_t k, unsigned int salt);
uint32_t biomcmc_hashint_mix_salted (uint32_t a, uint32_t b, unsigned int salt);
uint64_t biomcmc_hashint64_mix_salted (uint64_t a, uint64_t b, unsigned int salt);

uint32_t biomcmc_hashint_64to32_seed (uint64_t x, unsigned int seed);
uint32_t biomcmc_hashint_64to32 (uint64_t key);
/*! \brief 32bits hash value for bipartition */
uint32_t bipartition_hash (bipartition bip);

/*! \brief murmurhash3 using 64bits to return 128 bits (4 ints) of hash into out[] and also 64 bits as return value 
 * The 64bits is the format used internally (for speed), but the input can be a vector of any size (>1 byte) */
uint64_t biomcmc_murmurhash3_128bits ( const void *key, const size_t len, const uint32_t seed, void *out);
/*! \brief convenience function for calling mumurhash_128bits without an output vector */
uint64_t biomcmc_murmurhash3_64bits ( const void *key, const size_t len, const uint32_t seed);
/*! \brief murmurhash3 using 32bits to return 32 bits of hash as return value */
uint32_t biomcmc_murmurhash3_32bits (const void *data, size_t nbytes, const uint32_t seed);
/*! \brief xxhash function for 64 bits */
uint64_t biomcmc_xxh64 (const void *input, const size_t len, const uint32_t seed);

/**   Google HighwayHash for C https://github.com/google/highwayhash/ **/

/** \brief 256 bits for highwayhash */
typedef struct {
  uint64_t v0[4];
  uint64_t v1[4];
  uint64_t mul0[4];
  uint64_t mul1[4];
} HighwayHashState;

/* Initializes state with given key */
void HighwayHashReset (const uint64_t key[4], HighwayHashState* state);
/* Takes a packet of 32 bytes */
void HighwayHashUpdatePacket (const uint8_t* packet, HighwayHashState* state);
/* Adds the final 1..31 bytes, do not use if 0 remain */
void HighwayHashUpdateRemainder (const uint8_t* bytes, const size_t size_mod32, HighwayHashState* state);
/* Compute final hash value. Makes state invalid. */
uint64_t HighwayHashFinalize64 (HighwayHashState* state);
void HighwayHashFinalize128 (HighwayHashState* state, uint64_t hash[2]);
void HighwayHashFinalize256 (HighwayHashState* state, uint64_t hash[4]);

/* Non-cat API: single call on full data                                      */
uint64_t HighwayHash64 (const uint8_t* data, size_t size, const uint64_t key[4]);
void HighwayHash128 (const uint8_t* data, size_t size, const uint64_t key[4], uint64_t hash[2]);
void HighwayHash256 (const uint8_t* data, size_t size, const uint64_t key[4], uint64_t hash[4]);

/* Cat API: allows appending with multiple calls                              */

typedef struct {
  HighwayHashState state;
  uint8_t packet[32];
  int num;
} HighwayHashCat;

/* Allocates new state for a new streaming hash computation */
void HighwayHashCatStart (const uint64_t key[4], HighwayHashCat* state);
void HighwayHashCatAppend (const uint8_t* bytes, size_t num, HighwayHashCat* state);
/* Computes final hash value */
uint64_t HighwayHashCatFinish64 (const HighwayHashCat* state);
void HighwayHashCatFinish128 (const HighwayHashCat* state, uint64_t hash[2]);
void HighwayHashCatFinish256 (const HighwayHashCat* state, uint64_t hash[4]);

#endif
