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

/*! \file kmerhash.h 
 *  \brief k-mer handling of DNA sequences, with hash transformation  
 *
 */

#ifndef _biomcmc_kmerhash_h_
#define _biomcmc_kmerhash_h_

#include "alignment.h"

typedef struct kmer_params_struct* kmer_params;
typedef struct kmerhash_struct* kmerhash;

extern const char *biomcmc_kmer_class_string[];
// maybe string must be declared extern here and defined in .c 

struct kmer_params_struct
{
  uint64_t mask1[7], mask2[7]; 
  uint8_t n1, n2, shift1[7], shift2[7], size[14], nbytes[14]; // size = how many bases are stored (if dense, x2); nbytes = how many bytes (uint8_t) fit 
  uint32_t seed[14]; 
  uint64_t (*hashfunction) (const void *, const size_t, const uint32_t);
  int dense; /*! \brief 4bits per base or 2bits or 1 bit (GC content) */
  int kmer_class_mode;
  int ref_counter;
};

struct kmerhash_struct 
{
  kmer_params p;
  uint64_t *forward, *reverse; 
  uint64_t *hash, *kmer;  /*! \brief hash = 4mer, 8mer, etc. hashed ; kmer = original bitstring OR its complement, masked */
  int n_hash, n_f; /*! \brief n_f = 2 (128bits) */
  char *dna;
  size_t i, n_dna;
  int ref_counter;
};

kmer_params new_kmer_params (int mode);
void del_kmer_params (kmer_params p);
kmerhash new_kmerhash (int mode);
void link_kmerhash_to_dna_sequence (kmerhash kmer, char *dna, size_t dna_length);
void del_kmerhash (kmerhash kmer);
bool kmerhash_iterator (kmerhash kmer);

#endif
