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
typedef enum {kmer_class_fast, kmer_class_genome, kmer_class_full} kmer_class;
const char *kmer_class_string[] = {"fast", "genome", "full"};
// maybe string must be declared extern here and defined in .c 

struct kmer_params_struct
{
  uint64_t mask1[6] = {0xffffUL, 0xffffffUL, 0xffffffffUL, 0xffffffffffUL, 0xffffffffffffUL, 0xffffffffffffffffUL}; // 4, 6, 8, 10, 12, 16 4bits
  uint8_t shift1[6] = {      48,         40,           32,             24,               16,                    0};
  uint64_t mask2[4] = {0xffffUL, 0xffffffffUL, 0xffffffffffffUL, 0xffffffffffffffffUL};
  uint8_t  size[10] = {4, 6, 8, 10, 12, 16, 20, 24, 28, 32};
  bool dense=false;
  uint8_t n1 = 6, n2 = 4; // if n2 == 0 then do not use second element (128 bits)
  int ref_counter;
};

struct kmerhash_struct 
{
  kmer_params p;
  uint64_t forward[2], reverse[2]; 
  uint64_t *hash, *kmer;   /*! \brief hash = 4mer, 8mer, 16mer, and 32mer hashed ; kmer = 4mer, 8mer, 16mer original */
  int n_f, n_hash, n_kmer; /*! \brief n_f = 2 (128bits), n_hash = 4, n_kmer = 3 */
  char *dna;
  size_t i, n_dna;
  bool dense; /*! \brief 4bits per base (false) or 2bits (true) */
  int ref_counter;
};



kmerhash new_kmerhash (bool dense);
void link_kmerhash_to_dna_sequence (kmerhash kmer, char *dna, size_t dna_length);
void del_kmerhash (kmerhash kmer);
bool kmerhash_iterator (kmerhash kmer);

#endif
