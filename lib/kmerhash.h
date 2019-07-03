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

typedef struct kmerhash_struct* kmerhash;

struct kmerhash_struct 
{
  uint64_t *forward, *reverse; 
  uint64_t *hash, *kmer;   /*! \brief hash = 4mer, 8mer, 16mer, and 32mer hashed ; kmer = 4mer, 8mer, 16mer original */
  int n_f, n_hash, n_kmer; /*! \brief n_f = 2 (128bits), n_hash = 4, n_kmer = 3 */
  char *dna;
  size_t i, n_dna;
  bool dense; /*! \brief 4bits per base (false) or 2bits (true) */
};

kmerhash new_kmerhash_from_dna_sequence (char *dna, size_t dna_length, bool dense);
void del_kmerhash (kmerhash kmer);
bool kmerhash_iterator (kmerhash kmer);

#endif
