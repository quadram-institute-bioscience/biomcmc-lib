/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 */

/*! \file
 *  \brief Rolling hash representation of sequence k-mers
 *
 *  Inspired by the [linclust algorithm](https://github.com/soedinglab/MMseqs2), here we calculate both the forward and
 *  reverse-complement (not present in minclust) rolling kmers, storing the canonical as the minimum between both.
 *  The homopolymer compression does not take place here.
 *  
 */

#ifndef _biomcmc_rolling_hash_h_
#define _biomcmc_rolling_hash_h_

#include "alignment.h"

typedef struct rolling_hash_struct* rolling_hash;

/** \brief struct:: for one kmer, updated as sequence is scanned */ 
struct rolling_hash_struct 
{
  uint32_t kmer[2]; /**< forward and reverse rolling hashes */
  uint32_t **salted_hashcode; /**< translation between DNA bases and random numbers */
  char *dna;         /**< pointer to DNA sequence */
  uint8_t kmer_size, /**< kmer-size must be between 2 and 31 */
          rol,       /**< bits to roll */
          canonical; /**< boolean, indicating canonical direction: forward (zero) or reverse (one) */
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
