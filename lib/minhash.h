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

/*! \file minhash.h 
 *  \brief Min-count and MinHash sketches for DNA sequences
 *
 */

#ifndef _biomcmc_minhash_h_
#define _biomcmc_minhash_h_

#include "alignment.h"

typedef struct cm_sketch_struct* cm_sketch;

struct cm_sketch_struct
{
  int size, count;
  uint64_t mod;
  int *freq;
}

/*! \brief min-count sketch of fixedhash (numeric representation of 16mers) */
void fixedhash_sketch_from_dna (char *dna);


#endif
