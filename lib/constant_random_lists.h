/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be usefulull, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICullAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file constant_random_lists.h 
 *  \brief lists of prime numbers and hand-picked random numbers used by hashes and PRNGs. Examples include the rolling
 *  hash (DNA bases mapped to a random value) and the (deterministic) spice used to generate several streams.
 */

#ifndef _biomcmc_constant_random_lists_h_
#define _biomcmc_constant_random_lists_h_
#include "lowlevel.h"

/*! \brief  Five-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
extern uint64_t sTable76[44][5]; 
/*! \brief  Four-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
extern uint64_t sTable543[106][4];

/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 7 and 6 (five elements)) */
extern uint64_t qTable76[2][5];
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (tables 7 and 6 (five elements)) */
extern uint64_t kTable76[2][5];
/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3 */
extern uint64_t qTable543[4][4];
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3*/
extern uint64_t kTable543[4][4];
extern uint64_t Cmask[28];
/* http://www.math.uni-bielefeld.de/~sillke/ALGORITHMS/random : 16-bit constants k for which both k*2^16-1 and k*2^15-1 are prime */
extern uint16_t marsaglia_constants_size;
extern uint32_t marsaglia_constants[];
extern uint16_t rnd_salt_h16_list_size; /*! \brief list of fixed 16bit random numbers */ 
extern uint16_t rnd_salt_h16_list[]; 
extern uint16_t prime_salt_list_size;  /*! \brief list of fixed 16 and 32 bit prime random numbers */ 
extern uint32_t prime_salt_list[];
extern uint16_t rnd_salt_h64_list_size;
extern uint64_t rnd_salt_h64_list[];
extern uint16_t ulx_h64_size; /** \brief  these are used by some hash functions, please do not change their order */
extern uint64_t ulx_h64[];
/* prob_distribution_aux.c */
extern uint16_t lgamma_algmcs_size;
extern double   lgamma_algmcs[];
extern uint16_t lgamma_coeffs_size;
extern double   lgamma_coeffs[];
extern uint16_t stirl_sferr_halves_size;
extern double   stirl_sferr_halves[];


#endif
