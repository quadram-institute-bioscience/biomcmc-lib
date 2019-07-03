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

#include "kmerhash.h"

/** OBS: uint32_t d[2] --> d = [A B C D] [E F G H] (1 byte per letter) then 
 * uint8_t *x = d      --> x = [D] [C] [B] [A] [H] [G] [F] [E] (endianness)  **/

/* similar to char2bit[] from alignment[], but has forward and reverse */
static uint8_t dna_in_4_bits[256][2] = {{0xff}}; /* DNA base to bitpattern translation, with 1st element set to arbitrary value */
static uint8_t dna_in_2_bits[256][2] = {{0xff}}; /* no ambigous chars/indels, represented by 4 (0100 in bits) */

static void initialize_dna_to_bit_tables (void);

kmerhash
new_kmerhash_from_dna_sequence (char *dna, size_t dna_length, bool dense)
{
  int i;
  if (dna_in_4_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); // run once per program
  kmerhash kmer = (kmerhash) biomcmc_malloc (sizeof (struct kmerhash_struct));
  kmer->n_f = 2;
  kmer->forward = (uint64_t*) biomcmc_malloc (kmer->n_f * sizeof (uint64_t*));
  kmer->reverse = (uint64_t*) biomcmc_malloc (kmer->n_f * sizeof (uint64_t*));
  for (i=0; i < kmer->n_f; i++) kmer->forward[i] = kmer->reverse[i] = 0UL;

  kmer->dense = dense;
  kmer->n_hash = 4;
  kmer->n_kmer = 3;
  kmer->hash = (uint64_t*) biomcmc_malloc (kmer->n_hash * sizeof (uint64_t*));
  kmer->kmer = (uint64_t*) biomcmc_malloc (kmer->n_kmer * sizeof (uint64_t*));
  for (i=0; i < kmer->n_hash; i++) kmer->hash[i] = 0UL;
  for (i=0; i < kmer->n_kmer; i++) kmer->kmer[i] = 0UL;

  kmer->dna = dna;
  kmer->n_dna = dna_length; 
  kmer->i = 0;
  return kmer;
}

void
del_kmerhash (kmerhash kmer)
{
  if (!kmer) return;
  if (kmer->forward) free (kmer->forward);
  if (kmer->reverse) free (kmer->reverse);
  if (kmer->hash) free (kmer->hash);
  if (kmer->kmer) free (kmer->kmer);
  free (kmer);
}

bool
kmerhash_iterator (kmerhash kmer)
{
  unsigned int n_ops = 1, dnachar;
  uint64_t hf, hr;
  if (kmer->i == kmer->n_dna) return false;

  // assume for/rev are full, solve initial case later
  if (kmer->dense) {
    while ((kmer->i < kmer->n_dna) && (dna_in_2_bits[(int)(kmer->dna[kmer->i])][0] > 3)) kmer->i++;
    if (kmer->i == kmer->n_dna) return false;
    dnachar = (int) kmer->dna[kmer->i];
    n_ops = 2;
    //ABCD->BCDE: forward [C D] [A B] --> [D E] [B C] , reverse: [b a] [d c] -->  [c b] [e d]
    kmer->forward[1] = kmer->forward[1] << 2 | kmer->forward[0] >> 62UL;  // forward[1] at left of forward[0]
    kmer->forward[0] = kmer->forward[0] << 2 | dna_in_2_bits[dnachar][0]; // forward[0] receives new value

    kmer->reverse[0] = kmer->reverse[0] >> 2 | kmer->reverse[1] << 62UL; // reverse[0] is at left
    kmer->reverse[1] = kmer->reverse[1] >> 2 | ((uint64_t)(dna_in_2_bits[dnachar][1]) << 62UL); // reverse[1] receives new value
  }
  else { // 4 bits (dense==false) can handle non-orthodox DNA characters
    dnachar = (int) kmer->dna[kmer->i];
    kmer->forward[1] = kmer->forward[1] << 4 | kmer->forward[0] >> 60UL; 
    kmer->forward[0] = kmer->forward[0] << 4 | dna_in_4_bits[dnachar][0];

    kmer->reverse[0] = kmer->reverse[0] >> 4 | kmer->reverse[1] << 60UL;
    kmer->reverse[1] = kmer->reverse[1] >> 4 | ((uint64_t)(dna_in_4_bits[dnachar][1]) << 60UL);
  } // else if (dense)
  kmer->i++;
  if (kmer->i < (4 * n_ops)) kmerhash_iterator (kmer); //do not update hashes, just fill forward[] and reverse[] 

  // 16 bits: 4-mer (or 8-mer if dense)
  hf = kmer->forward[0] & 0xffffUL; hr = kmer->reverse[1] >> 48;
  if (hr < hf) hf = hr;
  kmer->kmer[0] = hf;
  kmer->hash[0] = biomcmc_xxh64 (&hf, 2, 171); // 171 is just a seed

  if (kmer->i >= (8 * n_ops)) { // 32 bits: 8-mer (or 16-mer if dense)
    hf = kmer->forward[0] & 0xffffffffUL; hr = kmer->reverse[1] >> 32;
    if (hr < hf) hf = hr;
    kmer->kmer[1] = hf;
    kmer->hash[1] = biomcmc_xxh64 (&hf, 4, 317);
  }

  if (kmer->i >= (16 * n_ops)) { // 64 bits: 16-mer (or 32-mer if dense)
    kmer->kmer[2] = kmer->forward[0]; 
    if (kmer->kmer[2] > kmer->reverse[1]) kmer->kmer[2] = kmer->reverse[1];
    kmer->hash[2] = biomcmc_xxh64 (&(kmer->kmer[2]), 8, 635);
  }

  if (kmer->i >= (32 * n_ops)) { // 128 bits: 32-mer (or 64-mer if dense) --> no equiv k-mer
    if (kmer->forward[0] < kmer->reverse[1]) kmer->hash[3] = biomcmc_xxh64 (kmer->forward, 16, 9); // 9 must be the same seed 
    else if ((kmer->forward[0] == kmer->reverse[1]) && (kmer->forward[1] < kmer->reverse[0])) kmer->hash[3] = biomcmc_xxh64 (kmer->forward, 16, 9);
    else kmer->hash[3] = biomcmc_xxh64 (kmer->reverse, 16, 9);
  } 
  return true;
}

static void
initialize_dna_to_bit_tables (void)
{
  int i;
  for (i = 0; i < 256; i++) dna_in_4_bits[i][0] = dna_in_4_bits[i][1] = 0;
  /* The ACGT is PAUP convention (and maybe DNAml, fastDNAml); PAML uses TCAG ordering */
  dna_in_4_bits['A'][0] = 1;   dna_in_4_bits['A'][1] = 8;  /* .   A */ /* 0001 */ /* reverse is 'T'    = 8  */
  dna_in_4_bits['B'][0] = 14;  dna_in_4_bits['B'][1] = 7;  /* .TGC  */ /* 1110 */ /* reverse is 'ACG'  = 7  */
  dna_in_4_bits['C'][0] = 2;   dna_in_4_bits['C'][1] = 4;  /* .  C  */ /* 0010 */ /* reverse is 'G'    = 4  */
  dna_in_4_bits['D'][0] = 13;  dna_in_4_bits['D'][1] = 11; /* .TG A */ /* 1101 */ /* reverse is 'TCA'  = 11 */
  dna_in_4_bits['G'][0] = 4;   dna_in_4_bits['G'][1] = 2;  /* . G   */ /* 0100 */ /* reverse is 'C'    = 2  */
  dna_in_4_bits['H'][0] = 11;  dna_in_4_bits['H'][1] = 13; /* .T CA */ /* 1011 */ /* reverse is 'TGA'  = 13 */
  dna_in_4_bits['K'][0] = 12;  dna_in_4_bits['K'][1] = 3;  /* .TG   */ /* 1100 */ /* reverse is 'AC'   = 3  */
  dna_in_4_bits['M'][0] = 3;   dna_in_4_bits['M'][1] = 12; /* .  CA */ /* 0011 */ /* reverse is 'TG'   = 12 */
  dna_in_4_bits['N'][0] = 15;  dna_in_4_bits['N'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna_in_4_bits['O'][0] = 15;  dna_in_4_bits['O'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna_in_4_bits['R'][0] = 5;   dna_in_4_bits['R'][1] = 10; /* . G A */ /* 0101 */ /* reverse is 'TC'   = 10 */
  dna_in_4_bits['S'][0] = 6;   dna_in_4_bits['S'][1] = 6;  /* . GC  */ /* 0110 */ /* reverse is 'GC'   = 6  */
  dna_in_4_bits['T'][0] = 8;   dna_in_4_bits['T'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna_in_4_bits['U'][0] = 8;   dna_in_4_bits['U'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna_in_4_bits['V'][0] = 7;   dna_in_4_bits['V'][1] = 14; /* . GCA */ /* 0111 */ /* reverse is 'TGC'  = 14 */
  dna_in_4_bits['W'][0] = 9;   dna_in_4_bits['W'][1] = 9;  /* .T  A */ /* 1001 */ /* reverse is 'TA'   = 9  */
  dna_in_4_bits['X'][0] = 15;  dna_in_4_bits['X'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna_in_4_bits['Y'][0] = 10;  dna_in_4_bits['Y'][1] = 5;  /* .T C  */ /* 1010 */ /* reverse is 'GA'   =  5 */
  dna_in_4_bits['?'][0] = 15;  dna_in_4_bits['?'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna_in_4_bits['-'][0] = 0;   dna_in_4_bits['-'][1] = 0;  /* .TGCA */ /* fifth state */

  for (i = 0; i < 256; i++) dna_in_2_bits[i][0] = dna_in_2_bits[i][1] = 4; // calling function must check if < 4
  dna_in_2_bits['A'][0] = 0;   dna_in_2_bits['A'][1] = 3;  /*  A  <-> T  */
  dna_in_2_bits['C'][0] = 1;   dna_in_2_bits['C'][1] = 2;  /*  C  <-> G  */
  dna_in_2_bits['G'][0] = 2;   dna_in_2_bits['G'][1] = 1;  /*  G  <-> C  */
  dna_in_2_bits['T'][0] = 3;   dna_in_2_bits['T'][1] = 0;  /*  T  <-> A  */
  dna_in_2_bits['U'][0] = 3;   dna_in_2_bits['U'][1] = 0;  /*  U  <-> A  */
}

//for (j = 0; j < 16; j++) printf ("%d ", (int)((hash_f >> 4*j) & 15LL)); for (j = 0; j < 16; j++) printf ("%c", dna[i-j]);
