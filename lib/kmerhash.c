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

static uint64_t _tbl_mask[] = {0xffffUL, 0xffffffUL, 0xffffffffUL, 0xffffffffffUL, 0xffffffffffffUL, 0xffffffffffffUL, 0xffffffffffffffffUL}; 
static uint8_t _tbl_shift[] = {      48,         40,           32,             24,               16,                8,                    0};
static uint8_t _tbl_nbyte[] = {       2,          3,            4,              5,                6,                7,                    8};
static uint32_t _tbl_seed[] = {0x9040a6, 0x10bea992,   0x50edd67d,     0xb05a4f09,       0xf07046c5,       0x9c9445ab,           0xb2500f29};

/* i1[] and i2[] (i.e. elements above to be used on first and second 64bits, respect. */
static uint8_t _idx_mode[][2][7] = {
   {{0,2,4,6,0,0,0}, {2,6,0,0,0,0,0}}, 
   {{0,1,2,4,6,0,0}, {0,2,6,0,0,0,0}},
   {{0,1,2,3,4,5,6}, {0,0,0,0,0,0,0}},
   {{0,1,2,3,4,5,6}, {0,1,2,6,0,0,0}}
};   


static void initialize_dna_to_bit_tables (void);

kmer_params
new_kmer_params (enum_kmer_class mode)
{
  int j, *i1, *i2, bases_per_byte = 2; // bases_per_byte is 4 if dense
  kmer_params p = (kmer_params) biomcmc_malloc (sizeof (struct kmer_params_struct));
  p->ref_counter = 1;
  p->hashfunction = &biomcmc_xxh64;

  if (dna_in_4_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); // run once per program

  p->kmer_class_mode = mode;
  switch (mode) {
    case kmer_class_fast:
      int local_i1[] = {0,2,4,6}, local_i2[] = {2,6}; // indexes are based on _tbl_ above, static allocation
      i1 = &local_i1; i2 = &local_i2;
      p->n1 = 4; p->n2 = 2; 
      p->dense = false; bases_per_byte = 2; 
      break;

    case kmer_class_genome:
      int local_i1[] = {0,1,2,4,6}, local_i2[] = {0,2,6};
      i1 = &local_i1; i2 = &local_i2;
      p->n1 = 5; p->n2 = 3; 
      p->dense = true; bases_per_byte = 4;
      break;

    case kmer_class_short:
      int local_i1[] = {0,1,2,3,4,5,6}, local_i2[] = {0}; 
      i1 = &local_i1; i2 = &local_i2;
      p->n1 = 7; p->n2 = 0; // zero makes sure it's not called
      p->dense = false; bases_per_byte = 2;
      break;

    default: 
      int local_i1[] = {0,1,2,3,4,5,6}, local_i2[] = {0,1,2,6};
      i1 = &local_i1; i2 = &local_i2;
      p->n1 = 7; p->n2 = 4; 
      p->dense = false; bases_per_byte = 2; 
      break;
  };

  for (j=0; j < p->n1; j++) {
    p->mask1[j]  =  _tbl_mask[ i1[j] ];
    p->shift1[j] = _tbl_shift[ i1[j] ];
    p->seed[j]   =  _tbl_seed[ i1[j] ];
    p->nbytes[j] = _tbl_nbyte[ i1[j] ];
    p->size[j]  =  _tbl_nbyte[ i1[j] ] * bases_per_byte;
  }
  for (j=0; j < p->n2; j++) {
    p->mask2[j]  =  _tbl_mask[ i2[j] ];
    p->shift2[j] = _tbl_shift[ i2[j] ];
    p->seed[j]  = (_tbl_seed[ i2[j] ] >> 2) + 0x420a1d; // very noise, so random
    p->nbytes[j+p->n1] = _tbl_nbyte[ i2[j] ] + 8;
    p->size[j+p->n1] =  (_tbl_nbyte[ i2[j] ] + 8) * bases_per_byte;
  }
  return p;
}

void
del_kmer_params (kmer_params p)
{
  if (!p) return;
  if (--p->ref_counter) return;
  free (p);
} 

kmerhash
new_kmerhash (enum_kmer_class mode)
{
  int i;
  kmerhash kmer = (kmerhash) biomcmc_malloc (sizeof (struct kmerhash_struct));
  kmer->n_f = 2;
  kmer->forward = (uint64_t*) biomcmc_malloc (kmer->n_f * sizeof (uint64_t*));
  kmer->reverse = (uint64_t*) biomcmc_malloc (kmer->n_f * sizeof (uint64_t*));
  for (i=0; i < kmer->n_f; i++) kmer->forward[i] = kmer->reverse[i] = 0UL;

  kmer->p = new_kmer_params (mode);

  kmer->n_hash = kmer->p->n1 + kmer->p->n2;
  // p->size :: kmer->nsites_kmer = (uint8_t*) biomcmc_malloc (kmer->n_hash * sizeof (uint8_t));
  // if (kmer->dense) for (i = 0; i < kmer->n_hash; i++) kmer->nsites_kmer[i] = 2 * _kmer_size[i]; 
  // else             for (i = 0; i < kmer->n_hash; i++) kmer->nsites_kmer[i] = _kmer_size[i]; 

  kmer->hash = (uint64_t*) biomcmc_malloc (kmer->n_hash * sizeof (uint64_t*));
  kmer->kmer = (uint64_t*) biomcmc_malloc (kmer->n1 * sizeof (uint64_t*));
  for (i=0; i < kmer->n_hash; i++) kmer->hash[i] = 0UL;
  for (i=0; i < kmer->n1; i++) kmer->kmer[i] = 0UL;

  kmer->dna = NULL;
  kmer->n_dna = 0; 
  kmer->i = 0;
  kmer->ref_counter = 1;
  return kmer;
}

void
link_kmerhash_to_dna_sequence (kmerhash kmer, char *dna, size_t dna_length)
{
  int i;
  kmer->dna = dna;
  kmer->n_dna = dna_length; 
  kmer->i = 0;
  kmer->kmer = 0UL;
  for (i=0; i < kmer->n_f; i++) kmer->forward[i] = kmer->reverse[i] = 0UL;
  for (i=0; i < kmer->n_hash; i++) kmer->hash[i] = 0UL;
  for (i=0; i < kmer->n1; i++) kmer->kmer[i] = 0UL;
}

void
del_kmerhash (kmerhash kmer)
{
  if (!kmer) return;
  if (--kmer->ref_counter) return;
  if (kmer->forward) free (kmer->forward);
  if (kmer->reverse) free (kmer->reverse);
  if (kmer->hash) free (kmer->hash);
  if (kmer->kmer) free (kmer->kmer);
  del_kmer_params (kmer->p)
  free (kmer);
}

bool
kmerhash_iterator (kmerhash kmer)
{
  unsigned int i, j, dnachar;
  uint64_t hf, hr;
  uint8_t rev8b;

  if (kmer->i == kmer->n_dna) return false;

  // assume for/rev are full, solve initial case later
  if (kmer->p->dense) {
    while ((kmer->i < kmer->n_dna) && (dna_in_2_bits[(int)(kmer->dna[kmer->i])][0] > 3)) kmer->i++;
    if (kmer->i == kmer->n_dna) return false;
    dnachar = (int) kmer->dna[kmer->i];
    //ABCD->BCDE: forward [C D] [A B] --> [D E] [B C] , reverse: [b a] [d c] -->  [c b] [e d]
    if (kmer->p->n2) { // only if we need 2 uint64_t elements forward[1] and reverse[0] are used
      kmer->forward[1] = kmer->forward[1] << 2 | kmer->forward[0] >> 62UL;  // forward[1] at left of forward[0]
      kmer->reverse[0] = kmer->reverse[0] >> 2 | kmer->reverse[1] << 62UL;  // reverse[0] is at left
    }
    kmer->forward[0] = kmer->forward[0] << 2 | dna_in_2_bits[dnachar][0]; // forward[0] receives new value
    kmer->reverse[1] = kmer->reverse[1] >> 2 | ((uint64_t)(dna_in_2_bits[dnachar][1]) << 62UL); // reverse[1] receives new value
  }
  else { // 4 bits (dense==false) can handle non-orthodox DNA characters
    dnachar = (int) kmer->dna[kmer->i];
    if (kmer->p->n2) { 
      kmer->forward[1] = kmer->forward[1] << 4 | kmer->forward[0] >> 60UL; 
      kmer->reverse[0] = kmer->reverse[0] >> 4 | kmer->reverse[1] << 60UL;
    }
    kmer->forward[0] = kmer->forward[0] << 4 | dna_in_4_bits[dnachar][0];
    kmer->reverse[1] = kmer->reverse[1] >> 4 | ((uint64_t)(dna_in_4_bits[dnachar][1]) << 60UL);
  } // else if (dense)
  kmer->i++;
  if (kmer->i < kmer->p->size[0]) kmerhash_iterator (kmer); //do not update hashes, just fill forward[] and reverse[] 

  for (i=0; i < kmer->p->n1; i++) { // fit into one uint64_t 
    if (kmer->i >= kmer->p->size[i]) { 
      hf = kmer->forward[0] & kmer->p->mask1[i]; hr = kmer->reverse[1] >> kmer->p->shift1[i];
      if (hr < hf) hf = hr;
      kmer->kmer[i] = hf;
      kmer->hash[i] = kmer->p->hashfunction (&hf, kmer->p->nbytes[i], kmer->p->seed[i]); // 313 is just a seed
    }
  }

  for (i=0; i < kmer->p->n2; i++) { // need to compare two uint64_t; no kmer->kmer[] possible 
    j = kmer->p->n1 + i;
    if (kmer->i >= kmer->p->size[j]) {
      //ABCDE : forward[0][1]=[DE][BC] reverse[0][1]=[cb][ed] ; 1. compare [DE] and [ed]
      if (kmer->forward[0] < kmer->reverse[1])
        kmer->hash[j] = kmer->p->hashfunction (kmer->forward, kmer->p->nbytes[j], kmer->p->seed[j]);
      else if ((kmer->forward[0] == kmer->reverse[1]) && // 2. compare [_C] and [c_] (where '_' means masked out)
               ((kmer->forward[1] & kmer->p->mask2[i]) < (kmer->reverse[0] >> kmer->p->shift2[i])) )
        kmer->hash[j] = kmer->p->hashfunction (kmer->forward, kmer->p->nbytes[j], kmer->p->seed[j]);
      else { // 3. skip b from reverse[0]
        rev8b = (uint8_t) kmer->reverse; // so that rev8b[3] is 3 bytes after start of reverse
        kmer->hash[j] = kmer->p->hashfunction (rev8b + kmer->shift2[i], kmer->p->nbytes[j], kmer->p->seed[j]);
      }
    } // if dna[i] > min size for this operation
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
