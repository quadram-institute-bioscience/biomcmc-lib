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

#include "minhash.h"

/* almost identical to dna4bits[] from alignment[], but has forward and reverse */
uint8_t  dna4bits[256][2] = {{0xff}}; /* DNA base to bitpattern translation, with 1st element set to arbitrary value */

void initialize_dna_to_bit_table (void);
void fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr);

cm_sketch
new_cm_sketch (int max_vector_size)
{ // this may not be a locally-sensitive hashing (LSH) since similar inputs go to distinct buckets
  int i;
  cm_sketch cm = (cm_sketch) biomcmc_malloc (sizeof (struct cm_sketch_struct));
  cm->size = max_vector_size;
  cm->mod = (0xffffffffffffffffUL / (uint64_t) (max_vector_size + 1)); // plusone for case hash == MAX 
  cm->count = 0; // FIXME: modulus is 32bits now
  cm->freq = (int*) biomcmc_malloc (cm->size * sizeof (int));
  for (i = 0; i < cm->size; i++) cm->freq[i] = 0;

  return cm;
}

void
del_cm_sketch (cm_sketch cm)
{
  if (!cm) return;
  if (cm->freq) free (cm->freq);
  free (cm);
}

void
update_cm_sketch_from_fixedhash (cm_sketch cm, uint64_t hash_f, uint64_t hash_r)
{
  uint32_t h32[4];
  int i; // FIXME: must hash 64bits before tovecotr() otherwise it's just 8-mer
  if (hash_f < hash_r) biomcmc_hashint64_to_vector (hash_f, h32); // one of them is redundant
  else                 biomcmc_hashint64_to_vector (hash_r, h32);
  for (i=0; i < 4; i++) cm->freq[i][ (int) h32[i]/cm->mod ]++;  
}

void 
fixedhash_sketch_from_dna (char *dna)
{
  int i;
  int length = strlen (dna);
  uint64_t hash_f = 0UL, hash_r = 0UL;
  if (dna4bits[0][0] == 0xff) initialize_dna_to_bit_table ();
  for (i = 0; i < 15; i++) fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r); 
  for (i = 15; i < length; i++) {
    fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r);
    // TODO: here we update the cm_sketch
    //for (j = 0; j < 16; j++) printf ("%d ", (int)((hash_f >> 4*j) & 15LL)); for (j = 0; j < 16; j++) printf ("%c", dna[i-j]);
  }
}

void
fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr)
{
  *hf = *hf << 4 | dna4bits[dnachar][0]; // forward
  *hr = *hf >> 4 | ((uint64_t)(dna4bits[dnachar][1]) << 60UL); // reverse
}


void
initialize_dna_to_bit_table (void)
{
  int i;
  for (i = 0; i < 256; i++) dna4bits[i][0] = dna4bits[i][1] = 0;
  /* The ACGT is PAUP convention (and maybe DNAml, fastDNAml); PAML uses TCAG ordering */
  dna4bits['A'][0] = 1;   dna4bits['A'][1] = 8;  /* .   A */ /* 0001 */ /* reverse is 'T'    = 8  */
  dna4bits['B'][0] = 14;  dna4bits['B'][1] = 7;  /* .TGC  */ /* 1110 */ /* reverse is 'ACG'  = 7  */
  dna4bits['C'][0] = 2;   dna4bits['C'][1] = 4;  /* .  C  */ /* 0010 */ /* reverse is 'G'    = 4  */
  dna4bits['D'][0] = 13;  dna4bits['D'][1] = 11; /* .TG A */ /* 1101 */ /* reverse is 'TCA'  = 11 */
  dna4bits['G'][0] = 4;   dna4bits['G'][1] = 2;  /* . G   */ /* 0100 */ /* reverse is 'C'    = 2  */
  dna4bits['H'][0] = 11;  dna4bits['H'][1] = 13; /* .T CA */ /* 1011 */ /* reverse is 'TGA'  = 13 */
  dna4bits['K'][0] = 12;  dna4bits['K'][1] = 3;  /* .TG   */ /* 1100 */ /* reverse is 'AC'   = 3  */
  dna4bits['M'][0] = 3;   dna4bits['M'][1] = 12; /* .  CA */ /* 0011 */ /* reverse is 'TG'   = 12 */
  dna4bits['N'][0] = 15;  dna4bits['N'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['O'][0] = 15;  dna4bits['O'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['R'][0] = 5;   dna4bits['R'][1] = 10; /* . G A */ /* 0101 */ /* reverse is 'TC'   = 10 */
  dna4bits['S'][0] = 6;   dna4bits['S'][1] = 6;  /* . GC  */ /* 0110 */ /* reverse is 'GC'   = 6  */
  dna4bits['T'][0] = 8;   dna4bits['T'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna4bits['U'][0] = 8;   dna4bits['U'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna4bits['V'][0] = 7;   dna4bits['V'][1] = 14; /* . GCA */ /* 0111 */ /* reverse is 'TGC'  = 14 */
  dna4bits['W'][0] = 9;   dna4bits['W'][1] = 9;  /* .T  A */ /* 1001 */ /* reverse is 'TA'   = 9  */
  dna4bits['X'][0] = 15;  dna4bits['X'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['Y'][0] = 10;  dna4bits['Y'][1] = 5;  /* .T C  */ /* 1010 */ /* reverse is 'GA'   =  5 */
  dna4bits['?'][0] = 15;  dna4bits['?'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['-'][0] = 0;   dna4bits['-'][1] = 0;  /* .TGCA */ /* fifth state */
}

