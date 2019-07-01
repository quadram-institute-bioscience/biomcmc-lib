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

#include "hashtable.h"

/* Aux functions for hashtable of strings */
uint32_t hash(hashtable ht, const char* str);
inline uint32_t hash1 (hashtable ht);
inline uint32_t hash2 (hashtable ht);
/* Aux functions for hashtable of bipartitions */
inline uint32_t biphash1 (bip_hashtable ht);
inline uint32_t biphash2 (bip_hashtable ht);


hashtable 
new_hashtable (int size)
{
  int i; 
  double x;
  hashtable ht;
  
  ht = (hashtable) biomcmc_malloc (sizeof (struct hashtable_struct));

  /* Setting hashtable to be size of power of 2. This is required so that proper open addressing can occur, i.e. 
   * eventually every entry in the table is considered.
   */
  x = ceil (log (size) / log (2));
  size = pow (2, x+1);

  ht->size = size;
  
  /* allocate memory for table */
  ht->table = (hashtable_item*) biomcmc_malloc(ht->size * sizeof (hashtable_item)); 

  /* initialize to NULLS */
  for(i = 0; i < ht->size; i++) ht->table[i] = 0; 
  ht->P = 2147483647; /* initialize P (large prime)*/
  ht->probelength = 0;
  
  /*initialize hash1 and hash2 variables*/
  srand (time(0)); /* the GSL library would be overkill... */
  ht->a1 = rand() * (ht->P - 1) + 1;
  ht->a2 = rand() * (ht->P - 1) + 1;
  ht->b1 = rand() * ht->P;
  ht->b2 = rand() * ht->P;

  ht->ref_counter = 1; /* at least one place (the calling structure/function) is using this hashtable */

  return ht;
}

void 
del_hashtable (hashtable ht)
{
  if (ht) {
    int i;
    if (--ht->ref_counter) return; /* some other place is using this hashtable, we cannot delete it yet */
    for(i=ht->size-1; i>=0; i--) if (ht->table[i]) { 
      if(ht->table[i]->key) free (ht->table[i]->key); 
      free (ht->table[i]); 
    }
    if (ht->table) free (ht->table);
    free (ht);
  }
}

/* This is called before hash1() and hash2(). It sets the  value of h to be used by hash1() and hash2() and saves 
 * excessive calls to the hash function. */
uint32_t 
hash (hashtable ht, const char* key) 
{
  uint32_t g;
  
  ht->h = 0;
  while (*key) {
    ht->h = (ht->h << 4) + *key++;
    g = ht->h & 0xF0000000L;
    if (g) ht->h ^= g >> 24;
    ht->h &= ~g;
  }
  return (ht->h % ht->P); /*P_ is a large prime*/
}

uint32_t 
hash1 (hashtable ht) 
{
  return ((((ht->a1 * ht->h) + ht->b1) % ht->P) % ht->size) ;
}

uint32_t 
hash2 (hashtable ht) 
{
  return ((((ht->a2 * ht->h + ht->b2) % ht->P) % (ht->size - 3)) | 1);
}

void 
insert_hashtable (hashtable ht, char* key, int value) 
{
  uint32_t i;
  int h1, h2;
  
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);

  ht->probelength = 0;
  
  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!strcmp (ht->table[i]->key, key)) /* key was already inserted */
      return; 
  }
  /* alloc space for new key */
  ht->table[i] = biomcmc_malloc (sizeof (struct hashtable_item_struct));
  ht->table[i]->key = (char*) biomcmc_malloc((strlen (key)+1) * sizeof (char));
  strcpy(ht->table[i]->key, key);
  ht->table[i]->value = value;
  return;
}

int 
lookup_hashtable (hashtable ht, char* key) 
{
  uint32_t i;
  int h1, h2;
  
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);
  
  ht->probelength = 0;

  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!(ht->table[i])) return -1;
    else if ( (ht->table[i]) && (!strcmp (ht->table[i]->key, key)) ) 
      return ht->table[i]->value;
  }
  return -2;
}

/* * Functions to work with hashtable of bipartitions * */

bip_hashtable 
new_bip_hashtable (int size)
{
  int i; 
  double x;
  bip_hashtable ht;
  
  ht = (bip_hashtable) biomcmc_malloc (sizeof (struct bip_hashtable_struct));

  /* Setting hashtable to be size of power of 2. This is required so that proper open addressing can occur, i.e. 
   * eventually every entry in the table is considered.  */
  x = ceil (log (size) / log (2));
  size = pow (2, x+1);
  ht->size = size;
  
  /* allocate memory for table */
  ht->table = (bip_hashitem*) biomcmc_malloc(ht->size * sizeof (bip_hashitem)); 
  /* initialize to NULLS */
  for(i = 0; i < ht->size; i++) ht->table[i] = NULL; 
  ht->P = 2147483647; /* initialize P (large prime)*/
  ht->probelength = 0;
  ht->maxfreq = 1; // will store count of most frequent bipartition, to normalize frequency 
  
  /*initialize hash1 and hash2 variables*/
  srand (time(0)); /* the GSL library would be overkill... */
  ht->a1 = rand() * (ht->P - 1) + 1;
  ht->a2 = rand() * (ht->P - 1) + 1;
  ht->b1 = rand() * ht->P;
  ht->b2 = rand() * ht->P;

  ht->ref_counter = 1; /* at least one place (the calling structure/function) is using this hashtable */

  return ht;
}

void 
del_bip_hashtable (bip_hashtable ht)
{
  if (ht) {
    int i;
    if (--ht->ref_counter) return; /* some other place is using this hashtable, we cannot delete it yet */
    for(i=ht->size-1; i>=0; i--) if (ht->table[i]) { 
      if(ht->table[i]->key) del_bipartition (ht->table[i]->key); 
      free (ht->table[i]); 
    }
    if (ht->table) free (ht->table);
    free (ht);
  }
}

uint32_t biphash1 (bip_hashtable ht) { return ((((ht->a1 * ht->h) + ht->b1) % ht->P) % ht->size) ; }

uint32_t biphash2 (bip_hashtable ht) { return ((((ht->a2 * ht->h + ht->b2) % ht->P) % (ht->size - 3)) | 1); }

void 
bip_hashtable_insert (bip_hashtable ht, bipartition key) 
{
  uint32_t i;
  int h1, h2;
  
  ht->h = bipartition_hash (key) % ht->P; /*P_ is a large prime*/
  h1 = biphash1 (ht);
  h2 = biphash2 (ht);

  ht->probelength = 0;
  
  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (bipartition_is_equal (ht->table[i]->key, key)) {  // already observed
      if (ht->maxfreq < ++ht->table[i]->count) ht->maxfreq = ht->table[i]->count; // notice the 'count++' doing main work 
      return; 
    }
  }
  /* alloc space for new key */
  ht->table[i] = biomcmc_malloc (sizeof (struct hashtable_item_struct));
  ht->table[i]->key = new_bipartition_copy_from (key); 
  ht->table[i]->count = 1;
  return;
}

double
bip_hashtable_get_frequency (bip_hashtable ht, bipartition key) 
{
  uint32_t i;
  int h1, h2;
  
  ht->h = bipartition_hash (key) % ht->P; /*P_ is a large prime*/
  h1 = biphash1 (ht);
  h2 = biphash2 (ht);
  
  ht->probelength = 0;

  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
    ht->probelength++;
    if (!(ht->table[i])) return 0.;
    else if ( (ht->table[i]) && (bipartition_is_equal (ht->table[i]->key, key)) ) return (double)(ht->table[i]->count)/(double)(ht->maxfreq);
  }
  return -1.;
}
