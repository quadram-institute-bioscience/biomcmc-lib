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

/*! \file hashtable.h
 *  \brief double hashing open-address hash table using strings as key -- also has distance matrix, that can be used in
 *  alignments and trees
 *
 * Hash tables allow us to search for the position of a key (taxa name) without scanning the whole vector (like in 
 * sequencial search). This code is derived from the software DCM3, released under the GPL license (Copyright (C) 
 * 2004 The University of Texas at Austin). 
 */

#ifndef _biomcmc_hashtable_h_ 
#define _biomcmc_hashtable_h_ 

#include "lowlevel.h"

typedef struct hashtable_struct* hashtable;
typedef struct hashtable_item_struct* hashtable_item;

/*! \brief key/value pair for hash table */
struct hashtable_item_struct {
  /*! \brief String (vector of char). */
  char* key; 
  /*! \brief Integer (position in vector where hashtable_item_struct::key can be found) */
  int value;
};

/*! \brief Hash table (vector indexed by strings). */
struct hashtable_struct 
{ 
  /*! \brief Table size. */
  int size;
  /*! \brief Number of collisions before empty slot is found. */
  int probelength;
  /*! \brief Value set by hash(). Used in hash1() and hash2() to avoid calling hash() again. */
  uint32_t h;
  uint32_t a1, /*!< \brief Random values used in hash functions. */
                a2, /*!< \brief Random values used in hash functions. */
                b1, /*!< \brief Random values used in hash functions. */
                b2, /*!< \brief Random values used in hash functions. */
                P;  /*!< \brief Random values used in hash functions. */
  /*! \brief Vector with key/value pairs. */
  hashtable_item* table;
  /*! \brief Counter of how many external references (structures sharing this hashtable) to avoid deletion */
  int ref_counter;
};

/*! \brief Insert key/value pair into hashtable. */
void insert_hashtable (hashtable ht, char* key, int value);
/*! \brief Return location (value) of corresponding key (string) or negative value if not found. */
int  lookup_hashtable (hashtable ht, char* key);
/*! \brief Create new hashtable of size elements. */
hashtable new_hashtable (int size);
/*! \brief Free hashtable space. */
void del_hashtable (hashtable ht);

/* extra hash functions */

uint32_t biomcmc_hashint_1 (uint32_t a);
uint32_t biomcmc_hashint_2 (uint32_t a);
uint32_t biomcmc_hashint_3 (uint32_t a);
uint32_t biomcmc_hashint_4 (uint32_t a);
uint32_t biomcmc_hashint_5 (uint32_t a);
uint32_t biomcmc_hashint_6 (uint32_t a);
uint32_t biomcmc_hashint_7 (uint32_t a);
uint32_t biomcmc_hashint_8 (uint32_t a);
uint32_t biomcmc_hashint_9 (uint32_t a);
uint64_t biomcmc_hashint_64bits (uint64_t key);
uint32_t biomcmc_hashint_64to32 (uint64_t key);
uint32_t biomcmc_hashint_mix (uint32_t a, uint32_t b, uint32_t c);
uint32_t biomcmc_hashstring_1 (unsigned char *str); /* CRC algortihm */
uint32_t biomcmc_hashstring_2 (unsigned char *str); /* PJW algorithm */
uint32_t biomcmc_hashstring_3 (unsigned char *str); /* dbj2 algorithm */
uint32_t biomcmc_hashstring_4 (unsigned char *str); /* dbj2 algorithm with XOR */
uint32_t biomcmc_hashstring_5 (unsigned char *str); /* sdbm algorithm */

#endif
