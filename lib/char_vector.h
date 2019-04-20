/* 
 *This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
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

/*! \file char_vector.h 
 *  \brief list of strings (each string is a vector of chars)
 */

#ifndef _biomcmc_char_vector_h_
#define _biomcmc_char_vector_h_

#include "lowlevel.h"

typedef struct char_vector_struct* char_vector;

/*! \brief vector of strings (char vectors) of variable length */
struct char_vector_struct
{
  char **string;  /*! \brief vector of strings */
  int nstrings;   /*! \brief how many strings */
  size_t *nchars; /*! \brief length of allocated memory for each string excluding the ending '\0' (the actual size in 
                    use needs strlen() or a call to char_vector_compress() over the structure )*/
  int ref_counter;/*! \brief how many times this char_vector_struct is being used */
  int next_avail; /*! \brief next available position (empty string) */
};

/*! \brief Create a vector of strings with initial size for each string of zero */
char_vector new_char_vector (int nstrings);
/*! \brief Create a vector of strings from subset of strings of another char_vector */
void new_char_vector_from_valid_strings_char_vector (char_vector vec, int *valid, int n_valid); 
/*! \brief Create a vector of strings where each string is assigned an initial value of nchars */
char_vector new_char_vector_fixed_length (int nstrings, int nchars);

/*! \brief Delete vector of strings only after nobody is using it */
void del_char_vector (char_vector vec);

/*! \brief Link a previously allocated string (to avoid copying all characters) */
void char_vector_link_string_at_position (char_vector vec, char *string, int position);

/*! \brief Add a new string (vector of characters) at specific location */
void char_vector_add_string_at_position (char_vector vec, char *string, int position);

/*! \brief Add a new string (vector of characters) at next available location */
void char_vector_add_string (char_vector vec, char *string);

/*! \brief Append string at the end of existing string at location */
void char_vector_append_string_at_position (char_vector vec, char *string, int position);

/*! \brief Append string at the end of existing string at most recently used location */
void char_vector_append_string (char_vector vec, char *string);

/*! \brief Increase size of vector of strings (called automatically by other functions) */
void char_vector_expand_nstrings (char_vector vec, int new_size);

/*! \brief update order of strings in vector based on a vector of new positions */
void char_vector_reorder_strings_from_external_order (char_vector vec, int *order);

/*! \brief Reduce size of vector of strings by removing empty strings (returns number of empty strings) */
int char_vector_remove_empty_strings (char_vector vec);
/*! \brief Remove identical strings and resizes char_vector_struct */
int char_vector_remove_duplicate_strings (char_vector vec);
/*! \brief reduce char_string_struct to only those elements indexed by valid[] */
void char_vector_reduce_to_valid_strings (char_vector vec, int *valid, int n_valid);

/*! \brief Order char_vector_struct elements from longer string to smaller, or lexicographically; can be used after calling
 * new_char_vector_from_file() but if topology etc. are associated to it, then order[] must be externally defined and
 * will have new locations, to keep track of changes */
void char_vector_reorder_by_size_or_lexicographically (char_vector vec, bool lexico, int *order);
/*! \brief If the two char_vectors are identical (same strings in same order), then delete one and make it point to the other one */
bool char_vector_link_address_if_identical (char_vector *v1, char_vector *v2);
/*! \brief find occurences of species->string[] inside gene->string[] filling indexes in sp_idx_in_gene.
 *
 *  The species are taxon names which may be associated with topologies or alignments, such that we can not reorder its
 *  elements here (without also modifing e.g. tree leaves). But ordering from longer to shorter is essential for pattern finding, 
 *  so it is assumed that the char_vector is already sorted UNLESS user provides the ordering. */ 
void index_species_gene_char_vectors (char_vector species, char_vector gene, int *sp_idx_in_gene, int *order_external);
void update_species_count_from_gene_char_vector (char_vector species, char_vector gene, int *sp_count);

#endif
