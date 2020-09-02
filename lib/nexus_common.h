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

/*! \file nexus_common.h 
 *  \brief File handling functions for nexus format in general.
 */

#ifndef _biomcmc_nexus_common_h_
#define _biomcmc_nexus_common_h_

#include "char_vector.h"
#include "empirical_frequency.h"

/*! \brief maximum name length for taxa (alignment and tree files). */
#define MAX_NAME_LENGTH 4096 

/*! \brief General function that stores file content into char_vector_struct, removing shell-type comments */
char_vector new_char_vector_from_file (char *filename);

/*! \brief Removes (possible nested/multiline) nexus comments of the form [] (brackets). */
char* remove_nexus_comments (char **string, size_t *stringsize, FILE *stream);
/*! \brief Changes uppercase characters by lowercase versions. */
char* lowercase_string (char *string);
/*! \brief Changes lowercase characters by uppercase versions. */
char* uppercase_string (char *string);
/*! \brief Removes spaces, tabs from string. */
char* remove_space_from_string (char *string);
/*! \brief returns bool::false if string is composed only of space characters (' ', '\n', '\r', '\t', etc). */
bool nonempty_string (char *string);
/*! \brief returns bool::false if first nonspace character of string is ';' (FASTA comment) or '#' (HUPO extension) */
bool nonempty_fasta_line (char *string);
/*! \brief returns bool::false if first nonspace character of string is '#' except for '##' which is pragma */
bool nonempty_gff3_line (char *string);

#endif
