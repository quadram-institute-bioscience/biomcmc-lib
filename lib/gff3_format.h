/* SPDX-License-Identifier: GPL-3.0-or-later 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.

 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file 
 *  \brief GFF3 format 
 */

#ifndef _biomcmc_gff3_format_h_
#define _biomcmc_gff3_format_h_

#include "alignment.h"

static char *feature_type_list[] = {"gene", "CDS", "mRNA", "region"}; // many more, but I dont use them

struct
{
  char *str;
  uint64_t hash;
  int id;
} gff3_string;

struct 
{
  int start, end;
  uint8_t pos_strand:1, phase:3, type:4;
  gff3_string seqid, source, attr_id, attr_parent;
} gff3_fields;

#endif
