/* SPDX-License-Identifier: GPL-3.0-or-later
 *
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

/*! \file suffix_tree.h 
 * Original code: https://github.com/mattporritt/ukkonen_suffix_tree ( GPL-3.0-or-later) 
 */ 

#ifndef _biomcmc_suffix_tree_h_
#define _biomcmc_suffix_tree_h_

#include "lowlevel.h"


#define size_of_char 256
typedef struct STNode_struct* STNode;
typedef struct suffix_tree_struct* suffix_tree;

struct STNode_struct 
{
  struct STNode children[size_of_char];
  struct STNode suffixLink;  //pointer to other node via suffix link
  int start, *end, suffixIndex;
};

struct suffix_tree_struct
{
  char *text; 
  bool text_allocated_here; // tells if text is link (make sure not freed upstream) or alloced (uses more memory)
  STNode root, lastNewSTNode, activeSTNode;
  int activeEdge, activeLength, remainingSuffixCount;
  int size, leafEnd,rootEnd, splitEnd; // rootend, splitend are pointers in original
};

#endif
