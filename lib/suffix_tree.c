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

/*! \file suffix_tree.c 
 *  \brief  Ukkonen's Suffix Tree Construction
 *
 *  Following ideas from git@github.com:Jodh/Ukkonen_Algorithm.git (see README_algorithms.md for full list) 
 */ 

#include "suffix_tree.h"

#define size_of_char 255 // max degree of each node (number of children)

st_matches new_st_matches (void);
void st_matches_insert (st_matches match, int id);
void freeSuffixTreeByPostOrder (STNode n);
STNode new_STNode (int start, int *end, suffix_tree suftre);
int edgeLength (STNode n);
int walkDown (STNode currSTNode, suffix_tree suftre);
void extendSuffixTree (suffix_tree suftre, int pos);
void setSuffixIndexByDFS (STNode n, int labelHeight, suffix_tree suftre);
int traverseEdge (STNode node, char* p, int p_length, int pos, suffix_tree suftre);
STNode findLocusSTNode (char* pattern, int pattern_length, suffix_tree suftre, st_matches match);
void subTreeDFS (STNode u, suffix_tree suftre, st_matches match);
int sizeof_suffix_tree_below_node (STNode u);

suffix_tree 
new_suffix_tree (char *input_text, size_t text_size, bool create_text_copy)
{
  int i;
  suffix_tree suftre = (suffix_tree) biomcmc_malloc (sizeof (struct suffix_tree_struct));
  suftre->text_allocated_here = create_text_copy; // used by del_suffix_tree() 
  if (suftre->text_allocated_here || (input_text[text_size] != '\0')) {
    suftre->text_allocated_here = true;
    suftre->text = (char*) biomcmc_malloc ((text_size+1) * sizeof (char));
    strncpy(suftre->text, input_text, text_size);
    suftre->text[text_size] = '\0'; // original code has '$' but not used... 
  }
  else suftre->text = input_text;

  suftre->splitEnd = NULL; 
  suftre->n_splitEnd = 0;

  suftre->size = (int) text_size;
  suftre->lastNewSTNode = suftre->activeSTNode = NULL;
  suftre->activeEdge = -1;
  suftre->activeLength = 0;
  suftre->remainingSuffixCount = 0;
  suftre->leafEnd = suftre->rootEnd = -1;
  suftre->root = new_STNode (-1, &(suftre->rootEnd), suftre);
  suftre->activeSTNode = suftre->root; 
  for (i = 0; i < suftre->size; i++) extendSuffixTree (suftre, i);
  setSuffixIndexByDFS (suftre->root, 0, suftre);
  return suftre;
}

void 
del_suffix_tree (suffix_tree suftre)
{ 
  if (!suftre) return;
  freeSuffixTreeByPostOrder(suftre->root);
  if ((suftre->text_allocated_here) && (suftre->text)) free (suftre->text);
  if (suftre->splitEnd) free (suftre->splitEnd);
  suftre->text = NULL;
  free (suftre);
}

st_matches 
new_st_matches (void)
{
  st_matches match = (st_matches) biomcmc_malloc (sizeof (struct st_matches_struct));
  match->idx = (int*) biomcmc_malloc (sizeof (int)); // could be NULL and let realloc take care  
  match->n_idx = match->length = 0;
  match->is_partial = false;
  return match;
}

void
del_st_matches (st_matches match)
{
  if (!match) return;
  if (match->idx) free (match->idx);
  free (match);
}

void 
st_matches_insert (st_matches match, int id)
{
  match->idx = (int*) biomcmc_realloc ((int*) match->idx, sizeof (int) * (match->n_idx + 1));
  match->idx[match->n_idx++] = id;
}

void 
freeSuffixTreeByPostOrder (STNode n)
{
  if (!n) return;
  int i;
  if (n->children) {
    for (i = size_of_char; i >= 0; i--) freeSuffixTreeByPostOrder (n->children[i]);
    free (n->children);
  }
  free (n);
  return;
}

STNode 
new_STNode (int start, int *end, suffix_tree suftre)
{
  int i;
  STNode node = (STNode) biomcmc_malloc (sizeof (struct STNode_struct));
  node->children = (STNode*) biomcmc_malloc (size_of_char * sizeof (STNode)); 
  for (i = 0; i < size_of_char; i++) node->children[i] = NULL;
  node->suffixLink = suftre->root;
  node->start = start;
  printf ("DBG::end::%p::%d\n", end, *end);
  node->end = end;
  node->suffixIndex = -1;
  return node;
}

int 
edgeLength (STNode n) {
  return *(n->end) - (n->start) + 1;
}

int 
walkDown (STNode currSTNode, suffix_tree suftre)
{
  int curr_edgeLength = edgeLength (currSTNode);
  if (suftre->activeLength >= curr_edgeLength) {
    suftre->activeEdge     += curr_edgeLength;
    suftre->activeLength   -= curr_edgeLength;
    suftre->activeSTNode = currSTNode;
    return 1;
  }
  return 0;
}

void 
extendSuffixTree (suffix_tree suftre, int pos)
{
  suftre->leafEnd = pos;
  suftre->remainingSuffixCount++;
  suftre->lastNewSTNode = NULL;

  while(suftre->remainingSuffixCount > 0) {
    if (suftre->activeLength == 0) suftre->activeEdge = pos; 
    if (!suftre->activeSTNode->children[(int)suftre->text[suftre->activeEdge]]) {  // first character of edge not found
      suftre->activeSTNode->children[(int)suftre->text[suftre->activeEdge]] = new_STNode (pos, &(suftre->leafEnd), suftre);
      if (suftre->lastNewSTNode) {
        suftre->lastNewSTNode->suffixLink = suftre->activeSTNode;
        suftre->lastNewSTNode = NULL;
      }
    }
    else {  // first character of edge found
      STNode next = suftre->activeSTNode->children[(int)suftre->text[suftre->activeEdge]];
      if (walkDown (next, suftre)) continue;
      if (suftre->text[next->start + suftre->activeLength] == suftre->text[pos]) {// Rule 3 found, end phase //
        if(suftre->lastNewSTNode && (suftre->activeSTNode != suftre->root)) {
          suftre->lastNewSTNode->suffixLink = suftre->activeSTNode;
          suftre->lastNewSTNode = NULL;
        }
        suftre->activeLength++;
        break;
      }
      // Rule 2 found, Create new internal node //
      suftre->splitEnd = (int*) biomcmc_realloc ((int*) suftre->splitEnd, (suftre->n_splitEnd + 1) * sizeof (int));
      suftre->splitEnd[suftre->n_splitEnd] = next->start + suftre->activeLength - 1;
      STNode split = new_STNode (next->start, suftre->splitEnd + suftre->n_splitEnd, suftre);
      suftre->activeSTNode->children[(int) suftre->text[suftre->activeEdge]] = split;
      suftre->n_splitEnd++;

      split->children[(int) suftre->text[pos]] = new_STNode (pos, &(suftre->leafEnd), suftre);
      next->start += suftre->activeLength;
      split->children[(int) suftre->text[next->start]] = next;

      if (suftre->lastNewSTNode) suftre->lastNewSTNode->suffixLink = split;
      suftre->lastNewSTNode = split;
    }

    suftre->remainingSuffixCount--;  // decrease remaining leaf nodes to be created 
    if ((suftre->activeSTNode == suftre->root) && (suftre->activeLength > 0)) { // update activeSTNode for next extension //
      suftre->activeLength--;
      suftre->activeEdge = pos - suftre->remainingSuffixCount + 1;
    }
    else if (suftre->activeSTNode != suftre->root)  suftre->activeSTNode = suftre->activeSTNode->suffixLink;
  }
}

void 
setSuffixIndexByDFS (STNode n, int labelHeight, suffix_tree suftre)
{
  if (!n)  return;
  int i;
  bool leaf = true;

  for (i = 0; i < size_of_char; i++) if (n->children[i]) {
    leaf = false;
    setSuffixIndexByDFS (n->children[i], labelHeight + edgeLength (n->children[i]), suftre);
  }
  if (leaf) n->suffixIndex = suftre->size - labelHeight;
}

int 
traverseEdge (STNode node, char* p, int p_length, int pos, suffix_tree suftre) {
  int i, flag=0;
  for(i = 0; (i < edgeLength (node)) || (p[pos] == '\0') || (i == p_length); i++, pos++) if (suftre->text[(node->start) + i] != p[pos]) {
    flag = -1;
    break;
  }
  if ((p[pos] == '\0') || (i == p_length)) return pos + 1; // perfect match
  else if (flag == -1) return - pos - 1; // mismatch
  return 0;
}

STNode 
findLocusSTNode (char* pattern, int pattern_length, suffix_tree suftre, st_matches match) {
  STNode u = suftre->root;
  int k, pos = 0;
  while ((pattern[pos] != '\0') && (pos < pattern_length)) {
    if (u->children[(int) pattern[pos]]) u = u->children[(int) pattern[pos]];
    else break;
    /*if len of p runs out returns 0, if whole edge matches returns 1, if missmatch occurs return -1*/
    k = traverseEdge (u, pattern, pattern_length, pos, suftre);
    printf ("DBG::k::%d, %d\n", k, pos);
    if (k == 0) { pos = pos + edgeLength (u); continue; }
    else if (k > 0) { match->is_partial = false; match->length =  k-1; return u; } 
    else if (k < 0) { match->is_partial = true;  match->length = -k-1; return u; } 
  
  } 
  match->is_partial = true; match->length = pos; return u; 
}

st_matches 
new_st_matches_from_pattern (char *pattern, suffix_tree suftre)
{
  st_matches match = new_st_matches();
  int pattern_length = (int) strlen (pattern);
  STNode u = findLocusSTNode (pattern, pattern_length, suftre, match);
  subTreeDFS (u, suftre, match);
  if (match->n_idx > 1) qsort (match->idx, match->n_idx, sizeof (int), compare_int_increasing);
  return match;
}

void 
subTreeDFS (STNode u, suffix_tree suftre, st_matches match)
{
  if(u && (u->suffixIndex != -1)) st_matches_insert (match, u->suffixIndex); // single hit
  if (u == NULL) return;
  int i;
  for(i = 0; i < size_of_char; i++) if (u->children[i]) {
    if(u->children[i]->suffixIndex == -1) subTreeDFS (u->children[i], suftre, match);
    else  st_matches_insert (match, u->children[i]->suffixIndex);
  }
}

// memory usage tracking function 
int 
sizeof_suffix_tree (suffix_tree suftre)
{
  return sizeof_suffix_tree_below_node (suftre->root) + suftre->n_splitEnd;
}
 
int 
sizeof_suffix_tree_below_node (STNode u) 
{		
  if (!u) return 0;
  if (u->suffixIndex != -1) return sizeof (u);
  int i, size = 0;
  for (i = 0; i < size_of_char; i++) if (u->children[i]) size += sizeof_suffix_tree_below_node (u->children[i]);
  return size;
}
