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

//void buildSuffixTree()
suffix_tree new_suffix_tree (char *input_text, bool create_text_copy)
{
  int i;
  suftre = (suffix_tree) biomcmc_malloc (sizeof (struct suffix_tree_struct));
  text_size = strlen(input_text);
  suftre->text_allocated_here = create_text_copy; // used by del_suffix_tree() 
  if (suftre->text_allocated_here) {
    suftre->text = (char*) biomcmc_malloc ((text_size+1) * sizeof (char));
    strncpy(suftre->text, input_text, text_size);
    suftre->text[text_size] = '\0'; // original code has '$' but not used... 
  }
  else suftre->text = input_text;

  suftre->lastNewSTNode = suftre->activeSTNode = NULL;
  suftre->activeEdge = -1;
  suftre->activeLength = 0;
  suftre->remainingSuffixCount = 0;
  suftre->splitEnd = suftre->leafEnd = suftre->rootEnd = suftre->size = -1;
  suftre->root = new_STNode (-1, &(suftre->rootEnd), suftre);
  suftre->activeSTNode = suftre->root; 
  for (i = 0; i < text_size; i++) extendSuffixTree (suftre, i);
  setSuffixIndexByDFS (suftre->root, 0, suftre);
}

void del_suffix_tree (suffix_tree suftre)
{ 
  if (!suftre) return;
  freeSuffixTreeByPostOrder(suftre->root);
  if ((suftre->text_allocated_here) && (suftre->text)) free (suftre->text);
  suftre->text = NULL;
  free (suftre);
}

void freeSuffixTreeByPostOrder (STNode n)
{
  if (!n) return;
  int i;
  for (i = 0; i < size_of_char; i++) if (n->children[i]) freeSuffixTreeByPostOrder (n->children[i]);
  if ((n->suffixIndex == -1) && (n->end)) free (n->end);
  free (n);
}

STNode new_STNode (int start, int *end, suffix_tree suftre)
{
  STNode node = (STNode*) biomcmc_malloc(sizeof(struct STNode_struct));
  int i;
  for (i = 0; i < size_of_char; i++) node->children[i] = NULL;
  node->suffixLink = suftre->root;
  node->start = start;
  node->end = end;
  node->suffixIndex = -1;
  return node;
}

int edgeLength (STNode n) {
  return *(n->end) - (n->start) + 1;
}

int walkDown (STNode currSTNode, suffix_tree suftre)
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

void extendSuffixTree (int pos, suffix_tree suftre)
{
  suftre->leafEnd = pos;
  suftre->remainingSuffixCount++;
  suftre->lastNewSTNode = NULL;

  while(suftre->remainingSuffixCount > 0) {
    if (suftre->activeLength == 0) suftre->activeEdge = pos; 
    if (!suftre->activeSTNode->children[suftre->text[suftre->activeEdge]]) {  // first character of edge not found
      suftre->activeSTNode->children[suftre->text[suftre->activeEdge]] = new_STNode (pos, &(suftre->leafEnd), suftre);
      if (suftre->lastNewSTNode) {
        suftre->lastNewSTNode->suffixLink = suftre->activeSTNode;
        suftre->lastNewSTNode = NULL;
      }
    }
    else {  // first character of edge found
      STNode next = suftre->activeSTNode->children[suftre->text[suftre->activeEdge]];
      if (walkDown (next)) continue;
      if (suftre->text[next->start + suftre->activeLength] == suftre->text[pos]) {// Rule 3 found, end phase //
        if(suftre->lastNewSTNode && (suftre->activeSTNode != suftre->root)) {
          suftre->lastNewSTNode->suffixLink = suftre->activeSTNode;
          suftre->lastNewSTNode = NULL;
        }
        suftre->activeLength++;
        break;
      }
      // Rule 2 found, Create new internal node //
      suftre->splitEnd = next->start + suftre->activeLength - 1;

      STNode split = new_STNode (next->start, &(suftre->splitEnd), suftre);
      suftre->activeSTNode->children[suftre->text[suftre->activeEdge]] = split;

      split->children[suftre->text[pos]] = new_STNode(pos, &(suftre->leafEnd));
      next->start += suftre->activeLength;
      split->children[suftre->text[next->start]] = next;

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

void setSuffixIndexByDFS (STNode n, int labelHeight, suffix_tree suftre)
{
  if (!n)  return;
  int i;
  bool leaf = true;

  for (i = 0; i < size_of_char; i++) if (n->children[i]) {
    leaf = false;
    setSuffixIndexByDFS(n->children[i], labelHeight + edgeLength (n->children[i]), suftre);
  }
  if (leaf) n->suffixIndex = suftre->size - labelHeight;
}

int traverseEdge (STNode node, char* p, int pos, suffix_tree suftre) {
  int i, flag=0;
  for(i = 0; i < edgeLength (node) || p[pos] == '\0'; i++, pos++) if(suftre->text[(node->start) + i] != p[pos]) {
    flag = -1;
    break;
  }
  if(p[pos] == '\0') return pos + 1;
  else if (flag == -1) return - pos - 1;
  return 0;
}

STNode findLocusSTNode (char* pattern, suffix_tree suftre) {
  STNodex u = suftre->root;
  int pos = 0;
  while(pattern[pos] != '\0') {
    if (u->children[pattern[pos]]) u = u->children[pattern[pos]];
    else break;
    /*if len of p runs out returns 0, if whole edge matches returns 1, if missmatch occurs return -1*/
    int k;
    k = traverseEdge (u, pattern, pos, suftre);
    if (k == 0) { pos = pos + edgeLength (u); continue; }
    else if (k > 0) { puts("locus node found\n"); return u; } // match size = pos -1
    else if (k < 0) return NULL;  //not found;  partial match of size = -k-1
  }
  return NULL; // partial match of size pos
}

void suffix_tree_find_match (char *pattern, suffix_tree subtre)
{
  STNode u = findLocusSTNode (pattern, suftre);
  subtreeDFS (u, subtre);
}

void subtreeDFS (Node* u, suffix_tree subtre)
{
  if(u && (u->suffixIndex != -1)) printf ("Index:%d\n", u->suffixIndex); // single hit
  if (u == NULL) return;
  int i;
  for(i = 0; i < size_of_char; i++) if (u->children[i]) {
    if(u->children[i]->suffixIndex == -1) subtreeDFS(u->children[i], suftre);
    else printf("Index:%d\n", u->children[i]->suffixIndex);//Leo:idx of match
  }
}

// memory usage tracking function 
int sizeof_suffix_tree (suffix_tree suftre)
{
  return sizeof_suffix_tree_below_node (suftre->root);
}
 
int sizeof_suffix_tree_below_node (STNode u) 
{		
  if (!u) return 0;
  if (u->suffixIndex != -1) return sizeof (u);
  int i, size = 0;
  for (i = 0; i < size_of_char; i++) if (u->children[i]) size += sizeof_suffix_tree_below_node (u->children[i]);
  return size;
}
