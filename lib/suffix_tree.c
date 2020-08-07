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
int size, leafEnd,rootEnd, splitEnd;
};


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

  suftre->lastNewSTNode = NULL;
  suftre->activeSTNode = NULL;
  suftre->activeEdge = -1;
  suftre->activeLength = 0;
  suftre->remainingSuffixCount = 0;
  suftre->splitEnd = -1;
  suftre->leafEnd = -1;
  suftre->rootEnd = -1;
  suftre->root = newSTNode (-1, &(suftre->rootEnd));
  suftre->activeSTNode = suftre->root; 
  for (i = 0; i < text_size; i++) extendSuffixTree (suftre, i);
  setSuffixIndexByDFS (suftre->root, 0);
}

void del_suffix_tree (suffix_tree suftre)
{ 
  if (!suftre) return;
  freeSuffixTreeByPostOrder(suftre->root);
  if ((suftre->text_allocated_here) && (suftre->text)) free (suftre->text);
  suftre->text = NULL;
  free (suftre);
}

void freeSuffixTreeByPostOrder(STNode *n)
{
  if (!n) return;
  int i;
  for (i = 0; i < size_of_char; i++) if (n->children[i]) freeSuffixTreeByPostOrder (n->children[i]);
  if ((n->suffixIndex == -1) && (n->end != NULL)) free (n->end);
  free (n);
}


// Structure and Functions for extending tree  
STNode *newSTNode(int start, int *end)
{
  STNode *node =(STNode*) malloc(sizeof(STNode));
  int i;
  for (i = 0; i < size_of_char; i++) node->children[i] = NULL;
  node->suffixLink = root;
  node->start = start;
  node->end = end;
  node->suffixIndex = -1;
  return node;
}

int edgeLength(STNode *n) {
  return *(n->end) - (n->start) + 1;
}

int walkDown(STNode *currSTNode)
{
  if (activeLength >= edgeLength(currSTNode)) {
    activeEdge += edgeLength(currSTNode);
    activeLength -= edgeLength(currSTNode);
    activeSTNode = currSTNode;
    return 1;
   }
  return 0;
}

void extendSuffixTree(int pos)
{

  leafEnd = pos;
  remainingSuffixCount++;
  lastNewSTNode = NULL;

  while(remainingSuffixCount > 0) {
    if (activeLength == 0) activeEdge = pos; 
  
    if (activeSTNode->children[text[activeEdge]] == NULL) {  // first character of edge not found
      activeSTNode->children[text[activeEdge]] = newSTNode(pos, &leafEnd);
      if (lastNewSTNode != NULL) {
        lastNewSTNode->suffixLink = activeSTNode;
        lastNewSTNode = NULL;
       }
     }
    else {  // first character of edge found
      STNode *next = activeSTNode->children[text[activeEdge]];
      if (walkDown(next)) continue;
      if (text[next->start + activeLength] == text[pos]) {// Rule 3 found, end phase //
        if(lastNewSTNode != NULL && activeSTNode != root) {
          lastNewSTNode->suffixLink = activeSTNode;
          lastNewSTNode = NULL;
         }
        activeLength++;
        break;
       }
      // Rule 2 found, Create new internal node //
      splitEnd = (int*) malloc(sizeof(int));
      *splitEnd = next->start + activeLength - 1;

      STNode *split = newSTNode(next->start, splitEnd);
      activeSTNode->children[text[activeEdge]] = split;

      split->children[text[pos]] = newSTNode(pos, &leafEnd);
      next->start += activeLength;
      split->children[text[next->start]] = next;

      if (lastNewSTNode != NULL)    lastNewSTNode->suffixLink = split;
      lastNewSTNode = split;
     }
  
    remainingSuffixCount--;  // decrease remaining leaf nodes to be created 
    if (activeSTNode == root && activeLength > 0) { // update activeSTNode for next extension //
      activeLength--;
      activeEdge = pos - remainingSuffixCount + 1;
     }
    else if (activeSTNode != root)  activeSTNode = activeSTNode->suffixLink;
  }
}
// End Structure and Functions for extending tree Definitions  

//Functions for printing Suffix Tree and Labeling Leaf nodes//
void setSuffixIndexByDFS(STNode *n, int labelHeight)
{
  if (n == NULL)  return;

  int leaf = 1;
  int i;
  for (i = 0; i < size_of_char; i++)  if (n->children[i] != NULL) {
      leaf = 0;
      setSuffixIndexByDFS(n->children[i], labelHeight + edgeLength(n->children[i]));
     }
  if (leaf == 1)  n->suffixIndex = size - labelHeight;
}

// free allocated memory where children are NULL// 
//End of Functions for printing Suffix Tree and Labeling Leaf nodes//

// pattern search  version2.0
STNode* pickEdge(STNode* node, char* p, int pos) {
  char c = p[pos];
  if(node->children[c] != NULL) return node->children[c];
  return NULL;
}

int traverseEdge(STNode* node, char* p, int pos) {
  int i, flag=0;
  for(i = 0; i < edgeLength(node) || p[pos] == '\0'; i++, pos++) if(text[(node->start) + i] != p[pos]) {
      flag = -1;
      break;
    }
  if(p[pos] == '\0') return 0;
  else if(flag == -1) return -1;
  return 1;
}

STNode* findLocusSTNode(char* p) {
  STNode* u = malloc(sizeof(STNode));
  u = root;
  int pos = 0;
  while(p[pos] != '\0') {
    u = pickEdge(u,p,pos); // give next child if exists, NULL other wise
    if (u == NULL) break;
    /*if len of p runs out returns 0, if whole edge matches returns 1, if missmatch occurs return -1*/
    int k;
    k = traverseEdge(u, p, pos);
    if (k == 1) { pos = pos + edgeLength(u); continue; }
    else if (k == 0) { puts("locus node found\n"); return u; }
    else if (k == -1) { puts("pattern not found\n"); return NULL; }	
  }
  return NULL;	
}

void subtreeDFS (STNode* u) 
{
  if(u != NULL && u->suffixIndex != -1) printf("Index:%d\n", u->suffixIndex);
  if (u == NULL) return;
  int i;
  for(i = 0; i < size_of_char; i++) if(u->children[i] != NULL) {
      if(u->children[i]->suffixIndex == -1) subtreeDFS(u->children[i]);
      else printf("Index:%d\n", u->children[i]->suffixIndex);//Leo:idx of match
    }
}
// end of pattern search v2.0 //

// memory usage tracking function 
int sizeofTree(STNode* u) 
{		
  if(u == NULL) return 0;
  if(u->suffixIndex != -1) return sizeof(u);
  int i, size = 0;
  for(i = 0; i < size_of_char; i++) if(u->children[i] != NULL) size += sizeofTree(u->children[i]);
  return size;
}
// end of memory usage tracking function //

int main()
{
  //puts("Enter string to build suffix tree on or press enter for default string: ");
  //gets(text);
  //if(strlen(text) == 0)
  //strcpy(text,mississippi$");
  char whole_text[22000] = "ggcttggatgaggctgtactgcgaccatttgccttttgtggtatcatatcctccctgatacaggcgccaatgttttattttaattgggtgattgacagtcaaaattgcctctgtggccggccggtcacgagtaagaatggtgattctcgaagcgaattgtttcggcgccggccggttcattattaaaacaagccgtcggtttaatttgagtaagaccagtggaaaccggaaattgccgcaggttatccaaccggaaattttatagccagtactatccagagcctcaatgtacgctgcggggacgttaccatagtcagtagagggcatgaaacccgattcggtggggattgcttggggataaaaatcccggatgaatattacatatccgcctaatttcaatttcaatcccgccgtaataaatacttcgcctttagatcgatctttgaaaggtgccccagtagtgatattcattagcaagagttcaggggcatattcttcgatccaaaattcggataactttagcgtaaagggcaaatcctgacgcataccattaagatctatgccgctggattctggtttgccctggtgaagaaccaggtttaacttgtgcttatcaccggctcccaaagccatggctactatgatcagccagaatcccaggtggctgaccaggaagcgaatattccggtagcttaggggaactgctcttttgaggattacaaaggttagagtagtcaggaggtatagccaggcaaagaaaaattgccaagagcgccagaggtgcgtaaaatccagtagagaatgttcatgtagattagaagattgaggtaccagtccaataattagtgtcatgacgccgaaatagcttatggcactgagggcggcaggcaaattagtcaaccactgaaccaggggatgatgtcttttaagacgataggcggtaattgccagcgtaagaaaagcgataacaattagtaggttaaaaggaaaaccaatcaggggtgggggaggatcaacaagcacctcaagcaaagtgccaaacaggcacaagcctattgctatcgctaaacttttacaataaccccacgggaattgccagatgatcatgtcatttaaaacttgcttttccattgtggagtttgggatgacatccgagcacactaagttgattcgtaaatattcggtgatgcgatctcacattaagaatatagaatgctctaagataattctggagataattcttgaggtacatataatttttattggatagcttgttcggagaaactgtcccaattatccgcgatctaggttggattgggatattagtagtgcatttcattaattaatcttgtcgttcaggtagtcccgcttgccctttttacgctttttggtgatttaaccgaataaataacaaggcacattaattatagagcaatcaatgccaaacgatttgtcacacctataggggcgataagttaacctttcaagagcgaatcatgactagcaacatcttgaaaggtgagactgcttttacagcgtgtggtttgtgggctggcatgatgagcatttcgccggattttagccggtataccttgcctgagatagtgacatcggcttcaccatctacgcaatacaccagagcatcaaagggggcggtgtgctcgcttagtccttctccctggctaaaggcaaacaaggttaccgatccggtctcacgttcgaagattgtgcggctgaccactgccccaggctgataattgattaaaccactcaaattaaatgcactgtccaagagtgtttccatagctctcctccatttggattagtatttaaaaaacacaattctgagtttattttcagcagattatgcttttatataatccgcgttgcgctttcgattggctttattttgattttttggattacttccaatgtacccaggctttcacctctggggtaaacagcaaattcgcggctataaaggatttgcccctagagcttggctcaaaatatcgcccgatcaagacaattattgccgtagttgagaatttccgcgctgtgacgcaaaccctgcttgaagtgaaggagggtgcttgacaattaataaccacgatagtaagggcagaggaccttttgccgatgcagaagagagaatttccgtttttggtattgacttgctggtatgtatggttgttgtgaatctaaaattgttcaaataatataaaccagcgcggattaaccgccgtatctatttcacgaggagcatttttttatattgggtaaagtcaccagaccggaggcgatataaataaattccgctggagagctgggcgccatcgaattgtagctgataatcaccgacaggcaaaattttattggtaagtacttgaatcaaccggccttgcaaatcaaaaacctcaattgttataggtgcttttgtgaggatgctgaagctaatccaggtgacagaattaaagggattaggataattctgattgagaacccatcgctgcgggatgccggtctcaggtgaattggcaatatctgcaagctgctcgatattggcatcgatccagccgattcgcaaagctagccagttttttaggtaattgatttcttgctggaaagtcgtggcggtgcccatttcaggccagcggaggtagttgcggtattgggcttcggttacaacgttttgcaaggaatcaatcatgcgataaagggtgtccaggtgaaaggtagtagagcgcagggcaaaccagcgctgagctgcatggtaggcgaatgtcgggcgtcggaccaacacttcccaccagaatggaatcttgggaatgtcgaaggttaaacgtttattatggtcgatttcccagtcgccataaaggttgcggtcttgctcgtaccaggcgtggccgaaagtcaggttcatgtcccagatagggccggcttgtagtttgccacctttactgtcgcgatttttatatagaaaggcgctcaggcgataggcgtcgatattacgggccagctcgtttagcaaaaaataatccacgaatgaagcttcgtctataatcccggcataaccactgataggatcgtccacttgcggacttttcatgatggattcaaaagcggcgatgaatttttggatgtactttttctgctgagcggtaatctcgctgggcttggggtaatcatactggtactgaatattgttaggagagtaccaatagcccacgttctcgccctcgatcttatccaacttgatgatgtagccgccagttaatgcatcaccgctcgtatcctgaggctcgcatttggcaatattgacacggtttttgtcccgtttgattttctccatcaggatataaatacctttgtaatcattattaatgatcaattcgaagtagcgggtacggctggcgtaacgtccgatcgcgcgggacaactcgaaaggaagcacattgcggataaaggttttgtcatcataaggaccatagaaaatccaatcggcttccgctggaaagcccagcaaagcaacatctaggtcttctccgttagcatcaatggtttcgacgcggtactgtttcttaggccataacattgagctgctgcctcgcaattcgatgacgatttttccatcgtaattgttatagggatcagtgatgtaattgcgattgcctttaccattgtcgatgacgcccatgtcggctacgatacggctttcatcgacgatggtctgcccgtgggtgtcaatgacgataattggcaggttagatgaggttaggttgacctgagccagtagccaccagggtattaggattagcagcgctactgacagcgagattgttgccagcggtttggtctggtttctgggtttgcttaaatcaatctggaaattactgagcagaccagttttggaatgcgagtctatcatcatgaattaatctatgcagtttaaatcctcgatgaaatgcttttttcctttggtaatgattgtatgcccctctttttcaaggaggtctttatggctggtggcgccgccggggaatttttcgttcagttcgccattgtttttgattatgcgccaatacggcaggtcacttccggcttctgcagcagcgttggcagctatcataacgaaaatagcggtggtcagcggacagcaaaattttgtgccgtattgttgggctaattgcgcgcaaatggtactgagcgtcatcagtttcccgaacggcacggtctggatgattgcgactacagcggaaggaggggccagtacgacactatcaccaggttgggcaccgctttttagtaaggctttcccacaggggaaccggggattaaaggacagaattttagggaagcctttgttatcggctaatttttcagcccaggttttctttttattagccatatcgacctgagttttttacctaacgatttttctcaatttgtgatacccgagccgctttggcgccgggcagtgaacgtgccttgtttagttggtggccagccgaaagttgtaaaagcaatgctccattcaccttgggcttgaccatcatgcgcagtgcctttgacgttaatcgtaaacggtgaggtattatagcctgccgggattgtaggattcggattggtagcgacgccttggccactgacctggatcgaatcgccatttacgatgagagttccctcttcgaccgggcattctactggatacccatcggatataaaaatccagatagccgtagtagtgatgctaccatcagatttcttatcgaaggtagcttcgccatgattggaagaatcattatccataatgaactcccagacttcggtagctacgatcgttggctcatcatcaaccgtggtacaacctaaaatgagcaggccaagccccatagcaattatcagactaattttggtcgacatagctcctccttaaatttatgttgtcagcgattagttctatgaaggacgcgtcatttctattagccctgaatgggtgaagtaatttcgtctccgagtggggagagtgaaaagccaaattcgagtggctcaaagtcgatccttgcgctgggttttttgataatttttcaccattttggagaaggatttttcgcgctttcgcctttcggcgagggtcgtttcaaaatgcgtcctttcccgctctagtttgagataattttcataagcagcccgatcaatcgctccatctgcaatcgctgccaggactgcacagccggcttcgtgggtgtgggtgcaatcgcgatacttacaattttgccccaaggcaagtatttggtcgaaggtcatttccaggccgttatcggtgtcggcaaggccaatttcacgcatcccgggggtgtcaatcaaaagacccccatttcccagcatgatcaactggcggtgagtagtaacatggcgccctttaagcgtgtgctggctgatggtgtcggttttcatcagttgtttgcctgccagattattgagtagagtagatttgcctacaccggaagaacccagcagacagtaggttttgcctttttcgataaagcggtcgatagcctcatagccttgacgagttgcattactaattgccaaaactggcacctggccaatccggactttgagattctcgatcaatcttgggacacaatcagcctcgatcaagtcgattttggtaagcacgattattggcataaccttggaagagtagcagattgtcaaataacgttccaggcggttcaggttgtagtcccgatcaatcgcttggagaataaaggcatagtcaatgttagtggcaatgatttgaacttcaccgtactggccgacggcttgccgtttgattacagaaaagcgtgggaggatgcggtggatcagcgcaaagtcgggttcgtaggtcgtcacgagtacccaatcaccgacgctggggaaatcttcctggctttgcgccgcatagcgcaagctccctgtaatttcggctttgaattcaccttgtcttgtgcggacaagatagagttcgcggtattcagctactacccgtcctacctcaaaagcactcgcatctatgccctgataggcgttttctaattctttattgtaacctaaatcaaataaatccataaatggtttgactaatctccttgtttattcctgattcttagatgctgtcgggtagctggtaggacagtcacgaagacaatcattgtcccgaggagggaatagataaagtggaatagatgaatagagcgaatgaacaggacctgaaaatcaatccagagaattaggattattccgatggcgaaagctagcgaccaagaccaatggtattctttgtctagattcagatactcaaaaaagcgacaatcgggacgcttgatcagggcataggcggtgtaccctggcagtagtcccaggatgcacaacaggaaaagtcccggaattagaaagttgggaaatggtgaatgctccagcagcgagagtggcatctgcagcaacgccccactgggactaatcaccagcgctaagccgccaccgatggcaccaagggattggaaaatgaggagacccaccaagacaaagagcgaaaatttttgaaatgaattcatgctcctattttgcctccaatttttcagttaccttcgtattactgaatcgataaataatgataattacctcgctttttttaccaattttcggtaataatataactcatttagattaagcggttaaaatgttgatgctacgtacggaggggagttttttacctaatattttttaatatcgcggtgatttgaatactggcggtaccatctccgaagatattggtttgcacggaagtaaaattagacgccttcgcgatttcttctggtatacgctgaaaatccgtgccgatgagctttgcccagccattttcaactaactcaacccattccgtctcggtacggagtacaagacacggcactccggcgaaataggcttctttctgaacaccgccggaatcggttagaatgatacgcgcacttttttccagcgccagcatatctagatagccgactggttcaatgaggcgcaggtaggaatttttgatggtcagctgctgttggcggatgaatttcatggttcgggggtggaccgggaagagaataggttcgtctaagtactctaggatgttgatcagtgctttcaaaatggtcggatcgtcagtgttttgagggcgatggatagttgtcaagatgaattttcttggggcaacaccatattttttaagaacggcctcgatgtcaactcgcgagagattgaacaagaacgaatcgtacatcacatcgcctaccaaatagacccccgtgttgatgccctcgttttgcaggtgtttgatcgccgtctgggtcggggcaaataagaactgggcaagggcatcagccacaatacgattaatctcttcaggcatgtcgcgattgaaactacgcagaccagcttcgacatgggcaagtggtagctgcaatttggcggcaactagtgcccccgccatagtagaattggtatccccatacaccatgaccaggtcgggtttgatatgcagcaatactttttcgagctcaatcagcatgcgaccggtttgctcgccatggttgccactaccgataccgaggttgtagtccggttctggaatggcgagctctttgaaaaaggcttgcgacatattatagtcgtaatgctgacctgtatggatcaggatttcgttgaagtgtttgcggatttcccgggaaactggtgcacatttgataaactgcggacgcgcaccgatgatggtggcgatggtgggcatcaattcaaactcctggttttagattggaaataaaatagcgcatcttaccggtggattccggcctaatttctggaactaaaataaaatcaattttgatatcggacgtctgaacggtgcgctcaaatattttgcggaatagttctaaatccgattcaatgaaatcacgactttttacaaaccgcaccgtaaaagcagacggcgcggtctgaacaatctgatactgctcaaatggatgacgccccaaatccaacagcgccaggttcaggtagtcaaagatttcactggaaaaagagcggccgtcactcataatcacaagatcggtggccttgcccgtgtgtaacttcatgattggcagattacatccacagggacatttttcagctacccaatcgccaatgtcgccgatcctataacggataaagggcatatattcattggttaacgatgtcccgatgatttctttggttccgtttttattttcgacaaattctagaaaaacattttcgatggaaacgtgccagttcccaaatgggcattcataggcaaagcccccgatctccgaacagccgtattcattggagatgggcactttgaagactgacgcaatgatttctcgctcgtgtggaaagagggtttcggcggtgacgaaaatcgccttcagccgggacggtactggttggttattatttctgaaataaagtgccagacgataaatgctgttggcgtaaccgtagaaatagcgtgggttaaagcgccgaatctggcggcaatattttgccagaatttcatcggagtattcgaacccggaaattagcagcgtgttgcggatacgacccttgaggacatcgataactctgcgataagtcgaatagagcggtcgcccccaaatggccacttccgggtcgccaatttgtattccgaaccaattccgtgcccgccaacgactagcccaatcccagctgttatattcgccagtaacataaaattctaacggaatgccagtcgagccgctggttttaccgggatatatttttccggcataatgcagattgagaaattgttcgaccttttcggcgattgtctttttcgttagtaacggcaggagtttaaaatcctctacggtgcggatgtcttccggtttaatgccgagccggtcgaaggtttcatgataatagggaatgttgcgataagccactttgagaatatgccagagtttggtcagctggattttctgcaattcttccaacggcagaggttcgttgtttttcaagattctgagatattgggtgattgattcacccttaagacgctgggtgagtggaaagtaaaaatgccggacgaattggctatacatcagcgaatttcctgtaaaaccctttcaatttttcgactagtgatctcccaggtgaaatattgcgtaacgattttgcgagccaattctcccaatgcgatccgccgtttttcgtctagcagcaattcaataagtgtctgggcaaatttttgctgatccttcacatcgacaatctcaacacacctaaaattagtaactaattctctaacggtttcaattggactagcaacgactgttttgccgcatgccaaataatcgaataacttcatgggcgaaccgctaattttattgaaaaattcggtatcccatggcgacacgcaaatatcgaatgaatttatataatacggcgccatctcgaacggcactttgccagtgaaggtgaacgatttactcacatttaattggcgactcagggctttccattcctctaattgggctccatcgccagcaactacgaaacggacctgcggcattttctcaagaatcaaaggcgcggcctgaatcaaatactgcacgccgtgataatgatagcaacttccgatgaaaccgatcgtcggagcatccggtagactggtttgggcaatagcttccgatttgggcattggccggaaaaggtttggattagtgccattcgcgatgggaataattttagcggctttgaccggataattggcgcaaatcaggtttttcaagccttggctgactggaataatgcgcgcgcagtagccgaagttgaggcgttgaatcaggcgcatataggcaatcttccatggtgagtaaccggtctggcgtaattcttccaacacccagccgttgatttcgatagtcagcggtacattgcagagtttcgagactgccagcggaataaatgacgaggtttcccggaaatataaatgatccggccggtgcttcaggatattgtatatcaggtaacagaacgaaaaaatctcgaatattatccagcgcacaagcggaaggttaatgaccggaatcgaatgaatggtgactttaatgtcgaaagtgggcttgacgatgtccggattgaaaaaatgcacctcatgcccaaattttaccaaccagccgatggtttcgcgaatgtgattgttaatgttacccagtcggctaaagccctcgaaacaaaagtaataaattctgaaaggtcgcttcaaagtgatttttaagcctggttaatattagagaaaggataaaggttcatgacttcaaatttaaaatatagttatggtccctatattcacattttgatatttcaggagggtaggtttttatggttttgaatatgatggagttagccattgccaagttggagtttctcgaaaattttaaggtactgttcggcatttacactccaggaagtggtgcgcggcgttatttttccagtattattgagagcatcgtgtatcgccttgattaaagcgttcaagtcattgggttgacatagaaaaccgcctgaataattttccagaatttcaggaataccaccgaccttagtagccacaaccggtttgcccatcgccatggcttccagaatcacattcgggacgccttccgcaagactgggcaaaacaatcacatcggaggctccgtaccagagcggcaactgcgaatgtggaatgttacccggaagaagaacctgcccttgcagtgagttttccgtgatgaaattctcgattaagttccggtcaatgccctgacctattaaaaccaatttcatattcggttggctcatggccaaatatgccgaaagcaacaggtcaatgcctttacggcggctcaaattaccgacaaaaagcagataggtggaatctggagacatatttaattcccggcgtgatcggagttgatcaattggatgaaaaagtgattcatcgacgccattgtagattacggcaattttcgaggcgggatagccacgtcgggttagctctgttttcagagcggctgaaactgttacgagctggtcagcttgagtcagcgttcgatgcatgaagcgcttcaagatcagccgggtgagataggtattgacgtccgagccgcgcagaccgactacggatttcttgccatatttttgcgctaaccggacggccgcatatccacccgggcagcccatttgtcctacaataacatcaaaagcatggatttctggtttcagatccgcttctgccaatctgtaaaattgcccggaacgccgccagatgcggtcgaagcgatcgatcactttaattcgctgaatccaggcaaatcgttgattagtgccgtaatttttacggaaactgactaaatctggatcacagaccaggcggttcctttccctgaaaatcacctggtaaaccgtaaggtcacattttaaggccaggtgaacagcccattgatagataaaaatcccgtggttgggattggttgggtgggggaatagatcggtaaggactaagactttcacgatgtcaattctaaatatacggcttctacatgattcacaaaatttgtaatagaatatttttcagccgttcttttagcgttgctagaaattttctgatataataccggatcgttatataactcagcaatagtttgggctagcgcttgggccgcatccaaccgattcccctggaaggtcttttcccatttttccggaaaagtcagtcgcccatcatacccattatggatgatttcagaaattacaccaacattgagtcccacgaccggtagtccgtaactcattgcttctaataaggcgatcccgaatccttcggtatgaatggagggaaacaggtaacaattcaggctggcgtacaccagactcaagtctttttgaaacggtaggatagtgactgattcttggagtttattctcagcgatcaggcgtttgtgctcggcctggctggcttcatcgcccatgacgacggcgtgtgccggaacgcccatggcgcgcagctttttcagggtcaagtagaagatgttaaagcctttgaaaggatcaaagcggctggtgttgccaatgagaaagcaattggagggtatatgccatcgttttttaaaggctttaataacggcggtgtcggatttgcgcggattttcgaccggattatatattagccggatgcgtttttgcggcgcatattcggtgttaattaaatgcgcagtaacggttttagagattgataaaatccgctcggccagaaagcggttgagcaggatccggccgagtttgaaaaaaagaccagttttaattcgaatattatgtttagtataaatgcaatggcgaacgccattcagcctggcgaaaaacgtcccaaagaattctagtgcgccataatgcgaatgagcgatatcgatgcgatgccgccggatgatttggcgcaattctctaaaagtctctaaatcatgccaggattgcagatgcaggtgatagagcgggatattcagagaacggaaatctttatagaagaaatcatttacaccgtcagaatacatggtaaccacaaccggattgaatttggtacgatcgatattttgcaataggagcaagagcgatttttcggctccgccgacctgtaaactgccgatgatatgcagaacattgatctttttcattttaggattcgggctacttctctaagcggttggtcacgcgtgattattttatcgatgaacagcaagtgccagagttccaaccagaagagcgaccagagccgatgtttgtgctcggcagtgtgggaacgttgttcgctccacagcttattgatgcagtcgtaattgaaataattacgttgccgggcaaatcggctgaaaacaaccgcctcgaaaattccactcaggtcattttgaatccagcggtcgactggaacgccgaatccacgtttcgcggtgtaaaggagctttttgggcagatatttttcggcgattttctttaggagatatttggtctggtggtttttgaatttcagtttaggggggatttgagcgacgtattccagtaaagccgtgtccagaaatggcgaacggcattcgatggagttggccatggttgcccggtcgattttgactaagtaatcgtagaccaagcgggtcattccatcggtgaaaagcacctgatcaattagattaattcctttgactttgttgaaattatcatgaaacgggtgaagtggattgtgttcggataaatgctgtaaatgttccggatgccatagctcatggcggtatttatgaaatcccatcgtattgtaaaagctttcgtacgggaaccctttcaggtagtcgatataaaacagagccttgccgaaaaccgaatcgctttcgggattgaccgacaaatttttcaaaatgttgccgactagcggttggagaaaaggtggaacaatccgggaaaataactgggcattgcgcgggactattgtccgaccataaccgccgaaagcttcatcgccaccatcgccggttagcatcaccgtaatgtattgacgcgcggctttggagatcaaataagtggggatcgctgaagaatcggcgaaaggctcaccgaattgccagacgagttcgagcagattatcagtaacattatgatcgacaatgacctcagtgtgatcggtgtcgtaaagtttggcaacttcccgggcggcaagcacatcttccatcgggtaatctttaaaacccatcgtaaaggtcttgatccgctgggtagagttacggcatagaatggcagtaatgacactggaatctaccccgccactgagcatcgtcccgagtggcacatcgctgatagtgcgtcgcatgactgcttcgacgatcaactgctcgcagtggtgcaagatctctggttcgctaccgaaagtttttaatgaataatcccaatgccaataggtagtggtggcggatgagttctcagaaaacactgtgaatgtacctggctggattttctgcacgccctcgaaaatcgtatgttgttccggaatagcgatgtggtgcaaatagcaatcgattgctgtcaggctgaccggtggttgtgacggcagcgccgcgatgatggctttaatgtccgaggcgaaataaactaggccggattgttcggcatagtagagaggctttttgccgacgcggtctcggccgagaattagttgatttttacgcgcatcgtaaatcgcgatggcaaacatgcccgtgatttttcggaataattctgttccccattcggcgtaaccatgcacgagcacttcggtatcagaattagattgaaattggtaaccagcagcgagtaaggcttgacgcagttcttgaaaattataaatctcaccgttgtaaactgtccaaacggtacggtctgcattagtcattggctgattggcggcgtctgacaggtcaatgatcttcaatcgtcggtgaccgagagctattttcggcagaattgcataaccggcggaatccggaccacgttcgatcatacaatcccgcatccgacagatttcagctggagtgacctcgattctaccctgatgaaaaatgccacaaatgccgcacatgttactcgcctaatttattaagcgtcaattgataaagattcactaggtgccggatattaggatggcccaattaaccacagaataattttttgataggcgttttcttccagaccagtgatttgaattgccctggtgcctcctatgaggcaaggattatggtcttactatcggcactgctcgcaaccaacctctacatgagcgatactattattgagacgctcacaacagtgttactcgcatccgctacaccccgcacagttattatcgattagcgaatattcatccttatgttttagataaaactccagtaaagccggtttgaatatatccttctcatcagctgaaaggtcgttatataaattcacctgctcccacatcatttcggcgatatcgtcttcgtccacggggaaagatatcaacgttcccagtacgcattcatccaaccccagcagacgacgctctgcggtgttgattttcaatataaattcctttttggtacagctacaacccatacaattctccctgttcgattaaaagatctggcggtcatcccggaaatatttgtttaccccaaaaatgataccaccgaagattagcaccgccaccggcacgacgatccagtgctggacccccaaggcctcgaatagattgccattaatggatctggcaagcagatgggctttgaaaaactggacattcagcgcatacaccaggattcccagcaccatgccgagtaatcccagccaggcgtcctttttgccttcgccaatgcgtgatagagcagtgccggggcaatagcccagtaaacccatgccgacgccaaaaatcagcccgccgactaatttgccccagtccagcggcttgggattgaactgaatcaagccaaggtcgctcataagcgtaaaacctaccagggaaatggtcacggcagtgagcatgaacttgagaattttgaaatcttttaacagcagtaagcccgtcacgcgcggatacttggtgatgccagtacgttgcagaatcaatccgaaagccaagcccatcaacagaccaatgactattttcaggatcatttcatacctccggtttttggatagagcagtttggcggtgatgatcccggtggcaaaaaacgaggctccggcgaccatactgcccagcgacagatgcgcccagcccatcaggatattaccggtcgtacagcctgatgccagaatcgccccgaagccaagacaaaatccgccgatcagcacaacgatcagccgtttccaaatattgggaccgtggtattttacccacaccgcgggaatggtttcaggtgaatactgacgtgatgatctagcaccgacatagccgccaatgatgattccgatcagcgcgaagaattccacaaaggcagctggtgaagagcttaggcttttgatgagcggattggttgccagggacgggaagatgtggcggaactgagcgctgatgtagccgacgccagcagtgatcccgaaaggataatccttgatgctcagcgtctgaaacagccagtaactgaacattcccagcagcgacgccagcaaaccgccggtgaaccagtgccattcgccgttttgattgaaaagttttttcatgtgctctccattttactaactatggtttgaatgtttcgtttcatctgtttttcagctattccggctcgatcatttttcatcgtttggtaatttaattgatgggttgttcgccgaggatgattcctggtgaagctttgatcccgaatttttcgatgatttcttcggcgcccggatcatccacatatttcaccacgtatggcaatcctgattccttcagaaagcgctcgacgtggttgctcttgcaacatttttgcgaacggattaagagaatttgcattttcactgagacttagctatttatgatctcaaccgatttaaaaaattaccatgcctattaaaccacagccaatgagcaaatagattggattcaactgtttaaactttaaaaaaacaaagaatgccactatggccagtagatacggacgaaaatccctgaaaccacccgctgccatgagaaaaacagtccgcataagcaggacgagtacgatgaccgtgatagcgcttaagaaatttttatagattttgttctctttgattttaactaaggcgatgatgacaatcagggctaaaaggaacggcaaaatatcgtaacctattgtagcaataaccgcccctaaaaatcccttgagcctgaaaccgataaatgtcgacagattcatcccgattggacccggagtggtttcggcgatcgagaccatatagagaaaatcctgttcggtgatccattgtagattttcaacgatttcacggtgtaaaatcgggaaaacggcattgccgccgccaaaagcaatgatgccaatcttgcagaattcaaggaaaattttaactaacatcgcttctcatcatcctactaatcggcatgttaatgagttttctgacatcaacggtaatatcaaatttgctgggaaagcaaagtggcattttcgggattctgctgacggctttacttggtgcggtcatgcatatcccttcgattgtcgcctttccgctggcagcgtcattgctggaaaaaggggcgtcggtcatggcaattgctactttcatcactaccctgacaatgattgggattgtcactttgccgttggaaatcaaagagttgggaaagaaattcgcgttctggcgcaacggactcagtttcattttggcaattctggttggcctcctgatgggagttgtcctgtgaaaccgccatcaaatcaaaaatcgaaactgaccaaacgtaccgactggctgatcctggcaattgcagttcttacgggcgttgtcctcttactaattttcccagagaaacaggaaatcgctgggataaccgcccggcgtattgcgacggagatgtttaccatttttccagcagtactaattttaatgggattattcgccgtctgggtgtcaaaagagacagtgatgcaatttctgggcgaaaactccgggttaaaaggcattgcgctggcgctgtttttcggggcattacccaccggaccgctctacatcgcttttccgttggccgccggtttgcaccaaaaaggagcaagtcttggtaatatcgccattttcctgactgcctgggcttgtatcaaactgccgcaggaattggtggagctgcagtttctgggagttcggtttatgcttgcccgtctgggcttgacaatcctggccgcagtgattatgggtatgattatcgataaaattgtccaggcggaccattttccggacaataaatcaagggagaatattctatgaaagtatttgatctcaaggccatgcaagccgctgactatgccgaacgcggtaaaaatgtgttctacagtacgccagaatttaaaacccggattatcgaactgccggctggaggtacgatgccacaatgtgaaatggcgtccgatgtcgtctttgtcgttatagatggcacagcaactgtcactgttaatcaacaggaagagcagttgaaagctggtcagtgcctgatcaccgaacccgccatgctctcaatgaaaacggaaaaaggcgttaaaattcttggtattcagattcaaaaacaaaaatagattagaattaaaaattataaggtacagaattatgataagaaaaatggctctgataactacgctcagtttatctgtcggattatttctcggttgccaggcaaaaaaggagaccgctccagtctctgaaaaagctgtagcaaaattagtgcaggttgaagtcgttaaaacgcgcaaaatggtggaaaccctggacctgactggcacattgcgggcggagaatgtcgccaatattctctcgaccgtcgagggtaaaatctcctgcctgttggtacgggaaggtgatcaggtagaacccaaccaagtcgtggcgatgatcagttcactggtgcgtgaagatattatcaacgcggcccggctccgcatggaagccgctaaaagacaattaaacgacaatcctgaaaatccccaattcaaagtcgcgtatgaacaggcgcaacaggattatgaatttgccgtgcaacaatacaaagaaattccggtaacatccccaatgcaagggctggtttcccggcgttgggtggatttgggggatatggtcccttccaaagctaaactctttgaaattcaaagcagcgccaaactgattgtggatgtaccggtatcagagcttgatcttagaaaactgaaaaatgagcagaaggctaaaatattcgttgatgcctgtccagaaaaacaattccaaggcgctattcagcgcattcatccccaggttgatgcccagacccgcaatggattggtggaaatcatccttcttgacccctgccctaatctgaaatccggcatgttcgtgcgggtcaccttcgatgtccgcaccattgaaaatgcaatcgcagttccggtgcaagccatcatcgaacgtcctcaatacaaaacctgttttatagtcaacgatgacaaggctcaagaaaaaataatcacaacgggtctggaaagcggcggctgggttgaaattctttccggtctttctgtcggagaaaaaatcgtgatcgaaggtcagcaacaactcaaaaccggttcgcctgttaaaataaaagggaagaagtaaccacaaaaacacaaagaacacagaaaaaaatatgttttaaaattttaaaataagtaatctgttttgtatcttctatgcctttgtggcaaaaagtaaagaaataatgtccaccgagacacttcagcaagagtctccaatgggacgaggaacccagagaattcattatgattactaaaaatttgaaaataaaaaaagtcttcgtgccctgcgtgtccaaaaaggaatcttaaatgaagctaacaactatttccattcatcgtccggtagcaaccattgtgttaattattaccgcggttgtacttggatttttcgggtacagccgcatgccgatcgatttcttgccggaaatcacttacccgatgatcaaggtttatgtctattggcgtggtgccacaccggaagaaattaacgacaacatcgccgacccgattgaaagaaccctcagcacggttgacgatctcgactacctggaatcgtcgtcaattgaaggcctctataccttattggtgaatttcgaatacggcaccgacgtagacgtcgcctatcaggatgttgtcgctaaaatgggactagcaacccgtaaattgccacctgatgtcgatccgccgataatttttaaagctgatccttctcaattgcccgttgttcagctaattgttcaatcggaaggccgtgatttggttaaactgcgcgaatggatcgaaaatgttgtgcaggatcagtttctcgccgtcaagggcgtcgccgggactgaaattctcggcggattgaaacgtgaaatccgcattcacctcgatcccaaacgaactacagcctacaatttaaccttaaatacggtgattaaacgcctccaagaagagaacatcgagcgtttagctggacgtgtcaccgaaagcgggcgtgaatttatcgtccgtactaatgcagaatttcgtaacctagatgacattcgcaatgtcattattttcaacgacggcaaatccatggtgcgcctgaaagatctggcaaccgttgaagacagtcacgaagaacaacgggtcattacccggtttaacagcaaaccagccatcaaactcaacattttaaaacaagccgatgctaacactgtcgatattgccgagcaggtcaacgcgcttattgcaaaactagatgatactgcgccgccggacatcacttttaatacagtcgaaaaccaggcggactatatcaagggcgccattgccggggtacgggattccgctattattgccgtcttactggtaatccttgtcatttatatttttctgggacattggcggcaggtagttatcatgctgattgccctcccggtaactattttattcaatttctttttgatgcaactgcttgggttttcaatcaatatcttttccttgggcggattagtggtcgcgatgggtgtggtgctagataattccatcgttgttattgagaatatcacccgcctgcatcttgaaaaaggcagccgtaccgatacagtcgcagacactgccaccagtgaggtcgcgacggcagtactggcttccactttaacttttatggcgatatttctgccgtttttgct";
  time_t start, end;
  double cpu_time[22];
  int size[22];
  int i  = 0;
  const char dollar[] = "$";
  for(i = 0; i < 22; i++) {
    strncpy(text, whole_text, i*1000);
    strcat(text, dollar);
    start = clock();
    buildSuffixTree();
    size[i] = sizeofTree(root);
    freeSuffixTreeByPostOrder(root);
    end = clock();
    cpu_time[i] = ((double) (end - start))/CLOCKS_PER_SEC;		
  }

  for(i = 0; i < 22; i++) {
    //printf("Time taken, #of characters:[%f sec, %d]\n", cpu_time[i], i*200);
    printf("%d,%f,%d\n", i*1000, cpu_time[i], size[i]);
  }

  buildSuffixTree();

  char pattern[1000];
  while(1) {
    puts("Enter pattern to search: ");
    gets(pattern);
    if(strlen(pattern) == 0)
      break;
    printf("Pattern: %s\n", pattern);
    subtreeDFS(findLocusSTNode(pattern));	

  }
  freeSuffixTreeByPostOrder(root);
  return 0;
}
