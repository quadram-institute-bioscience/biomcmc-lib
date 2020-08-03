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

#define MAX_CHAR 256
unsigned char *text; //Input string 1 (for suffix tree)
unsigned char *text2; //Input string 2 (for generalized suffix tree)
extern char *tree_string; //String to hold text version of tree
unsigned char print_enabled;

struct SuffixTreeNode {
    struct SuffixTreeNode *children[MAX_CHAR];
    struct SuffixTreeNode *suffixLink;    // pointer to other node via suffix link
    /*(start, end) interval specifies the edge, which will be stored in the child node. Lets say there are two nodes A and B
     connected by an edge with indices (5, 8) then this indices (5, 8) will be stored in node B. */
    int start;
    int *end;
    int suffixIndex;    /*for leaf nodes, it stores the index of suffix for the path  from root to leaf */
};

typedef struct SuffixTreeNode Node;

extern Node *root; //Pointer to root node

/*lastNewNode will point to newly created internal node,
  waiting for it's suffix link to be set, which might get
  a new suffix link (other than root) in next extension of
  same phase. lastNewNode will be set to NULL when last
  newly created internal node (if there is any) got it's
  suffix link reset to new internal node created in next
  extension of same phase. */
extern Node *lastNewNode;
extern Node *activeNode;

struct substring {
        int stringstart;
        int stringend;
};

struct substring subtrings[MAX_CHAR];

/**
 * Create new Suffix Tree Node
 */
Node *newNode(int start, int *end);

/*
 * Return the edge length of given node
 */
int edgeLength(Node *n);

/*activePoint change for walk down (APCFWD) using
 Skip/Count Trick  (Trick 1). If activeLength is greater
 than current edge length, set next  internal node as
 activeNode and adjust activeEdge and activeLength
 accordingly to represent same activePoint*/
int walkDown(Node *currNode);

void extendSuffixTree(int pos);

void print(int i, int j);

void setSuffixIndexByDFS(Node *n, int labelHeight);

void freeSuffixTreeByPostOrder(Node *n);

/*Build the suffix tree and print the edge labels along with
suffixIndex. suffixIndex for leaf edges will be >= 0 and
for non-leaf edges will be -1*/
int buildSuffixTree(unsigned char *string1, unsigned char *string2, unsigned char print_tree);

int doTraversal(Node *n, int labelHeight, int* maxHeight, int* substringStartIndex);
int allCommonSubstringsTraversal(Node *n, int labelHeight);
int substringAllOccurenceTraversal(Node *n, unsigned char* str, int idx);
void longestRepeatedSubstringTraversal(Node *n, int labelHeight,
                                       int* maxHeight,
                                       int* substringStartIndex);

char *getLongestCommonSubstring(unsigned char *string1, unsigned char *string2, unsigned char print_tree);
char *getAllCommonSubstrings(unsigned char *string1, unsigned char *string2, unsigned char print_tree);

int checkForSubString(unsigned char* search_string, unsigned char* source_string);
int *checkAllSubStringOccurences(unsigned char* search_string, unsigned char* source_string);

char * getLongestRepeatedSubstring(unsigned char* source_string);

int doTraversalToCountLeaf(Node *n);
int countLeaf(Node *n);

int substringTraversal(Node *n, unsigned char* str, int idx);

int traverseEdge(unsigned char *str, int idx, int start, int end);

/*
 * builds string that contains text version of suffix tree
 * for displaying on console
 */
void buildString(char** current_text, const char *new_text);


#endif
