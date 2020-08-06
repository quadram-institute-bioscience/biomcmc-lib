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

/** Create new Suffix Tree Node  */
Node *newNode(int start, int *end);

/* Return the edge length of given node  */
int edgeLength(Node *n);

/* activePoint change for walk down (APCFWD) using Skip/Count Trick  (Trick 1). If activeLength is greater
   than current edge length, set next  internal node asactiveNode and adjust activeEdge and activeLength
   accordingly to represent same activePoint*/
int walkDown(Node *currNode);

void extendSuffixTree(int pos);

void print(int i, int j);

void setSuffixIndexByDFS(Node *n, int labelHeight);

void freeSuffixTreeByPostOrder(Node *n);

/* Build the suffix tree and print the edge labels along with suffixIndex. suffixIndex for leaf edges will be >= 0 and
  for non-leaf edges will be -1*/
int buildSuffixTree(unsigned char *string1, unsigned char *string2, unsigned char print_tree);

int doTraversal(Node *n, int labelHeight, int* maxHeight, int* substringStartIndex);
int allCommonSubstringsTraversal(Node *n, int labelHeight);
int substringAllOccurenceTraversal(Node *n, unsigned char* str, int idx);
void longestRepeatedSubstringTraversal (Node *n, int labelHeight, int* maxHeight, int* substringStartIndex);

char *getLongestCommonSubstring(unsigned char *string1, unsigned char *string2, unsigned char print_tree);
char *getAllCommonSubstrings(unsigned char *string1, unsigned char *string2, unsigned char print_tree);

int checkForSubString(unsigned char* search_string, unsigned char* source_string);
int *checkAllSubStringOccurences(unsigned char* search_string, unsigned char* source_string);
char * getLongestRepeatedSubstring(unsigned char* source_string);
int doTraversalToCountLeaf(Node *n);
int countLeaf(Node *n);
int substringTraversal(Node *n, unsigned char* str, int idx);
int traverseEdge(unsigned char *str, int idx, int start, int end);

/* builds string that contains text version of suffix tree for displaying on console */
void buildString(char** current_text, const char *new_text);

#define MAX_CHAR 256

struct SuffixTreeNode {
  struct SuffixTreeNode *children[MAX_CHAR];
  //pointer to other node via suffix link
  struct SuffixTreeNode *suffixLink;
  int start;
  int *end;
  int suffixIndex;
};

typedef struct SuffixTreeNode Node;
// Global variable declaration 
char text[22000]; 
Node *root = NULL; 

Node *lastNewNode = NULL;
Node *activeNode = NULL;

<<<<<<< HEAD
int activeEdge = -1;
int activeLength = 0;

=======
/*activeEdge is represented as input string character index (not the character itself)*/
int activeEdge = -1;
int activeLength = 0;

// remainingSuffixCount tells how many suffixes yet to be added in tree
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
int remainingSuffixCount = 0;
int leafEnd = -1;
int *rootEnd = NULL;
int *splitEnd = NULL;
int size = -1; 
// End GLobal Variable Declaration

// Structure and Functions for extending tree  
Node *newNode(int start, int *end)
{
  Node *node =(Node*) malloc(sizeof(Node));
  int i;
<<<<<<< HEAD
  for (i = 0; i < MAX_CHAR; i++) node->children[i] = NULL;
  node->suffixLink = root;
  node->start = start;
  node->end = end;
=======
  for (i = 0; i < MAX_CHAR; i++)
    node->children[i] = NULL;
  /*For root node, suffixLink will be set to NULL For internal nodes, suffixLink will be set to root
    by default in  current extension and may change in next extension*/
  node->suffixLink = root;
  node->start = start;
  node->end = end;
  /*suffixIndex will be set to -1 by default and actual suffix index will be set later for leaves at the end of all phases*/
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
  node->suffixIndex = -1;
  return node;
}

<<<<<<< HEAD
int edgeLength(Node *n) {
=======
int edgeLength (Node *n) {
  if(n == root) return 0;
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
  return *(n->end) - (n->start) + 1;
}

int walkDown(Node *currNode)
{
<<<<<<< HEAD
=======
  /* activePoint change for walk down (APCFWD) using Skip/Count Trick  (Trick 1). If activeLength is greater
     than current edge length, set next  internal node as activeNode and adjust activeEdge and activeLength
     accordingly to represent same activePoint*/
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
  if (activeLength >= edgeLength(currNode)) {
    activeEdge += edgeLength(currNode);
    activeLength -= edgeLength(currNode);
    activeNode = currNode;
    return 1;
<<<<<<< HEAD
   }
=======
  }
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
  return 0;
}

void extendSuffixTree(int pos)
{
<<<<<<< HEAD

  leafEnd = pos;
  remainingSuffixCount++;
  lastNewNode = NULL;

  while(remainingSuffixCount > 0) {
    if (activeLength == 0) activeEdge = pos; 
  
    if (activeNode->children[text[activeEdge]] == NULL) {  // first character of edge not found
      activeNode->children[text[activeEdge]] = newNode(pos, &leafEnd);
=======
  /*Extension Rule 1, this takes care of extending all leaves created so far in tree*/
  leafEnd = pos;
  /*Increment remainingSuffixCount indicating that a new suffix added to the list of suffixes yet to be added in tree*/
  remainingSuffixCount++;
  /*set lastNewNode to NULL while starting a new phase, indicating there is no internal node waiting for it's suffix link reset in current phase*/
  lastNewNode = NULL;
  /* Add all suffixes (yet to be added) one by one in tree */
  while(remainingSuffixCount > 0) {
    if (activeLength == 0) activeEdge = pos; //APCFALZ
    // There is no outgoing edge starting with activeEdge from activeNode
    if (activeNode->children[text[activeEdge]] == NULL) {
      //Extension Rule 2 (A new leaf edge gets created)
      activeNode->children[text[activeEdge]] = newNode(pos, &leafEnd);
      /* A new leaf edge is created in above line starting from  an existng node (the current activeNode), and
        if there is any internal node waiting for it's suffix link get reset, point the suffix link from that last
        internal node to current activeNode. Then set lastNewNode to NULL indicating no more node waiting for suffix link reset.*/
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
      if (lastNewNode != NULL) {
        lastNewNode->suffixLink = activeNode;
        lastNewNode = NULL;
       }
     }
<<<<<<< HEAD
    else {  // first character of edge found
      Node *next = activeNode->children[text[activeEdge]];
      if (walkDown(next)) continue;
      if (text[next->start + activeLength] == text[pos]) {// Rule 3 found, end phase //
=======
    // There is an outgoing edge starting with activeEdge from activeNode
    else {  // Get the next node at the end of edge starting with activeEdge
      Node *next = activeNode->children[text[activeEdge]];
      if (walkDown(next)) continue; //Start from next node (the new activeNode)
      /*Extension Rule 3 (current character being processed is already on the edge)*/
      if (text[next->start + activeLength] == text[pos]) {
        //If a newly created node waiting for it's suffix link to be set, then set suffix link of that waiting node to curent active node
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
        if(lastNewNode != NULL && activeNode != root) {
          lastNewNode->suffixLink = activeNode;
          lastNewNode = NULL;
         }
<<<<<<< HEAD
        activeLength++;
        break;
       }
      // Rule 2 found, Create new internal node //
      splitEnd = (int*) malloc(sizeof(int));
      *splitEnd = next->start + activeLength - 1;

      Node *split = newNode(next->start, splitEnd);
      activeNode->children[text[activeEdge]] = split;

=======
        activeLength++; // APCFER3
        break;   /*STOP all further processing in this phase and move on to next phase*/
       }
      /*We will be here when activePoint is in middle of the edge being traversed and current character
        being processed is not  on the edge (we fall off the tree). In this case, we add a new internal node
        and a new leaf edge going out of that new node. This is Extension Rule 2, where a new leaf edge and a new internal node get created*/
      splitEnd = (int*) malloc(sizeof(int));
      *splitEnd = next->start + activeLength - 1;
      //New internal node
      Node *split = newNode(next->start, splitEnd);
      activeNode->children[text[activeEdge]] = split;
      //New leaf coming out of new internal node
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
      split->children[text[pos]] = newNode(pos, &leafEnd);
      next->start += activeLength;
      split->children[text[next->start]] = next;

<<<<<<< HEAD
      if (lastNewNode != NULL)    lastNewNode->suffixLink = split;
      lastNewNode = split;
     }
  
    remainingSuffixCount--;  // decrease remaining leaf nodes to be created 
    if (activeNode == root && activeLength > 0) { // update activeNode for next extension //
      activeLength--;
      activeEdge = pos - remainingSuffixCount + 1;
     }
    else if (activeNode != root)  activeNode = activeNode->suffixLink;
=======
      /*We got a new internal node here. If there is any internal node created in last extensions of same
        phase which is still waiting for it's suffix link reset, do it now.*/
      if (lastNewNode != NULL) lastNewNode->suffixLink = split;   /*suffixLink of lastNewNode points to current newly created internal node*/

      /*Make the current newly created internal node waiting for it's suffix link reset (which is pointing to root
        at present). If we come across any other internal node (existing or newly created) in next extension of same
        phase, when a new leaf edge gets added (i.e. when Extension Rule 2 applies is any of the next extension
        of same phase) at that point, suffixLink of this node will point to that internal node.*/
      lastNewNode = split;
     }

    /* One suffix got added in tree, decrement the count of suffixes yet to be added.*/
    remainingSuffixCount--;
    if (activeNode == root && activeLength > 0) { //APCFER2C1
      activeLength--;
      activeEdge = pos - remainingSuffixCount + 1;
    }
    else if (activeNode != root)   activeNode = activeNode->suffixLink; // APCFER2C2
  }
}

void print(int i, int j)
{
  int k;
  char *printer;
  for (k=i; k<=j && text[k] != '+'; k++){
    asprintf(&printer, "%c", text[k]);
    buildString(&tree_string, printer);
    free(printer);
  }
  if(k<=j){
    buildString(&tree_string, "+");
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
  }
}
// End Structure and Functions for extending tree Definitions  

//Functions for printing Suffix Tree and Labeling Leaf nodes//
void setSuffixIndexByDFS(Node *n, int labelHeight)
{
  if (n == NULL)  return;
<<<<<<< HEAD

  int leaf = 1;
  int i;
  for (i = 0; i < MAX_CHAR; i++)  if (n->children[i] != NULL) {
      leaf = 0;
      setSuffixIndexByDFS(n->children[i], labelHeight + edgeLength(n->children[i]));
     }
  if (leaf == 1)  n->suffixIndex = size - labelHeight;
=======
  //A non-root node
  if (n->start != -1) if (print_enabled == 1) print(n->start, *(n->end));  //Print the label on edge from parent to current node
  int leaf = 1;
  int i;
  char *printer;
  for (i = 0; i < MAX_CHAR; i++) if (n->children[i] != NULL) {
    if (print_enabled == 1) if (leaf == 1 && n->start != -1) {
      asprintf(&printer, " [%d]\n", n->suffixIndex);
      buildString(&tree_string, printer);
      free(printer);
    }
    leaf = 0;   // Current node is not a leaf as it has outgoing edges from it.
    setSuffixIndexByDFS(n->children[i], labelHeight + edgeLength(n->children[i]));
  }
  if (leaf == 1) {
    for(i= n->start; i<= *(n->end); i++) if(text[i] == '+') {
      n->end = (int*) malloc(sizeof(int));
      *(n->end) = i;
    }
    n->suffixIndex = totalStringLength - labelHeight;
    if (print_enabled == 1){
      asprintf(&printer, " [%d]\n", n->suffixIndex);
      buildString(&tree_string,printer);
      free(printer);
    }
  }
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
}

// free allocated memory where children are NULL// 
void freeSuffixTreeByPostOrder(Node *n)
{
  if (n == NULL) return;
  int i;
  for (i = 0; i < MAX_CHAR; i++) if (n->children[i] != NULL) freeSuffixTreeByPostOrder(n->children[i]);
  if (n->suffixIndex == -1) free(n->end);
  free(n);
<<<<<<< HEAD
=======
}

/*Build the suffix tree and print the edge labels along with suffixIndex. suffixIndex for leaf edges will be >= 0 and for non-leaf edges will be -1*/
int buildSuffixTree (unsigned char *string1, unsigned char *string2, unsigned char print_tree) {
  if(string1 == NULL) return 1; // First check if we have at least one valid string
  // Check if we are printing suffix tree
  print_enabled = print_tree;

  //Check if we have two strings
  //Two strings mean we are building a generalized suffix tree
  if(string2 != NULL){
    asprintf (&text, "%s%c%s%c", string1, 43,  string2, 36);
    string1Length = strlen((const char*) string1);
  }
  else asprintf (&text, "%s%c", string1, 36);

  totalStringLength = strlen((const char*) text);
  int i;
  rootEnd = (int*) malloc(sizeof(int));
  *rootEnd = - 1;

  root = newNode (-1, rootEnd); /* Root is a special node with start and end indices as -1 */

  activeNode = root; //First activeNode will be root
  for (i=0; i<totalStringLength; i++) extendSuffixTree(i);
  int labelHeight = 0;
  setSuffixIndexByDFS (root, labelHeight);

  return 0;
}

int doTraversal(Node *n, int labelHeight, int* maxHeight, int* substringStartIndex)
{
  if(n == NULL) return 0;
  int i=0;
  int ret = -1;
// If it is internal node 
  if(n->suffixIndex < 0) for (i = 0; i < MAX_CHAR; i++) if(n->children[i] != NULL) {
    ret = doTraversal(n->children[i], labelHeight + edgeLength(n->children[i]), maxHeight, substringStartIndex);
    if(n->suffixIndex == -1) n->suffixIndex = ret;
    else if((n->suffixIndex == -2 && ret == -3) || (n->suffixIndex == -3 && ret == -2) || n->suffixIndex == -4) {
      n->suffixIndex = -4; //Mark node as XY

      if(*maxHeight < labelHeight) { // Keep track of deepest node
        *maxHeight = labelHeight;
        *substringStartIndex = *(n->end) - labelHeight + 1;
      }
    }
  }
  // suffix of X
  else if(n->suffixIndex > -1 && n->suffixIndex < string1Length) return -2; //Mark node as X
  // suffix of Y
  else if(n->suffixIndex >= string1Length) return -3; //Mark node as Y

  return n->suffixIndex;
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
}

void buildSuffixTree()
{
<<<<<<< HEAD
  size = strlen(text);
  int i;
  rootEnd = (int*) malloc(sizeof(int));
  *rootEnd = - 1;
  root = newNode(-1, rootEnd);
  activeNode = root; 
  for (i=0; i<size; i++) extendSuffixTree(i);
  int labelHeight = 0;
  setSuffixIndexByDFS(root, labelHeight);
=======
  if(n == NULL)   return 0;
  int i=0;
  int ret = -1;
  // If it is internal node
  if(n->suffixIndex < 0)  for (i = 0; i < MAX_CHAR; i++) if(n->children[i] != NULL) {
    ret = allCommonSubstringsTraversal(n->children[i], labelHeight + edgeLength(n->children[i]));

    if(n->suffixIndex == -1) n->suffixIndex = ret;
    else if((n->suffixIndex == -2 && ret == -3) || (n->suffixIndex == -3 && ret == -2) || n->suffixIndex == -4) {
      n->suffixIndex = -4; // Mark node as XY
      if(labelHeight > 2) if (subtrings[i].stringend == 0 || n->start < subtrings[i].stringstart) {
          subtrings[i].stringstart = n->start;
          subtrings[i].stringend = *(n->end);
        }
    }
  }
  //suffix of X
  else if(n->suffixIndex > -1 && n->suffixIndex < string1Length) return -2;//Mark node as X
  //suffix of Y
  else if(n->suffixIndex >= string1Length) return -3;//Mark node as Y
  return n->suffixIndex;
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
}
//End of Functions for printing Suffix Tree and Labeling Leaf nodes//

<<<<<<< HEAD
// pattern search  version2.0
Node* pickEdge(Node* node, char* p, int pos) {
  char c = p[pos];
  if(node->children[c] != NULL) return node->children[c];
  return NULL;
}

int traverseEdge(Node* node, char* p, int pos) {
  int i, flag=0;
  for(i = 0; i < edgeLength(node) || p[pos] == '\0'; i++, pos++) if(text[(node->start) + i] != p[pos]) {
      flag = -1;
      break;
    }
  if(p[pos] == '\0') return 0;
  else if(flag == -1) return -1;
  return 1;
}

Node* findLocusNode(char* p) {
  Node* u = malloc(sizeof(Node));
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

void subtreeDFS (Node* u) 
{
  if(u != NULL && u->suffixIndex != -1) printf("Index:%d\n", u->suffixIndex);
  if (u == NULL) return;
  int i;
  for(i = 0; i < MAX_CHAR; i++) if(u->children[i] != NULL) {
      if(u->children[i]->suffixIndex == -1) subtreeDFS(u->children[i]);
      else printf("Index:%d\n", u->children[i]->suffixIndex);//Leo:idx of match
    }
}
// end of pattern search v2.0 //

// memory usage tracking function 
int sizeofTree(Node* u) 
{		
  if(u == NULL) return 0;
  if(u->suffixIndex != -1) return sizeof(u);
  int i, size = 0;
  for(i = 0; i < MAX_CHAR; i++) if(u->children[i] != NULL) size += sizeofTree(u->children[i]);
  return size;
=======
char *getLongestCommonSubstring(unsigned char *string1, unsigned char *string2, unsigned char print_tree)
{
  int maxHeight = 0;
  int substringStartIndex = 0;

  //First check if we have two strings
  if((string1 == NULL) || (string2 == NULL)) return NULL;
  buildSuffixTree(string1, string2, print_tree);
  doTraversal(root, 0, &maxHeight, &substringStartIndex);
  static char *longest_substring = NULL;
  char *printer;
  int k;
  buildString(&longest_substring,"");
  for (k=0; k<maxHeight; k++) {
    asprintf(&printer, "%c", text[k + substringStartIndex]);
    buildString(&longest_substring,printer);
    free(printer);
  }
  if(k == 0) buildString(&longest_substring,"No common substring");
  else {
    asprintf(&printer, ", of length: %d",maxHeight);
    buildString(&longest_substring,printer);
    free(printer);
  }
  buildString(&longest_substring,"\n");

  //Free the dynamically allocated memory
  freeSuffixTreeByPostOrder(root);
  return longest_substring;
}

char *getAllCommonSubstrings(unsigned char *string1, unsigned char *string2, unsigned char print_tree){
  if((string1 == NULL) || (string2 == NULL)) return NULL;
  buildSuffixTree(string1, string2, print_tree);

  int i = 0;
  allCommonSubstringsTraversal(root, 0);

  for (i = 0; i < MAX_CHAR; i++) if(subtrings[i].stringend != 0) {
    printf("string start: %d \n",subtrings[i].stringstart);
    printf("string end: %d \n",subtrings[i].stringend);
  }

  return NULL;
}

int traverseEdge(unsigned char *str, int idx, int start, int end) {
  int k = 0;
  //Traverse the edge with character by character matching
  for(k=start; k<=end && str[idx] != '\0'; k++, idx++) if(text[k] != str[idx]) return -1;  // mo match
  if(str[idx] == '\0') return 1;  // match
  return 0;  // more characters yet to match
}

int substringTraversal(Node *n,unsigned char* str, int idx) {
  if(n == NULL) return -1; // no match
  int res = -1;
  //If node n is not root node, then traverse edge from node n's parent to node n.
  if(n->start != -1) {
    res = traverseEdge(str, idx, n->start, *(n->end));
    if(res != 0) return res;  // match (res = 1) or no match (res = -1)
  }
  idx = idx + edgeLength(n);  // Get the character index to search
  //If there is an edge from node n going out with current character str[idx], travrse that edge
  if(n->children[str[idx]] != NULL) return substringTraversal(n->children[str[idx]], str, idx);
  else return -1;  // no match
}

int doTraversalToCountLeaf(Node *n){
  if(n == NULL) return 0;
  if(n->suffixIndex > -1) {
    vector_append(&vector, n->suffixIndex);
    return 1;
  }
  int count = 0;
  int i = 0;
  for (i = 0; i < MAX_CHAR; i++) if(n->children[i] != NULL) count += doTraversalToCountLeaf(n->children[i]);
  return count;
}

int countLeaf(Node *n) {
  if(n == NULL) return 0;
  return doTraversalToCountLeaf(n);
}

int substringAllOccurenceTraversal(Node *n, unsigned char* str, int idx){
  if(n == NULL) return -1; // no match
  int res = -1;
  //If node n is not root node, then traverse edge from node n's parent to node n.
  if(n->start != -1) {
    res = traverseEdge(str, idx, n->start, *(n->end));
    if(res == -1) return -1;//no match
    if(res == 1) {//match
      if(n->suffixIndex > -1) {
        vector_set(&vector, 0, 1);
        vector_append(&vector, n->suffixIndex);
      }
      else {
        int leaf_count = countLeaf(n);
        vector_set(&vector, 0, leaf_count);
      }
      return 1;
    }
  }
  // Get the character index to search
  idx = idx + edgeLength(n);
  // If there is an edge from node n going out with current character str[idx], travrse that edge
  if(n->children[str[idx]] != NULL) return substringAllOccurenceTraversal(n->children[str[idx]], str, idx);
  else return -1;  // no match
}

void longestRepeatedSubstringTraversal(Node *n, int labelHeight, int* maxHeight, int* substringStartIndex) {
  if(n == NULL) return;
  int i=0;
  // If it is internal node
  if(n->suffixIndex == -1) for (i = 0; i < MAX_CHAR; i++) if(n->children[i] != NULL)
    longestRepeatedSubstringTraversal(n->children[i], labelHeight + edgeLength(n->children[i]), maxHeight, substringStartIndex);
  else if (n->suffixIndex > -1 && (*maxHeight < labelHeight - edgeLength(n))) {
    *maxHeight = labelHeight - edgeLength(n);
    *substringStartIndex = n->suffixIndex;
  }
}

int checkForSubString(unsigned char* search_string, unsigned char* source_string) {
  buildSuffixTree(source_string, NULL, 0);
  int res = substringTraversal(root, search_string, 0);
  freeSuffixTreeByPostOrder(root); // Free the dynamically allocated memory
  if(res == 1) return 1;
  else return 0;
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
}
// end of memory usage tracking function //

<<<<<<< HEAD
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
    subtreeDFS(findLocusNode(pattern));	

  }
  freeSuffixTreeByPostOrder(root);
  //checkForSubString("abc");
  return 0;
=======
int *checkAllSubStringOccurences(unsigned char* search_string, unsigned char* source_string){
  vector_init(&vector); //Initiate the dynamic array
  vector_set(&vector, 0, 0); //Make the first element 0

  buildSuffixTree(source_string, NULL, 0);
  int res = substringAllOccurenceTraversal(root, search_string, 0);
  freeSuffixTreeByPostOrder(root);  //Free the dynamically allocated memory

  const int len = vector.size * sizeof(int);
  int * result_array = malloc(len);

  if(res == 1) memcpy(result_array, vector.data, len);
  else result_array[0] = 0;
  vector_free(&vector);
  return result_array;
}

char * getLongestRepeatedSubstring(unsigned char* source_string){
  buildSuffixTree(source_string, NULL, 0);
  int maxHeight = 0;
  int substringStartIndex = 0;
  static char *longest_repeat;
  char *printer;

  longestRepeatedSubstringTraversal(root, 0, &maxHeight, &substringStartIndex);
  freeSuffixTreeByPostOrder(root);

  longest_repeat = malloc(maxHeight * 4);
  longest_repeat[0] = '\0';

  int k;
  for (k=0; k<maxHeight; k++){
    asprintf(&printer, "%c", text[k + substringStartIndex]);
    strcat(longest_repeat, printer);
    free(printer);
  }
  if(k == 0) buildString(&longest_repeat,"No repeated substring");

  return longest_repeat;
}

void buildString(char** current_text, const char *new_text){
  size_t new_len = strlen(new_text) + 1; /* + 1 for terminating NULL */
  if (*current_text == NULL) *current_text = malloc(new_len);
  size_t current_len = strlen(*current_text);
  *current_text = realloc(*current_text, (new_len + current_len));
  strncat(*current_text, new_text, new_len);
}

/*  Vector (Dynamic Array) helper functions  */

void vector_init(Vector *vector) {
  vector->size = 0;  // initialize size and capacity
  vector->capacity = 1;
  vector->data = malloc(sizeof(int) * vector->capacity);// allocate memory for vector->data
}

void vector_append(Vector *vector, int value) {
  vector_double_capacity_if_full(vector);// make sure there's room to expand into
  vector->data[vector->size++] = value;  // append the value and increment vector->size
}

int vector_get(Vector *vector, int index) {
  if (index >= vector->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, vector->size);
    exit(1);
  }
  return vector->data[index];
}

void vector_set(Vector *vector, int index, int value) {
  // zero fill the vector up to the desired index
  while (index >= vector->size) vector_append(vector, 0);
  vector->data[index] = value;// set the value at the desired index
}

void vector_double_capacity_if_full(Vector *vector) {
  if (vector->size >= vector->capacity) {
    // double vector->capacity and resize the allocated memory accordingly
    vector->capacity *= 2;
    vector->data = realloc(vector->data, sizeof(int) * vector->capacity);
  }
}

void vector_free(Vector *vector) {
  free(vector->data);
}

void self_test(){/**  run self tests  */
  printf("Running Self Tests... \n");

  printf("Build suffix tree test: \n");
  char *tree_output =     "$ [3]\n"
    "abc$ [0]\n"
    "bc$ [1]\n"
    "c$ [2]\n"; //Expected output
  printf("Building suffix tree for string: abc \n");
  buildSuffixTree((unsigned char *)"abc", NULL, 1); //Build Suffix tree for this string
  freeSuffixTreeByPostOrder(root);
  printf(tree_string, 's');
  int compare_result = strcmp(tree_string, tree_output);
  assert(compare_result == 0);
  memset(tree_string,0,strlen(tree_string)); //clear string from previous test
  char *tree_output2 =    "$ [10]\n"
    "ab [-1]\n"
    "c [-1]\n"
    "abxabcd$ [0]\n"
    "d$ [6]\n"
    "xabcd$ [3]\n"
    "b [-1]\n"
    "c [-1]\n"
    "abxabcd$ [1]\n"
    "d$ [7]\n"
    "xabcd$ [4]\n"
    "c [-1]\n"
    "abxabcd$ [2]\n"
    "d$ [8]\n"
    "d$ [9]\n"
    "xabcd$ [5]\n"; //Expected output
  printf("Building suffix tree for string: abcabxabcd \n");
  buildSuffixTree((unsigned char *)"abcabxabcd", NULL, 1); //Build Suffix tree for this string
  freeSuffixTreeByPostOrder(root);
  printf(tree_string, 's');
  int compare_result2 = strcmp(tree_string, tree_output2);
  assert(compare_result2 == 0);
  memset(tree_string,0,strlen(tree_string)); //clear string from previous test
  printf("Suffix tree build test: Passed\n\n");

  printf("Longest substrings test: \n");
  char *lcs;

  printf("Longest Common Substring in xabxac and abcabxabcd is: ");
  lcs = getLongestCommonSubstring((unsigned char *)"xabxac",
                                  (unsigned char *)"abcabxabcd",
                                  0);
  printf(lcs, 's');
  int compare_result3 = strcmp(lcs, "abxa, of length: 4\n");
  assert(compare_result3 == 0);
  memset(lcs,0,strlen(lcs)); //clear string from previous test; //clear string from previous test

  printf("Longest Common Substring in xabxaabxa and babxba is: ");
  lcs = getLongestCommonSubstring((unsigned char *)"xabxaabxa",
                                  (unsigned char *)"babxba",
                                  0);

  printf(lcs, 's');
  int compare_result4 = strcmp(lcs, "abx, of length: 3\n");
  assert(compare_result4 == 0);
  memset(lcs,0,strlen(lcs)); //clear string from previous test;
  printf("Longest Common Substring test: Passed\n\n");

  printf("Substrings test: \n");
  int is_substring = 0;

  printf("\"test\" is a substring of \"this is a test\" \n");
  is_substring = checkForSubString((unsigned char *)"test",
                                   (unsigned char *)"this is a test");
  assert(is_substring == 1);

  printf("\"foo\" is a substring of \"this is a test\" \n");
  is_substring = checkForSubString((unsigned char *)"foo",
                                   (unsigned char *)"this is a test");
  assert(is_substring == 0);
  printf("Substrings test: Passed\n\n");

  printf("All occurrences of Substring test: \n");
  int *all_substrings;
  size_t i = 0;

  printf("Text: AABAACAADAABAAABAA, Pattern to search: AABA\n");
  all_substrings = checkAllSubStringOccurences((unsigned char *)"AABA",
                                               (unsigned char *)"AABAACAADAABAAABAA");
  assert(all_substrings[0] == 3);
  assert(all_substrings[1] == 13);
  assert(all_substrings[2] == 9);
  assert(all_substrings[3] == 0);

  printf("Substring found count: %d\n", all_substrings[0]);
  printf("Substrings found at positions: ");
  for( i = 1; i <all_substrings[0] + 1; i++)
   {
    printf ("%d,",all_substrings[i]);
   }
  free(all_substrings);

  printf("\n\nText: AABAACAADAABAAABAA, Pattern to search: AABAACAAD\n");
  all_substrings = checkAllSubStringOccurences((unsigned char *)"AABAACAAD",
                                               (unsigned char *)"AABAACAADAABAAABAA");
  assert(all_substrings[0] == 1);
  assert(all_substrings[1] == 0);

  printf("Substring found count: %d\n", all_substrings[0]);
  printf("Substrings found at position: %d", all_substrings[1]);
  free(all_substrings);

  printf("\n\nText: AABAACAADAABAAABAA, Pattern to search: AA\n");
  all_substrings = checkAllSubStringOccurences((unsigned char *)"AA",
                                               (unsigned char *)"AABAACAADAABAAABAA");
  assert(all_substrings[0] == 7);
  assert(all_substrings[1] == 16);
  assert(all_substrings[2] == 12);
  assert(all_substrings[3] == 13);
  assert(all_substrings[4] == 9);
  assert(all_substrings[5] == 0);
  assert(all_substrings[6] == 3);
  assert(all_substrings[7] == 6);

  printf("Substring found count: %d\n", all_substrings[0]);
  printf("Substrings found at positions: ");
  for( i = 1; i <all_substrings[0] + 1; i++)
   {
    printf ("%d,",all_substrings[i]);
   }
  free(all_substrings);


  printf("\n\nText: AABAACAADAABAAABAA, Pattern to search: ZZ\n");
  all_substrings = checkAllSubStringOccurences((unsigned char *)"ZZ",
                                               (unsigned char *)"AABAACAADAABAAABAA");
  assert(all_substrings[0] == 0);
  printf("No Substrings found\n");
  printf("\nAll occurences of Substring test: Passed\n\n");

  printf("Longest Repeated Substring test: \n");
  char * lrs;
  int lcs_compare;

  printf("Longest Repeated Substring in AAAAAAAAAA is: ");
  lrs = getLongestRepeatedSubstring((unsigned char *)"AAAAAAAAAA");
  lcs_compare = strcmp(lrs, "AAAAAAAAA");
  assert(lcs_compare == 0);
  printf(lrs, 's');
  memset(lrs,0,strlen(lrs)); //clear string from previous test;

  printf("\nLongest Repeated Substring in ABABABA is: ");
  lrs = getLongestRepeatedSubstring((unsigned char *)"ABABABA");
  lcs_compare = strcmp(lrs, "ABABA");
  assert(lcs_compare == 0);
  printf(lrs, 's');
  memset(lrs,0,strlen(lrs)); //clear string from previous test;

  printf("\nLongest Repeated Substring in ABCDEFG is: ");
  lrs = getLongestRepeatedSubstring((unsigned char *)"ABCDEFG");
  lcs_compare = strcmp(lrs, "No repeated substring");
  assert(lcs_compare == 0);
  printf(lrs, 's');
  memset(lrs,0,strlen(lrs)); //clear string from previous test;

  printf("\nLongest Repeated Substring in pqrpqpqabab is: ");
  lrs = getLongestRepeatedSubstring((unsigned char *)"pqrpqpqabab");
  lcs_compare = strcmp(lrs, "ab");
  assert(lcs_compare == 0);
  printf(lrs, 's');
  memset(lrs,0,strlen(lrs)); //clear string from previous test;

  printf("\nLongest Repeated Substring in abcpqrabpqpq is: ");
  lrs = getLongestRepeatedSubstring((unsigned char *)"abcpqrabpqpq");
  lcs_compare = strcmp(lrs, "ab");
  assert(lcs_compare == 0);
  printf(lrs, 's');
  memset(lrs,0,strlen(lrs)); //clear string from previous test;

  printf("\nLongest Repeated Substring test: Passed\n\n");

  printf("done");
>>>>>>> 05a2c908e43f06141fb8df26956e7b77980cd34f
}
