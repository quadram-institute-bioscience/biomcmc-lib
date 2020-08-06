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
 * Original code: https://github.com/mattporritt/ukkonen_suffix_tree ( GPL-3.0-or-later) 
 */ 

#include "suffix_tree.h"


// Define a vector type
typedef struct {
  int size;      // slots used so far
  int capacity;  // total available slots
  int *data;     // array of integers we're storing
} Vector;

unsigned char print_enabled = 0;
unsigned char *text = NULL; //Input string 1 (for suffix tree)
unsigned char *text2 = NULL; //Input string 2 (for generalized suffix tree)
char *tree_string = NULL; //String to hold text version of tree
Node *root = NULL; //Pointer to root node
Node *lastNewNode = NULL;
Node *activeNode = NULL;
Vector vector;

/*activeEdge is represented as input string character
  index (not the character itself)*/
int activeEdge = -1;
int activeLength = 0;

// remainingSuffixCount tells how many suffixes yet to
// be added in tree
int remainingSuffixCount = 0;
int leafEnd = -1;
int *rootEnd = NULL;
int *splitEnd = NULL;
int totalStringLength = -1; //Length of input string
int string1Length = 0; //totalStringLength of 1st string

Node *newNode(int start, int *end)
{
    Node *node =(Node*) malloc(sizeof(Node));
    int i;
    for (i = 0; i < MAX_CHAR; i++)
          node->children[i] = NULL;

    /*For root node, suffixLink will be set to NULL
    For internal nodes, suffixLink will be set to root
    by default in  current extension and may change in
    next extension*/
    node->suffixLink = root;
    node->start = start;
    node->end = end;

    /*suffixIndex will be set to -1 by default and
      actual suffix index will be set later for leaves
      at the end of all phases*/
    node->suffixIndex = -1;
    return node;
}

int edgeLength(Node *n) {
    if(n == root)
        return 0;
    return *(n->end) - (n->start) + 1;
}

int walkDown(Node *currNode)
{
    /*activePoint change for walk down (APCFWD) using
     Skip/Count Trick  (Trick 1). If activeLength is greater
     than current edge length, set next  internal node as
     activeNode and adjust activeEdge and activeLength
     accordingly to represent same activePoint*/
    if (activeLength >= edgeLength(currNode))
    {
        activeEdge += edgeLength(currNode);
        activeLength -= edgeLength(currNode);
        activeNode = currNode;
        return 1;
    }
    return 0;
}

void extendSuffixTree(int pos)
{
    /*Extension Rule 1, this takes care of extending all
    leaves created so far in tree*/
    leafEnd = pos;

    /*Increment remainingSuffixCount indicating that a
    new suffix added to the list of suffixes yet to be
    added in tree*/
    remainingSuffixCount++;

    /*set lastNewNode to NULL while starting a new phase,
     indicating there is no internal node waiting for
     it's suffix link reset in current phase*/
    lastNewNode = NULL;

    //Add all suffixes (yet to be added) one by one in tree
    while(remainingSuffixCount > 0) {

        if (activeLength == 0)
            activeEdge = pos; //APCFALZ

        // There is no outgoing edge starting with
        // activeEdge from activeNode
        if (activeNode->children[text[activeEdge]] == NULL)
        {
            //Extension Rule 2 (A new leaf edge gets created)
            activeNode->children[text[activeEdge]] =
                                          newNode(pos, &leafEnd);

            /*A new leaf edge is created in above line starting
             from  an existng node (the current activeNode), and
             if there is any internal node waiting for it's suffix
             link get reset, point the suffix link from that last
             internal node to current activeNode. Then set lastNewNode
             to NULL indicating no more node waiting for suffix link
             reset.*/
            if (lastNewNode != NULL)
            {
                lastNewNode->suffixLink = activeNode;
                lastNewNode = NULL;
            }
        }
        // There is an outgoing edge starting with activeEdge
        // from activeNode
        else
        {
            // Get the next node at the end of edge starting
            // with activeEdge
            Node *next = activeNode->children[text[activeEdge]];
            if (walkDown(next))//Do walkdown
            {
                //Start from next node (the new activeNode)
                continue;
            }
            /*Extension Rule 3 (current character being processed
              is already on the edge)*/
            if (text[next->start + activeLength] == text[pos])
            {
                //If a newly created node waiting for it's
                //suffix link to be set, then set suffix link
                //of that waiting node to curent active node
                if(lastNewNode != NULL && activeNode != root)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = NULL;
                }

                //APCFER3
                activeLength++;
                /*STOP all further processing in this phase
                and move on to next phase*/
                break;
            }

            /*We will be here when activePoint is in middle of
              the edge being traversed and current character
              being processed is not  on the edge (we fall off
              the tree). In this case, we add a new internal node
              and a new leaf edge going out of that new node. This
              is Extension Rule 2, where a new leaf edge and a new
            internal node get created*/
            splitEnd = (int*) malloc(sizeof(int));
            *splitEnd = next->start + activeLength - 1;

            //New internal node
            Node *split = newNode(next->start, splitEnd);
            activeNode->children[text[activeEdge]] = split;

            //New leaf coming out of new internal node
            split->children[text[pos]] = newNode(pos, &leafEnd);
            next->start += activeLength;
            split->children[text[next->start]] = next;

            /*We got a new internal node here. If there is any
              internal node created in last extensions of same
              phase which is still waiting for it's suffix link
              reset, do it now.*/
            if (lastNewNode != NULL)
            {
            /*suffixLink of lastNewNode points to current newly
              created internal node*/
                lastNewNode->suffixLink = split;
            }

            /*Make the current newly created internal node waiting
              for it's suffix link reset (which is pointing to root
              at present). If we come across any other internal node
              (existing or newly created) in next extension of same
              phase, when a new leaf edge gets added (i.e. when
              Extension Rule 2 applies is any of the next extension
              of same phase) at that point, suffixLink of this node
              will point to that internal node.*/
            lastNewNode = split;
        }

        /* One suffix got added in tree, decrement the count of
          suffixes yet to be added.*/
        remainingSuffixCount--;
        if (activeNode == root && activeLength > 0) //APCFER2C1
        {
            activeLength--;
            activeEdge = pos - remainingSuffixCount + 1;
        }
        else if (activeNode != root) //APCFER2C2
        {
            activeNode = activeNode->suffixLink;
        }
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
    }
}

//Print the suffix tree as well along with setting suffix index
//So tree will be printed in DFS manner
//Each edge along with it's suffix index will be printed
void setSuffixIndexByDFS(Node *n, int labelHeight)
{
    if (n == NULL)  return;

    if (n->start != -1) //A non-root node
    {
            if (print_enabled == 1){
                    //Print the label on edge from parent to current node
                    print(n->start, *(n->end));
            }
    }
    int leaf = 1;
    int i;
    char *printer;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
                if (print_enabled == 1){
                        //print suffix index
                        if (leaf == 1 && n->start != -1){
                                asprintf(&printer, " [%d]\n", n->suffixIndex);
                                buildString(&tree_string, printer);
                                free(printer);
                        }
                }

            //Current node is not a leaf as it has outgoing
            //edges from it.
            leaf = 0;
            setSuffixIndexByDFS(n->children[i], labelHeight +
                                  edgeLength(n->children[i]));
        }
    }
    if (leaf == 1)
    {
        for(i= n->start; i<= *(n->end); i++)
        {
            if(text[i] == '+')
            {
                n->end = (int*) malloc(sizeof(int));
                *(n->end) = i;
            }
        }
        n->suffixIndex = totalStringLength - labelHeight;
        if (print_enabled == 1){
                //print suffix index
                asprintf(&printer, " [%d]\n", n->suffixIndex);
                buildString(&tree_string,printer);
                free(printer);
        }
    }
}

void freeSuffixTreeByPostOrder(Node *n)
{
    if (n == NULL)
        return;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
            freeSuffixTreeByPostOrder(n->children[i]);
        }
    }
    if (n->suffixIndex == -1)
        free(n->end);
    free(n);
}

/*Build the suffix tree and print the edge labels along with
suffixIndex. suffixIndex for leaf edges will be >= 0 and
for non-leaf edges will be -1*/
int buildSuffixTree(unsigned char *string1, unsigned char *string2, unsigned char print_tree){
        //First check if we have at least one valid string
        if(string1 == NULL){
                return 1;
        }
        // Check if we are printing suffix tree
        print_enabled = print_tree;

        //Check if we have two strings
        //Two strings mean we are building a generalized suffix tree
        if(string2 != NULL){
                asprintf(&text, "%s%c%s%c", string1, 43,  string2, 36);
                string1Length = strlen((const char*) string1);
        }
        else{
                asprintf(&text, "%s%c", string1, 36);
        }

        totalStringLength = strlen((const char*) text);
    int i;
    rootEnd = (int*) malloc(sizeof(int));
    *rootEnd = - 1;

    /*Root is a special node with start and end indices as -1,
    as it has no parent from where an edge comes to root*/
    root = newNode(-1, rootEnd);

    activeNode = root; //First activeNode will be root
    for (i=0; i<totalStringLength; i++)
        extendSuffixTree(i);
    int labelHeight = 0;
    setSuffixIndexByDFS(root, labelHeight);

    return 0;
}

int doTraversal(Node *n, int labelHeight, int* maxHeight,
int* substringStartIndex)
{
    if(n == NULL)
    {
        return 0;
    }
    int i=0;
    int ret = -1;
    if(n->suffixIndex < 0) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != NULL)
            {
                ret = doTraversal(n->children[i], labelHeight +
                    edgeLength(n->children[i]),
                    maxHeight, substringStartIndex);

                if(n->suffixIndex == -1)
                    n->suffixIndex = ret;
                else if((n->suffixIndex == -2 && ret == -3) ||
                    (n->suffixIndex == -3 && ret == -2) ||
                    n->suffixIndex == -4)
                {
                    n->suffixIndex = -4; //Mark node as XY

                    if(*maxHeight < labelHeight) //Keep track of deepest node
                    {
                        *maxHeight = labelHeight;
                        *substringStartIndex = *(n->end) -
                            labelHeight + 1;
                    }
                }
            }
        }
    }
    else if(n->suffixIndex > -1 && n->suffixIndex < string1Length)//suffix of X
        return -2;//Mark node as X
    else if(n->suffixIndex >= string1Length)//suffix of Y
        return -3;//Mark node as Y
    return n->suffixIndex;
}

int allCommonSubstringsTraversal(Node *n, int labelHeight)
{
    if(n == NULL)
    {
        return 0;
    }
    int i=0;
    int ret = -1;
    if(n->suffixIndex < 0) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != NULL)
            {
                ret = allCommonSubstringsTraversal(n->children[i], labelHeight + edgeLength(n->children[i]));

                if(n->suffixIndex == -1)
                    n->suffixIndex = ret;
                else if((n->suffixIndex == -2 && ret == -3) ||
                    (n->suffixIndex == -3 && ret == -2) ||
                    n->suffixIndex == -4)
                {
                    n->suffixIndex = -4;//Mark node as XY
                    if(labelHeight > 2){
                            if (subtrings[i].stringend == 0 || n->start < subtrings[i].stringstart){
                                    subtrings[i].stringstart = n->start;
                                    subtrings[i].stringend = *(n->end);
                            }
                    }
                }
            }
        }
    }
    else if(n->suffixIndex > -1 && n->suffixIndex < string1Length)//suffix of X
        return -2;//Mark node as X
    else if(n->suffixIndex >= string1Length)//suffix of Y
        return -3;//Mark node as Y
    return n->suffixIndex;
}

char *getLongestCommonSubstring(unsigned char *string1, unsigned char *string2, unsigned char print_tree)
{
    int maxHeight = 0;
    int substringStartIndex = 0;


    //First check if we have two strings
    if((string1 == NULL) || (string2 == NULL)){
            return NULL;
    }
    buildSuffixTree(string1, string2, print_tree);
    doTraversal(root, 0, &maxHeight, &substringStartIndex);
    static char *longest_substring = NULL;
    char *printer;
    int k;
    buildString(&longest_substring,"");
    for (k=0; k<maxHeight; k++){
            asprintf(&printer, "%c", text[k + substringStartIndex]);
            buildString(&longest_substring,printer);
            free(printer);

    }
    if(k == 0){
            buildString(&longest_substring,"No common substring");
    }
    else{
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
        //First check if we have two strings
        if((string1 == NULL) || (string2 == NULL)){
                return NULL;
        }
        buildSuffixTree(string1, string2, print_tree);

        int i = 0;
        allCommonSubstringsTraversal(root, 0);

        for (i = 0; i < MAX_CHAR; i++){
                if(subtrings[i].stringend != 0){
                        printf("string start: %d \n",subtrings[i].stringstart);
                        printf("string end: %d \n",subtrings[i].stringend);
                }
        }

        return NULL;
}

int traverseEdge(unsigned char *str, int idx, int start, int end){
    int k = 0;
    //Traverse the edge with character by character matching
    for(k=start; k<=end && str[idx] != '\0'; k++, idx++)
    {
        if(text[k] != str[idx])
            return -1;  // mo match
    }
    if(str[idx] == '\0')
        return 1;  // match
    return 0;  // more characters yet to match
}

int substringTraversal(Node *n,unsigned char* str, int idx){
    if(n == NULL)
    {
        return -1; // no match
    }
    int res = -1;
    //If node n is not root node, then traverse edge
    //from node n's parent to node n.
    if(n->start != -1)
    {
        res = traverseEdge(str, idx, n->start, *(n->end));
        if(res != 0)
            return res;  // match (res = 1) or no match (res = -1)
    }
    //Get the character index to search
    idx = idx + edgeLength(n);
    //If there is an edge from node n going out
    //with current character str[idx], travrse that edge
    if(n->children[str[idx]] != NULL)
        return substringTraversal(n->children[str[idx]], str, idx);
    else
        return -1;  // no match
}

int doTraversalToCountLeaf(Node *n){
    if(n == NULL)
        return 0;
    if(n->suffixIndex > -1)
    {
        vector_append(&vector, n->suffixIndex);
        return 1;
    }
    int count = 0;
    int i = 0;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if(n->children[i] != NULL)
        {
            count += doTraversalToCountLeaf(n->children[i]);
        }
    }
    return count;
}

int countLeaf(Node *n){
    if(n == NULL)
        return 0;
    return doTraversalToCountLeaf(n);
}

int substringAllOccurenceTraversal(Node *n, unsigned char* str, int idx){
    if(n == NULL)
    {
        return -1; // no match
    }
    int res = -1;
    //If node n is not root node, then traverse edge
    //from node n's parent to node n.
    if(n->start != -1)
    {
        res = traverseEdge(str, idx, n->start, *(n->end));
        if(res == -1)  //no match
            return -1;
        if(res == 1) //match
        {
            if(n->suffixIndex > -1){
                vector_set(&vector, 0, 1);
                vector_append(&vector, n->suffixIndex);
            }
            else{
                int leaf_count = countLeaf(n);
                vector_set(&vector, 0, leaf_count);
            }
            return 1;
        }
    }
    //Get the character index to search
    idx = idx + edgeLength(n);
    //If there is an edge from node n going out
    //with current character str[idx], travrse that edge
    if(n->children[str[idx]] != NULL)
        return substringAllOccurenceTraversal(n->children[str[idx]], str, idx);
    else
        return -1;  // no match
}

void longestRepeatedSubstringTraversal(Node *n, int labelHeight,
                                       int* maxHeight,
                                       int* substringStartIndex){
    if(n == NULL)
    {
        return;
    }
    int i=0;
    if(n->suffixIndex == -1) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != NULL)
            {
                    longestRepeatedSubstringTraversal(n->children[i], labelHeight +
                                edgeLength(n->children[i]), maxHeight,
                                 substringStartIndex);
            }
        }
    }
    else if(n->suffixIndex > -1 &&
                (*maxHeight < labelHeight - edgeLength(n)))
    {
        *maxHeight = labelHeight - edgeLength(n);
        *substringStartIndex = n->suffixIndex;
    }
}

int checkForSubString(unsigned char* search_string, unsigned char* source_string){
        buildSuffixTree(source_string, NULL, 0);
        int res = substringTraversal(root, search_string, 0);
        //Free the dynamically allocated memory
        freeSuffixTreeByPostOrder(root);

        if(res == 1)
                return 1;
        else
                return 0;
}

int *checkAllSubStringOccurences(unsigned char* search_string, unsigned char* source_string){
        vector_init(&vector); //Initiate the dynamic array
        vector_set(&vector, 0, 0); //Make the first element 0

        buildSuffixTree(source_string, NULL, 0);
        int res = substringAllOccurenceTraversal(root, search_string, 0);
        //Free the dynamically allocated memory
        freeSuffixTreeByPostOrder(root);

        const int len = vector.size * sizeof(int);
        int * result_array = malloc(len);

        if(res == 1){
                memcpy(result_array, vector.data, len);
        }
        else{
                result_array[0] = 0;
        }
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
        if(k == 0){
                buildString(&longest_repeat,"No repeated substring");
        }

        return longest_repeat;
}

void buildString(char** current_text, const char *new_text){
        size_t new_len = strlen(new_text) + 1; /* + 1 for terminating NULL */
        if (*current_text == NULL){
                *current_text = malloc(new_len);
        }
        size_t current_len = strlen(*current_text);
        *current_text = realloc(*current_text, (new_len + current_len));
        strncat(*current_text, new_text, new_len);
}

/*
 * Vector (Dynamic Array) helper functions
 */

void vector_init(Vector *vector) {
  // initialize size and capacity
  vector->size = 0;
  vector->capacity = 1;

  // allocate memory for vector->data
  vector->data = malloc(sizeof(int) * vector->capacity);
}

void vector_append(Vector *vector, int value) {
  // make sure there's room to expand into
  vector_double_capacity_if_full(vector);

  // append the value and increment vector->size
  vector->data[vector->size++] = value;
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
  while (index >= vector->size) {
    vector_append(vector, 0);
  }

  // set the value at the desired index
  vector->data[index] = value;
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


//// start of main.c
/*
 * Command line interface for suffix tree functions
 */
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <argp.h>
#include "suffix_tree.h"

void self_test();
unsigned char *input_text; //Input string 1 (for suffix tree)
unsigned char *input_text2; //Input string 2 (for generalized suffix tree)
int opt_print_tree = 0;
int opt_lcs = 0;
int opt_acs = 0;

// Process command line options and arguments
static int
parse_opt (int key, char *arg, struct argp_state *state)
{
        int *arg_count = state->input;

        switch (key)
        {
                case 't': //Run the self test suite
                {
                        self_test();
                        exit(0);
                }
                case 'p': //Print the generated tree
                {
                        opt_print_tree = 1;
                        break;
                }
                case 'l': //Find Longest Common Substring (LCS)
                {
                        opt_lcs = 1;
                        break;
                }
                case 'a': //Find Longest Common Substring (LCS)
                {
                        opt_acs = 1;
                        break;
                }
                case ARGP_KEY_ARG: //Process the command line arguments
                {
                        (*arg_count)--;
                        if (*arg_count == 3){
                                input_text = (unsigned char*)arg;
                        }
                        else if (*arg_count == 2){
                                input_text2 = (unsigned char*)arg;
                        }

                }
                break;
                case ARGP_KEY_END:
                {
                        printf ("\n");
                        if (*arg_count >= 4){
                                argp_failure (state, 1, 0, "too few arguments");
                        }
                        else if (*arg_count < 0){
                                argp_failure (state, 1, 0, "too many arguments");
                        }
                        else {
                                if (opt_print_tree){
                                        // Construct the tree and process based on supplied options
                                        buildSuffixTree(input_text, input_text2, opt_print_tree);
                                        freeSuffixTreeByPostOrder(root);
                                        printf(tree_string, 's');
                                }
                                if (opt_lcs){
                                        char *lcs;
                                        if(!input_text2){
                                                argp_failure (state, 1, 0, "missing comparison string");
                                        }
                                        printf("Longest Common Substring in %s and %s is: ",
                                               input_text,
                                               input_text2);
                                        lcs = getLongestCommonSubstring(input_text,
                                                                        input_text2,
                                                                        opt_print_tree);

                                        printf(lcs, 's');

                                }
                                if (opt_acs){
                                        char *acs;
                                        if(!input_text2){
                                                argp_failure (state, 1, 0, "missing comparison string");
                                        }
                                        printf("All Common Substrings in %s and %s are: ",
                                               input_text,
                                               input_text2);
                                        acs = getAllCommonSubstrings(input_text,
                                                                     input_text2,
                                                                     opt_print_tree);

                                        printf(acs, 's');

                                }
                        }
                }
                break;
        }
        return 0;
}

/**
 * run self tests
 */
void self_test(){
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

//        printf("All Common substrings test: \n");
//
//        printf("All Common Substrings in orangeisatypeoffruit and fruitsomestugfruitgoeshereorange are: \n");
//        getAllCommonSubstrings((unsigned char *)"orangeisatypeoffruit",
//                                        (unsigned char *)"fruitsomestugfruitgoeshereorange",
//                                        0);

        printf("done");

}

int main(int argc, char **argv)
{
        struct argp_option options[] =
        {
                { "print", 'p', 0, 0, "Print the tree for the supplied string"},
                { "test", 't', 0, 0, "Run self tests"},
                { "lcs", 'l', 0, 0, "Find longest common substring"},
                { "acs", 'a', 0, 0, "Find all common substring"},
                { 0 }
        };

        int arg_count = 4;
        struct argp argp = { options, parse_opt, "STRING [STRING]" };
        return argp_parse (&argp, argc, argv, 0, 0, &arg_count);

//        strcpy(text, "abcabxabcd$"); buildSuffixTree(); freeSuffixTreeByPostOrder(root);
//        size1 = 21;
//        printf("All Common Substrings in orangeisatypeoffruit and fruitsomestugfruitgoeshereorange are: \n");
//        strcpy((char*)text, "orangeisatypeoffruit#fruitsomestugfruitgoeshereorange$");
//        buildSuffixTree();
//        getLongestCommonSubstring();
//        //Free the dynamically allocated memory
//        freeSuffixTreeByPostOrder(root);

//    return 0;
}
