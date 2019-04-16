/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.c 
 *  \brief Lowest level basic functions, that should be available to all other modules. 
 */

#include "lowlevel.h"
void  hungarian_solve_integer (hungarian p, int this_size);
void  hungarian_solve_double (hungarian p, int this_size);

/* error-safe memory allocation functions */
void *
biomcmc_malloc (size_t size)
{
  void *value = malloc (size);
  if (value == NULL) biomcmc_error ( "biomcmc_malloc error allocating %d bites", size);
  return value;
}

void *
biomcmc_realloc (void *ptr, size_t size)
{
  void *value = (void *) realloc ((void *)ptr, size);
  if (value == NULL) biomcmc_error ( "biomcmc_realloc error on pointer 0x%08X of %d bites\n", ptr, size);
  return value;
}

FILE *
biomcmc_fopen (const char *path, const char *mode)
{
  FILE *fp = fopen (path, mode);
  if (fp == NULL) {
    fprintf (stderr, "Please check if path is correct, if there are non-ASCII characters in file name,\n");
    fprintf (stderr, "if you have enough permissions (to read/write). Remember that paths are relative to\n");
    fprintf (stderr, "where this program is being called\n");
    biomcmc_error ( "problem opening file \"%s\" with mode \"%s\"", path, mode);
  }
  return fp;
}

void
biomcmc_error (const char *template, ...)
{
  va_list ap;

  fprintf (stderr, "%s error: ", PACKAGE_STRING);
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
  fprintf (stderr, "\n");
  fprintf (stderr, "[note to developers] If you want to debug, set a breakpoint on function biomcmc_error()\n");
  fflush (stderr);
  exit (EXIT_FAILURE);
}

int
compare_int_increasing (const void *a, const void *b)
{
  return (*(int *) a - *(int *) b);
}

int
compare_int_decreasing (const void *a, const void *b)
{
  return (*(int *) b - *(int *) a);
}

int
compare_double_increasing (const void *a, const void *b)
{
  if ((double *) a > (double *) b) return 1;
  if ((double *) a < (double *) b) return -1;
  return 0;
}

int
compare_double_decreasing (const void *a, const void *b)
{
  if ((double *) a < (double *) b) return 1;
  if ((double *) a > (double *) b) return -1;
  return 0;
}

/* \brief size, in bytes, when extending the buffer of biomcmc_getline() */
#define MIN_CHUNK 64
int
biomcmc_getline (char **lineptr, size_t *n, FILE *stream)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline()");
  if (!stream)  biomcmc_error ("lack of input file in biomcmc_getline()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = getc (stream);

    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the 
     * first char) NUL-terminate the line buffer.  
     */
    if ((*lineptr + *n) != (read_pos + nchars_avail)) 
      biomcmc_error ("problem setting string size in biomcmc_getline()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) 
        biomcmc_error ("problem setting string size in biomcmc_getline()");
    }

    if (ferror (stream)) return -1;

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}

/* The hungarian method below is copied from http://www.informatik.uni-freiburg.de/~stachnis/misc.html
 * The (edited) original message follows:
 *
 ** libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the 
 ** Hungarian Method.         ** This file may be freely copied and distributed! **
 **
 ** Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 ** As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 ** "Stanford GraphGase" distrubition! */

void
hungarian_reset (hungarian p)
{
  int i, j;

  if (p->is_double) for (i = 0; i < p->size; i++) {
    p->col_mate[i] = p->unchosen_row[i] = p->slack_row[i] = p->row_mate[i] = p->parent_row[i] = 0;
    p->row_dec_d[i] = p->col_inc_d[i] = p->slack_d[i] = 0.;
    for (j = 0; j < p->size; j++) p->dcost[i][j] = 0;
  }
  else for (i = 0; i < p->size; i++) {
    p->col_mate[i] = p->unchosen_row[i] = p->row_dec[i] = p->slack_row[i] = p->row_mate[i] = p->parent_row[i] = p->col_inc[i] = p->slack[i] = 0;
    for (j = 0; j < p->size; j++) p->cost[i][j] = 0;
  }
  p->final_cost = 0; p->final_dcost = 0.;
}

hungarian
new_hungarian (int size, bool is_double)
{
  int i;
  hungarian p;

  p = (hungarian) biomcmc_malloc (sizeof (struct hungarian_struct)); 
  p->size = size; /* n_rows = n_columns; if it's not, fill with zeroes (no cost) */
  p->is_double = is_double;
  if (is_double) {
    p->dcost = (double**) biomcmc_malloc (size * sizeof (double*));
    for (i = 0; i < p->size; i++) p->dcost[i] = (double*) biomcmc_malloc (size * sizeof (double));
    p->cost = NULL;
    p->row_dec = p->col_inc = p->slack = NULL;
    p->row_dec_d = (double*) biomcmc_malloc (size * sizeof (double));
    p->col_inc_d = (double*) biomcmc_malloc (size * sizeof (double));
    p->slack_d   = (double*) biomcmc_malloc (size * sizeof (double));
  }
  else {
    p->cost = (int**) biomcmc_malloc (size * sizeof (int*));
    for (i = 0; i < p->size; i++) p->cost[i] = (int*) biomcmc_malloc (size * sizeof (int));
    p->dcost = NULL;
    p->row_dec_d = p->col_inc_d = p->slack_d = NULL;
    p->row_dec = (int*) biomcmc_malloc (size * sizeof (int));
    p->col_inc = (int*) biomcmc_malloc (size * sizeof (int));
    p->slack   = (int*) biomcmc_malloc (size * sizeof (int));
  }
  /* edges would be assignment_matrix[ i * ncols + col_mate[i] ] = true; and other elems "false" (but we don't use the matrix notation) */
  p->col_mate     = (int*) biomcmc_malloc (size * sizeof (int)); /* for a given row node, col_mate[row] is the assigned col node */
  p->unchosen_row = (int*) biomcmc_malloc (size * sizeof (int));
  p->slack_row    = (int*) biomcmc_malloc (size * sizeof (int));
  p->row_mate     = (int*) biomcmc_malloc (size * sizeof (int));
  p->parent_row   = (int*) biomcmc_malloc (size * sizeof (int));

  hungarian_reset (p);
  return p;
}

void
hungarian_update_cost (hungarian p, int row, int col, void *cost)
{
  if (row >= p->size) return;
  if (col >= p->size) return;
  if (p->is_double) p->dcost[row][col] = (double)(*(double*)cost);
  else p->cost[row][col] =  (int)(*(int*)cost);
}

void 
del_hungarian (hungarian p)
{
  int i;
  if (!p) return;
  if (p->cost) {
    for (i = p->size - 1; i >= 0; i--) if (p->cost[i]) free (p->cost[i]);
    free (p->cost);
  }
  if (p->dcost) {
    for (i = p->size - 1; i >= 0; i--) if (p->dcost[i]) free (p->dcost[i]);
    free (p->dcost);
  }
  if (p->col_mate) free (p->col_mate); /* this is the important one, with i assigned to col_mate[i] */
  if (p->slack)      free (p->slack);
  if (p->slack_d)    free (p->slack);
  if (p->col_inc)    free (p->col_inc);
  if (p->col_inc_d)  free (p->col_inc);
  if (p->row_dec)    free (p->row_dec);
  if (p->row_dec_d)  free (p->row_dec);
  if (p->parent_row) free (p->parent_row);
  if (p->row_mate)   free (p->row_mate);
  if (p->slack_row)  free (p->slack_row);
  if (p->unchosen_row) free (p->unchosen_row);
  free (p);
}

void 
hungarian_solve (hungarian p, int this_size)
{
  if (p->is_double) hungarian_solve_double (p, this_size);
  else              hungarian_solve_integer (p, this_size);
}

void 
hungarian_solve_integer (hungarian p, int this_size)
{
  int i, j, nrows = this_size, ncols = this_size, k, l, s, t, q, unmatched;
  p->final_cost = p->initial_cost = 0;

  if (this_size > p->size) { p->final_cost = -1; return; } /* we don't call biomcmc_error(), but it *is* an error! */

  for (l = 0; l < ncols; l++) { // Begin subtract column minima in order to start with lots of zeroes 12
    s = p->cost[0][l];
    for (k = 1; k < nrows; k++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->initial_cost += s; /* this should be added to final_cost to have classical assignment cost; here we distinguish them */
    if (s!=0) for (k = 0; k < nrows; k++) p->cost[k][l] -= s;
  } // End subtract column minima in order to start with lots of zeroes 12
 
  /*for (i=0;i<nrows; i++) {
    for (j=0;j<ncols; j++) printf ("%4d ", p->cost[i][j]);
    printf (" DEBUGcost\n");
  }  // DEBUG */ 

  // Begin initial state 16
  t=0;
  for (l = 0; l < ncols; l++)  { // n => num_cols
    p->row_mate[l]= -1;
    p->parent_row[l]= -1;
    p->col_inc[l]=0;
    p->slack[l]= 0x7FFFFFFF;
  }
  for (k = 0; k < nrows; k++) { // m => num_rows
    s = p->cost[k][0];
    for (l = 1; l < ncols; l++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->row_dec[k]=s;
    for (l = 0; l < ncols; l++) if ((s==p->cost[k][l]) && (p->row_mate[l] < 0)) {
      p->col_mate[k] = l;
      p->row_mate[l] = k;  // fprintf(stderr, "matching col %d==row %d\n",l,k);
      goto row_done;
    }
    p->col_mate[k] = -1;  // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    p->unchosen_row[t++] = k;
row_done:
    ;
  }
  // End initial state 16

  // Begin Hungarian algorithm 18
  if (t==0)    goto done;
  unmatched=t;
  while (1) {
    q=0; // fprintf(stderr, "Matched %d rows.\n",m-t);
    while (1)	{
      while (q<t) {
         { // Begin explore node q of the forest 19
          k = p->unchosen_row[q];
          s=p->row_dec[k];
          for (l=0;l<ncols;l++) if (p->slack[l]) {
            int del;
            del = p->cost[k][l] - s + p->col_inc[l];
            if (del < p->slack[l]) {
              if (del==0) {
                if (p->row_mate[l]<0)  goto breakthru;
                p->slack[l]=0;
                p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n", t,row_mate[l],l,k);
                p->unchosen_row[t++]=p->row_mate[l];
              }
              else { p->slack[l]=del; p->slack_row[l]=k; }
            }
          }
         } // End explore node q of the forest 19
        q++;
      }

      // Begin introduce a new zero into the matrix 21
      s = 0x7FFFFFFF;
      for (l = 0;l < ncols; l++) if (p->slack[l] && p->slack[l] < s) s = p->slack[l];
      for (q = 0; q < t; q++) p->row_dec[ p->unchosen_row[q] ] += s;
      for (l = 0; l < ncols; l++) if (p->slack[l]) {
        p->slack[l]-=s;
        if (p->slack[l]==0) {  // Begin look at a new zero 22
          k = p->slack_row[l]; // fprintf(stderr, "Decreasing uncovered elements by %d produces zero at [%d,%d]\n", s,k,l);
          if (p->row_mate[l]<0)  {
            for (j=l+1;j<ncols;j++)  if (p->slack[j]==0) p->col_inc[j]+=s;
            goto breakthru;
          }
          else {
            p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
            p->unchosen_row[t++]=p->row_mate[l];
          }
        } // End look at a new zero 22

      }
      else  p->col_inc[l]+=s;
      // End introduce a new zero into the matrix 21
    }
breakthru:
    // fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
    while (1)	{    // Begin update the matching 20
      j=p->col_mate[k];
      p->col_mate[k]=l;
      p->row_mate[l]=k; // fprintf(stderr, "rematching col %d==row %d\n",l,k);
      if (j<0)   break;
      k=p->parent_row[j];
      l=j;
    }    // End update the matching 20
    if (--unmatched==0)	goto done;
    // Begin get ready for another stage 17
    t=0;
    for (l=0;l<ncols;l++) {
      p->parent_row[l]= -1;
      p->slack[l]=0x7FFFFFFF;
    }
    for (k=0;k<nrows;k++) if (p->col_mate[k]<0) p->unchosen_row[t++]=k; // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    // End get ready for another stage 17
  }
done:
  // Begin doublecheck the solution 23
  for (k = 0; k < nrows; k++) for (l=0;l<ncols;l++) if (p->cost[k][l] < p->row_dec[k] - p->col_inc[l]) { p->final_cost = -1; printf ("\n**\n"); return;}
  for (k = 0; k < nrows; k++) {
    l=p->col_mate[k];
    if ((l < 0) || (p->cost[k][l] != p->row_dec[k] - p->col_inc[l])) { p->final_cost = -1; return; }
  }
  k=0;
  for (l=0;l<ncols;l++) if (p->col_inc[l])  k++;
  if (k>nrows) { p->final_cost = -1; return; }
  // End doublecheck the solution 23
  // End Hungarian algorithm 18

  for (k = 0; k < nrows; ++k) for (l = 0; l < ncols; ++l) p->cost[k][l] = p->cost[k][l] - p->row_dec[k] + p->col_inc[l];
  for (i = 0; i < nrows; i++) p->final_cost += p->row_dec[i];
  for (i = 0; i < ncols; i++) p->final_cost -= p->col_inc[i]; // fprintf(stderr, "Cost is %d\n",cost);
}

void 
hungarian_solve_double (hungarian p, int this_size)
{
  int i, j, nrows = this_size, ncols = this_size, k, l, t, q, unmatched;
  double s;
  p->final_cost = p->initial_cost = 0;

  if (this_size > p->size) { p->final_dcost = -1; return; } /* we don't call biomcmc_error(), but it *is* an error! */

  for (l = 0; l < ncols; l++) { // <double> Begin subtract column minima in order to start with lots of zeroes 12
    s = p->dcost[0][l];
    for (k = 1; k < nrows; k++) if (p->dcost[k][l] < s) s = p->dcost[k][l];
    p->initial_dcost += s; /* this should be added to final_cost to have classical assignment cost; here we distinguish them */
    if (s > DBL_MIN) for (k = 0; k < nrows; k++) p->dcost[k][l] -= s;
  } // End subtract column minima in order to start with lots of zeroes 12
 
  /*for (i=0;i<nrows; i++) {
    for (j=0;j<ncols; j++) printf ("%4d ", p->cost[i][j]);
    printf (" DEBUGcost\n");
  }  // DEBUG */ 

  // Begin initial state 16
  t=0;
  for (l = 0; l < ncols; l++)  { // n => num_cols
    p->row_mate[l]= -1;
    p->parent_row[l]= -1;
    p->col_inc_d[l]=0.;
    p->slack_d[l]= DBL_MAX;
  }
  for (k = 0; k < nrows; k++) { // m => num_rows
    s = p->dcost[k][0];
    for (l = 1; l < ncols; l++) if (p->dcost[k][l] < s) s = p->dcost[k][l];
    p->row_dec_d[k]=s;
    for (l = 0; l < ncols; l++) if ((s==p->dcost[k][l]) && (p->row_mate[l] < 0)) {
      p->col_mate[k] = l;
      p->row_mate[l] = k;  // fprintf(stderr, "matching col %d==row %d\n",l,k);
      goto row_done;
    }
    p->col_mate[k] = -1;  // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    p->unchosen_row[t++] = k;
row_done:
    ;
  }
  // End initial state 16

  // Begin Hungarian algorithm 18
  if (t==0)    goto done;
  unmatched=t;
  while (1) {
    q=0; // fprintf(stderr, "Matched %d rows.\n",m-t);
    while (1)	{
      while (q<t) {
         { // Begin explore node q of the forest 19
          k = p->unchosen_row[q];
          s=p->row_dec_d[k];
          for (l=0;l<ncols;l++) if (p->slack_d[l] > DBL_MIN) {
            double del = p->dcost[k][l] - s + p->col_inc_d[l];
            if (del < p->slack_d[l]) {
              if (del < DBL_MIN) {
                if (p->row_mate[l]<0)  goto breakthru;
                p->slack_d[l]=del;
                p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n", t,row_mate[l],l,k);
                p->unchosen_row[t++]=p->row_mate[l];
              }
              else { p->slack_d[l]=del; p->slack_row[l]=k; }
            }
          }
         } // End explore node q of the forest 19
        q++;
      }

      // Begin introduce a new zero into the matrix 21
      s = DBL_MAX;
      for (l = 0;l < ncols; l++) if ((p->slack_d[l] > DBL_MIN) && p->slack_d[l] < s) s = p->slack_d[l];
      for (q = 0; q < t; q++) p->row_dec_d[ p->unchosen_row[q] ] += s;
      for (l = 0; l < ncols; l++) if (p->slack_d[l] > DBL_MIN) {
        p->slack_d[l]-=s;
        if (p->slack_d[l]==0) {  // Begin look at a new zero 22
          k = p->slack_row[l]; // fprintf(stderr, "Decreasing uncovered elements by %d produces zero at [%d,%d]\n", s,k,l);
          if (p->row_mate[l]<0)  {
            for (j=l+1;j<ncols;j++)  if (p->slack_d[j] < DBL_MIN) p->col_inc_d[j] += s;
            goto breakthru;
          }
          else {
            p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
            p->unchosen_row[t++]=p->row_mate[l];
          }
        } // End look at a new zero 22

      }
      else  p->col_inc_d[l]+=s;
      // End introduce a new zero into the matrix 21
    }
breakthru:
    // fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
    while (1)	{    // Begin update the matching 20
      j=p->col_mate[k];
      p->col_mate[k]=l;
      p->row_mate[l]=k; // fprintf(stderr, "rematching col %d==row %d\n",l,k);
      if (j<0)   break;
      k=p->parent_row[j];
      l=j;
    }    // End update the matching 20
    if (--unmatched==0)	goto done;
    // Begin get ready for another stage 17
    t=0;
    for (l=0;l<ncols;l++) {
      p->parent_row[l]= -1;
      p->slack_d[l]=DBL_MAX;
    }
    for (k=0;k<nrows;k++) if (p->col_mate[k]<0) p->unchosen_row[t++]=k; // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    // End get ready for another stage 17
  }
done:
  // Begin doublecheck the solution 23
  for (k = 0; k < nrows; k++) for (l=0;l<ncols;l++) if (p->dcost[k][l] < p->row_dec_d[k] - p->col_inc_d[l]) { p->final_dcost = -1; printf ("\n**\n"); return;}
  for (k = 0; k < nrows; k++) {
    l=p->col_mate[k];
    if ((l < 0) || fabs(p->dcost[k][l] - p->row_dec_d[k] + p->col_inc_d[l]) > 2. * DBL_MIN) { p->final_cost = -1; return; }
  }
  k=0;
  for (l=0;l<ncols;l++) if (p->col_inc_d[l] > DBL_MIN)  k++;
  if (k>nrows) { p->final_cost = -1; return; }
  // End doublecheck the solution 23
  // End Hungarian algorithm 18

  for (k = 0; k < nrows; ++k) for (l = 0; l < ncols; ++l) p->dcost[k][l] = p->dcost[k][l] - p->row_dec_d[k] + p->col_inc_d[l];
  for (i = 0; i < nrows; i++) p->final_dcost += p->row_dec_d[i];
  for (i = 0; i < ncols; i++) p->final_dcost -= p->col_inc_d[i]; 
}


