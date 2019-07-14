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

/*  Quickselect routine based on algorithm described in "Numerical recipes in C" Section 8.5 (ISBN 0-521-43108-5)
 *  This code by Nicolas Devillard - 1998. Public domain.  http://ndevilla.free.fr/median
 */
#include "quickselect_quantile.h"

#define ELEM_SWAP_DOUBLE(a,b) { register double t=(a);(a)=(b);(b)=t; }

double 
biomcmc_quantile_double (double *original_vector, int n, double quantile) 
{
  int low, high, q, middle, ll, hh;
  double *v = (double*) biomcmc_malloc (n * sizeof (double));
  double result;
  low = 0 ; high = n-1 ; q = (int) ((double)(n) * quantile);
  if (quantile <=0) q = 0;
  if (quantile > n-1) q = n-1;
  for (hh = 0; hh < n; hh++) v[hh] = original_vector[hh];
  for (;;) {
    if (high <= low) {result = v[q]; free (v); return result; }  /* One element only */
    if (high == low + 1) {  /* Two elements only */
      if (v[low] > v[high]) ELEM_SWAP_DOUBLE(v[low], v[high]);
      result = v[q]; free (v); return result;
    }
    /* Find quantile of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (v[middle] > v[high]) ELEM_SWAP_DOUBLE(v[middle], v[high]);
    if (v[low] > v[high])    ELEM_SWAP_DOUBLE(v[low], v[high]);
    if (v[middle] > v[low])  ELEM_SWAP_DOUBLE(v[middle], v[low]);
    ELEM_SWAP_DOUBLE(v[middle], v[low+1]); /* Swap low item (now in position middle) into position (low+1) */

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1; hh = high;
    for (;;) {
      do ll++; while (v[low] > v[ll]) ;
      do hh--; while (v[hh]  > v[low]) ;
      if (hh < ll) break;
      ELEM_SWAP_DOUBLE(v[ll], v[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP_DOUBLE(v[low], v[hh]);
    /* Re-set active partition */
    if (hh <= q) low = ll;
    if (hh >= q) high = hh - 1;
  }
}

void
biomcmc_quantile_vector_double (double *original_vector, int n, double *quantile, int n_quantile, double *result) 
{ // for several values, we create temp vector (that's modified by wirth algo) only once
  double *v = (double*) biomcmc_malloc (n * sizeof (double));
  int i, q;
  for (q = 0; q < n; q++) v[q] = original_vector[q];
  for (i = 0; i < n_quantile; i++) {
    q = (int) ((double)(n) * quantile[i]);
    if (q < 0) q = 0;
    if (q > (n-2)) q = n - 1;
    result[i] = biomcmc_wirth_algorithm (v, n, q); // a bit slower than quickselect for median, I don't know general case
  }
  free (v);
}

/* Author: Wirth, Niklaus ( implementation by N. Devillard)
"Algorithms + data structures = programs". Englewood Cliffs: Prentice-Hall, 1976 */ 

double 
biomcmc_wirth_algorithm (double *a, int n, int k)
{
  int i,j,l,m;
  double x;

  l = 0; m = n-1;
  while (l < m) {
    x = a[k]; i = l; j = m;
    do {
      while (a[i] < x) i++;
      while (x < a[j]) j--;
      if (i <= j) {
        ELEM_SWAP_DOUBLE(a[i],a[j]);
        i++; j--;
      }
    } while (i <= j);
    if (j < k) l = i;
    if (k < i) m = j;
  }
  return a[k];
}

#undef ELEM_SWAP_DOUBLE
