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

/*! \file quickselect_quantile.h
 *  \brief find k smallest element in vector  
 *
 */

#ifndef _biomcmc_quickselect_quantile_h_
#define _biomcmc_quickselect_quantile_h_

#include "lowlevel.h"

double biomcmc_quantile_double (double *original_vector, int n, double quantile);
void   biomcmc_quantile_vector_double (double *original_vector, int n, double *quantile, int n_quantile, double *result);
double biomcmc_wirth_algorithm (double *a, int n, int k); /*!< \brief find k-smallest element, changing vector a[] */

#endif
