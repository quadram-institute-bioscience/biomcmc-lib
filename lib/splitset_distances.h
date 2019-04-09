/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.

 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file splitset_distances.h
 *  \brief Low-level functions that use only the split bipartitions of topologies -- treating them as unrooted usually. 
 *
 *  Use a "splitset" structure that copies the bipartition information of all nodes (so that original trees are
 *  untouched) and then modifies this splitset. These functions assume a gene tree (mul-tree) and a species tree (NOT
 *  mul-tree). Compared to guenomu and genefam-dist, I removed the simpler 'orthologous' functions since they assumed
 *  _same_leaves on both trees, which is not usual even without multrees. 
 */

#ifndef _biomcmc_splitset_distances_h_
#define _biomcmc_splitset_distances_h_

#include "reconciliation.h" 

/*! \brief Splitset structure for dSPR calculation (also allocates aux vectors) */
splitset new_splitset_genespecies (topology gene, topology species, reconciliation rec);
/*! \brief free memory allocated for splitset structure */
void del_splitset (splitset split);
/*! \brief approximate dSPR between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species (topology gene, topology species, splitset split);
/*! \brief RF distance between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species_rf (topology gene, topology species, splitset split);
/*! \brief h distance (edge disagreement assignment cost) between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species_hdist (topology gene, topology species, splitset split);

#endif
