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

/*! \file upgma.h 
 *  \brief UPGMA and bioNJ from (onedimensional representation of) distance matrices 
 *
 */

#ifndef _biomcmc_upgma_h_
#define _biomcmc_upgma_h_

#include "topology_randomise.h" 

/*! \brief lowlevel UPGMA (or single-linkage) function that depends on a topology and a matrix_distance */
void upgma_from_distance_matrix (topology tree, distance_matrix dist, bool single_linkage);
/*! \brief lowlevel bioNJ function (Gascuel and Cuong implementation) that depends on a topology and a matrix_distance */
void bionj_from_distance_matrix (topology tree, distance_matrix dist) ;
/*! \brief updates distances between species based on genes and gene-to-species mapping, with min on upper and mean on lower diagonal  */
void fill_species_dists_from_gene_dists (distance_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
/*! \brief update global (over loci) species distances besed on local (within locus) species distances */
void update_species_dists_from_spdist (distance_matrix global, distance_matrix local, int *spexist);

int prepare_spdistmatrix_from_gene_species_map (spdist_matrix spdist, int *sp_id, int n_sp_id);
void fill_spdistmatrix_from_gene_dists (spdist_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
void update_spdistmatrix_from_spdistmatrix (spdist_matrix global, spdist_matrix local);

#endif
