/* 
 * This file is part of biomcmc-lib, , a low-level library for phylogenomic analysis. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file optics_clustering.h 
 *  \brief distance calculation between generic objects,without generating full matrix beforehand 
 */

#ifndef _biomcmc_optics_clustering_h_
#define _biomcmc_optics_clustering_h_

#include "distance_generator.h"

typedef struct optics_cluster_struct* optics_cluster;

struct optics_cluster_struct 
{
  int n_samples; 
  int *cluster, *order;
  bool *core;
  double *core_distance, *reach_distance;
};

optics_cluster new_optics_cluster (int n_samples);
void del_optics_cluster (optics_cluster oc);
void optics_cluster_reset (optics_cluster oc);
void optics_cluster_run (optics_cluster oc, distance_generator dg, int minPoints, double minDist, double clustDist);

#endif
