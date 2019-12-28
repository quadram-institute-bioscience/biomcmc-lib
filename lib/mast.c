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

#include "mast.h"

void
biomcmc_mast(topology t1, topology t2)
{
  int i1, i2;
  bipartition **m;
  bip_size bip;
  trees_with_common_leaves(t1, t2);
  bip = new_bipsize (t1->nleaves);

  m = (bipartition**) biomcmc_malloc ((2*t1->nleaves - 2) * sizeof (bipartition*));
  for (i1 = 0; i1 < 2*t1->nleaves-2; i1++) { // btw nleaves should be the same for t1 and t2 
    mi[i1] = (bipartition*) biomcmc_malloc ((2*t2->nleaves - 2) * sizeof (bipartition));
    for (i2 = 0; i2 < 2*t2->nleaves-2; i2++) m[i1][i2] = new_bipartition_from_bipsize (bip->n);
  }

  fill_mast_matrix_from_trees (m, t1, t2);

  if (m) {
    for (i1 = 2*t1->nleaves-3; i1 >=0; i1--) if (m[i1]) {
      for (i2 = 2*t2->nleaves-3; i2 >=0; i2--) del_bipartition (m[i1][i2]);
      free (m[i1]);
    }
    free (m);
  }
  del_bipsize (bip);
}

void  
fill_mast_matrix_from_trees (bipartition **m, topology t1, topology t2)
{
  int i1, i2;
  bipartition bx, by, bz;
  bx = new_bipartition_from_bipsize (m[0][0]->n);
  by = new_bipartition_from_bipsize (m[0][0]->n);
  bz = new_bipartition_from_bipsize (m[0][0]->n);

  for (i1 = 0; i1 < t1->nleaves; i1++) for (i2 = 0; i2 < t2->nleaves; i2++) {
    bipartition_AND (m[i1][i2], t1->nodelist[i1]->split, t2->nodelist[i2]->split, false);
  }
  for (i1 = 0; i1 < t1->nleaves; i1++) for (i2 = 0; i2 < t2->nleaves; i2++) { // STOPHERE
    bipartition_AND (m[i1][i2], t1->nodelist[i1]->split, t2->nodelist[i2]->split, false);
  }
  for (i1 = 0; i1 < t1->nleaves-2; i1++) for (i2 = 0; i2 < t2->nleaves-2; i2++) {
    t1->postorder[i1]->id
    t2->postorder[i2]->id
    // etc
    m[t1->nleaves + i1][t2->nleaves + i2] = bx;
  }

  del_bipartition (bx);
  del_bipartition (by);
  del_bipartition (bz);
}
