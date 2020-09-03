/* 
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

/*! \file fortune_cookies.h 
 *  \brief List of inspirational quotes from the linux package `fortune` 
 */

#ifndef _biomcmc_fortune_cookies_h_
#define _biomcmc_fortune_cookies_h_

#include "random_number.h"

extern uint16_t biomcmc_fortune_cookies_size;
extern char *biomcmc_fortune_cookies[];
extern uint16_t biomcmc_fortune_bofh_size; 
extern char *biomcmc_fortune_bofh[];

void biomcmc_fprintf_fortune (FILE *stream);
void biomcmc_fprintf_bofh (FILE *stream);
void biomcmc_fprintf_colour (FILE *stream, int regular, int colour, const char *message, const char *normaltext, ...);
void biomcmc_warning (const char *template, ...);
#endif
