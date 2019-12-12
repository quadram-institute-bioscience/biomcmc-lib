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

/*! \file random_number_lists.h 
 *  \brief lists of prime numbers and hand-picked random numbers used by hashes and PRNGs. Examples include the rolling
 *  hash (DNA bases mapped to a random value) and the (deterministic) spice used to generate several streams.
 */

#ifndef _biomcmc_random_number_lists_h_
#define _biomcmc_random_number_lists_h_

/*! \brief  Five-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
uint64_t sTable76[44][5] = {
  /* table 7 of MathsComput(1999)p261 (from elements 0 to 19 here) - safe for up to \f$10^{35}\f$ samples */
   {9ull, 34ull, 5ull, 26ull, 18ull},  {9ull, 32ull, 5ull, 31ull, 6ull},   {9ull, 25ull, 5ull, 37ull, 22ull},
   {10ull, 24ull, 5ull, 7ull, 12ull},  {12ull, 17ull, 5ull, 14ull, 8ull},  {12ull, 40ull, 5ull, 16ull, 22ull}, 
   {12ull, 26ull, 5ull, 34ull, 23ull}, {17ull, 27ull, 5ull, 13ull, 9ull},  {17ull, 8ull, 5ull, 37ull, 19ull},
   {20ull, 41ull, 5ull, 14ull, 6ull},  {22ull, 40ull, 5ull, 4ull, 18ull},  {22ull, 19ull, 5ull, 14ull, 19ull},
   {22ull, 41ull, 5ull, 16ull, 6ull},  {22ull, 16ull, 5ull, 32ull, 4ull},  {26ull, 9ull, 5ull, 11ull, 14ull},
   {26ull, 19ull, 5ull, 29ull, 3ull},  {44ull, 20ull, 5ull, 8ull, 6ull},   {44ull, 31ull, 5ull, 22ull, 14ull},
   {53ull, 8ull, 5ull, 23ull, 17ull},  {53ull, 12ull, 5ull, 31ull, 18ull},
   /* table 6 (from elements 20 to 43 here) - safe for up to \f$10^{31}\f$ samples) */
   {10ull, 5ull, 29ull, 23ull, 8ull},  {12ull, 5ull, 11ull, 16ull, 15ull}, {17ull, 5ull, 16ull, 6ull, 7ull},
   {17ull, 5ull, 19ull, 16ull, 14ull}, {18ull, 5ull, 37ull, 7ull, 3ull},   {19ull, 5ull, 31ull, 15ull, 13ull},
   {20ull, 5ull, 11ull, 13ull, 6ull},  {22ull, 5ull, 17ull, 10ull, 11ull}, {23ull, 5ull, 37ull, 13ull, 7ull},
   {24ull, 5ull, 7ull, 16ull, 8ull},   {26ull, 5ull, 22ull, 4ull, 9ull},   {26ull, 5ull, 26ull, 13ull, 12ull},
   {26ull, 5ull, 31ull, 14ull, 13ull}, {36ull, 5ull, 32ull, 16ull, 8ull},  {36ull, 5ull, 32ull, 21ull, 8ull},
   {39ull, 5ull, 19ull, 6ull, 8ull},   {43ull, 5ull, 14ull, 20ull, 15ull}, {44ull, 5ull, 14ull, 15ull, 15ull},
   {44ull, 5ull, 29ull, 6ull, 13ull},  {44ull, 5ull, 34ull, 25ull, 9ull},  {45ull, 5ull, 16ull, 21ull, 8ull},
   {51ull, 5ull, 28ull, 3ull, 12ull},  {53ull, 5ull, 26ull, 16ull, 8ull},  {54ull, 5ull, 28ull, 13ull, 3ull}
};

/*! \brief  Four-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
uint64_t sTable543[106][4] = {
   /* table 5, first 4 lines of MathsComput(1999)p261 (from elements 0 to 3 here) - safe for up to \f$10^{21}\f$ samples */
   {30ull, 23ull, 17ull, 18ull}, {13ull, 23ull, 26ull, 5ull}, {17ull, 38ull, 23ull, 24ull}, {26ull, 47ull, 17ull, 19ull},
   /* table 5, last 2 lines (from elements 4 to 5 here) - safe for up to \f$10^{21}\f$ samples */
   {26ull, 34ull, 20ull, 17ull}, {29ull, 38ull, 28ull, 18ull},
   /* table 4 of (from elements 6 to 97 here) - safe for up to \f$10^{17}\f$ samples */
   {18ull, 10ull, 23ull, 11ull}, {26ull, 10ull, 13ull, 11ull}, {48ull, 17ull, 30ull, 11ull}, {27ull, 20ull, 9ull, 11ull},  
   {46ull, 22ull, 9ull, 11ull},  {23ull, 29ull, 24ull, 11ull}, {25ull, 29ull, 13ull, 11ull}, {34ull, 29ull, 9ull, 11ull},  
   {50ull, 7ull, 38ull, 12ull},  {15ull, 8ull, 19ull, 12ull},  {44ull, 22ull, 16ull, 12ull}, {6ull, 23ull, 29ull, 12ull},
   {16ull, 5ull, 22ull, 13ull},  {11ull, 10ull, 25ull, 13ull}, {18ull, 11ull, 40ull, 13ull}, {19ull, 16ull, 30ull, 13ull}, 
   {45ull, 23ull, 24ull, 13ull}, {17ull, 7ull, 9ull, 14ull},   {52ull, 11ull, 20ull, 14ull}, {52ull, 22ull, 30ull, 14ull}, 
   {25ull, 23ull, 26ull, 14ull}, {27ull, 7ull, 19ull, 15ull},  {25ull, 11ull, 13ull, 15ull}, {6ull, 26ull, 31ull, 15ull},
   {19ull, 28ull, 25ull, 15ull}, {38ull, 28ull, 37ull, 15ull}, {53ull, 28ull, 18ull, 15ull}, {50ull, 29ull, 32ull, 15ull}, 
   {17ull, 32ull, 41ull, 15ull}, {39ull, 8ull, 12ull, 16ull},  {53ull, 13ull, 33ull, 16ull}, {12ull, 5ull, 13ull, 17ull},  
   {16ull, 5ull, 11ull, 17ull},  {25ull, 7ull, 32ull, 17ull},  {54ull, 10ull, 36ull, 17ull}, {45ull, 11ull, 29ull, 17ull},
   {30ull, 20ull, 18ull, 17ull}, {39ull, 20ull, 43ull, 17ull}, {19ull, 22ull, 22ull, 17ull}, {50ull, 23ull, 25ull, 17ull}, 
   {11ull, 26ull, 19ull, 17ull}, {19ull, 26ull, 11ull, 17ull}, {13ull, 29ull, 40ull, 17ull}, {46ull, 32ull, 29ull, 17ull}, 
   {20ull, 4ull, 31ull, 18ull},  {5ull, 10ull, 33ull, 18ull},  {43ull, 16ull, 31ull, 18ull}, {38ull, 23ull, 37ull, 18ull},
   {46ull, 25ull, 39ull, 18ull}, {47ull, 4ull, 26ull, 19ull},  {33ull, 7ull, 27ull, 19ull},  {18ull, 11ull, 17ull, 19ull}, 
   {43ull, 11ull, 37ull, 19ull}, {5ull, 14ull, 13ull, 19ull},  {53ull, 20ull, 27ull, 19ull}, {24ull, 25ull, 25ull, 19ull}, 
   {30ull, 25ull, 27ull, 19ull}, {34ull, 29ull, 41ull, 19ull}, {18ull, 5ull, 36ull, 20ull},  {15ull, 11ull, 18ull, 20ull},
   {52ull, 11ull, 34ull, 20ull}, {5ull, 22ull, 10ull, 20ull},  {9ull, 22ull, 10ull, 20ull},  {16ull, 23ull, 38ull, 20ull}, 
   {17ull, 23ull, 26ull, 20ull}, {40ull, 23ull, 37ull, 20ull}, {46ull, 23ull, 5ull, 20ull},  {6ull, 28ull, 27ull, 20ull},  
   {25ull, 28ull, 33ull, 20ull}, {5ull, 32ull, 26ull, 20ull},  {13ull, 7ull, 37ull, 21ull},  {26ull, 8ull, 41ull, 21ull},
   {37ull, 10ull, 43ull, 21ull}, {38ull, 10ull, 11ull, 21ull}, {30ull, 13ull, 39ull, 21ull}, {38ull, 16ull, 43ull, 21ull}, 
   {9ull, 17ull, 32ull, 21ull},  {34ull, 25ull, 17ull, 21ull}, {38ull, 26ull, 41ull, 21ull}, {8ull, 28ull, 31ull, 21ull},  
   {19ull, 29ull, 12ull, 21ull}, {37ull, 32ull, 27ull, 21ull}, {27ull, 8ull, 5ull, 22ull},   {8ull, 10ull, 29ull, 22ull},
   {41ull, 10ull, 25ull, 22ull}, {50ull, 13ull, 4ull, 22ull},  {55ull, 13ull, 37ull, 22ull}, {50ull, 17ull, 36ull, 22ull}, 
   {39ull, 26ull, 29ull, 22ull}, {55ull, 26ull, 23ull, 22ull}, {13ull, 28ull, 16ull, 22ull}, {51ull, 32ull, 10ull, 22ull},
   /* table 3 (from elements 98 to 105 here) - safe for up to \f$10^{14}\f$ samples */
   {18ull, 28ull, 7ull, 8ull},   {26ull, 20ull, 11ull, 7ull},  {19ull, 25ull, 12ull, 9ull},  {18ull, 31ull, 13ull, 6ull},
   {18ull, 22ull, 16ull, 6ull},  {30ull, 28ull, 17ull, 9ull},  {17ull, 28ull, 18ull, 6ull},  {12ull, 8ull, 22ull, 9ull}
};

/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 7 and 6 (five elements)) */
uint64_t qTable76[2][5] = { {1ull, 7ull, 24ull, 3ull, 5ull}, {1ull, 24ull, 3ull, 5ull, 3ull} };
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (tables 7 and 6 (five elements)) */
uint64_t kTable76[2][5] = { {63ull, 57ull, 55ull, 52ull, 47ull}, {63ull, 55ull, 52ull, 47ull, 41ull} };

/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3 */
uint64_t qTable543[4][4] = {
   {31ull, 1ull, 19ull, 22ull}, {31ull, 11ull, 19ull, 22ull}, {1ull, 19ull, 7ull, 24ull}, {31ull, 19ull, 24ull, 21ull}
};
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3*/
uint64_t kTable543[4][4] = {
   {63ull, 60ull, 58ull, 57ull}, {63ull, 60ull, 58ull, 57ull}, {63ull, 58ull, 57ull, 55ull}, {63ull, 58ull, 55ull, 47ull}
};

uint64_t Cmask[28] = {
  0xfffffffff0000000ull, 0xfffffffff8000000ull, 0xfffffffffc000000ull, 0xfffffffffe000000ull, /* 36-39 */
  0xffffffffff000000ull, 0xffffffffff800000ull, 0xffffffffffc00000ull, 0xffffffffffe00000ull, /* 40-43 */
  0xfffffffffff00000ull, 0xfffffffffff80000ull, 0xfffffffffffc0000ull, 0xfffffffffffe0000ull, /* 44-47 */
  0xffffffffffff0000ull, 0xffffffffffff8000ull, 0xffffffffffffc000ull, 0xffffffffffffe000ull, /* 48-51 */
  0xfffffffffffff000ull, 0xfffffffffffff800ull, 0xfffffffffffffc00ull, 0xfffffffffffffe00ull, /* 52-55 */
  0xffffffffffffff00ull, 0xffffffffffffff80ull, 0xffffffffffffffc0ull, 0xffffffffffffffe0ull, /* 56-59 */
  0xfffffffffffffff0ull, 0xfffffffffffffff8ull, 0xfffffffffffffffcull, 0xfffffffffffffffeull  /* 60-63 */
};

/* http://www.math.uni-bielefeld.de/~sillke/ALGORITHMS/random 
 * some possible 16-bit constants k for which both k*2^16-1 and k*2^15-1 are prime */
uint16_t marsaglia_constants_size = 81;
uint32_t marsaglia_constants[] = {
  18000, 18030, 18273, 18513, 18879, 19074, 19098, 19164, 19215, 19584, 19599, 19950, 20088, 20508, 20544, 20664, 20814,
  20970, 21153, 21243, 21423, 21723, 21954, 22125, 22188, 22293, 22860, 22938, 22965, 22974, 23109, 23124, 23163, 23208,
  23508, 23520, 23553, 23658, 23865, 24114, 24219, 24660, 24699, 24864, 24948, 25023, 25308, 25443, 26004, 26088, 26154,
  26550, 26679, 26838, 27183, 27258, 27753, 27795, 27810, 27834, 27960, 28320, 28380, 28689, 28710, 28794, 28854, 28959,
  28980, 29013, 29379, 29889, 30135, 30345, 30459, 30714, 30903, 30963, 31059, 31083, 36969
};

static uint16_t rnd_salt_h16_list_size = 256; /*! \brief hardcoded table size (*must* be power of 4) */ 
static uint16_t rnd_salt_h16_list[] = { // 16 bits 
  0xe02f, 0x2076, 0xea7b, 0x8547, 0x1c49, 0x211b, 0x7af3, 0x5460, 0x3e49, 0xc657, 0xa0e7, 0x169c, 0x26c2, 0x04e9, 0xcaa4, 0x88d0,
  0x5ce8, 0xa00c, 0x5c21, 0xcdf2, 0x024a, 0xbac6, 0xc8ac, 0x0a76, 0x973d, 0x5fd7, 0x79aa, 0x99cf, 0xbe46, 0x2e28, 0x4ff0, 0x4c33,
  0xbee5, 0x7ef3, 0xd911, 0x7b59, 0xa574, 0xe5d7, 0xecf4, 0xcada, 0x79ac, 0xea92, 0xbcb8, 0x19b3, 0x0998, 0xab8f, 0xfdac, 0x399b,
  0x3dea, 0x30c4, 0x00a9, 0xbe33, 0x32b8, 0x46df, 0x2931, 0x99e5, 0x5dc6, 0xe750, 0xc8cf, 0x2a15, 0xffdb, 0x2f1f, 0x99a1, 0x8c84, 
  0xcc11, 0xb0d5, 0x2123, 0x981f, 0x50c2, 0x7afd, 0x6300, 0xa28e, 0x9c25, 0xbfc7, 0x3f26, 0x32b5, 0xc3f1, 0x95f9, 0x8acd, 0x800e, 
  0x8c64, 0xd492, 0xeba4, 0xda95, 0x8e6c, 0x56f3, 0x74d8, 0x75b8, 0x9b7a, 0xf584, 0x4a39, 0x86ca, 0x879a, 0x3bea, 0x92f6, 0xdd71, 
  0xd93c, 0x8a10, 0x0166, 0xb109, 0x269e, 0xe172, 0xc726, 0x4bcf, 0xc317, 0xa53c, 0x6a31, 0x616c, 0x9acf, 0x4ba0, 0x5519, 0xc128, 
  0x2e5a, 0x128e, 0x55d8, 0xe7c5, 0x81b3, 0xc49d, 0xdbf5, 0x1713, 0xc247, 0xb34a, 0x47d5, 0xebfe, 0x2841, 0x06da, 0x2086, 0x26f5, 
  0xa973, 0x03a4, 0x6ecb, 0x9e4c, 0xc62f, 0x543d, 0xff54, 0xb68f, 0x3753, 0xde90, 0x7287, 0x7018, 0xf2e4, 0x6eea, 0x384f, 0x2a21, 
  0x40ee, 0x147b, 0xd95a, 0xa19e, 0x24bc, 0xfb51, 0x4d04, 0x2c21, 0xf7cb, 0xf77b, 0x55fd, 0x7729, 0xc0fe, 0x1573, 0x2a3a, 0x9f26, 
  0x6dfd, 0xd8c5, 0x7d13, 0xc4fd, 0xfd30, 0xdc86, 0x8585, 0xe05c, 0x622a, 0xcc15, 0xfa04, 0xb3eb, 0xe5df, 0xd387, 0xa310, 0xd2ba, 
  0x07df, 0x21e3, 0xa3bd, 0xd8d4, 0x6de1, 0xd8c2, 0x9b8f, 0x926b, 0x6688, 0x7a7c, 0x7fc8, 0xb2ab, 0x792a, 0xc832, 0x7b49, 0x6401, 
  0x9901, 0xc972, 0xdaf2, 0x4d96, 0x0f76, 0xad0e, 0xd390, 0xda01, 0x81eb, 0x214c, 0xc055, 0xf215, 0xbb8b, 0x9176, 0x8af9, 0x1688, 
  0xb442, 0xbcba, 0x82bf, 0xfea7, 0xf55b, 0x856c, 0x4b3f, 0x28d9, 0x169a, 0x04e2, 0x985e, 0xedd9, 0x19b8, 0xfffb, 0xd0f8, 0x323c, 
  0x8447, 0xa395, 0x9782, 0x0ea8, 0xa19d, 0xe59f, 0x3e68, 0xea44, 0x41a2, 0xd45a, 0x37bd, 0x898b, 0x4201, 0xea92, 0x2d5b, 0x5541, 
  0x2f8b, 0xc0d5, 0x2721, 0x1335, 0xdd9e, 0x1593, 0x37d1, 0xd22f, 0x82e5, 0x612f, 0x2a5c, 0xc027, 0x3b54, 0xb025, 0x7d9b, 0xf47e
};
  
// random tip: "printf '0x%08x\n' $(< o1) | pr -ts", " --columns 16 " will transform from decimal to hex and create multicolumn 
static uint16_t prime_salt_list_size = 512; /*! \brief hardcoded table size */
static uint32_t prime_salt_list[] = { // 256 x 16.6 bits (<100k) + 256 x 32 bits ==> not very random
  0x0,     0x17e87, 0x14369, 0x1245b, 0x0addb, 0x14d85, 0x1479d, 0x0d3a5, 0x14aef, 0x1737d, 0x04879, 0x02851, 0x10715, 0x0909d, 0x15b57, 0x0d685,
  0x0b8b7, 0x17dfd, 0x00d1f, 0x0349d, 0x09f7d, 0x0a2af, 0x159e3, 0x112f9, 0x02e25, 0x14d37, 0x00ffb, 0x07edf, 0x110bd, 0x035bd, 0x12c95, 0x10f79,
  0x03ebf, 0x13b95, 0x0dcf9, 0x044d7, 0x024cd, 0x07e41, 0x063ef, 0x059b1, 0x0327b, 0x040f9, 0x0cf79, 0x0fb4d, 0x0b6ef, 0x12713, 0x08425, 0x10505,
  0x11129, 0x16e2b, 0x182e1, 0x133f7, 0x06d0d, 0x13091, 0x149cf, 0x00d03, 0x13481, 0x038e1, 0x13a17, 0x02e39, 0x04735, 0x05245, 0x15f0b, 0x184eb,
  0x09103, 0x06475, 0x0312d, 0x168ef, 0x0f91d, 0x069e3, 0x0f293, 0x0c6cb, 0x0836f, 0x03497, 0x11333, 0x07bd3, 0x15677, 0x0135d, 0x016eb, 0x004e1,
  0x007c3, 0x0e90b, 0x11525, 0x0e9b9, 0x00907, 0x0e027, 0x0c047, 0x17fa5, 0x034f1, 0x1316f, 0x080cb, 0x0c457, 0x04a6b, 0x148fd, 0x14a35, 0x119e1,
  0x13e9b, 0x066b5, 0x1673f, 0x00bcb, 0x08617, 0x17fcb, 0x16933, 0x068e1, 0x0715d, 0x0c6f1, 0x1387d, 0x02a8b, 0x17d73, 0x11a55, 0x0e83d, 0x147e3,
  0x0f5cf, 0x11fff, 0x079ab, 0x122e7, 0x0c893, 0x0655f, 0x11117, 0x12c71, 0x15197, 0x15a5b, 0x0944b, 0x035a1, 0x00e9b, 0x146cf, 0x1417d, 0x0f047,
  0x0a9b1, 0x142b5, 0x080a5, 0x0a475, 0x06ebd, 0x0a6ab, 0x0d3eb, 0x152b7, 0x082f1, 0x156d7, 0x003d1, 0x0ba25, 0x0481d, 0x11861, 0x02807, 0x03517,
  0x06c9b, 0x0e5a3, 0x12689, 0x12e69, 0x0cea5, 0x132ad, 0x11255, 0x02a55, 0x13b8f, 0x0b599, 0x10795, 0x0abc5, 0x08f1d, 0x02cb3, 0x05dad, 0x17141,
  0x0b377, 0x09793, 0x128d5, 0x03b21, 0x0a3a9, 0x064db, 0x0dbb9, 0x0d9f7, 0x055ed, 0x07e79, 0x07721, 0x07207, 0x01357, 0x046e5, 0x00b47, 0x109db,
  0x08551, 0x0594d, 0x11e2b, 0x0e7c5, 0x06bff, 0x08ec7, 0x08a49, 0x0653f, 0x10409, 0x0c3a7, 0x06f65, 0x16f1b, 0x069d3, 0x07f6d, 0x13aad, 0x10061,
  0x173d5, 0x07991, 0x13cbb, 0x07cd5, 0x07efb, 0x0b6a5, 0x06641, 0x15901, 0x159cb, 0x0fe9b, 0x11ce3, 0x0d159, 0x0649f, 0x0961f, 0x0c3ef, 0x07a9f,
  0x0d459, 0x05047, 0x12323, 0x01741, 0x08225, 0x064ab, 0x0cabb, 0x06e29, 0x15f53, 0x11213, 0x0f763, 0x114dd, 0x023e3, 0x178df, 0x0321d, 0x00dcd,
  0x0bee7, 0x0a351, 0x0a111, 0x14a4f, 0x0a7b1, 0x111bf, 0x01cb5, 0x15967, 0x04c93, 0x00e87, 0x0ae4d, 0x159c1, 0x17abf, 0x17a0b, 0x16fc7, 0x155b9,
  0x16f9d, 0x0993b, 0x160f7, 0x0f277, 0x0e1bb, 0x0a1c5, 0x0f8ab, 0x0570b, 0x17113, 0x13af3, 0x0fe8f, 0x03c6b, 0x0f001, 0x02b4f, 0x17153, 0x1285d,
  0x1a3f8511, 0x2459568d, 0x24230afb, 0x2e1b0e4d, 0x23b8680f, 0x0310bd13, 0x23a645d9, 0x20d18c47, 0x353fc16d, 0x2e338595, 0x1aefa2c9, 0x38457ea1, // 32 bits 
  0x2da4b0ed, 0x2d7383f5, 0x385c112b, 0x242fe7cb, 0x2d738407, 0x2df900cb, 0x3475c6bb, 0x0321af5d, 0x2d619825, 0x0766148b, 0x24034277, 0x031dd179, 
  0x2dede32f, 0x38b14529, 0x386e9417, 0x19f61a83, 0x0177b927, 0x03eb850f, 0x0401330b, 0x38c2c29b, 0x1a3f8447, 0x2d2e7b1b, 0x20c57881, 0x037d2669, 
  0x19e9fe1b, 0x34708d75, 0x2e3a9d5f, 0x01416379, 0x2deb4089, 0x2d5060ad, 0x037668d1, 0x01b03d31, 0x0343eaa9, 0x20f0d9b3, 0x23b82395, 0x38af1321, 
  0x34c7a8ff, 0x38457e99, 0x349d798b, 0x017a8a49, 0x20819aa3, 0x2dd83c1f, 0x01b1b66d, 0x03a3ce91, 0x385c1135, 0x03ab06f1, 0x01ca7913, 0x01d78b49, 
  0x3531b37b, 0x1ab17089, 0x2dd93b67, 0x03bd0e15, 0x2d335bb5, 0x38c2fe43, 0x03ebc7db, 0x39148711, 0x23ab7a7b, 0x2016d181, 0x2df9f3f1, 0x2457f99d, 
  0x01516d29, 0x01579895, 0x20ea16f1, 0x01d6233f, 0x2052dc49, 0x1a77c6f1, 0x031ee243, 0x3820db73, 0x201f2c47, 0x1fe76613, 0x2dd83bb3, 0x0138f0d9, 
  0x03338ea3, 0x02ffd1f3, 0x2dcad921, 0x1a802fcb, 0x019c9a05, 0x355d9191, 0x2d49028f, 0x0340b23b, 0x03ea51f5, 0x2011f29b, 0x3834c5a7, 0x1ad2ebad, 
  0x0847a069, 0x1a298b57, 0x03cdcf03, 0x03623c5f, 0x3512fee7, 0x01df3581, 0x2367352d, 0x01021ff3, 0x2d38c785, 0x34ca3dbf, 0x3477e071, 0x07b8fb7f, 
  0x03f54d8d, 0x080831e7, 0x355ac9e9, 0x2449c9f9, 0x1a73d453, 0x030f3c41, 0x0137a1a5, 0x389aee09, 0x2df900df, 0x1a4c6131, 0x355d91ad, 0x0856bf6b, 
  0x033a6e1f, 0x38c0e241, 0x2d33a5bf, 0x1a836a83, 0x01262e39, 0x23ab7b17, 0x2002722d, 0x1a483759, 0x34a8a125, 0x34bb246b, 0x3945aaad, 0x0786c463, 
  0x2e273af5, 0x011c7be3, 0x390e5c8d, 0x03553d21, 0x01954e77, 0x0385b7df, 0x205a8ea3, 0x24204899, 0x2e47fbdb, 0x0343eaa3, 0x34ca745b, 0x1fe02ed3, 
  0x0105038f, 0x35439fd7, 0x20d7841b, 0x23dcaca5, 0x2440a065, 0x2367867d, 0x0366cb17, 0x34bc6e59, 0x0831e2f1, 0x34d810f9, 0x03d44d57, 0x38225e21, 
  0x38f34dd3, 0x2020a771, 0x2091aec5, 0x01a7c155, 0x2df988b7, 0x01615307, 0x03c5c58d, 0x34bc6e21, 0x386e255f, 0x237953c5, 0x239b564d, 0x2d11f66f, 
  0x03a11923, 0x03885fa3, 0x03095fb9, 0x2d9d4e7b, 0x247d33fb, 0x23f93083, 0x2055aa09, 0x20b7250b, 0x35144aa7, 0x24455f81, 0x34e487df, 0x07ac8d2f, 
  0x0760eeb9, 0x2d2e244b, 0x079f0a7b, 0x23f2c505, 0x34a8a10d, 0x34d0387b, 0x34ca3e33, 0x20a91403, 0x07b86221, 0x07fe3e8d, 0x23d24859, 0x23b867db, 
  0x07afec2d, 0x238b33c7, 0x2df988df, 0x07bfb063, 0x0105ead9, 0x07d0bac3, 0x24932f59, 0x1fe69c73, 0x38bd8805, 0x081d0423, 0x03a3cea1, 0x38783c93, 
  0x01e00191, 0x393ff301, 0x2452b13d, 0x0183185f, 0x20f0d9fd, 0x07f2f02b, 0x1ad07415, 0x2d9886c9, 0x2459da69, 0x3878588b, 0x0195c499, 0x2e419667, 
  0x24325b05, 0x19eed341, 0x202a13e9, 0x078eae8d, 0x393ff2c5, 0x0168c329, 0x1a43faf3, 0x3834c58d, 0x3820c583, 0x2002138d, 0x0847a093, 0x39027fd9, 
  0x2dffd6b1, 0x010f92d9, 0x07ef8753, 0x1abef6bb, 0x19d425eb, 0x07a14c55, 0x390e5c51, 0x2d33a5cd, 0x01262ebd, 0x011c7bed, 0x1a4c60ef, 0x2df3e7bb, 
  0x1ab19101, 0x03637d81, 0x037cf541, 0x01817a43, 0x0117315b, 0x1fe81199, 0x01796767, 0x353aaba7, 0x2e1b0ed9, 0x34f34d5b, 0x3903032b, 0x07c779d7, 
  0x387def0d, 0x19e2937b, 0x2381a6c7, 0x1a48fd73
};

static uint16_t rnd_salt_h64_list_size = 256;
static uint64_t rnd_salt_h64_list[] = { 
  0x247f608720f0da15ull, 0x207de5d735679355ull, 0x2df63f7434982647ull, 0x0868a4712dda9fefull, 0x076e0f9538d30781ull, 0x0771a5d31a2fafcbull, 0x390912b320ea16ddull, 0x1a48fd0b38ac7547ull,
  0x38658e0908525ba5ull, 0x2d27db5f018a9139ull, 0x1ac977a3390e5c3full, 0x33356675350209d5ull, 0x0760eec91aa3566bull, 0x01b0f28f07fccf25ull, 0x02fab3892d4902bfull, 0x0352e06b389698f9ull,
  0x2d2e7b4519dfb76bull, 0x3844516d1a5e4279ull, 0x2df48ce92d94bdf9ull, 0x03f95fb33544acddull, 0x035b309d11ded3a3ull, 0x38432e2b1a75de81ull, 0x1aefa26d34cc6eebull, 0x3943bfcb03521465ull,
  0x01a6d15101e87251ull, 0x38fec2110310f9bdull, 0x11b4368c1a3a0845ull, 0x0a0dd3642407210dull, 0x01c0f7b32e3ab38bull, 0x24156bb93497a6efull, 0x045550152408c14dull, 0x247f60df20e1c9d5ull,
  0x19e293412df341c3ull, 0x081158472029eff1ull, 0x39376a4501d55647ull, 0x2deb40e11a76ba07ull, 0x2d6fc96b03571151ull, 0x07c9e92f34aa8fa1ull, 0x011c6eff34824e51ull, 0x2df29629017a8a07ull,
  0x2048904107e5e873ull, 0x0306cc7d01534a51ull, 0x386ee0fb20935f93ull, 0x38e68d55204b150bull, 0x2e419e05387dee89ull, 0x1fef027d2e1bb211ull, 0x086cb82538182e1full, 0x2dedb0c92d27db3bull,
  0x080192092d335b97ull, 0x1a0e3867204979bbull, 0x23b8679309389998ull, 0x38991f5b19df6d59ull, 0x04049bc101ca795full, 0x209c07a3394621d3ull, 0x1a913bd51fef0229ull, 0x018974f32029ef83ull,
  0x19df5b6d0355f7e7ull, 0x0924bd34387970b1ull, 0x1ab2349900f2be07ull, 0x07bfb041209bc0c3ull, 0x0974c7161ffb2c0bull, 0x01761c2f20935eefull, 0x0391d2f123869c4full, 0x247d33dd38d839e9ull,
  0x2391ad2d0114d8e3ull, 0x07d04a0f03bd0de7ull, 0x1a19a67f03c5c563ull, 0x3550346303d44d43ull, 0x1fc768e124451b89ull, 0x39362b4907ef870dull, 0x01ab365b0849a3afull, 0x1ab234b103549f65ull,
  0x071c1a263517caa1ull, 0x0771a5fb11bd9fa0ull, 0x207db8bf23b86829ull, 0x38280817355c26d7ull, 0x2d5060ab011078dbull, 0x38c2fea316bef4d3ull, 0x039303470795cbf9ull, 0x38499d512d94bdbbull,
  0x1018c562016bf22full, 0x248b26310ba1df06ull, 0x1a8fe19d351064a5ull, 0x01a6d14920126197ull, 0x2460df9f086cb865ull, 0x1281f9512ea249ceull, 0x081d048d23935735ull, 0x38e68d091fe76627ull,
  0x2457f92903637d47ull, 0x244931ef3537bca7ull, 0x34c1f5432db1651dull, 0x38a9698501ead6e3ull, 0x244bf47f2d2cf2e7ull, 0x3855b965080e097bull, 0x25b65e1e38225e39ull, 0x351925970123a189ull,
  0x2407210124379601ull, 0x3869e25f388fa3abull, 0x085f85cb2df2963bull, 0x20cc0a23204b154bull, 0x38e26adb2d1f89f5ull, 0x2d94bdd91abd955dull, 0x2440a07723b41e63ull, 0x01087e7124230b1full,
  0x031dd1692e47fb81ull, 0x20b0f4cf07d44b85ull, 0x012f87d3011c6eedull, 0x20b0f4e30342d59bull, 0x01021fb71fef0259ull, 0x23e1109b1a24275dull, 0x34e0cb2f2d6a68a3ull, 0x2e389661209fbf93ull,
  0x0764673b38747ef9ull, 0x2a67b5df1abfedbbull, 0x2dffd69f34e71b07ull, 0x011c6f1d018a1837ull, 0x389b770d01673fefull, 0x39315b512d2e7afbull, 0x07b12c180cc9193aull, 0x032f8001355ac9dbull,
  0x01db6b1e0df1e688ull, 0x24451bb3242fe787ull, 0x086c514723fa4979ull, 0x1aaa5e8f23b41e71ull, 0x20e93c6701ddf1a3ull, 0x1fd40c5b353133d7ull, 0x39148669208255e1ull, 0x03f367911af8e6a7ull,
  0x2e3ab3c323fe9393ull, 0x085f85c303d7da69ull, 0x3820c55b3533e8ebull, 0x3547521508083201ull, 0x0137a20938d306fdull, 0x03a1196d3946210dull, 0x19e2938735524e37ull, 0x2d38c7c520b411afull,
  0x0331b5d52d7617bfull, 0x07ab12873857ee65ull, 0x018a17d12e3ab355ull, 0x014161a938499d7full, 0x38a5f4a92df9892full, 0x0d1e4eb134e0cb8full, 0x03553c7b381d5cd1ull, 0x3867c7ff38c63569ull,
  0x2e3896652f8182a2ull, 0x086deafb011c6eb1ull, 0x1aa7a33f0865156dull, 0x19f12cda1aee5a05ull, 0x2dcad935237299a5ull, 0x238b58cd1aca628dull, 0x07f8a1ab12ce4837ull, 0x1aa356750352148bull,
  0x07764a7b24932f1dull, 0x35144afb1febed35ull, 0x1a5afe3f24707e85ull, 0x017a89ef2c731ea2ull, 0x23b998af0856bf8full, 0x208d12031a500a89ull, 0x19e0397915975ef7ull, 0x20cf242920c578bdull,
  0x010f094526863126ull, 0x202f7b9703f95fa1ull, 0x388fa38123ea6dddull, 0x2271b5172df9f3fbull, 0x23de48e92402b4fbull, 0x20e270752d33a5e9ull, 0x242fe7cf2d2e247full, 0x0372263f348db8e7ull,
  0x17a19039037d2651ull, 0x03b7547f03676447ull, 0x2459564b24690835ull, 0x1aee6e5f080e09abull, 0x2069fa5534b22c93ull, 0x1aee5a133529fa57ull, 0x2d9da603076b3307ull, 0x0755a367209021b1ull,
  0x35442e5f242e8167ull, 0x34fcd7a31a3b24d3ull, 0x033f7cd3153c1176ull, 0x38cfedc92e3ab403ull, 0x1a26423f34e0a0cdull, 0x3853c8bb23ae4ebfull, 0x011731af2e1d7a4dull, 0x2d990c292dd93b65ull,
  0x2eb53eae01bbb51dull, 0x0186fab1010d2a17ull, 0x19d1c75534a5a155ull, 0x03ebc80538658e63ull, 0x02fe779d10e84f9dull, 0x20598b2f0308ce1bull, 0x34aa8f2707f9708full, 0x3471154319e039bfull,
  0x033f7cab20beaca5ull, 0x20ac4b4f0153e1c5ull, 0x204bc0a703b27659ull, 0x2dd993cd20e43855ull, 0x03ff83d723869badull, 0x35356a6338d83a31ull, 0x011e551b20e4391bull, 0x2deb40f31ab19113ull,
  0x0392ecb538b77d0full, 0x0112b7191aa8f359ull, 0x086206332dedb16bull, 0x07bc0dd30139132dull, 0x2d95e00923ea6dcbull, 0x01ca78f3019c99c9ull, 0x1fc8df2d1a7f79c7ull, 0x2d99b21535442e57ull,
  0x033f7cb70b8844d9ull, 0x01e8728938d839a7ull, 0x01c2c35320beac91ull, 0x03ff83ef34c1f527ull, 0x35356a51241b113full, 0x037225f93855b9f3ull, 0x018a180903866d8full, 0x38a04a392386c075ull,
  0x201f2c393823a3f5ull, 0x2474eded38182e5bull, 0x34e71b1b2e374e27ull, 0x08250a8d38910f01ull, 0x3896c93908a22fb3ull, 0x34aa8f8d23bcc12dull, 0x207db8c903206263ull, 0x32cead6423749a15ull,
  0x03ffcdff356bbadfull, 0x01ddf1e520a60f01ull, 0x388a3a4f1aabf761ull, 0x03af669f355586e3ull, 0x20d48cb334711571ull, 0x03338ecd2a9773a2ull, 0x0182a3912dd9933full, 0x2df3e7cb2d5cb469ull,
  0x38d2ec4520677c7bull, 0x07db8d7334e68debull, 0x01b258a1384ef0b7ull, 0x1a4cea41036b945bull, 0x207de5a134c7a95bull, 0x03ea52030406d189ull, 0x2d2cf2fd1adce2b5ull, 0x0862985f34982609ull,
  0x178b91782d38c7e5ull, 0x01a7c1571ac4ed3bull, 0x3820dbeb24379639ull, 0x2365a32d03e1d0dfull, 0x1ffb2c4d204979cdull, 0x1ab95bc323cb5b89ull, 0x240a64dd23601185ull, 0x349d79d1245cc0bfull,
  0x031ee2571fe81161ull, 0x34bc6dcd1b1d6baaull, 0x35974f6f1aed083dull, 0x39238d3b204d4c75ull, 0x2d79b4330814831dull, 0x23ae4ee3037668ddull, 0x20677c850aad9847ull, 0x2dd2efef353cc79bull
};

static uint16_t ulx_h64_size = 12; // these are used by some hash functions, please do not change their order 
static uint64_t ulx_h64[] = {
  0x65d200ce55b19ad8UL, 0x4f2162926e40c299UL, 0x162dd799029970f8UL, // 0...2
  0x68b665e6872bd1f4UL, 0xb6cfcf9d79b51db2UL, 0x7a2b92ae912898c2UL, // 3...5
  0xff51afd7ed558ccdUL, 0xc4ceb9fe1a85ec53UL, 0x87c37b91114253d5UL, 0x4cf5ad432745937fUL, // 6...9 (murmurhash) 
  0x52dce729UL, 0x38495ab5UL, // 10...11 
  11400714785074694791ull, 14029467366897019727ull, 1609587929392839161ull, // 12..14 used by xxhash
  9650029242287828579ull, 2870177450012600261ull // 15..16 used by xxhash
};

#endif
