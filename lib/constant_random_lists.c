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
 * 
 * Collection of prime numbers: http://www.primos.mat.br/ 
 * to subsample and print as hex:
 *  ```
 *   for i in 2T_part*; do shuf $i | head -16 | cut -f 1; done | sort -n > o
 *   for i in `cat o`; do echo "obase=16; $i" | bc; done 
 *  ```
 * Or `printf '0x%08xULL\n' $(< o1) | pr -ts", " --columns 6 ` will transform from decimal to hex and create multicolumn 
 */

#include "constant_random_lists.h"

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

uint16_t rnd_salt_h16_list_size = 256; /*! \brief hardcoded table size (*must* be power of 4) */ 
uint16_t rnd_salt_h16_list[] = { // 16 bits 
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
uint16_t prime_salt_list_size = 1024; /*! \brief hardcoded table size */
uint32_t prime_salt_list[] = { // 256 x 16.6 bits (<100k) + 768 x 32 bits ==> not very random
  0x1f,    0x17e87, 0x14369, 0x1245b, 0x0addb, 0x14d85, 0x1479d, 0x0d3a5, 0x14aef, 0x1737d, 0x04879, 0x02851, 0x10715, 0x0909d, 0x15b57, 0x0d685,
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
  0xc7b455db, 0xcdd2c1fb, 0x44175e99, 0x60d11091, 0xc1bfde0d, 0xce4ce165, 0x67f8e975, 0xb4e6a25f, 0x07dd6bdd, 0x1ca99d5f, 0x6291c78d, 0x5fc3f0b7,
  0xb32eec23, 0x677477cb, 0x09693281, 0xb6c68dc5, 0x6b021905, 0x54c65cdb, 0x04c93c85, 0x1a9948b7, 0x0510b34d, 0x3e0e919b, 0xcc3df6bb, 0x58016b99,
  0x06dfabf1, 0x3f780571, 0xb0c57b07, 0xad7645df, 0x200646a3, 0xaf638967, 0x06a5be5b, 0xaee80609, 0xb66bdfff, 0x019d2eaf, 0xb0e14cd7, 0xcb0d19cf,
  0x1be2347b, 0x5df966df, 0x545b2a4f, 0xad81f9cf, 0xc1018217, 0x03118cfd, 0xd0d0d223, 0xc93574f1, 0x05eacb31, 0xb6b3258d, 0x467ad2e5, 0xb616d9c1,
  0x6b7779ef, 0xb3a84f55, 0x3e282e99, 0xb4408769, 0x02d26b73, 0x6847ed19, 0x580e770b, 0x1842bf4b, 0xad497dc9, 0x202d21ab, 0xb34431c3, 0x005ebf1d,
  0x0a17f353, 0xb7547e1f, 0x04320e8f, 0x20e5d801, 0x53a74dd7, 0x5fe6ffa9, 0x6a52e40b, 0xd194adfb, 0xc4e403b1, 0xbb66eea5, 0xc7e5e8a7, 0xb1c0610d,
  0x59e91509, 0xc39edf85, 0x22238b3f, 0x69ebc80d, 0x039095ad, 0xc474786d, 0xce0684c3, 0x02e935ad, 0xae52287b, 0x09cdde7b, 0x550dd163, 0x46dea5d3,
  0x58685ec3, 0x541db5e3, 0xb0b29f6d, 0x088e438d, 0xb334e515, 0xbb8f636b, 0x18347ba3, 0x3c2a3da9, 0x5ab1f719, 0xb6abd46f, 0x5766e053, 0x43136a0b,
  0x1f9023cb, 0x64df1415, 0x69e6cf63, 0x0225e5f9, 0x059c81f9, 0x3c62b505, 0xb3f8b833, 0x5d2d8ddf, 0xbcea4983, 0xad94d531, 0x406c25b9, 0x63e406ad,
  0xc1a21b27, 0xd1524e13, 0x00d76fd5, 0x65b71d0f, 0xb17db5a1, 0x01f84d0f, 0x569d4595, 0x695be625, 0x0472d1a1, 0x21540a93, 0xbe088377, 0xbef7d8c9,
  0xc1e1cf25, 0xba01e7ab, 0xb6c74291, 0xd36e3cdd, 0x208f33b5, 0xb1d6bc61, 0x553f5e77, 0x0404be1f, 0x57b5582d, 0x683b3175, 0x099dec9d, 0x04a7506f,
  0xcfd79d89, 0x5ca8b913, 0xaef8d2f9, 0xb3dbbdbd, 0x5860b001, 0x174b4817, 0xb16a9cc7, 0xb497390d, 0x09106823, 0x421bf1c7, 0x3af0bc19, 0x69a45d65,
  0xb30eb769, 0x083227c7, 0x5aed1ba7, 0x3b584691, 0x4558c9db, 0xc862ed6d, 0xcb6d4369, 0x61019a6d, 0x1c2e035f, 0x59dd7935, 0xb5a3247b, 0x5c91cf13,
  0x086bbba5, 0x090a552f, 0x53da8541, 0x08bdc2e1, 0xc941c93b, 0xb297626d, 0x0479ddc3, 0x60f67d05, 0xc0bb1a0d, 0xb23352b9, 0x0a11e335, 0x3fbf0cab,
  0x1b10b2c9, 0x1e18fff5, 0x41c5487f, 0xbbf2f1f9, 0x6137188f, 0xb0cf0fcb, 0xb8c27973, 0x5dd12379, 0xaef00c17, 0x02e5fe81, 0x5fca0eab, 0x65b9636d,
  0x45331709, 0xb4e648c5, 0xb9a400a9, 0xb7083ce7, 0xccb23fa3, 0x07998f73, 0x1f067903, 0x1d32d94d, 0x1eb0957b, 0x5c655b0f, 0xc2b75dcb, 0x07ac968f,
  0x1b0f8491, 0x5683d61b, 0x6878a2a1, 0x09e94901, 0x5aa1e7a5, 0xca726a89, 0xb8a20e67, 0x097f0f25, 0xca614b41, 0x1a3b616f, 0x608fa215, 0x1f2e2aa5,
  0xb68fd783, 0x45156bbf, 0x5d99a6e7, 0x5dca2145, 0x06a4346b, 0xbb8c3ba7, 0xada03b09, 0x408c741b, 0x097048a9, 0x3f613319, 0x3d8b6b97, 0x4393bcf1,
  0x6cb9a6d9, 0x07df6a5f, 0x0a680e83, 0xb261b249, 0x65e59b6f, 0x678baf9b, 0x1c009699, 0xcb3ae4e9, 0xbcb3f95b, 0x045bcc41, 0xc7111a93, 0x6045c88d,
  0x0aab78e7, 0x0a53a261, 0xb4c42bcb, 0x6aeb7e1d, 0x617d35cf, 0x0159aca5, 0x1fac2ec7, 0xd1533f3d, 0x44e91e3f, 0xcbae5df9, 0x08006afd, 0x199ffbfb,
  0x04028003, 0x57313477, 0x41f71393, 0x53fdb537, 0x6775e52d, 0xad1babcb, 0xb3cea29d, 0xbb847b3d, 0xd3f7a6e5, 0xc1bf4e37, 0xd1e55837, 0x6275120f,
  0x65de7d5d, 0xb9451c59, 0xc2739019, 0xcda8d4d3, 0x681c7e4f, 0xbe194cad, 0x6c792e09, 0xc3f188c9, 0xc58254e3, 0xbe3a33e9, 0x032cafc7, 0xb7379255,
  0xd0d53ed1, 0xb0dc12fd, 0xc485723b, 0x0573a685, 0x07a51d3d, 0x04be9687, 0x41ec8e09, 0x65e042e7, 0x09ada4bb, 0x65956e4d, 0x1a18466b, 0xafda49df,
  0x5f80fde7, 0xc612131b, 0x6b88ef6f, 0x3d6240eb, 0xb849bf7b, 0xc63c05a1, 0xb210230f, 0x41fe1edf, 0xca7d82d9, 0x0a79048f, 0xb032a35d, 0x3be78b13,
  0xb0fa54b3, 0x1810c5af, 0xb621b615, 0x021a2d21, 0x6537b4b9, 0x43fe766d, 0x628869d9, 0xb412e31d, 0x6cae7e1d, 0xd1d1f727, 0xaed645db, 0x073f9049,
  0x019b1907, 0x0176c44d, 0xc7a319d9, 0xafd73a63, 0x673cffb5, 0xc97b0427, 0x014f500d, 0xb4508abf, 0x5e4768f3, 0xce1be721, 0xc1f9a159, 0x03478177,
  0x60be7969, 0xb5667323, 0xbd15f623, 0x21ab1fb1, 0xae476e67, 0x007f91a1, 0xc91126b1, 0x3e6e2e53, 0x09ba9ccb, 0x0014f1d9, 0xcd671107, 0x3c5ea64d,
  0x460c0961, 0x5533ab9f, 0x5d3d7c7b, 0xc7fe4513, 0x686aae1d, 0xaf66c715, 0xc921f0c7, 0xbef38043, 0xae2ad459, 0xaefc94db, 0x3e092873, 0xae4db17b,
  0xafcb0d1b, 0xb0dedc17, 0xba4d8241, 0x2117d89f, 0x539e418f, 0xd22bd1fb, 0xb2270055, 0xafc30585, 0xb03d8f63, 0xc895741f, 0x62618b81, 0xb2eca055,
  0xb3ead151, 0x1a8f458b, 0x56441837, 0x06fc0079, 0xc373b41b, 0x034df60b, 0xbfce03db, 0xc5dc30d9, 0x5ed608eb, 0xae5975b3, 0xcd4a52f7, 0x063798c7,
  0x04211293, 0x04b57d8d, 0xcafff669, 0xb1503605, 0xba74a88f, 0x1e6282e5, 0xb9f7a4d7, 0x5ebc6899, 0xcb865e35, 0xb6c760ef, 0x45c2424b, 0x02c3f699,
  0xc18af38b, 0xbce2f6a3, 0x69a44ff1, 0x427c1017, 0xae80ac4d, 0xcc6fd5fb, 0x21feb217, 0xb86b2cad, 0x626a64db, 0x5edbdb6d, 0x0553ee87, 0x02a3d2d9,
  0x04a9fb39, 0xaf403419, 0xb48a4f1b, 0x60ebbccd, 0xc31ee497, 0x6416e655, 0xb40d91d5, 0xced31d7d, 0x64b5af0d, 0xb0f0d28d, 0x00718ee9, 0x04891c2b,
  0xcaaa2a23, 0x3b0ee01f, 0xb66e5947, 0x3fdfacb9, 0x3c8478cd, 0x41b2fc71, 0xb3e770c5, 0x03975893, 0x0696a30d, 0x5a1cec49, 0x3ec55013, 0x44c92e6b,
  0xbce9b0eb, 0x5d74ced5, 0x40d9ef93, 0x072dc6d5, 0xcfdaae6b, 0x1bec9ca7, 0x5cd8a803, 0xcfdafbaf, 0x563c4359, 0xb1f530fb, 0xaf8c29e7, 0x1721487d,
  0x09bca32f, 0xd2152747, 0x5886ccd7, 0x4561214f, 0x0852aee9, 0xb72b74c7, 0xc028e87d, 0xc3b3d133, 0x19a96a79, 0xcdf81511, 0x5c6d6083, 0x08c364a9,
  0x640e11d7, 0xbb7d3fbb, 0x0a6561c9, 0x0aa04a81, 0xc5d248f5, 0x58a34155, 0x69ec27b7, 0x00e0857f, 0x630ad5e7, 0xbbb6425b, 0x6bd34f69, 0xb014dacb,
  0x007231db, 0xb3e2550d, 0xb1b5995d, 0xad7832e7, 0xc2b583b9, 0x6948b5b5, 0xada2cb81, 0x692fea0f, 0xb9149681, 0x5d00eb4b, 0xb4264f55, 0x41451185,
  0xcd4a441d, 0xb65f7a8f, 0x0449963b, 0x69f401db, 0xb4e346bb, 0xb4111085, 0x1f2e5267, 0xb8602787, 0xb9a3debf, 0x3bf28391, 0x4176cfa9, 0xbb3e6817,
  0x18fbf7a1, 0x3cd52275, 0xc1595231, 0xcc9c7e49, 0xb754e54b, 0xc6217e83, 0xc2b21141, 0xb01e6443, 0x0030782b, 0xb18e14cb, 0xb840a239, 0x044f30c9,
  0x09fe3683, 0x5c003fef, 0xb7cee2c7, 0x08a3d5fd, 0xc0396f43, 0xcae7dc9b, 0x0752707d, 0x1c32e92f, 0x0463e731, 0x0638c273, 0xadcdcc6d, 0x17b066bb,
  0x559f33df, 0x5b156631, 0xca26c1ff, 0x0a32f713, 0x41eb9adb, 0xbe222d6b, 0x6c6d9ce7, 0x052be533, 0x45fbebdf, 0xb8aec2f5, 0x42f36845, 0x09d6de93,
  0xbe99a051, 0xb58cbaef, 0x4086abfd, 0x1f112da9, 0x4144d101, 0xaeda69bd, 0xb85aaefd, 0x6403b78d, 0x387def0d, 0x19e2937b, 0x2381a6c7, 0x1a48fd73
};

uint16_t rnd_salt_h64_list_size = 512; 
uint64_t rnd_salt_h64_list[] = {  // 256 "regular" random numbers, and 256 prime numbers (smaller than e^13)
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
  0x031ee2571fe81161ull, 0x34bc6dcd1b1d6baaull, 0x35974f6f1aed083dull, 0x39238d3b204d4c75ull, 0x2d79b4330814831dull, 0x23ae4ee3037668ddull, 0x20677c850aad9847ull, 0x2dd2efef353cc79bull,
  // random prime numbers from http://www.primos.mat.br/ (64bits, but smaller than 10^13); remember that they're all odd numbers 
  0x03b18f4e21ull, 0x5bdb030315ull, 0x062cd9d891ull, 0x18f71bda41ull, 0x333d260b05ull, 0x1b02e7d3d3ull, 0x382ac3403bull, 0xe0050dd397ull, 0x96c2a19d33ull, 0xb250018151ull,
  0xb662abc811ull, 0xd99ac79dcdull, 0xa73dd57d93ull, 0x00c6f3bea3ull, 0xd3088b3205ull, 0x96b5085357ull, 0x897d67317dull, 0x0134a77999ull, 0x3962a757cfull, 0x977266eaa3ull,
  0x56daf29089ull, 0x05fc8d5145ull, 0x1f25c219b9ull, 0xe854622ef9ull, 0x8975f33255ull, 0x01b26256a9ull, 0xb66219ca4dull, 0xe001bcb685ull, 0x6dbea981f5ull, 0x7c5da9136dull,
  0x24ca1e99b7ull, 0x06a834a927ull, 0x173e1ea6a7ull, 0x7eb6902dddull, 0x00c8db55cbull, 0x473376b1d7ull, 0xa740417b95ull, 0x0504b91619ull, 0x12469103c3ull, 0x15d89edc61ull,
  0x174076faf3ull, 0x02a1b1fd8full, 0x96f4e2b147ull, 0x01ab34cca1ull, 0x4730b832a1ull, 0xbc3aaf4395ull, 0x54cbf0aff5ull, 0x09ef97ea79ull, 0x02aa7d4a71ull, 0xe0042ca541ull,
  0xa741036da7ull, 0x23f9a1be83ull, 0x006a81966dull, 0x5158c1b65dull, 0x686224f7e3ull, 0x2400238ccfull, 0x85bea4984dull, 0x77e29acfa5ull, 0x1ecebba231ull, 0xe895eee73dull,
  0xa96a6f052bull, 0x08d4129087ull, 0x00b2dac023ull, 0x81d6f061e3ull, 0x4f1d4b8b51ull, 0x1b06cf717dull, 0x2a49bd5d75ull, 0x2c66e5060dull, 0x0c9bcdd921ull, 0x6da4fbfe6full,
  0x0134fa7e55ull, 0x182c4dfb05ull, 0x472c126981ull, 0x92a43f49c3ull, 0x092a58d1fdull, 0xc4a631a97bull, 0x48cbb66af5ull, 0x1566277f0full, 0x2f17b6bb6dull, 0x00194c0a65ull,
  0x97699fb445ull, 0xcb8240b529ull, 0x2e21454443ull, 0xbc77c3f183ull, 0x5559065abbull, 0x81d76d8c0full, 0x06ce261555ull, 0x74f912e11dull, 0x049cf6572bull, 0x74eeafe835ull,
  0x2c68b51de1ull, 0x6172522fabull, 0x02a15c3c49ull, 0x45419de4bbull, 0x0254a786efull, 0x2e53f47c45ull, 0x73de24775bull, 0x6dc41a0555ull, 0x45466b74efull, 0x74f7844f39ull,
  0x04fce22b61ull, 0x0421487645ull, 0x005704be03ull, 0x616efe7d7bull, 0x05fcf8ebefull, 0x96fffafe9bull, 0x96fa8345a9ull, 0x1564bd37c5ull, 0x2f1c042c49ull, 0x2a4fe2952bull,
  0x0606041c59ull, 0xe287251eefull, 0xe84f6c149dull, 0x515a24cbe7ull, 0x00cfcb4febull, 0xcc16894e79ull, 0x2ac0f8fb5bull, 0x74ec50fbffull, 0x7eb26bc5a1ull, 0x09265656c7ull,
  0x00b1128329ull, 0x381d89279bull, 0x032ab87573ull, 0x0041b1d3f1ull, 0x08cf626ccdull, 0x18f2648ba9ull, 0xb65b1a6d29ull, 0x07a74e5e4bull, 0xbe6d7dba4full, 0x6175d2eb0full,
  0x48ce06fc23ull, 0x024f862f33ull, 0x85c131a5dfull, 0xa96ec61a0full, 0x3f8fe94bc7ull, 0x016cd6e70full, 0x0a8c455ff3ull, 0xd306104467ull, 0x3f90f19357ull, 0xcb8238bd0bull,
  0xa2b76f327dull, 0x6536870653ull, 0x396282b429ull, 0x454843e7e5ull, 0x4235bb3667ull, 0x60e376f233ull, 0x2c6e96c2d1ull, 0x1f24177241ull, 0x24047be08bull, 0x0215a6479full,
  0x024eb8a735ull, 0x000a5ee733ull, 0x092d56be07ull, 0x16893242f3ull, 0x02a8c22611ull, 0x01409208b3ull, 0x89781441b3ull, 0x77df4b9c59ull, 0x020c237f5dull, 0x60ea5dba1dull,
  0xcc1653a7f3ull, 0x0f3801ab0bull, 0x7c5581e673ull, 0x06050f52fdull, 0x00052be533ull, 0x09ef18591bull, 0xae0696f815ull, 0x0a8ddf37dfull, 0x15aa543ce9ull, 0xe85056975full,
  0x96f55cdd73ull, 0xb6569998d7ull, 0x61948dc5cdull, 0xb656e126f1ull, 0x5c91dfd58bull, 0x0507f38927ull, 0x3331329863ull, 0xb6c4c41b1bull, 0x168d87c787ull, 0x15d80dc831ull,
  0x60e774d715ull, 0x000338ab6bull, 0xe283453ac1ull, 0x03328ef347ull, 0x032e4dd145ull, 0x24c7610527ull, 0x381da1bbb7ull, 0xa2b46c8135ull, 0x5c94561b2dull, 0x0420d62c05ull,
  0x23f7a306b5ull, 0x15dc156cf3ull, 0x653a44c6cfull, 0xbc799b8c57ull, 0x617d147305ull, 0xa04c71a27full, 0x079f948fe5ull, 0x74f2fabe4full, 0x155a5ca45bull, 0xc53236c09full,
  0xe84fe09e83ull, 0x1242f87d4full, 0xbc3c1e04d1ull, 0x7b173bbb6bull, 0x00613f2b81ull, 0x07a37805b7ull, 0x24cbd46859ull, 0xa966fa938dull, 0xc49e708ca7ull, 0xe28c6eb937ull,
  0x24c5d97d1dull, 0x126b4ab6c1ull, 0xbc3d0b7ee9ull, 0x515c1fdc37ull, 0x1742e9b9f7ull, 0x515a2ab2a5ull, 0xe89390b725ull, 0xbc378f1f8bull, 0x00cfcaab03ull, 0xba085d666bull,
  0x00b745bf7dull, 0x32557d5bdfull, 0x1b0522ba4dull, 0xd99bf1d999ull, 0xbe75cbc385ull, 0x8246d92b1bull, 0x422ce594f1ull, 0x9768472283ull, 0xcb7ee4d2a5ull, 0x652b51ebb5ull,
  0x0008526bf9ull, 0x4e032dc271ull, 0x125f099609ull, 0x92a08df5bdull, 0x9ebaac2e39ull, 0x049a87d245ull, 0x062dbe654full, 0x24c3cf8afbull, 0x06c668faf5ull, 0x032ae2039dull,
  0x041f102735ull, 0x0a8a889429ull, 0x049d690e3full, 0x15df3da437ull, 0x2f11cd4d3full, 0x0500ba44dfull, 0x014b710fc5ull, 0x041bc2f991ull, 0x08d74bce73ull, 0x0630611a6bull,
  0x15e0e739f5ull, 0x1564864e5bull, 0x85c205c5ffull, 0xb65c00e8f3ull, 0x978fc169c9ull, 0x155e1ac961ull
};

uint16_t ulx_h64_size = 256; // contains constants used by xxhash etc, mixed with other random numbers
uint64_t ulx_h64[] = {
  0x65d200ce55b19ad8ULL, 0xa9582618e03fc9aaULL, 0x08fe9de839e32377ULL, 0x1ba81b5424bcb2a8ULL, 0x0d8d7b7e013c01bdULL, 0x395537df01e0994fULL, 0x32511b0928203f96ULL, 0x2d0bc01f27b91db0ULL,
  0x4f2162926e40c299ULL, 0x39abdc4529b1661cULL, 0x3934374303a56c91ULL, 0x2224f3d51801d004ULL, 0x11765d9128047766ULL, 0x12ed6f8b236ad874ULL, 0x17effb4d119e2c01ULL, 0x05a4ec980f86d273ULL,
  0x162dd799029970f8ULL, 0x76e15d3efefdcbbfULL, 0x2f63130f255d6ec8ULL, 0x32f614f40c83a2efULL, 0x30fecd1334b92ffaULL, 0x2864ee591b5f632fULL, 0x1fb601f132ffbacfULL, 0x3aad6e7f2dc86536ULL,
  0x68b665e6872bd1f4ULL, 0xc5004e441c522fb3ULL, 0x020510f734dc102cULL, 0x200b51e436ada8d3ULL, 0x2ebb48542d4594c1ULL, 0x3b47d17a0c990d08ULL, 0x00224bec1cffbb9dULL, 0x31464e25314190bfULL,
  0xb6cfcf9d79b51db2ULL, 0x77710069854ee241ULL, 0x08cc0d2f39a5e791ULL, 0x31ae8f3a0df2e857ULL, 0x0670ac731af10527ULL, 0x171dd77c033bdf99ULL, 0x18a299f01aff8d72ULL, 0x0262de3211dcabcdULL,
  0x7a2b92ae912898c2ULL, 0x6d0c71D67aeb5b9dULL, 0x2628d32b15d8fb0cULL, 0x29d5f28e2b434440ULL, 0x1cd88e5b35d0938eULL, 0x0ffa773c1c7da99eULL, 0x37ef4ac915fbbe4fULL, 0x17287dcf0c2ee8f5ULL,
  0xff51afd7ed558ccdULL, 0x42b6b9b6e2274c79ULL, 0x162fa3db26d526fcULL, 0x2a4120651e473bc7ULL, 0x0054b2202680f0f1ULL, 0x1245d2c432bb2376ULL, 0x18360ba81eff64b3ULL, 0x21c4b1413aed1726ULL,
  0xc4ceb9fe1a85ec53ULL, 0x00B502aF529770f8ULL, 0x2a0c6d6f16a57d79ULL, 0x38f7c83d27b73d02ULL, 0x025802e0151a34efULL, 0x23a0aab233240811ULL, 0x30220fcc227fe6d8ULL, 0x2627cffb3b498a14ULL,
  0x87c37b91114253d5ULL, 0xa13c8e4680b583ebULL, 0x1f42f4e60f297da2ULL, 0x1ae50ec910926d6bULL, 0x24d0621e1ed53f7dULL, 0x2c86cccd3070950cULL, 0x1af2a5b7034e6466ULL, 0x1e9b44802883ace7ULL,
  0x4cf5ad432745937fULL, 0x61c8864680b583ebULL, 0x111f40df1eec279dULL, 0x0b4c804d0e9b725fULL, 0x181264d328907bc9ULL, 0x268fe438370c29eaULL, 0x31525baa075a9c16ULL, 0x0090e9940ea5f115ULL,
  0x9e3779b185ebca87ULL, 0x033f7cd3153caa76ULL, 0x3955d141111bfef4ULL, 0x3330cfd826887730ULL, 0x360f07823b94702dULL, 0x29462fb130f9e10cULL, 0x313ed639167a8928ULL, 0x035dca9b05f08e88ULL,
  0xc2b2ae3d27d4eb4fULL, 0x3837edc92e67b403ULL, 0x365ec72e3780197aULL, 0x16d4b48423a9d533ULL, 0x0134da85094a8ef4ULL, 0x23cc668200352433ULL, 0x1cab9eca0d980fe9ULL, 0x1696c58b14d5976fULL,
  0x165667b19e3779f9ULL, 0x1aa1923f34e0a0cdULL, 0x2e26ca942dfe50cbULL, 0x0acbe5d6329ff53bULL, 0x36080a3d32e28bf5ULL, 0x2d8475d2143dc72aULL, 0x1fe4815205cce088ULL, 0x24b4a50233eae1e2ULL,
  0x85ebca77c2b2ae63ULL, 0x12a60c064e612766ULL, 0x1ea960e23547ef4cULL, 0x21d1a6753289bf4fULL, 0x145c5df425ba4beeULL, 0x0f0dd7a104f29dd3ULL, 0x2cf5100126c7ea9fULL, 0x01468b203a5fcfadULL,
  0x27d4eb2f165667c5ULL, 0xB5026F5AA96619E9ULL, 0x16ede4ac2bf5e8b0ULL, 0x16bfdb2239ee38acULL, 0x0ad3117903f614c2ULL, 0x337c4648264d7b39ULL, 0x0de1567038d77254ULL, 0x1e867f430dae63d8ULL,
  0xdf900294d8f554a5ULL, 0x54891a96c5c02018ULL, 0x3b6286f016686582ULL, 0x040afbb426d170d7ULL, 0x365d62482490ee41ULL, 0x03e98904197f7b1dULL, 0x24df6c552da4963bULL, 0x07f718d6231becb3ULL,
  0x170865df4b3201fcULL, 0x3432f3f4a0793005ULL, 0x2e7b325f0c14beedULL, 0x1ed8855907592707ULL, 0x2ab6d1270ba9a899ULL, 0x2a835b9c35604503ULL, 0x0630eb69175d4b95ULL, 0x300d0caf1c179d72ULL,
  0xd2a98b26625eee7bULL, 0x9d2c5545ebe6e8a0ULL, 0x0b3d81dc1df2237fULL, 0x116468e13a90979bULL, 0x12adfe2d07c60749ULL, 0x2db810f42bada6cfULL, 0x306422432624c7fdULL, 0x2dad3c2b39e3e1ceULL,
  0xdddf9b1090aa7ac1ULL, 0x1812a4b3ce30a253ULL, 0x09fe3f1c14bc7c90ULL, 0x2a5c20d32b92d804ULL, 0x34b6f4b739ce9b58ULL, 0x2b093d02216ab186ULL, 0x1291cc6634ca3628ULL, 0x2a132db41d99d8e9ULL,
  0x2bd7a6a6e99c2ddcULL, 0x15dd913930b9ef26ULL, 0x0f1c8aa805ea9871ULL, 0x0d0e0de228f205d7ULL, 0x35f263cc1e0f9433ULL, 0x20ca32e01d10aa4eULL, 0x080b36ea02b57850ULL, 0x0fbe74f6023b18b1ULL,
  0x0992ccaf6a6fca05ULL, 0xe118a5f414798669ULL, 0x0232c419184bdf0cULL, 0x2e450b743172e15eULL, 0x36dfea420f5ca2d5ULL, 0x1d334a2e3a610ad1ULL, 0x1f15639432f17875ULL, 0x02e85ad226213f8bULL,
  0x360fd5f2cf8d5d99ULL, 0x30085cb207003845ULL, 0x03f654082d3a5b7aULL, 0x2db499af0c397c7eULL, 0x32d646c30847d74bULL, 0x11e02c6a376b677bULL, 0x0a780b94248c7976ULL, 0x16a151f233403148ULL,
  0x9c6e6877736c46e3ULL, 0x286293c355aa1757ULL, 0x0b1cb9b139e68eb3ULL, 0x24dea4a90cc9e0d6ULL, 0x2eeea8410c3b1375ULL, 0x0bb206dd25139aeaULL, 0x39cc78ac30db8601ULL, 0x2ab63e051196bbb7ULL,
  0x180ec6d33cfd0abaULL, 0x29288647153e7672ULL, 0x303f21642d7453e3ULL, 0x08a451840733dc98ULL, 0x1d08c7972a31a258ULL, 0x0a58df7a321c215bULL, 0x237ffb610f7439c8ULL, 0x1e44e5872f4871bbULL,
  0xd5a61266f0c9392cULL, 0x343121cb09abcd8bULL, 0x008b8f31239581f5ULL, 0x298eb6b7189f7b04ULL, 0x1d0be2a82b553d3dULL, 0x2fbeb318125c938eULL, 0x405b5d623ffa6f3bULL, 0x0146ef3efaa584bdULL,
  0x0060cea59c375f71ULL, 0x009193fc97aee695ULL, 0x05f113fed3d1f427ULL, 0x31036d9ca7f322f5ULL, 0x03668392fbb5507bULL, 0x028572fbbe600d27ULL, 0x1177117c5ad6b1a3ULL, 0x29bae3b69cbcd763ULL,
  0x1bd718de3c39a4ffULL, 0x20d150ac1ec85289ULL, 0x0fc5dc00cdb2ccb9ULL, 0x0a5a24b5141ac261ULL, 0x2a3ab43beb6eed39ULL, 0x00e7e0b1d82c947bULL, 0x0b8674e7ad18580bULL, 0x02516ecb791e0fc7ULL,
  0x000e7fa08e236fa7ULL, 0x03451bcbfe7a217bULL, 0x30efd0db104d8dafULL, 0x3e7e2acc8ecc3e6fULL, 0x4c8d2732afaac517ULL, 0x05e88282095e664dULL, 0x31573560d50a9309ULL, 0x79feb89b28279fc3ULL,
  0x39dfb711290073d5ULL, 0x733fdab3e5d66f7dULL, 0x2067fab1dd0d0759ULL, 0x002aa3eab047dc19ULL, 0x09aa2775c5d7c643ULL, 0x42866edcd0c408c3ULL, 0x07325e893427f599ULL, 0x08e1235af97f32d9ULL,
  0x234401cae6140d25ULL, 0xe380feea11accb75ULL, 0x67fbdef7659bfca1ULL, 0x02a8df697fba79d7ULL, 0x1273b8f434aa097dULL, 0x8a4a05147b4b0223ULL, 0x0cf9cce3aed12db5ULL, 0x4469d75ff0b73103ULL,
  0x28d382bacbb82e65ULL, 0x070b55847165cdbbULL, 0x2cc015d8000b145dULL, 0x0e27a5c6f2e3275dULL, 0x000c290d7063c621ULL, 0x0c8286001f685655ULL, 0x7d781f3a7e2f277bULL, 0x0d083c3601293725ULL,
  0x07251ac2f2de3b45ULL, 0x40e7cccd929517e7ULL, 0x0d5c095bafe9d09dULL, 0x0f8e473e5b88f8b1ULL, 0x944ed79994c25157ULL, 0x00113748774a681fULL, 0x0d1fe98512579261ULL, 0x68f2e70f4d541e2dULL
};

/* * * from prob_distribution_aux.c * * */

uint16_t lgamma_algmcs_size = 15;
double   lgamma_algmcs[] = {
  +.1666389480451863247205729650822e+0,  -.1384948176067563840732986059135e-4,  +.9810825646924729426157171547487e-8, 
  -.1809129475572494194263306266719e-10, +.6221098041892605227126015543416e-13, -.3399615005417721944303330599666e-15, 
  +.2683181998482698748957538846666e-17, -.2868042435334643284144622399999e-19, +.3962837061046434803679306666666e-21, 
  -.6831888753985766870111999999999e-23, +.1429227355942498147573333333333e-24, -.3547598158101070547199999999999e-26,
  +.1025680058010470912000000000000e-27, -.3401102254316748799999999999999e-29, +.1276642195630062933333333333333e-30
};

uint16_t lgamma_coeffs_size = 40;
double   lgamma_coeffs[] = { /* (zeta(2)-1)/2, (zeta(3)-1)/3 ... (zeta(40+1)-1)/(40+1) */
  0.3224670334241132182362075833230126e-0, 0.6735230105319809513324605383715000e-1,
  0.2058080842778454787900092413529198e-1, 0.7385551028673985266273097291406834e-2,
  0.2890510330741523285752988298486755e-2, 0.1192753911703260977113935692828109e-2,
  0.5096695247430424223356548135815582e-3, 0.2231547584535793797614188036013401e-3,
  0.9945751278180853371459589003190170e-4, 0.4492623673813314170020750240635786e-4,
  0.2050721277567069155316650397830591e-4, 0.9439488275268395903987425104415055e-5,
  0.4374866789907487804181793223952411e-5, 0.2039215753801366236781900709670839e-5,
  0.9551412130407419832857179772951265e-6, 0.4492469198764566043294290331193655e-6,
  0.2120718480555466586923135901077628e-6, 0.1004322482396809960872083050053344e-6,
  0.4769810169363980565760193417246730e-7, 0.2271109460894316491031998116062124e-7,
  0.1083865921489695409107491757968159e-7, 0.5183475041970046655121248647057669e-8,
  0.2483674543802478317185008663991718e-8, 0.1192140140586091207442548202774640e-8,
  0.5731367241678862013330194857961011e-9, 0.2759522885124233145178149692816341e-9,
  0.1330476437424448948149715720858008e-9, 0.6422964563838100022082448087644648e-10,
  0.3104424774732227276239215783404066e-10, 0.1502138408075414217093301048780668e-10,
  0.7275974480239079662504549924814047e-11, 0.3527742476575915083615072228655483e-11,
  0.1711991790559617908601084114443031e-11, 0.8315385841420284819798357793954418e-12,
  0.4042200525289440065536008957032895e-12, 0.1966475631096616490411045679010286e-12,
  0.9573630387838555763782200936508615e-13, 0.4664076026428374224576492565974577e-13,
  0.2273736960065972320633279596737272e-13, 0.1109139947083452201658320007192334e-13
};

uint16_t stirl_sferr_halves_size = 31;
double   stirl_sferr_halves[] = { /* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0. */
  0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
  0.0413406959554092940938221,  0.03316287351993628748511048, 0.02767792568499833914878929,
  0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
  0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
  0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
  0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
  0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
  0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
  0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
  0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690 
};
