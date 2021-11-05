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

/*! \file 
 *  \brief file compression/decompression handler 
 *  lzma functions adapated from https://github.com/mpimd-csc/flexiblas and https://github.com/kobolabs/liblzma/blob/master/doc/examples/
 */

#ifndef _biomcmc_file_compression_h_
#define _biomcmc_file_compression_h_

#include "lowlevel.h"
#include <errno.h>
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef HAVE_LZMA
#include <lzma.h>
#endif
#ifdef HAVE_BZIP2 // #ifdef HAVE_LIBBZ2  
#include <bzlib.h>
#endif

typedef struct file_compress_struct* file_compress_t;

/* lower level function (to directly access a specific compression library) */
#ifdef HAVE_ZLIB
gzFile biomcmc_gzopen (const char *path, const char *mode);
/*! \beief zlib library version of getline */
int biomcmc_getline_gz (char **lineptr, size_t *n, gzFile zstream);
#endif

#ifdef HAVE_LZMA
/* Internal structure to represent an uncompressed file. */
typedef struct {
  FILE *fp;
  char *path, mode;
  lzma_stream strm;
  uint8_t *inbuf, *outbuf, *readbuf;
  size_t buffer_size, getc_avail, getc_pos;
  uint8_t eof;
  lzma_action action;
} xz_file_t;

xz_file_t * biomcmc_xz_open (const char *path, const char *mode, size_t buffer_size);
void biomcmc_xz_close (xz_file_t *f);
size_t biomcmc_xz_read (xz_file_t *f);
size_t biomcmc_xz_write (xz_file_t *f, char *cbuf, size_t len);
int biomcmc_xz_getc (xz_file_t *f);
int biomcmc_getline_xz (char **lineptr, size_t *n, xz_file_t *f);
#endif // HAVE_LZMA

#ifdef HAVE_BZIP2 // #ifdef HAVE_LIBBZ2  
// https://www.sourceware.org/bzip2/manual/manual.html#libprog
#include <bzlib.h>
typedef struct {
  BZFILE *fp;
  char *path, mode;
  uint8_t *readbuf;
  int buffer_size, getc_avail, getc_pos; // cannot be size_t since bz2 returns negative values 
} bz2_file_t;

bz2_file_t * biomcmc_bz2_open (const char *path, const char *mode, size_t buffer_size);
void biomcmc_bz2_close (bz2_file_t *f);
size_t biomcmc_bz2_read (bz2_file_t *f);
int biomcmc_bz2_getc (bz2_file_t *f);
int biomcmc_getline_bz2 (char **lineptr, size_t *n, bz2_file_t *f);
#endif // HAVE_LIBBZ2

/*! \brief Memory-safe fopen() function.
 *
 * Opens the file whose name is the string pointed to by path and associates a stream with it. An error message is 
 * thrown in case of failure.
 * \param[in] path file name 
 * \param[in] mode opening mode ("r" for reading, "w" for writing, etc)
 * \result pointer to file stream */
FILE *biomcmc_fopen (const char *path, const char *mode);

/*! \brief read file line-by-line (like homonymous function from GNU C library)
 *
 * This implementation is originally from the CvsGui project (http://www.wincvs.org/). The explanation from the 
 * original file adapted to our system  follows:
 * \verbatim 
   Read up to (and including) a newline ("\n") from STREAM into *LINEPTR and null-terminate it. *LINEPTR is a pointer 
   returned from malloc (or NULL), pointing to *N characters of space.  It is realloc'd as necessary.  Return the 
   number of characters read (not including the null terminator), or -1 on error or EOF. \endverbatim */
int biomcmc_getline (char **lineptr, size_t *n, FILE *stream);

struct file_compress_struct
{
  int8_t format;
#ifdef HAVE_LZMA
  xz_file_t *xz;
#endif
#ifdef HAVE_BZIP2 // #ifdef HAVE_LIBBZ2 is set if using AC_CHECK_LIB 
  bz2_file_t *bz2;
#endif
#ifdef HAVE_ZLIB
  gzFile gz;
#endif
  FILE *raw;
  char *filename;
};

file_compress_t biomcmc_open_compress (const char *path, const char *mode); // new_file_compress_t()
/*! \brief if suffix is .xz, .bz, or .gz then opens respective file for writting; otherwise assume raw (uncompressed) output */
file_compress_t biomcmc_create_compress_from_suffix (const char *path);
int biomcmc_getline_compress (char **lineptr, size_t *n, file_compress_t fc);
void biomcmc_close_compress (file_compress_t fc); // del_file_compress_t()
int biomcmc_write_compress (file_compress_t fc, char *string);

#endif
