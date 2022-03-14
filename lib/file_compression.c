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
 */

#include "file_compression.h"

enum {FORMAT_XZ, FORMAT_BZ2, FORMAT_GZ, FORMAT_RAW};

file_compress_t
biomcmc_open_compress (const char *path, const char *mode)
{
  if (!path) biomcmc_error ("No file name was given to biomcmc_open_compress() (null pointer)\n");
  file_compress_t fc = (file_compress_t) biomcmc_malloc (sizeof (struct file_compress_struct));

  fc->filename = (char*) biomcmc_malloc ((strlen (path) + 1) * sizeof(char));
  strncpy (fc->filename, path, strlen (path));
  fc->filename[strlen(path)] = '\0';

#ifdef HAVE_LZMA
  fc->xz = biomcmc_xz_open (path, mode, 1<<16); // return NULL if file is not lzma
  if (fc->xz != NULL) { fc->format = FORMAT_XZ; return fc; }
#endif
#ifdef HAVE_BZIP2
  fc->bz2 = biomcmc_bz2_open (path, mode, 1<<16); // return NULL if file is not bzip2
  if (fc->bz2 != NULL) { fc->format = FORMAT_BZ2; return fc; }
#endif
#ifdef HAVE_ZLIB
  fc->gz = biomcmc_gzopen (path, mode);
  if (fc->gz != NULL) { fc->format = FORMAT_GZ; return fc; }
#endif
  fc->raw = biomcmc_fopen (path, mode);
  fc->format = FORMAT_RAW; 
  return fc;
}

file_compress_t
biomcmc_create_compress_from_suffix (const char *path)
{
  if (!path) biomcmc_error ("No file name was given to biomcmc_create_compress_from_suffix() (null pointer)\n");
  size_t last = strlen (path);
  file_compress_t fc = (file_compress_t) biomcmc_malloc (sizeof (struct file_compress_struct));

  fc->filename = (char*) biomcmc_malloc ((strlen (path) + 1) * sizeof(char));
  strncpy (fc->filename, path, strlen (path));
  fc->filename[strlen(path)] = '\0';

  if ((path[last-3] == '.') && (path[last-1] == 'z')) {
#ifdef HAVE_LZMA
    if (path[last-2] == 'x') {
      fc->xz = biomcmc_xz_open (path, "w", 1<<16); // return NULL if file is not lzma
      if (fc->xz != NULL) { fc->format = FORMAT_XZ; return fc; }
    }
#endif
#ifdef HAVE_BZIP2
    if (path[last-2] == 'b') {
      fc->bz2 = biomcmc_bz2_open (path, "w", 1<<16); // return NULL if file is not bzip2
      if (fc->bz2 != NULL) { fc->format = FORMAT_BZ2; return fc; }
    }
#endif
#ifdef HAVE_ZLIB
    if (path[last-2] == 'g') {
      fc->gz = biomcmc_gzopen (path, "w");
      if (fc->gz != NULL) { fc->format = FORMAT_GZ; return fc; }
    }
#endif
    //if arrived here, then required library was not available; remove suffix and save in raw format
    fc->filename = (char*) biomcmc_realloc ((char*)fc->filename, (last - 2) * sizeof(char)); // one extra for null char
    fc->filename[last-3] = '\0'; // replace "." by null char 
  }
  fc->raw = biomcmc_fopen (fc->filename, "w");
  fc->format = FORMAT_RAW; 
  return fc;
}

int 
biomcmc_getline_compress (char **lineptr, size_t *n, file_compress_t fc)
{
#ifdef HAVE_LZMA
  if (fc->format == FORMAT_XZ) return biomcmc_getline_xz (lineptr, n, fc->xz);
#endif
#ifdef HAVE_BZIP2
  if (fc->format == FORMAT_BZ2) return biomcmc_getline_bz2 (lineptr, n, fc->bz2);
#endif
#ifdef HAVE_ZLIB
  if (fc->format == FORMAT_GZ) return biomcmc_getline_gz (lineptr, n, fc->gz);
#endif
  return biomcmc_getline (lineptr, n, fc->raw);
}

void
biomcmc_close_compress (file_compress_t fc)
{
  if (!fc) return;
#ifdef HAVE_LZMA
  if (fc->format == FORMAT_XZ) biomcmc_xz_close (fc->xz);
#endif
#ifdef HAVE_BZIP2
  if (fc->format == FORMAT_BZ2) biomcmc_bz2_close (fc->bz2);
#endif
#ifdef HAVE_ZLIB
  if (fc->format == FORMAT_GZ) gzclose (fc->gz);
#endif
  if (fc->format == FORMAT_RAW) fclose (fc->raw);
  if (fc->filename) free (fc->filename);
  free (fc);
  return;
}

int
biomcmc_write_compress (file_compress_t fc, char *string)
{
  if (!fc) return 0;
#ifdef HAVE_LZMA
  if (fc->format == FORMAT_XZ) return biomcmc_xz_write (fc->xz, string, strlen (string));
#endif
#ifdef HAVE_BZIP2
  if (fc->format == FORMAT_BZ2) return BZ2_bzwrite (fc->bz2, string, strlen (string));
#endif
#ifdef HAVE_ZLIB
  if (fc->format == FORMAT_GZ) return gzprintf (fc->gz, "%s", string);
#endif
  return fprintf (fc->raw, "%s", string);
}

/**   lowlevel functions (can be used independently)   **/

FILE *
biomcmc_fopen (const char *path, const char *mode)
{
  FILE *fp = fopen (path, mode);
  if (fp == NULL) {
    fprintf (stderr, "Please check if path is correct, if there are non-ASCII characters in file name,\n");
    fprintf (stderr, "if you have enough permissions (to read/write). Remember that paths are relative to\n");
    fprintf (stderr, "where this program is being called\n");
    biomcmc_error ( "problem opening file \"%s\" with mode \"%s\"", path, mode);
  }
  return fp;
}

#ifdef HAVE_ZLIB
gzFile
biomcmc_gzopen (const char *path, const char *mode)
{
  gzFile zfp = gzopen (path, mode); // can also process raw (uncompressed); o.w. we cannot exit with failure
  if (zfp == NULL) { // gzopen handles raw data naturally 
    fprintf (stderr, "Please check if path is correct, if there are non-ASCII characters in file name,\n");
    fprintf (stderr, "if you have enough permissions (to read/write). Remember that paths are relative to\n");
    fprintf (stderr, "where this program is being called\n");
    biomcmc_error ( "problem opening file \"%s\" with mode \"%s\" with zlib", path, mode);
  }
  const char *err_mesg;
  int err_number;
  err_mesg = gzerror (zfp, &err_number);
  if (err_number != Z_OK) biomcmc_error ("zlib error when trying to open \"%s\". zlib error message:\n%s", path, err_mesg);
  return zfp;
}
#endif

/* \brief size, in bytes, when extending the buffer of biomcmc_getline() */
#define MIN_CHUNK 256
int
biomcmc_getline (char **lineptr, size_t *n, FILE *stream)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline()");
  if (!stream)  biomcmc_error ("lack of input file in biomcmc_getline()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = getc (stream);

    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the 
     * first char) NUL-terminate the line buffer.  
     */
    if ((*lineptr + *n) != (read_pos + nchars_avail))  biomcmc_error ("problem_1 setting string size in biomcmc_getline()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) biomcmc_error ("problem_2 setting string size in biomcmc_getline()");
    }

    if (ferror (stream)) return -1;

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}

#ifdef HAVE_ZLIB
int
biomcmc_getline_gz (char **lineptr, size_t *n, gzFile zstream)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */
  const char *err_mesg;
  int err_number;

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline_gz() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline_gz()");
  if (!zstream) biomcmc_error ("lack of input file in biomcmc_getline_gz()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = gzgetc (zstream);
    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the 
     * first char) NUL-terminate the line buffer.  
     */
    if ((*lineptr + *n) != (read_pos + nchars_avail))  biomcmc_error ("problem_1 setting string size in biomcmc_getline_gz()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) biomcmc_error ("problem_2 setting string size in biomcmc_getline_gz()");
    }

    err_mesg = gzerror (zstream, &err_number);
    if (err_number != Z_OK) biomcmc_error ("error in biomcmc_getline_gz()::  %s", err_mesg);

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }
  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}
#endif

#ifdef HAVE_LZMA
bool init_xz_encoder (lzma_stream *strm, uint32_t preset); // preset should be 3~7 (default 6)
bool init_xz_decoder (lzma_stream *strm);
void del_xz_file_t (xz_file_t *f);

bool 
init_xz_encoder (lzma_stream *strm, uint32_t preset) // preset should be 3~7 (default 6)
{
#ifdef _OPENMP 
	lzma_mt mt_options = {
		.flags = 0,
//		.threads = sysconf(_SC_NPROCESSORS_ONLN),
		.block_size = 0,
		.timeout = 0,
		.filters = NULL,
    .preset = preset,
		.check = LZMA_CHECK_CRC64,
	};
  mt_options.threads = lzma_cputhreads();
  lzma_ret ret = lzma_stream_encoder_mt (strm, &mt_options);
#else
  lzma_ret ret = lzma_easy_encoder (strm, preset, LZMA_CHECK_CRC64);
#endif

  if (ret == LZMA_OK) return true;    // Return successfully if the initialization went fine.
  const char *msg;  /* error */
  switch (ret) {
    case LZMA_MEM_ERROR: msg = "LZMA:: Memory allocation failed"; break;
    case LZMA_OPTIONS_ERROR: msg = "LZMA:: Specified preset is not supported"; break;
    case LZMA_UNSUPPORTED_CHECK: msg = "LZMA:: Specified integrity check is not supported"; break;
    default: msg = "LZMA:: Unknown error, possibly a bug"; break;
  }
  fprintf (stderr, "LZMA:: Error initializing the encoder: %s (error code %u)\n", msg, ret);
  return false;
}

bool
init_xz_decoder (lzma_stream *strm)
{
  lzma_ret ret = lzma_stream_decoder( strm, UINT64_MAX, 0x00 | LZMA_CONCATENATED );
  if (ret == LZMA_OK) return true;  // Return successfully if the initialization went fine.
  const char *msg;  /* error */
  switch (ret) {
    case LZMA_MEM_ERROR: msg = "LZMA:: Memory allocation failed"; break;
    case LZMA_OPTIONS_ERROR: msg = "LZMA:: Unsupported decompressor flags"; break;
    default: msg = "LZMA:: Unknown error, possibly a bug"; break;
  }
  fprintf(stderr, "LZMA:: Error initializing the decoder: %s (error code %u)\n",msg, ret);
  return false;
}

xz_file_t * 
biomcmc_xz_open (const char *path, const char *mode, size_t buffer_size) 
{
  xz_file_t *f;
  int err;

  if ((*mode != 'w') && (*mode != 'r')) {fprintf (stderr, "xz_open():: unrecognised mode %c\n", *mode); return NULL; }
  /* Initlize the data structure  */
  f = (xz_file_t *) malloc (sizeof (xz_file_t));
  memset(&(f->strm), 0, sizeof (f->strm));
  f->path = strdup(path);
  f->mode = *mode;
  if (buffer_size < 1024) buffer_size = 1024;
  f->buffer_size = buffer_size;
  f->getc_pos = f->getc_avail = 0;
  f->inbuf   = (uint8_t*) malloc (f->buffer_size * sizeof (uint8_t));
  f->outbuf  = (uint8_t*) malloc (f->buffer_size * sizeof (uint8_t));
  f->readbuf = (uint8_t*) malloc (f->buffer_size * sizeof (uint8_t)); // this interfaces with ext functions
  memset (f->readbuf, 0, f->buffer_size * sizeof (uint8_t));

  if ( f->mode == 'w' ) {  /* Open for writing */
    f->fp = fopen(f->path, "w");
    if ( !(f->fp)) {
      err = errno;
      fprintf (stderr, " Opening %s as XZ for writing failed, errno: %03d - %s\n", path, err, strerror(err));
      del_xz_file_t (f);
      return NULL;
    }
    if (! init_xz_encoder (&(f->strm), 6)) {
      fprintf (stderr, "Can not initialize lzma (XZ) encoder for writing on file %s\n", path);
      del_xz_file_t (f);
      return NULL;
    }
    f->strm.next_in = NULL;
    f->strm.avail_in = 0;
    f->strm.next_out = f->outbuf;
    f->strm.avail_out = f->buffer_size;
    f->action = LZMA_RUN;
  } else { /* Open for READING */
    f->fp =fopen(f->path, "r");
    if ( !(f->fp)) {
      err = errno;
      fprintf (stderr, " Opening %s as XZ for reading failed, errno: %03d - %s\n", path, err, strerror(err));
      del_xz_file_t (f);
      return NULL;
    }
    if (!init_xz_decoder(&(f->strm))) {
      fprintf (stderr, "Can not initialize lzma (XZ) decoder for reading from file %s\n", path);
      del_xz_file_t (f);
      return NULL;
    }
    f->strm.next_in = NULL;
    f->strm.avail_in = 0;
    f->strm.next_out = f->outbuf;
    f->strm.avail_out = f->buffer_size;
    f->eof = 0;
    f->action = LZMA_RUN;
    // reads one block to check if file is actually XZ
    f->getc_avail = biomcmc_xz_read (f);
    if (!f->getc_avail) { del_xz_file_t (f); return NULL; }

    if (biomcmc_xz_getc (f) == EOF) return NULL;
    f->getc_pos--; // if file is valid xz then we already read one
  }
  return f;
}

void
biomcmc_xz_close (xz_file_t *f)
{
  lzma_ret ret;
  size_t write_size = 0;

  if (!f) { fprintf(stderr, "LZMA:: Error: data == NULL\n"); return; }
  if (f->mode == 'w') {
    f->action = LZMA_FINISH;
    ret = LZMA_OK;
    while ( ret != LZMA_STREAM_END ) {  /* Finish the encoding   */
      ret = lzma_code(&(f->strm), f->action);
      // If the output buffer is full or if the compression finished successfully, write the data from the output buffer to the output file.
      if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) {
        // When lzma_code() has returned LZMA_STREAM_END, the output buffer is likely to be only partially full. 
        // Calculate how much new data there is to be written to the output file.
        write_size = f->buffer_size - f->strm.avail_out; /*originally:: write_size = sizeof (f->outbuf) - f->strm.avail_out;*/
        if (fwrite (f->outbuf, 1, write_size, f->fp) != write_size) { fprintf (stderr, "LZMA:: Write error on closing: %s\n", strerror(errno)); return; }
        // Reset next_out and avail_out.
        f->strm.next_out = f->outbuf;
        f->strm.avail_out = f->buffer_size; /* originally  f->strm.avail_out = sizeof (f->outbuf); */
      }
      if ( ret != LZMA_STREAM_END && ret != LZMA_OK) { fprintf(stderr, "LZMA Encode error\n"); break; }
    }
  }
  del_xz_file_t (f);
  return;
}

void
del_xz_file_t (xz_file_t *f)
{
  if (!f) return;
  if (&(f->strm)) lzma_end (&(f->strm)); 
  if (f->fp) { fflush (f->fp); fclose(f->fp); }
  f->fp = NULL;
  if (f->path)    free (f->path);
  if (f->inbuf)   free (f->inbuf);
  if (f->outbuf)  free (f->outbuf);
  if (f->readbuf) free (f->readbuf);
  free (f);
}

/* Reads a block from an uncompressed file and stores in readbuf */ 
size_t 
biomcmc_xz_read (xz_file_t *f) 
{
  lzma_ret ret;
  size_t write_size = 0;

  if ((f==NULL) || (f->eof)) return 0;
  f->strm.next_out = f->outbuf;
  f->strm.avail_out = f->buffer_size;

  while (!f->eof) {
    if (f->strm.avail_in == 0 && !feof(f->fp)) {
      f->strm.next_in = f->inbuf;
      f->strm.avail_in = fread (f->inbuf, 1, sizeof (f->inbuf), f->fp);
      if (ferror(f->fp)) { fprintf(stderr, "LZMA:: Read error: %s\n",strerror(errno)); return 0; }
      // Once the end of the input file has been reached, we need to tell lzma_code() that no more input will be coming.
      if (feof(f->fp)) f->action = LZMA_FINISH;
    }

    ret = lzma_code (&(f->strm), f->action);

    if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) { // if output buffer is full or decompression finished
      write_size = f->buffer_size - f->strm.avail_out;
      memcpy (f->readbuf, f->outbuf, sizeof (uint8_t) * write_size);
      f->strm.next_out = f->outbuf;
      f->strm.avail_out = f->buffer_size;
      if (ret == LZMA_STREAM_END) f->eof = 1;
      break;
    } // if output buffer full or decompression finished

    if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
      const char *msg;
      switch (ret) {
        case LZMA_MEM_ERROR: msg = "LZMA read:: Memory allocation failed"; break;
        case LZMA_FORMAT_ERROR: msg = "LZMA read:: The input is not in the .xz format"; break; // FIXME: xz_open must check for this code
        case LZMA_OPTIONS_ERROR: msg = "LZMA read:: Unsupported compression options"; break;
        case LZMA_DATA_ERROR: msg = "LZMA read::Compressed file is corrupt"; break;
        case LZMA_BUF_ERROR: msg = "LZMA read::Compressed file is truncated or otherwise corrupt"; break;
        case LZMA_MEMLIMIT_ERROR: msg = "LZMA read::The memory limit for decompression is too small"; break;
        default: msg = "LZMA read:: Unknown error, possibly a bug"; break;
      }
      if (ret != LZMA_FORMAT_ERROR) fprintf (stderr, "LZMA read:: Decoder error: %s (error code %u)\n",msg, ret);
      return 0; // if LZMA_FORMAT_ERROR then getc will return EOF (getc is called to check if really is xz)
    }
  } // while not eof
  return write_size;
}

size_t 
biomcmc_xz_write (xz_file_t *f, char *cbuf, size_t len) 
{
  lzma_ret ret;
  size_t pos = 0, write_size = 0;

  if ((f == NULL) || ( cbuf == NULL)) return 0;

  while (true) { // Fill the input buffer if it is empty.
    if (f->strm.avail_in == 0 && pos != len) {
      size_t tocopy = BIOMCMC_MIN (f->buffer_size, (len-pos));
      memcpy (f->inbuf, cbuf+pos, tocopy);
      f->strm.next_in = f->inbuf;
      f->strm.avail_in = tocopy;
      pos = pos + tocopy;
    } else if ( pos >= len && f->strm.avail_in == 0 ) return len;

    // Up to strm->avail_out bytes of compressed output will be written starting from strm->next_out. avail_out and next_out
    // will be incremented by an equal amount to match the number of output bytes written.
    ret = lzma_code (&(f->strm), f->action);

    if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) {
      write_size = f->buffer_size - f->strm.avail_out; /*originally:: write_size = sizeof (f->outbuf) - f->strm.avail_out;*/
      if (fwrite(f->outbuf, 1, write_size, f->fp) != write_size) { fprintf(stderr, "LZMA:: Write error: %s\n", strerror(errno)); return 0; }
      f->strm.next_out = f->outbuf;      // Reset next_out and avail_out.
      f->strm.avail_out = f->buffer_size; /* originally  f->strm.avail_out = sizeof (f->outbuf); */
    }

    if (ret != LZMA_OK) {  // Normally the return value of lzma_code() will be LZMA_OK until everything has been encoded.
      const char *msg = "";
      if (ret == LZMA_STREAM_END) return true;
      switch (ret) {
        case LZMA_MEM_ERROR: msg = "LZMA write:: Memory allocation failed"; break;
        case LZMA_DATA_ERROR: msg = "LZMA write:: File size limits exceeded"; break;
        default: break;
      }
      fprintf (stderr, "LZMA write:: Encoder error: %s (error code %u)\n",msg, ret);
      return 0;
    }
  }
  return 0;
}

int
biomcmc_xz_getc (xz_file_t *f)
{
  if (!f) return EOF;
  if (f->getc_pos == f->getc_avail) {
    f->getc_avail = biomcmc_xz_read (f);
    if (f->getc_avail <= 0) { f->eof = 1; return EOF; }
    f->getc_pos = 0;
  }
  if (f->getc_pos < f->getc_avail) return (int) f->readbuf[f->getc_pos++];
  return EOF;
}

int
biomcmc_getline_xz (char **lineptr, size_t *n, xz_file_t *f)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline_xz() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline_xz()");
  if (!f)       biomcmc_error ("lack of input file in biomcmc_getline_xz()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = biomcmc_xz_getc (f);
    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the first char) NUL-terminate the line buffer.*/
    if ((*lineptr + *n) != (read_pos + nchars_avail))  biomcmc_error ("problem_1 setting string size in biomcmc_getline_xz()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) biomcmc_error ("problem_2 setting string size in biomcmc_getline_xz()");
    }

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}

#endif  //HAVE_LZMA

#ifdef HAVE_BZIP2
bz2_file_t * 
biomcmc_bz2_open (const char *path, const char *mode, size_t buffer_size) 
{
  bz2_file_t *f;

  if ((*mode != 'w') && (*mode != 'r')) {fprintf (stderr, "unrecognised mode %c for bzip2\n", *mode); return NULL; }
  /* Initialize the data structure  */
  f = (bz2_file_t *) malloc (sizeof (bz2_file_t));
  f->path = strdup (path);
  f->mode = *mode;
  if (buffer_size < 1024) buffer_size = 1024;
  f->buffer_size = buffer_size;
  f->getc_pos = f->getc_avail = 0;
  f->readbuf = (uint8_t*) malloc (f->buffer_size * sizeof (uint8_t)); // this interfaces with ext functions
  memset (f->readbuf, 0, f->buffer_size * sizeof (uint8_t));
  if ( !(f->fp = BZ2_bzopen(path, mode))) {
    int err = errno;
    fprintf (stderr, "Opening bzip2 file %s failed errno: %03d - %s \n", path, err, strerror (err));
    biomcmc_bz2_close (f);
    return NULL;
  }
  // reads one block to check if file is actually bzip2
  f->getc_avail = biomcmc_bz2_read (f);
  if (f->getc_avail < 1) { biomcmc_bz2_close (f); return NULL; }
  return f;
}

void
biomcmc_bz2_close (bz2_file_t *f)
{
  if (!f) return;
  if (f->fp) BZ2_bzclose (f->fp);
  f->fp = NULL;
  if (f->path)    free (f->path);
  if (f->readbuf) free (f->readbuf);
  free (f);
  return;
}

size_t 
biomcmc_bz2_read (bz2_file_t *f) 
{
  if (!f) return 0; // default is int BZ2_bzRead ( int *bzerror, BZFILE *b, void *buf, int len ) which gives error code
  return BZ2_bzread (f->fp, (void*) f->readbuf, (int) f->buffer_size);
}

int
biomcmc_bz2_getc (bz2_file_t *f)
{
  if (!f) return EOF;
  if (f->getc_pos == f->getc_avail) {
    f->getc_avail = biomcmc_bz2_read (f);
    if (f->getc_avail <= 0) return EOF; 
    f->getc_pos = 0;
  }
  if (f->getc_pos < f->getc_avail) return (int) f->readbuf[f->getc_pos++];
  return EOF;
}

int
biomcmc_getline_bz2 (char **lineptr, size_t *n, bz2_file_t *f)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */
  const char *err_mesg;
  int err_number;

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline_bz2() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline_bz2()");
  if (!f)       biomcmc_error ("lack of input file in biomcmc_getline_bz2()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = biomcmc_bz2_getc (f);
    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the first char) NUL-terminate the line buffer.*/
    if ((*lineptr + *n) != (read_pos + nchars_avail))  biomcmc_error ("problem_1 setting string size in biomcmc_getline_xz()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) biomcmc_error ("problem_2 setting string size in biomcmc_getline_xz()");
    }
    
    err_mesg = BZ2_bzerror (f->fp, &err_number);
    if (err_number != Z_OK) biomcmc_error ("error %d in biomcmc_getline_bz2 ()::  %s", err_number, err_mesg);

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}

#endif  //HAVE_BZIP2
