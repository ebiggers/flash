/*
 * util.c: miscellaneous useful functions for FLASH
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013 Eric Biggers
 *
 * This file is part of FLASH, a fast tool to merge overlapping paired-end
 * reads.
 *
 * FLASH is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FLASH; if not, see http://www.gnu.org/licenses/.
 */

#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "util.h"

#ifdef __WIN32__
    /* Get the pthread mutex declarations as a replacement for flockfile() and
     * funlockfile(). */
#  include <pthread.h>
    /* Get the GetSystemInfo() declaration as replacement for
     * sysconf(_SC_NPROCESSORS_ONLN). */
#  include <windows.h>
#endif

#ifndef O_BINARY
#  define O_BINARY 0
#endif

const char canonical_ascii_tab[256] = {
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', /* A */ 'A', 'N', /* C */ 'C', 'N', 'N', 'N', /* G */ 'G',
                'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', /* T */ 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', /* a */ 'A', 'N', /* c */ 'C', 'N', 'N', 'N', /* g */ 'G',
                'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', /* t */ 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
};

const char complement_tab[] = {
	['A'] = 'T',
	['C'] = 'G',
	['G'] = 'C',
	['T'] = 'A',
	['N'] = 'N',
};

#define PROGRAM_TAG "[FLASH] "

/* Prints an error message and exits the program with failure status. */
void fatal_error(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);
	fflush(stdout);
	fputs(PROGRAM_TAG "ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	exit(1);
}

void fatal_error_with_errno(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);
	fflush(stdout);
	fputs(PROGRAM_TAG "ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	fprintf(stderr, ": %s\n", strerror(errno));
	va_end(va);
	exit(1);
}

/* Prints a warning message. */
void warning(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);
	fputs(PROGRAM_TAG "WARNING: ", stderr);
	fprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
}

#ifdef __WIN32__
static pthread_mutex_t stdout_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

/* Prints an informational message. */
void info(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);

#ifdef __WIN32__
	pthread_mutex_lock(&stdout_lock);
#else
	flockfile(stdout);
#endif
	fputs(PROGRAM_TAG, stdout);
	vprintf(msg, va);
	putchar('\n');
	fflush(stdout);
	va_end(va);
#ifdef __WIN32__
	pthread_mutex_unlock(&stdout_lock);
#else
	funlockfile(stdout);
#endif
}

/* Returns the number of processors (if it can be determined), otherwise returns
 * 1. */
unsigned get_default_num_threads()
{
#ifdef __WIN32__
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#else
	long nthreads = sysconf(_SC_NPROCESSORS_ONLN);
	if (nthreads == -1) {
		warning("Could not deteremine number of processors! Assuming 1");
		return 1;
	} else {
		return (unsigned)nthreads;
	}
#endif
}

/* malloc(), exiting the program with failure status if memory allocation fails.
 * */
void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		fatal_error("Out of memory: tried to allocate %zu bytes", size);
	return p;
}

void *xrealloc(void *ptr, size_t size)
{
	void *p = realloc(ptr, size);
	if (!p)
		fatal_error("Out of memory: tried to reallocate %zu bytes", size);
	return p;
}

/* fopen(), exiting the program with failure status if file open fails, and
 * interpreting "-" as standard output.  */
void *xfopen(const char *filename, const char *mode)
{
	if (strcmp(filename, "-") == 0)
		return stdout;

	FILE *fp = fopen(filename, mode);
	if (!fp) {
		fatal_error_with_errno("Could not open the file \"%s\" for %s",
			    filename,
			    strchr(mode, 'w') ? "writing" : "reading");
	}
	return fp;
}

/* gzopen(), exiting the program with failure status if file open fails, and
 * interpreting "-" as standard output or standard input depending on the
 * requested mode.  */
void *xgzopen(const char *filename, const char *mode)
{
	gzFile f;
	if (strcmp(filename, "-") == 0)  {
		int fd;
		if (strchr(mode, 'w'))
			fd = STDOUT_FILENO;
		else
			fd = STDIN_FILENO;
		f = gzdopen(fd, mode);
	} else {
		f = gzopen(filename, mode);
	}
	if (!f) {
		fatal_error_with_errno("Failed to open the file \"%s\"",
				       filename);
	}
	return f;
}

/* Open a pipe to the compression command @compress_prog (global variable),
 * redirecting the output to the file @filename.  If @filename is "-", the
 * compressed data is sent to stdout. */
void *xpopen(const char *filename, const char *mode)
{
	size_t len = strlen(compress_prog) + 100 +
		     strlen(filename) + strlen(compress_prog_args);

	char command[len + 1];

	if (filename[0] == '-' && filename[1] == '\0') {
		/* write to stdout */
		sprintf(command, "%s %s -c -", compress_prog,
			compress_prog_args);
	} else {
		/* redirect to a file */
		sprintf(command, "%s %s -c - > '%s'", compress_prog,
			compress_prog_args, filename);
	}

	FILE *fp = popen(command, mode);
	if (!fp)
		fatal_error_with_errno("Could not launch the command \"%s\"", command);
	return fp;
}

void xfclose(void *fp)
{
	if (fp && fclose((FILE*)fp) != 0)
		fatal_error_with_errno("Failed to close output file");
}

void xpclose(void *fp)
{
	if (fp && pclose((FILE*)fp) == -1)
		fatal_error_with_errno("Failed to close pipe to output file");
}

void xgzclose(void *fp)
{
	if (fp && gzclose((gzFile)fp) != Z_OK)
		fatal_error_with_errno("Failed to close output file");
}

#ifdef __WIN32__
#  define mkdir(path, mode) mkdir(path)
#endif

/* Like `mkdir -p': create directory, and all parent directories, as needed,
 * failing only if a needed directory cannot be created. */
void mkdir_p(const char *dir)
{
	size_t len = strlen(dir);
	char dir_copy[len + 1];
	char *p = dir_copy;
	for (size_t i = 0; i < len; i++) {
		/* Copy the directory name to the @dir_copy array, squashing
		 * together consecutive path separators. */
		if ((dir[i] != '/' && dir[i] != '\\')
			|| (dir[i + 1] != '/' &&
			    dir[i + 1] != '\\'))
		{
			*p++ = dir[i];
		}
	}
	*p = '\0';

	p = dir_copy;
	do {
		if (p != dir_copy && (*p == '/' || *p == '\\' || *p == '\0')) {
			char orig_char = *p;
			*p = '\0';

			if (mkdir(dir_copy, 0755) != 0 && errno != EEXIST) {
				fatal_error_with_errno("Cannot create "
						       "directory \"%s\"",
						       dir_copy);
			}
			*p = orig_char;
		}
	} while (*p++ != '\0');
}

/* Read data from gzFile, with error checking.  Returns number of bytes
 * successfully read; 0 implies stream is at end-of-file.  Aborts on error.  */
static size_t xgzread(void *_fp, void *buf, size_t count)
{
	gzFile gzfp = (gzFile)_fp;
	int ret;

retry:
	ret = gzread(gzfp, buf, count);
	if (ret >= 0) {
		return ret;
	} else if (gzeof(gzfp)) {
		return 0;
	} else {
		int errnum;
		const char *error_str;

		error_str = gzerror(gzfp, &errnum);
		if (errnum == Z_ERRNO) {
			if (errno == EINTR) {
				goto retry;
			} else {
				fatal_error_with_errno("Error reading "
						       "input file");
			}
		} else {
			fatal_error("zlib error while reading input "
				    "file: %s", error_str);
		}
	}
}

/* Read data from file descriptor, with error checking.  Returns number of bytes
 * successfully read; 0 implies stream is at end-of-file.  Aborts on error.  */
static size_t xread(void *_fp, void *buf, size_t count)
{
	int fd = (int)(intptr_t)_fp;
	ssize_t ret;

retry:
	ret = read(fd, buf, count);
	if (ret >= 0)
		return ret;
	else if (errno == EINTR)
		goto retry;
	else
		fatal_error_with_errno("Error reading input file");
}

static void xclose(void *_fp)
{
	close((int)(intptr_t)_fp);
}

/* Initializes an input stream to read lines from the file specified by
 * @filename.
 *
 * Gzip files are auto-detected.
 */
void init_input_stream(struct input_stream *in, const char *filename)
{
	unsigned char magic[2] = {0, 0};
	const size_t bufsize = 32768;
	FILE *tmp_fp;

	/* Test for gzip magic bytes { 0x1f, 0x8b}  */

	tmp_fp = xfopen(filename, "rb");
	fread(magic, sizeof(magic[0]), ARRAY_LEN(magic), tmp_fp);
	fclose(tmp_fp);
	if (magic[0] == 0x1f && magic[1] == 0x8b) {
		/* GZIP file  */
		in->fp = xgzopen(filename, "rb");
		in->read = xgzread;
		in->close = (void(*)(void*))gzclose;
	} else {
		/* Other file (hopefully uncompressed...)  */
		in->fp = (void*)(intptr_t)open(filename, O_RDONLY | O_BINARY);
		in->read = xread;
		in->close = xclose;
	}

	/* Allocate internal buffer  */
	in->buf_begin     = xmalloc(bufsize);
	in->buf_end       = in->buf_begin + bufsize;
	in->buf_cur_begin = in->buf_begin;
	in->buf_cur_end   = in->buf_begin;
}

/* Close an input stream initialized with init_input_stream()  */
void destroy_input_stream(struct input_stream *in)
{
	(*in->close)(in->fp);
	free(in->buf_begin);
#ifndef NDEBUG
	memset(in, 0xfd, sizeof(*in));
#endif
}

/* Read a line from an input stream.  Semantics are like getline(), but aborts
 * on read error.  */
ssize_t input_stream_getline(struct input_stream *in, char **lineptr, size_t *n)
{
	/* offset = number of bytes copied to *lineptr buffer, excluding
	 *          terminating null byte  */
	size_t offset = 0;

	for (;;) {
		size_t navail;
		char *nl_ptr;
		size_t copysize;

		navail = in->buf_cur_end - in->buf_cur_begin;

		if (navail == 0) {
			/* No more data in internal buffer; try to fill it  */

			in->buf_cur_begin = in->buf_begin;
			navail = (*in->read)(in->fp,
					     in->buf_cur_begin,
					     in->buf_end - in->buf_cur_begin);
			in->buf_cur_end = in->buf_cur_begin + navail;

			if (navail == 0)  /* At end-of-file  */
				break;
		}

		/* Find the first newline in the internal buffer.  If found,
		 * copy up to and including the newline, then return (break
		 * loop).  If not found, copy all the data and try to read more
		 * (continue loop).  */
		nl_ptr = memchr(in->buf_cur_begin, '\n', navail);
		if (nl_ptr)
			copysize = nl_ptr - in->buf_cur_begin + 1;
		else
			copysize = navail;

		if (*n < offset + copysize + 1) {
			*n = max(*n * 3 / 2, offset + copysize + 1);
			*n = max(*n, 128);
			*lineptr = xrealloc(*lineptr, *n);
		}

		memcpy(*lineptr + offset, in->buf_cur_begin, copysize);
		offset += copysize;
		in->buf_cur_begin += copysize;
		if (in->buf_cur_begin[-1] == '\n')
			break;
	}

	if (offset == 0)
		return -1;

	(*lineptr)[offset] = '\0';
	return offset;
}
