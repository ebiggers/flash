/*
 * util.c: miscellaneous useful functions for FLASH
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
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
#include <limits.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "util.h"

#ifdef __WIN32__
    /* Get the pthread mutex declarations as a replacement for flockfile() and
     * funlockfile(). */
#  include <pthread.h>
    /* Get the GetSystemInfo() declaration as replacement for
     * sysconf(_SC_NPROCESSORS_ONLN). */
#  include <windows.h>
#endif

/*
 * A table mapping ASCII characters (actually, 8-bit bytes) to "canonical" form
 * for base processing:
 *
 * a, A => A
 * c, C => C
 * g, G => G
 * t, T => T
 * everything else => N
 */
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

/* A table mapping an ASCII base character in canonical form to its complement.
 * Unknown bases (N's) stay unknown.  */
const char complement_tab[] = {
	['A'] = 'T',
	['C'] = 'G',
	['G'] = 'C',
	['T'] = 'A',
	['N'] = 'N',
};


#ifdef __WIN32__
static pthread_mutex_t stdout_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t stderr_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

static void
lock_stdout(void)
{
#ifdef __WIN32__
	pthread_mutex_lock(&stdout_lock);
#else
	flockfile(stdout);
#endif
}

static void
lock_stderr(void)
{
#ifdef __WIN32__
	pthread_mutex_lock(&stderr_lock);
#else
	flockfile(stderr);
#endif
}

static void
unlock_stdout(void)
{
#ifdef __WIN32__
	pthread_mutex_unlock(&stdout_lock);
#else
	funlockfile(stdout);
#endif
}

static void
unlock_stderr(void)
{
#ifdef __WIN32__
	pthread_mutex_unlock(&stderr_lock);
#else
	funlockfile(stderr);
#endif
}

#define PROGRAM_TAG "[FLASH] "

/* Prints an error message and exits the program with failure status.  */
void
fatal_error(const char *msg, ...)
{
	va_list va;

	lock_stderr();

	va_start(va, msg);
	fflush(stdout);
	fputs(PROGRAM_TAG "ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	exit(1);
}

/* Prints an error message, with added text for errno if it is nonzero, and
 * exits the program with failure status.  */
void
fatal_error_with_errno(const char *msg, ...)
{
	va_list va;

	lock_stderr();

	va_start(va, msg);
	fflush(stdout);
	fputs(PROGRAM_TAG "ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	if (errno)
		fprintf(stderr, ": %s\n", strerror(errno));
	else
		putc('\n', stderr);
	va_end(va);
	exit(1);
}

/* Prints a warning message.  */
void
warning(const char *msg, ...)
{
	va_list va;

	lock_stderr();

	va_start(va, msg);
	fputs(PROGRAM_TAG "WARNING: ", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);

	unlock_stderr();
}

/* Prints an informational message.  */
void
info(const char *msg, ...)
{
	va_list va;

	lock_stdout();

	va_start(va, msg);
	fputs(PROGRAM_TAG, stdout);
	vprintf(msg, va);
	putchar('\n');
	fflush(stdout);
	va_end(va);

	unlock_stdout();
}

/* Like malloc(), but aborts if out of memory, and always returns non-NULL, even
 * if 0 bytes were requested.  */
void *
xmalloc(size_t size)
{
	void *p = malloc(size);
	if (p)
		return p;
	if (!size) {
		p = malloc(1);
		if (p)
			return p;
	}
	fatal_error("Out of memory: tried to allocate %zu bytes", size);
}

void *
xzalloc(size_t size)
{
	return memset(xmalloc(size), 0, size);
}

/* Like strdup(), but aborts if out of memory.  */
char *
xstrdup(const char *str)
{
	return strcpy(xmalloc(strlen(str) + 1), str);
}

/* Like realloc(), but aborts if out of memory, and always returns non-NULL,
 * even if 0 bytes were requested.  */
void *
xrealloc(void *ptr, size_t size)
{
	void *p = realloc(ptr, size);
	if (p)
		return p;
	if (!size) {
		p = malloc(1);
		if (p)
			return p;
	}
	fatal_error("Out of memory: tried to reallocate %zu bytes", size);
}

#ifndef NDEBUG
void
xfree(void *p, size_t size)
{
	if (p) {
		memset(p, 0xfd, size);
		free(p);
	}
}
#endif

/* Returns the number of available processors if it can be determined.
 * Otherwise returns 1.  */
unsigned
get_default_num_threads(void)
{
#ifdef __WIN32__
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	if (si.dwNumberOfProcessors > 0 && si.dwNumberOfProcessors <= UINT_MAX)
		return si.dwNumberOfProcessors;
#else
	long nproc = sysconf(_SC_NPROCESSORS_ONLN);
	if (nproc > 0 && nproc <= UINT_MAX)
		return nproc;
#endif
	warning("Could not determine number of processors! Assuming 1");
	return 1;
}


/* Returns true if the specified character is a path separator on the current
 * platform.  */
static bool
is_path_separator(char c)
{
#ifdef __WIN32__
	return (c == '/') || (c == '\\');
#else
	return (c == '/');
#endif
}


/* mkdir() on Windows doesn't take a mode argument.  */
#ifdef __WIN32__
#  define mkdir(path, mode) mkdir(path)
#endif

/* Like `mkdir -p': create the specified directory, and all parent directories,
 * as needed, failing only if a needed directory cannot be created.  */
void
mkdir_p(const char *dir)
{
	size_t len = strlen(dir);
	char dir_copy[len + 1];
	char *p = dir_copy;
	/* Copy the directory name to the @dir_copy array, squashing together
	 * consecutive path separators.  */
	for (size_t i = 0; i < len; i++) {
		if (!is_path_separator(dir[i]) ||
		    !is_path_separator(dir[i + 1]))
			*p++ = dir[i];
	}
	*p = '\0';

	p = dir_copy;
	do {
		if (p != dir_copy && (*p == '\0' || is_path_separator(*p))) {
			char orig_char = *p;
			*p = '\0';

			if (mkdir(dir_copy, 0755) != 0 && errno != EEXIST) {
				fatal_error_with_errno("Failed to create "
						       "directory \"%s\"",
						       dir_copy);
			}
			*p = orig_char;
		}
	} while (*p++ != '\0');
}

pthread_t
create_thread(void *(*proc)(void *), void *params)
{
	int result;
	pthread_t t;

	result = pthread_create(&t, NULL, proc, params);
	if (result) {
		errno = result;
		fatal_error_with_errno("Failed to create new thread");
	}
	return t;
}

void
join_thread(pthread_t t)
{
	int result = pthread_join(t, NULL);
	if (result) {
		errno = result;
		fatal_error_with_errno("Failed to join thread");
	}
}
