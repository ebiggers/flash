#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <zlib.h>
#include "util.h"

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

/* Prints an error message and exits the program with failure status. */
void fatal_error(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);
	fflush(stdout);
	fputs("ERROR: ", stderr);
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
	fputs("ERROR: ", stderr);
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
	fputs("WARNING: ", stderr);
	fprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
}

#ifdef MULTITHREADED
/* Returns the number of processors (if it can be determined), otherwise returns
 * 1. */
unsigned get_default_num_threads()
{
	long nthreads = sysconf(_SC_NPROCESSORS_ONLN);
	if (nthreads == -1) {
		warning("Could not deteremine number of processors! Assuming 1");
		return 1;
	} else {
		return (unsigned)nthreads;
	}
}
#endif

/* malloc(), exiting the program with failure status if memory allocation fails.
 * */
void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p) {
		fatal_error("Out of memory: tried to allocate %lu bytes",
			     size);
	}
	return p;
}

/* fopen(), exiting the program with failure status if file open fails. */
void *xfopen(const char *filename, const char *mode)
{ 
	if (strcmp(filename, "-") == 0)
		return (void*)stdout;

	FILE *fp = fopen(filename, mode);
	if (!fp) {
		fatal_error_with_errno("Could not open the file \"%s\" for %s",
			    filename, 
			    strchr(mode, 'w') ? "writing" : "reading");
	}
	return (void*)fp;
}

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
		const char *err_str;
		int errnum;
		err_str = gzerror(f, &errnum);
		if (errnum == Z_ERRNO)
			fatal_error_with_errno("Failed to open the file \"%s\"",
					       filename);
		else
			fatal_error("zlib error opening \"%s\": %s",
				    filename, err_str);
	}
	return (void*)f;
}

void xfclose(void *fp)
{
	if (fp && fclose((FILE*)fp) != 0)
		fatal_error_with_errno("Failed to close output file");
}

void xgzclose(void *fp)
{
	if (fp && gzclose((gzFile)fp) != Z_OK)
		fatal_error_with_errno("Failed to close output file");
}

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
