#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
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
FILE *xfopen(const char *filename, const char *mode)
{ 
	FILE *fp = fopen(filename, mode);
	if (!fp) {
		fatal_error("Could not open the file \"%s\" for %s: %s",
			    filename, 
			    strchr(mode, 'w') ? "writing" : "reading",
			    strerror(errno));
	}
	return fp;
}

void xfclose(FILE *fp)
{
	if (fclose(fp) != 0)
		fatal_error("Failed to close output file");
}

/* Like `mkdir -p': create directory, and all parent directories, as needed,
 * failing only if a needed directory cannot be created. */
void mkdir_p(const char *dir)
{
	size_t len = strlen(dir);
	char dir_copy[len + 1];
	char *p = dir_copy;
	for (size_t i = 0; i < len; i++) {
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

			if (mkdir(dir_copy, 0755) != 0 && errno != EEXIST)
				fatal_error("Cannot create directory \"%s\": "
					    "%s", dir_copy, strerror(errno));
			*p = orig_char;
		}
	} while (*p++ != '\0');
}
