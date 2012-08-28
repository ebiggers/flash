#ifndef _UTIL_H
#define _UTIL_H

#include <stddef.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define ZERO_ARRAY(A) (memset((A), 0, sizeof(A)))
#define ARRAY_LEN(A) (sizeof(A) / sizeof(A[0]))

extern void fatal_error(const char *msg, ...) __attribute__((noreturn,cold));
extern void warning(const char *msg, ...) __attribute__((cold));
extern void *xmalloc(size_t size);
extern FILE *xfopen(const char *filename, const char *mode);
extern void xfclose(FILE *fp);
extern void mkdir_p(const char *dir);

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
static inline size_t trim(char *s)
{
	size_t len = strlen(s);
	while (len != 0 && isspace(s[len - 1]))
		s[--len] = '\0';
	return len;
}

extern const char canonical_ascii_tab[];
extern const char complement_tab[];

/* Turns lowercase a, c, g, t into uppercase; 
 * uppercase A, C, G, T stay the same;
 * everything else turns into 'N'.  */
static inline char canonical_ascii_char(char c)
{
	return canonical_ascii_tab[(unsigned char)c];
}

/* Complements a canonical ASCII base (A, C, G, T, N). */
static inline char complement(char c)
{
	return complement_tab[(unsigned char)c];
}

/* Reverse a string. */
static inline void reverse(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p < pp) {
		char tmp = *p;
		*p++ = *pp;
		*pp-- = tmp;
	}
}

/* Reverse complement an ASCII DNA sequence. */
static inline void reverse_complement(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p < pp) {
		char tmp = *p;
		*p++ = complement(*pp);
		*pp-- = complement(tmp);
	}
	if (p == pp) /* Odd sequence length. */
		*p = complement(*p);
}

#endif
