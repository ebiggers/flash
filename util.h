#ifndef _FLASH_UTIL_H_
#define _FLASH_UTIL_H_

#include <ctype.h>
#include <stddef.h>
#include <pthread.h>
#include <stdio.h>

#define ARRAY_LEN(A) (sizeof(A) / sizeof((A)[0]))

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
# 		define __cold __attribute__((cold))
#	else
#		define __cold
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#	define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#	define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })
#else
#	define __noreturn
#	define __cold
#	define __format
#	define max(a,b) (((a) > (b)) ? (a) : (b))
#	define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

extern void
fatal_error(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);

extern void
fatal_error_with_errno(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);

extern void
warning(const char *msg, ...) __cold __format(printf, 1, 2);

extern FILE *infofile;

extern void
info(const char *msg, ...) __format(printf, 1, 2);

extern void *
xmalloc(size_t size);

#ifdef NDEBUG
#  define xfree(p, size) free(p)
#else
extern void
xfree(void *p, size_t size);
#endif

extern void *
xzalloc(size_t size);

extern char *
xstrdup(const char *str);

extern void *
xrealloc(void *ptr, size_t size);

extern unsigned
get_default_num_threads(void);

extern void
mkdir_p(const char *dir);

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
static inline size_t
trim(char *s, size_t len)
{
	while (len != 0 && isspace(s[len - 1]))
		s[--len] = '\0';
	return len;
}

extern const char canonical_ascii_tab[];
extern const char complement_tab[];

/* Turns lowercase a, c, g, t into uppercase;
 * uppercase A, C, G, T stay the same;
 * everything else turns into 'N'.  */
static inline char
canonical_ascii_char(char c)
{
	return canonical_ascii_tab[(unsigned char)c];
}

/* Complements a canonical ASCII base (A, C, G, T, N). */
static inline char
complement(char c)
{
	return complement_tab[(unsigned char)c];
}

/* Reverse a string. */
static inline void
reverse(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p < pp) {
		char tmp = *p;
		*p = *pp;
		*pp = tmp;
		p++;
		pp--;
	}
}

/* Reverse complement an ASCII DNA sequence. */
static inline void
reverse_complement(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p <= pp) {
		char tmp = *p;
		*p = complement(*pp);
		*pp = complement(tmp);
		p++;
		pp--;
	}
}

extern pthread_t
create_thread(void *(*proc)(void *), void *params);

extern void
join_thread(pthread_t t);

#endif /* _FLASH_UTIL_H_ */
