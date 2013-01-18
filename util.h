#ifndef _UTIL_H
#define _UTIL_H

#include <ctype.h>
#include <string.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define ARRAY_LEN(A) (sizeof(A) / sizeof((A)[0]))
#define ZERO_ARRAY(A) (memset((A), 0, sizeof(A)))
#define ZERO_OBJECT(obj) (memset((&obj), 0, sizeof(obj)))

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
# 		define __cold __attribute__((cold))
#	else
#		define __cold
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#else
#	define __noreturn
#	define __cold
#	define __format
#endif

struct read;

/* These are only for writing files; reading files currently always uses a
 * gzFile and its associated functions. */
typedef void* (*open_file_t)(const char *filename, const char *mode);
typedef void (*close_file_t)(void *);
typedef void (*write_read_t)(struct read *read, void *fp, int phred_offset);

struct file_operations {
	char *name;
	char *suffix;
	open_file_t  open_file;
	close_file_t close_file;
	write_read_t write_read;
};

extern struct file_operations gzip_fops;
extern struct file_operations normal_fops;
extern struct file_operations pipe_fops;
extern char *compress_prog;
extern char *compress_prog_args;

extern void fatal_error(const char *msg, ...) \
			__noreturn __cold __format(printf, 1, 2);
extern void fatal_error_with_errno(const char *msg, ...) \
			__noreturn __cold __format(printf, 1, 2);
extern void warning(const char *msg, ...) __cold __format(printf, 1, 2);
extern void info(const char *msg, ...) __cold __format(printf, 1, 2);
extern void *xmalloc(size_t size) __cold;
extern void *xrealloc(void *ptr, size_t size) __cold;
extern unsigned get_default_num_threads() __cold;
extern void mkdir_p(const char *dir) __cold;


extern void *xfopen(const char *filename, const char *mode) __cold;
extern void *xpopen(const char *filename, const char *mode) __cold;
extern void *xgzopen(const char *filename, const char *mode) __cold;
extern void xfclose(void *fp) __cold;
extern void xgzclose(void *fp) __cold;
extern void xpclose(void *fp) __cold;

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
static inline size_t trim(char *s, size_t len)
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
		*p = *pp;
		*pp = tmp;
		p++;
		pp--;
	}
}

/* Reverse complement an ASCII DNA sequence. */
static inline void reverse_complement(char *p, size_t len)
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

#endif
