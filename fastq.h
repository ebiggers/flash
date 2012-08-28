#ifndef _FASTQ_H
#define _FASTQ_H

#include <zlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* Read from a FASTQ file. */
struct read {
	int tag_bufsz;
	int read_bufsz;
	int seq_len;
	char *tag;
	char *seq;
	char *separator;
	char *qual;
};

/* Mate pair: two reads on opposite strands in opposite directions, oriented
 * towards each other. */
struct mate_pair {
	struct read read_1;
	struct read read_2;
};


extern void init_mate_pair(struct mate_pair *pair, int max_tag_len, 
			   int max_read_len);

extern bool next_mate_pair(struct mate_pair *pair, gzFile mates1_g, 
			   gzFile mates2_g, int phred_offset);

extern void write_read(struct read *read, FILE *fp, int phred_offset);

#endif
