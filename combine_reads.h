#ifndef _COMBINE_H
#define _COMBINE_H

struct read;

extern int combine_reads(const struct read *read_1, const struct read *read_2,
			 char combined_seq[], char combined_qual[],
			 int min_overlap, int max_overlap,
			 float max_mismatch_density);
#endif
