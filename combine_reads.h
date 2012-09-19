#ifndef _COMBINE_H
#define _COMBINE_H

#include <stdbool.h>

struct read;

extern bool combine_reads(const struct read *read_1, const struct read *read_2,
			  struct read *combined_read, int min_overlap,
			  int max_overlap, float max_mismatch_density);
#endif
