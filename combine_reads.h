#ifndef _FLASH_COMBINE_READS_H_
#define _FLASH_COMBINE_READS_H_

#include <stdbool.h>

struct read;

struct combine_params {
	int min_overlap;
	int max_overlap;
	float max_mismatch_density;
	bool cap_mismatch_quals;
};

extern bool
combine_reads(const struct read *read_1, const struct read *read_2,
	      struct read *combined_read,
	      const struct combine_params *params);

#endif /* _FLASH_COMBINE_READS_H_  */
