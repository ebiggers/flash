#ifndef _COMBINE_H
#define _COMBINE_H

#include "fastq.h"

extern int combine_reads(const struct mate_pair *pair, char combined_seq[], 
			 char combined_qual[], int min_overlap, 
			 int max_overlap, float max_mismatch_density);
#endif
