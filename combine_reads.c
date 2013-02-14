/*
 * combine_reads.c:  This file contains the code implementing the core algorithm
 * to combine reads in FLASH.
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013 Eric Biggers
 *
 * This file is part of FLASH, a fast tool to merge overlapping paired-end
 * reads.
 *
 * FLASH is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FLASH; if not, see http://www.gnu.org/licenses/.
 */

#include "combine_reads.h"
#include "fastq.h"
#include "util.h"


/*
 * Align the read @seq_1 to the read @seq_2, where @overlap_len is the length of
 * the attempted overlap and the length of @seq_1, and @seq_2_len is the length
 * of @seq_2.
 *
 * In the output paramater @mismatch_density_ret, return the density of base
 * pair mismatches in the overlapped region.  In the output parameter
 * @qual_score_ret, return the average of the quality scores of all mismatches
 * in the overlapped region.
 *
 * The maximum overlap length considered in the scoring is @max_overlap, while
 * the minimum is @min_overlap.
 */
static void align_position(const char * restrict seq_1,
			   const char * restrict seq_2,
			   const char * restrict qual_1,
			   const char * restrict qual_2,
			   int overlap_len,
			   int min_overlap,
			   int max_overlap,
			   float *mismatch_density_ret,
			   float *qual_score_ret)
{
	int num_non_dna_chars = 0;
	int num_mismatches = 0;
	int mismatch_qual_total = 0;

	/* Calculate the number of mismatches, the number of 'N' characters, and
	 * the sum of the minimum quality values in the mismatched bases over
	 * the overlap region. */
	for (int i = 0; i < overlap_len; i++) {
		if (seq_1[i] == 'N' || seq_2[i] == 'N') {
			num_non_dna_chars++;
		} else {
			if (seq_1[i] != seq_2[i])  {
				num_mismatches++;
				mismatch_qual_total += min(qual_1[i],
							   qual_2[i]);
			}
		}
	}

	/* Reduce the length of overlap by the number of N's. */
	overlap_len -= num_non_dna_chars;

	/* Set the qual_score and mismatch_density.  qual_score is for
	 * mismatches only, so lower qual_score is preferable since lower
	 * qual_score means less confidence in base calls at mismatches. */
	if (overlap_len > max_overlap) {
		*qual_score_ret = (float)mismatch_qual_total / max_overlap;
		*mismatch_density_ret = (float)num_mismatches / max_overlap;
	} else  {
		if (overlap_len >= min_overlap) {
			*qual_score_ret = (float)mismatch_qual_total / overlap_len;
			*mismatch_density_ret = (float)num_mismatches / overlap_len;
		} else {
			*qual_score_ret = 0.0f;
			*mismatch_density_ret = 10000.0f;
		}
	}
}


/*
 * Return the position of the first overlapping character in the first read
 * (0-based) of the best overlap between two reads, or -1 if the best overlap
 * does not satisfy the required threshold.
 */
static int pair_align(const struct read *read_1, const struct read *read_2,
		      int min_overlap, int max_overlap,
		      float max_mismatch_density)
{
	/* Best (smallest) mismatch density that has been found so far in an
	 * overlap. */
	float best_mismatch_density = max_mismatch_density + 1.0f;
	float best_qual_score = 0.0f;
	int best_position = -1;

	/* Require at least min_overlap bases overlap, and require that the
	 * second read is not overlapped such that it is completely contained in
	 * the first read.  */
	int start = max(0, read_1->seq_len - read_2->seq_len);
	int end = read_1->seq_len - min_overlap + 1;
	for (int i = start; i < end; i++) {
		float mismatch_density;
		float qual_score;

		align_position(read_1->seq + i,
			       read_2->seq,
			       read_1->qual + i,
			       read_2->qual,
			       read_1->seq_len - i,
			       min_overlap, max_overlap,
			       &mismatch_density, &qual_score);

		if (mismatch_density < best_mismatch_density) {
			best_qual_score       = qual_score;
			best_mismatch_density = mismatch_density;
			best_position         = i;
		} else if (mismatch_density == best_mismatch_density) {
			if (qual_score < best_qual_score) {
				best_qual_score = qual_score;
				best_position = i;
			}
		}
	}

	if (best_mismatch_density > max_mismatch_density)
		return -1;
	else
		return best_position;
}

/* This is the entry point for the core algorithm of FLASH.  The following
 * function attempts to combine @read_1 with @read_2, and writes the result into
 * @combined_read.  %true is returned iff combination was successful.
 *
 * Note: @read_2 is provided to this function after having been
 * reverse-complemented.  Hence, the code just aligns the reads in the forward
 * orientation, which is equivalent to aligning the original reads in the
 * desired reverse-complement orientation.
 *
 * Please see the help output of FLASH for the description of the min_overlap,
 * max_overlap, and max_mismatch_density parameters.  (--min-overlap,
 * --max-overlap, and --max-mismatch-density on the command line).
 *
 * You may also want to read the original FLASH publication for a description of
 * the algorithm used here:
 *
 *  Title:   FLASH: fast length adjustment of short reads to improve genome assemblies
 *  Authors: Tanja Magoč and Steven L. Salzberg
 *  URL:     http://bioinformatics.oxfordjournals.org/content/27/21/2957.full
 *
 */
bool combine_reads(const struct read *read_1, const struct read *read_2,
		   struct read *combined_read, int min_overlap,
		   int max_overlap, float max_mismatch_density)
{
	/* Starting position of the alignment in the first read, 0-based. */
	int overlap_begin;

	/* Length of the overlapping part of two reads. */
	int overlap_len;

	/* Length of the part of the second read not overlapped with the first
	 * read. */
	int remaining_len;

	/* Length of the combined read. */
	int combined_seq_len;

	const char * restrict seq_1 = read_1->seq;
	const char * restrict seq_2 = read_2->seq;
	const char * restrict qual_1 = read_1->qual;
	const char * restrict qual_2 = read_2->qual;
	char * restrict combined_seq;
	char * restrict combined_qual;

	/* Do the alignment. */
	overlap_begin = pair_align(read_1, read_2, min_overlap, max_overlap,
				   max_mismatch_density);

	/* If no alignment found, return false */
	if (overlap_begin < 0)
		return false;

	/* Fill in the combined read. */

	overlap_len = read_1->seq_len - overlap_begin;
	remaining_len = read_2->seq_len - overlap_len;
	combined_seq_len = read_1->seq_len + remaining_len;

	if (combined_read->seq_bufsz < combined_seq_len) {
		combined_read->seq = xrealloc(combined_read->seq,
					      combined_seq_len);
		combined_read->qual = xrealloc(combined_read->qual,
					       combined_seq_len);
		combined_read->seq_bufsz = combined_seq_len;
	}

	combined_seq = combined_read->seq;
	combined_qual = combined_read->qual;
	combined_read->seq_len = combined_seq_len;

	/* Copy the beginning of the first read. */
	while (overlap_begin--) {
		*combined_seq++ = *seq_1++;
		*combined_qual++ = *qual_1++;
	}

	/* Copy the overlapping part. */
	while (overlap_len--) {
		if (*seq_1 == *seq_2) {
			/* Same base in both reads.  Take the higher quality
			 * value. */
			*combined_seq = *seq_1;
			*combined_qual = max(*qual_1, *qual_2);
		} else {
			/* Different bases in the two reads.  Take the base from
			 * the read that has the higher quality value, but use
			 * the lower quality value, and use a quality value of
			 * at most 2 (+ phred_offset in the final output--- here
			 * the quality values are all scaled to start at 0.) */

			*combined_qual = min(*qual_1, *qual_2);
			*combined_qual = min(*combined_qual, 2);

			if (*qual_1 > *qual_2) {
				*combined_seq = *seq_1;
			} else if (*qual_1 < *qual_2) {
				*combined_seq = *seq_2;
			} else {
				/* Same quality value; take the base from the
				 * first read if the base from the second read
				 * is an 'N'; otherwise take the base from the
				 * second read. */
				if (*seq_2 == 'N')
					*combined_seq = *seq_1;
				else
					*combined_seq = *seq_2;
			}
		}
		combined_seq++;
		combined_qual++;
		seq_1++;
		seq_2++;
		qual_1++;
		qual_2++;
	}
	/* Copy the part of the second read in the mate pair that is not in the
	 * overlapped region. */
	while (remaining_len--) {
		*combined_seq++ = *seq_2++;
		*combined_qual++ = *qual_2++;
	}
	return true;
}
