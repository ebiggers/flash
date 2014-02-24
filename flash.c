/*
 * flash.c:  parse arguments and set up and run the FLASH pipeline.
 *
 * Please see combine_reads.c if you are looking for the core algorithm used to
 * combine reads in FLASH.
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
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

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "combine_reads.h"
#include "iostream.h"
#include "read.h"
#include "read_io.h"
#include "read_queue.h"
#include "util.h"

#define VERSION_STR "v1.2.9"

#ifdef __WIN32__
#  define PAGER "more"
#else
#  define PAGER "less"
#endif

static void usage(const char *argv0)
{
	const char *usage_str =
"Usage: flash [OPTIONS] MATES_1.FASTQ MATES_2.FASTQ\n"
"       flash [OPTIONS] --interleaved-input (MATES.FASTQ | -)\n"
"       flash [OPTIONS] --tab-delimited-input (MATES.TAB | -)\n"
"\n"
"----------------------------------------------------------------------------\n"
"                                 DESCRIPTION                                \n"
"----------------------------------------------------------------------------\n"
"\n"
"FLASH (Fast Length Adjustment of SHort reads) is an accurate and fast tool\n"
"to merge paired-end reads that were generated from DNA fragments whose\n"
"lengths are shorter than twice the length of reads.  Merged read pairs result\n"
"in unpaired longer reads, which are generally more desired in genome\n"
"assembly and genome analysis processes.\n"
"\n"
"Briefly, the FLASH algorithm considers all possible overlaps at or above a\n"
"minimum length between the reads in a pair and chooses the overlap that\n"
"results in the lowest mismatch density (proportion of mismatched bases in\n"
"the overlapped region).  Ties between multiple overlaps are broken by\n"
"considering quality scores at mismatch sites.  When building the merged\n"
"sequence, FLASH computes a consensus sequence in the overlapped region.\n"
"More details can be found in the original publication\n"
"(http://bioinformatics.oxfordjournals.org/content/27/21/2957.full).\n"
"\n"
"Limitations of FLASH include:\n"
"   - FLASH cannot merge paired-end reads that do not overlap.\n"
"   - FLASH cannot merge read pairs that have an outward orientation, either\n"
"     due to being \"jumping\" reads or due to excessive trimming.\n"
"   - FLASH is not designed for data that has a significant amount of indel\n"
"     errors (such as Sanger sequencing data).  It is best suited for Illumina\n"
"     data.\n"
"\n"
"----------------------------------------------------------------------------\n"
"                               MANDATORY INPUT\n"
"----------------------------------------------------------------------------\n"
"\n"
"The most common input to FLASH is two FASTQ files containing read 1 and read 2\n"
"of each mate pair, respectively, in the same order.\n"
"\n"
"Alternatively, you may provide one FASTQ file, which may be standard input,\n"
"containing paired-end reads in either interleaved FASTQ (see the\n"
"--interleaved-input option) or tab-delimited (see the --tab-delimited-input\n"
"option) format.  In all cases, gzip compressed input is autodetected.  Also,\n"
"in all cases, the PHRED offset is, by default, assumed to be 33; use the\n"
"--phred-offset option to change it.\n"
"\n"
"----------------------------------------------------------------------------\n"
"                                   OUTPUT\n"
"----------------------------------------------------------------------------\n"
"\n"
"The default output of FLASH consists of the following files:\n"
"\n"
"   - out.extendedFrags.fastq      The merged reads.\n"
"   - out.notCombined_1.fastq      Read 1 of mate pairs that were not merged.\n"
"   - out.notCombined_2.fastq      Read 2 of mate pairs that were not merged.\n"
"   - out.hist                     Numeric histogram of merged read lengths.\n"
"   - out.histogram                Visual histogram of merged read lengths.\n"
"\n"
"FLASH also logs informational messages to standard output.  These can be\n"
"redirected to a file, as in the following example:\n"
"\n"
"  $ flash reads_1.fq reads_2.fq | tee flash.log\n"
"\n"
"In addition, FLASH supports several features affecting the output:\n"
"\n"
"   - Writing the merged reads directly to standard output (--to-stdout)\n"
"   - Writing gzip compressed output files (-z) or using an external\n"
"     compression program (--compress-prog)\n"
"   - Writing the uncombined read pairs in interleaved FASTQ format\n"
"     (--interleaved-output)\n"
"   - Writing all output reads to a single file in tab-delimited format\n"
"     (--tab-delimited-output)\n"
"\n"
"----------------------------------------------------------------------------\n"
"                                   OPTIONS\n"
"----------------------------------------------------------------------------\n"
"\n"
"  -m, --min-overlap=NUM   The minimum required overlap length between two\n"
"                          reads to provide a confident overlap.  Default:\n"
"                          10bp.\n"
"\n"
"  -M, --max-overlap=NUM   Maximum overlap length expected in approximately\n"
"                          90% of read pairs.  It is by default set to 70bp,\n"
"                          which works well for 100bp reads generated from a\n"
"                          180bp library, assuming a normal distribution of\n"
"                          fragment lengths.  Overlaps longer than the maximum\n"
"                          overlap parameter are still considered as good\n"
"                          overlaps, but the mismatch density (explained below)\n"
"                          is calculated over the first max_overlap bases in\n"
"                          the overlapped region rather than the entire\n"
"                          overlap.  Default: 70bp, or calculated from the\n"
"                          specified read length, fragment length, and fragment\n"
"                          length standard deviation.\n"
"\n"
"  -x, --max-mismatch-density=NUM\n"
"                          Maximum allowed ratio between the number of\n"
"                          mismatched base pairs and the overlap length.\n"
"                          Two reads will not be combined with a given overlap\n"
"                          if that overlap results in a mismatched base density\n"
"                          higher than this value.  Note: Any occurence of an\n"
"                          'N' in either read is ignored and not counted\n"
"                          towards the mismatches or overlap length.  Our\n"
"                          experimental results suggest that higher values of\n"
"                          the maximum mismatch density yield larger\n"
"                          numbers of correctly merged read pairs but at\n"
"                          the expense of higher numbers of incorrectly\n"
"                          merged read pairs.  Default: 0.25.\n"
"\n"
"  -p, --phred-offset=OFFSET\n"
"                          The smallest ASCII value of the characters used to\n"
"                          represent quality values of bases in FASTQ files.\n"
"                          It should be set to either 33, which corresponds\n"
"                          to the later Illumina platforms and Sanger\n"
"                          platforms, or 64, which corresponds to the\n"
"                          earlier Illumina platforms.  Default: 33.\n"
"\n"
"  -r, --read-len=LEN\n"
"  -f, --fragment-len=LEN\n"
"  -s, --fragment-len-stddev=LEN\n"
"                          Average read length, fragment length, and fragment\n"
"                          standard deviation.  These are convenience parameters\n"
"                          only, as they are only used for calculating the\n"
"                          maximum overlap (--max-overlap) parameter.\n"
"                          The maximum overlap is calculated as the overlap of\n"
"                          average-length reads from an average-size fragment\n"
"                          plus 2.5 times the fragment length standard\n"
"                          deviation.  The default values are -r 100, -f 180,\n"
"                          and -s 18, so this works out to a maximum overlap of\n"
"                          65 bp.  If --max-overlap is specified, then the\n"
"                          specified value overrides the calculated value.\n"
"\n"
"                          If you do not know the standard deviation of the\n"
"                          fragment library, you can probably assume that the\n"
"                          standard deviation is 10% of the average fragment\n"
"                          length.\n"
"\n"
"  --cap-mismatch-quals    Cap quality scores assigned at mismatch locations\n"
"                          to 2.  This was the default behavior in FLASH v1.2.7\n"
"                          and earlier.  Later versions will instead calculate\n"
"                          such scores as max(|q1 - q2|, 2); that is, the\n"
"                          absolute value of the difference in quality scores,\n"
"                          but at least 2.  Essentially, the new behavior\n"
"                          prevents a low quality base call that is likely a\n"
"                          sequencing error from significantly bringing down\n"
"                          the quality of a high quality, likely correct base\n"
"                          call.\n"
"\n"
"  --interleaved-input     Instead of requiring files MATES_1.FASTQ and\n"
"                          MATES_2.FASTQ, allow a single file MATES.FASTQ that\n"
"                          has the paired-end reads interleaved.  Specify \"-\"\n"
"                          to read from standard input.\n"
"\n"
"  --interleaved-output    Write the uncombined pairs in interleaved FASTQ\n"
"                          format.\n"
"\n"
"  -I, --interleaved       Equivalent to specifying both --interleaved-input\n"
"                          and --interleaved-output.\n"
"\n"
"  -Ti, --tab-delimited-input\n"
"                          Assume the input is in tab-delimited format\n"
"                          rather than FASTQ, in the format described below in\n"
"                          '--tab-delimited-output'.  In this mode you should\n"
"                          provide a single input file, each line of which must\n"
"                          contain either a read pair (5 fields) or a single\n"
"                          read (3 fields).  FLASH will try to combine the read\n"
"                          pairs.  Single reads will be written to the output\n"
"                          file as-is if also using --tab-delimited-output;\n"
"                          otherwise they will be ignored.  Note that you may\n"
"                          specify \"-\" as the input file to read the\n"
"                          tab-delimited data from standard input.\n"
"\n"
"  -To, --tab-delimited-output\n"
"                          Write output in tab-delimited format (not FASTQ).\n"
"                          Each line will contain either a combined pair in the\n"
"                          format 'tag <tab> seq <tab> qual' or an uncombined\n"
"                          pair in the format 'tag <tab> seq_1 <tab> qual_1\n"
"                          <tab> seq_2 <tab> qual_2'.\n"
"\n"
"  -o, --output-prefix=PREFIX\n"
"                          Prefix of output files.  Default: \"out\".\n"
"\n"
"  -d, --output-directory=DIR\n"
"                          Path to directory for output files.  Default:\n"
"                          current working directory.\n"
"\n"
"  -c, --to-stdout         Write the combined reads to standard output.  In\n"
"                          this mode, with FASTQ output (the default) the\n"
"                          uncombined reads are discarded.  With tab-delimited\n"
"                          output, uncombined reads are included in the\n"
"                          tab-delimited data written to standard output.\n"
"                          In both cases, histogram files are not written,\n"
"                          and informational messages are sent to standard\n"
"                          error rather than to standard output.\n"
"\n"
"  -z, --compress          Compress the output files directly with zlib,\n"
"                          using the gzip container format.  Similar to\n"
"                          specifying --compress-prog=gzip and --suffix=gz,\n"
"                          but may be slightly faster.\n"
"\n"
"  --compress-prog=PROG    Pipe the output through the compression program\n"
"                          PROG, which will be called as `PROG -c -',\n"
"                          plus any arguments specified by --compress-prog-args.\n"
"                          PROG must read uncompressed data from standard input\n"
"                          and write compressed data to standard output when\n"
"                          invoked as noted above.\n"
"                          Examples: gzip, bzip2, xz, pigz.\n"
"\n"
"  --compress-prog-args=ARGS\n"
"                          A string of additional arguments that will be passed\n"
"                          to the compression program if one is specified with\n"
"                          --compress-prog=PROG.  (The arguments '-c -' are\n"
"                          still passed in addition to explicitly specified\n"
"                          arguments.)\n"
"\n"
"  --suffix=SUFFIX, --output-suffix=SUFFIX\n"
"                          Use SUFFIX as the suffix of the output files\n"
"                          after \".fastq\".  A dot before the suffix is assumed,\n"
"                          unless an empty suffix is provided.  Default:\n"
"                          nothing; or 'gz' if -z is specified; or PROG if\n"
"                          --compress-prog=PROG is specified.\n"
"\n"
"  -t, --threads=NTHREADS  Set the number of worker threads.  This is in\n"
"                          addition to the I/O threads.  Default: number of\n"
"                          processors.  Note: if you need FLASH's output to\n"
"                          appear deterministically or in the same order as\n"
"                          the original reads, you must specify -t 1\n"
"                          (--threads=1).\n"
"\n"
"  -q, --quiet             Do not print informational messages.\n"
"\n"
"  -h, --help              Display this help and exit.\n"
"\n"
"  -v, --version           Display version.\n"
	;
	fputs(usage_str, stdout);
	if (isatty(STDOUT_FILENO)) {
		/* Just to be extra user-friendly... */
		printf("\nRun `%s --help | "PAGER"' to "
		       "prevent this text from scrolling by.\n", argv0);
	}
}

static void usage_short(const char *argv0)
{
	fprintf(stderr,
		"Usage: flash [OPTIONS] MATES_1.FASTQ MATES_2.FASTQ\n"
		"Run `%s --help | "PAGER"' for more information.\n", argv0);
}

static void version(void)
{
	fputs(
"FLASH "VERSION_STR"\n"
"Copyright (C) 2012 Tanja Magoc\n"
"Copyright (C) 2012, 2013, 2014 Eric Biggers\n"
"License GPLv3+; GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n"
"This is free software: you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n"
"\n"
"Report bugs to flash.comment@gmail.com or https://sourceforge.net/p/flashpage/bugs\n"
	, stdout);
}

enum {
	INTERLEAVED_INPUT_OPTION = 257,
	INTERLEAVED_OUTPUT_OPTION,
	CAP_MISMATCH_QUALS_OPTION,
	COMPRESS_PROG_OPTION,
	COMPRESS_PROG_ARGS_OPTION,
	SUFFIX_OPTION,
	TAB_DELIMITED_INPUT_OPTION,
	TAB_DELIMITED_OUTPUT_OPTION,
};

static const char *optstring = "m:M:x:p:r:f:s:IT:o:d:czt:qhv";
static const struct option longopts[] = {
	{"min-overlap",          required_argument,  NULL, 'm'},
	{"max-overlap",          required_argument,  NULL, 'M'},
	{"max-mismatch-density", required_argument,  NULL, 'x'},
	{"phred-offset",         required_argument,  NULL, 'p'},
	{"read-len",             required_argument,  NULL, 'r'},
	{"fragment-len",         required_argument,  NULL, 'f'},
	{"fragment-len-stddev",  required_argument,  NULL, 's'},
	{"cap-mismatch-quals",   no_argument,        NULL, CAP_MISMATCH_QUALS_OPTION},
	{"interleaved",          no_argument,        NULL, 'I'},
	{"interleaved-input",    no_argument,        NULL, INTERLEAVED_INPUT_OPTION},
	{"interleaved-output",   no_argument,        NULL, INTERLEAVED_OUTPUT_OPTION},
	{"tab-delimited-input",  no_argument,        NULL, TAB_DELIMITED_INPUT_OPTION},
	{"tab-delimited-output", no_argument,        NULL, TAB_DELIMITED_OUTPUT_OPTION},
	{"output-prefix",        required_argument,  NULL, 'o'},
	{"output-directory",     required_argument,  NULL, 'd'},
	{"to-stdout",            no_argument,        NULL, 'c'},
	{"compress",             no_argument,        NULL, 'z'},
	{"compress-prog",        required_argument,  NULL, COMPRESS_PROG_OPTION},
	{"compress-prog-args",   required_argument,  NULL, COMPRESS_PROG_ARGS_OPTION},
	{"suffix",               required_argument,  NULL, SUFFIX_OPTION},
	{"output-suffix",        required_argument,  NULL, SUFFIX_OPTION},
	{"threads",              required_argument,  NULL, 't'},
	{"quiet",                no_argument,        NULL, 'q'},
	{"help",                 no_argument,        NULL, 'h'},
	{"version",              no_argument,        NULL, 'v'},
	{NULL, 0, NULL, 0}
};

static char *
input_format_str(char *buf, size_t bufsize,
		 const struct read_format_params *iparams,
		 bool interleaved)
{
	switch (iparams->fmt) {
	case READ_FORMAT_FASTQ:
		snprintf(buf, bufsize, "FASTQ, phred_offset=%d%s",
			 iparams->phred_offset,
			 (interleaved ? ", interleaved" : ""));
		break;
	case READ_FORMAT_TAB_DELIMITED:
		snprintf(buf, bufsize, "Tab-delimited, phred_offset=%d",
			 iparams->phred_offset);
		break;
	default:
		assert(0);
		break;
	}
	return buf;
}

static char *
output_format_str(char *buf, size_t bufsize,
		  const struct read_format_params *oparams,
		  bool interleaved_output,
		  enum out_compression_type out_ctype,
		  const char *compress_prog,
		  const char *compress_prog_args)
{
	input_format_str(buf, bufsize, oparams, interleaved_output);
	switch (out_ctype) {
	case OUT_COMPRESSION_NONE:
		if (compress_prog) {
			char *p = strchr(buf, '\0');
			snprintf(p, &buf[bufsize] - p,
				", filtered through '%s %s'",
				compress_prog, compress_prog_args);
		}
		break;
	case OUT_COMPRESSION_GZIP:
		{
			char *p = strchr(buf, '\0');
			snprintf(p, &buf[bufsize] - p, ", gzip");
		}
		break;
	}
	return buf;
}

/* This is just a dynamic array used as a histogram.  It's needed to count the
 * frequencies of the lengths of the combined reads. */
struct histogram {
	unsigned long *array;
	size_t len;
};

static void hist_init(struct histogram *hist)
{
	hist->array = NULL;
	hist->len = 0;
}

static void hist_destroy(struct histogram *hist)
{
	xfree(hist->array, hist->len * sizeof(hist->array[0]));
}

static void hist_add(struct histogram *hist, unsigned long idx,
		     unsigned long amount)
{
	unsigned long *array = hist->array;
	size_t old_len = hist->len;
	if (idx >= old_len) {
		size_t new_len = idx + 1;
		array = xrealloc(array, new_len * sizeof(array[0]));
		memset(&array[old_len], 0,
		       (new_len - old_len) * sizeof(array[0]));
		hist->len = new_len;
		hist->array = array;
	}
	array[idx] += amount;
}

static void hist_inc(struct histogram *hist, unsigned long idx)
{
	hist_add(hist, idx, 1);
}

static void hist_combine(struct histogram *hist, const struct histogram *other)
{
	for (size_t i = 0; i < other->len; i++)
		hist_add(hist, i, other->array[i]);
}

static unsigned long hist_count_at(const struct histogram *hist,
				   unsigned long idx)
{
	assert(idx < hist->len);
	return hist->array[idx];
}

static unsigned long hist_total_zero(const struct histogram *hist)
{
	if (hist->len == 0)
		return 0;
	else
		return hist->array[0];
}

static unsigned long hist_total(const struct histogram *hist)
{
	unsigned long total = 0;
	for (size_t i = 0; i < hist->len; i++)
		total += hist->array[i];
	return total;
}

static void hist_stats(const struct histogram *hist,
		       unsigned long *max_freq_ret,
		       long *first_nonzero_idx_ret,
		       long *last_nonzero_idx_ret)
{
	*max_freq_ret = 0;
	*first_nonzero_idx_ret = -1;
	*last_nonzero_idx_ret = -2;
	for (size_t i = 1; i < hist->len; i++) {
		unsigned long freq = hist->array[i];
		if (freq != 0) {
			if (*first_nonzero_idx_ret == -1)
				*first_nonzero_idx_ret = i;
			*last_nonzero_idx_ret = i;
			if (freq > *max_freq_ret)
				*max_freq_ret = freq;
		}
	}
}


static void write_hist_file(const char *hist_file,
			    const struct histogram *hist,
			    long first_nonzero_idx, long last_nonzero_idx)
{
	FILE *fp = xfopen(hist_file, "w");
	for (long i = first_nonzero_idx; i <= last_nonzero_idx; i++) {
		unsigned long count = hist_count_at(hist, i);
		if (count != 0)
			if (fprintf(fp, "%ld\t%lu\n", i, count) < 0)
				goto write_error;
	}
	xfclose(fp, hist_file);
	return;

write_error:
	fatal_error_with_errno("Error writing to \"%s\"", hist_file);
}

static void write_histogram_file(const char *histogram_file,
				 const struct histogram *hist,
				 long first_nonzero_idx,
				 long last_nonzero_idx,
				 unsigned long max_freq)
{
	const double max_num_asterisks = 72;
	double scale = max_num_asterisks / (double)max_freq;

	FILE *fp = xfopen(histogram_file, "w");

	for (long i = first_nonzero_idx; i <= last_nonzero_idx; i++) {
		if (fprintf(fp, "%ld\t", i) < 0)
			goto write_error;
		size_t num_asterisks = (size_t)(scale * (double)hist_count_at(hist, i));
		while (num_asterisks--)
			if (fputc('*', fp) == EOF)
				goto write_error;
		if (fputc('\n', fp) == EOF)
			goto write_error;
	}
	xfclose(fp, histogram_file);
	return;
write_error:
	fatal_error_with_errno("Error writing to \"%s\"", histogram_file);
}

struct common_combiner_thread_params {
	struct read_io_handle *iohandle;
	struct combine_params alg_params;
};

struct combiner_thread_params {
	struct common_combiner_thread_params *common;
	struct histogram *combined_read_len_hist;
};

/* Buffer for read_sets for which all the read pointers have been invalidated.
 */
struct empty_sets {
	struct read_set *q1[2];
	struct read_set *q2[2];
};

static void
hold_empty_set(struct read_set *q[2], struct read_set *s)
{
	if (!q[0]) {
		q[0] = s;
	} else {
		assert(!q[1]);
		q[1] = s;
	}
}

static void
hold_empty_sets(struct empty_sets *e, struct read_set *s1, struct read_set *s2)
{
	hold_empty_set(e->q1, s1);
	hold_empty_set(e->q2, s2);
}

static struct read_set *
get_empty_set(struct read_set *q[2])
{
	struct read_set *s;

	assert(q[0]);
	s = q[0];
	q[0] = q[1];
	q[1] = NULL;

	return s;
}

static void
get_empty_sets(struct empty_sets *e,
	       struct read_set **s1_ret, struct read_set **s2_ret)
{

	*s1_ret = get_empty_set(e->q1);
	*s2_ret = get_empty_set(e->q2);
}

static void
free_empty_sets(struct empty_sets *e)
{
	free_read_set(e->q1[0]);
	free_read_set(e->q1[1]);

	free_read_set(e->q2[0]);
	free_read_set(e->q2[1]);
}

/* This procedure is executed in parallel by all the combiner threads. */
static void *combiner_thread_proc(void *_params)
{
	struct combiner_thread_params *params = _params;

	struct histogram *combined_read_len_hist = params->combined_read_len_hist;
	struct read_io_handle *iohandle = params->common->iohandle;
	const struct combine_params *alg_params = &params->common->alg_params;

	struct read_set *s_avail_1 = new_empty_read_set(iohandle);
	struct read_set *s_avail_2 = new_empty_read_set(iohandle);
	struct read_set *s_uncombined_1 = new_empty_read_set(iohandle);
	struct read_set *s_uncombined_2 = new_empty_read_set(iohandle);
	struct read_set *s_combined = get_avail_read_set(iohandle);
	struct empty_sets empty = {};

	struct read_set *s1;
	struct read_set *s2;

	/* While there are read pairs to process ...  */
	while (get_unprocessed_read_pairs(iohandle, &s1, &s2)) {

		/* ... process each read pair.  */
		for (size_t i = 0; i < s1->filled; i++) {
			struct read *r1 = s1->reads[i];
			struct read *r2 = s2->reads[i];
			struct read *r_combined;

			s1->reads[i] = NULL;
			s2->reads[i] = NULL;

			reverse_complement(r2);

			/* Get available read in which to try the combination.
			 */
			r_combined = s_combined->reads[s_combined->filled];

			/* Try combining the reads.  */
			if (combine_reads(r1, r2, r_combined, alg_params)) {

				/* Combination was successful.  */

				/* Uncombined read structures are unneeded; mark
				 * them as available.  */

				if (s_avail_1->filled == s_avail_1->num_reads) {
					put_avail_read_pairs(iohandle,
							     s_avail_1,
							     s_avail_2);
					get_empty_sets(&empty,
						       &s_avail_1,
						       &s_avail_2);
				}
				s_avail_1->reads[s_avail_1->filled++] = r1;
				s_avail_2->reads[s_avail_2->filled++] = r2;

				/* Tally length of combined read.  */
				hist_inc(combined_read_len_hist,
					 r_combined->seq_len);

				/* Compute tag for combined read.  */
				get_combined_tag(r1, r2, r_combined);

				/* Send combined read.  */
				if (++s_combined->filled == s_combined->num_reads) {
					put_combined_reads(iohandle, s_combined);
					s_combined = get_avail_read_set(iohandle);
				}
			} else {
				/* Tally 0 length for pair that could not be
				 * combined.  */
				hist_inc(combined_read_len_hist, 0);

				/* Send uncombined reads.  */
				if (s_uncombined_1->filled == s_uncombined_1->num_reads) {
					put_uncombined_read_pairs(iohandle,
								  s_uncombined_1,
								  s_uncombined_2);
					get_empty_sets(&empty,
						       &s_uncombined_1,
						       &s_uncombined_2);
				}

				s_uncombined_1->reads[s_uncombined_1->filled++] = r1;
				s_uncombined_2->reads[s_uncombined_2->filled++] = r2;
				reverse_complement(r2);
			}
		}

		s1->filled = 0;
		s2->filled = 0;
		hold_empty_sets(&empty, s1, s2);
	}

	/* No more reads to combine.  */

	/* Free read sets owned by this thread  */
	free_read_set(s_avail_1);
	free_read_set(s_avail_2);
	free_empty_sets(&empty);

	/* Send out any remaining uncombined and combined reads.
	 * If there are none, free the corresponding read sets.  */

	if (s_uncombined_1->filled) {
		put_uncombined_read_pairs(iohandle, s_uncombined_1, s_uncombined_2);
	} else {
		free_read_set(s_uncombined_1);
		free_read_set(s_uncombined_2);
	}
	if (s_combined->filled)
		put_combined_reads(iohandle, s_combined);
	else
		free_read_set(s_combined);

	notify_combiner_terminated(iohandle);

	xfree(params, sizeof(*params));
	return NULL;
}

int main(int argc, char **argv)
{
	infofile = stdout;

	const char *argv0 = argv[0];

	struct combine_params alg_params = {
		.max_overlap = 0,
		.min_overlap = 10,
		.max_mismatch_density = 0.25,
		.cap_mismatch_quals = false,
	};
	struct read_format_params iparams = {
		.fmt = READ_FORMAT_FASTQ,
		.phred_offset = 33,
	};
	struct read_format_params oparams = {
		.fmt = READ_FORMAT_FASTQ,
		.phred_offset = 33,
	};
	int read_len               = 100;
	int fragment_len           = 180;
	int fragment_len_stddev    = 18;
	const char *prefix         = "out";
	const char *output_dir     = ".";
	bool to_stdout             = false;
	bool verbose               = true;
	bool interleaved_input     = false;
	bool interleaved_output    = false;
	struct input_stream *mates1_in = NULL;
	struct input_stream *mates2_in = NULL;
	enum out_compression_type out_ctype = OUT_COMPRESSION_NONE;
	const char *compress_prog = NULL;
	char *compress_prog_args = "-c -";
	bool compress_prog_args_allocated = false;
	struct output_stream *out_combined = NULL;
	struct output_stream *out_notcombined_1 = NULL;
	struct output_stream *out_notcombined_2 = NULL;
	const char *out_filetype   = "fastq";
	char *out_suffix           = "";
	bool out_suffix_allocated  = false;
	unsigned long num_combiner_threads = 0;
	int c;
	char *tmp;
	struct timeval start_time;
	gettimeofday(&start_time, NULL);

	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
		switch (c) {
		case 'm':
			alg_params.min_overlap = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || alg_params.min_overlap < 1)
				fatal_error("Minimum overlap must be a "
					    "positive integer!  Please check "
					    "option -m.");
			break;
		case 'M':
			alg_params.max_overlap = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || alg_params.max_overlap < 1)
				fatal_error("Maximum overlap must be "
					    "a positive integer!  Please check "
					    "option -M.");
			break;
		case 'x':
			alg_params.max_mismatch_density = strtod(optarg, &tmp);
			if (tmp == optarg || *tmp || alg_params.max_mismatch_density < 0.0 ||
			    alg_params.max_mismatch_density > 1.0)
			{
				fatal_error("Max mismatch density must be a "
					    "number in the interval [0, 1]! "
					    "Please check option -x.");
			}
			break;
		case 'p':
			oparams.phred_offset =
				iparams.phred_offset =
					strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp ||
			    iparams.phred_offset < 0 ||
			    iparams.phred_offset > 127)
			{
				fatal_error("Phred offset must be an integer "
					    "in the range [0, 127]!  Please "
					    "check option -p.");
			}
			if (iparams.phred_offset != 33 &&
			    iparams.phred_offset != 64)
			{
				warning("Phred offset is usually either "
				        "64 (for earlier Illumina data) or 33 "
				        "(for Sanger and later Illumina data).");
			}
			break;
		case 'f':
			fragment_len = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || fragment_len <= 0)
				fatal_error("Fragment length must be a "
					    "positive integer!  Please check "
					    "option -f.");
			break;
		case 's':
			fragment_len_stddev = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || fragment_len_stddev <= 0)
				fatal_error("Fragment length standard deviation "
					    "must be a positive integer!  "
					    "Please check option -s.");
			break;
		case 'r':
			read_len = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || read_len <= 0)
				fatal_error("Read length must be a "
					    "positive integer!  Please check "
					    "option -r.");
			break;
		case CAP_MISMATCH_QUALS_OPTION:
			alg_params.cap_mismatch_quals = true;
			break;
		case 'I':
			interleaved_input = true;
			interleaved_output = true;
			break;
		case INTERLEAVED_INPUT_OPTION:
			interleaved_input = true;
			break;
		case INTERLEAVED_OUTPUT_OPTION:
			interleaved_output = true;
			break;
		case 'T':
			if ((*optarg != 'i' && *optarg != 'o') ||
			    *(optarg + 1))
			{
				fatal_error("Invalid option -T%s!  Use -Ti "
					    "(short for --tab-delimited-input)\n\t"
					    "or -To (short for "
					    "--tab-delimited-output)", optarg);
			}
			if (*optarg == 'i') {
		case TAB_DELIMITED_INPUT_OPTION:
				iparams.fmt = READ_FORMAT_TAB_DELIMITED;
			} else {
		case TAB_DELIMITED_OUTPUT_OPTION:
				oparams.fmt = READ_FORMAT_TAB_DELIMITED;
				out_filetype = "tab";
			}
			break;
		case 'o':
			prefix = optarg;
			break;
		case 'd':
			output_dir = optarg;
			break;
		case 'c':
			to_stdout = true;
			infofile = stderr;
			break;
		case 'z':
			out_ctype = OUT_COMPRESSION_GZIP;
			if (out_suffix_allocated)
				free(out_suffix);
			out_suffix = ".gz";
			out_suffix_allocated = false;
			compress_prog = NULL;
			break;
		case COMPRESS_PROG_OPTION:
			if (out_suffix_allocated)
				free(out_suffix);
			out_suffix = xmalloc(strlen(optarg) + 2);
			sprintf(out_suffix, ".%s", optarg);
			out_suffix_allocated = true;
			compress_prog = optarg;
			out_ctype = OUT_COMPRESSION_NONE;
			break;
		case COMPRESS_PROG_ARGS_OPTION:
			if (compress_prog_args_allocated)
				free(compress_prog_args);
			compress_prog_args = xmalloc(strlen(optarg) + 6);
			sprintf(compress_prog_args, "%s -c -", optarg);
			compress_prog_args_allocated = true;
			break;
		case SUFFIX_OPTION:
			if (out_suffix_allocated)
				free(out_suffix);
			if (*optarg) {
				out_suffix = xmalloc(strlen(optarg) + 2);
				sprintf(out_suffix, ".%s", optarg);
				out_suffix_allocated = true;
			} else {
				out_suffix = optarg;
				out_suffix_allocated = false;
			}
			break;
		case 't':
			num_combiner_threads = strtoul(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || num_combiner_threads < 1 ||
			    num_combiner_threads > UINT_MAX) {
				fatal_error("Number of threads must be "
					    "a positive integer!  Please "
					    "check option -t.");
			}
			break;
		case 'q':
			verbose = false;
			break;
		case 'v':
			version();
			return 0;
		case 'h':
			usage(argv0);
			return 0;
		default:
			usage_short(argv0);
			return 2;
		}
	}

	if (alg_params.max_overlap == 0)
		alg_params.max_overlap = (int)(2 * read_len - fragment_len +
					       2.5 * fragment_len_stddev);

	if (alg_params.max_overlap < alg_params.min_overlap) {
		fatal_error(
"Maximum overlap (%d) cannot be less than the minimum overlap (%d).\n"
"Please make sure you have provided the read length and fragment length\n"
"correctly.  Or, alternatively, specify the minimum and maximum overlap\n"
"manually with the --min-overlap and --max-overlap options.",
			alg_params.max_overlap, alg_params.min_overlap);
	}

	if (num_combiner_threads == 0)
		num_combiner_threads = get_default_num_threads();

	argc -= optind;
	argv += optind;

	if (argc == 0 || argc > 2) {
		usage_short(argv0);
		return 2;
	}
	if (interleaved_input && argc != 1)
		fatal_error("With --interleaved-input, only 1 input "
			    "file is allowed!");

	if (interleaved_input && iparams.fmt != READ_FORMAT_FASTQ)
		fatal_error("--interleaved-input is only relevant for FASTQ input!");

	if (argc == 1 && !interleaved_input && iparams.fmt == READ_FORMAT_FASTQ)
		fatal_error("Only 1 input file was specified!  Specify "
			    "--interleaved-input\n"
			    "\tif you're providing an interleaved FASTQ file, "
			    "or --tab-delimited-input\n"
			    "\tif you're providing a tab-delimited input file.  "
			    "Or specify two input\n"
			    "\tfiles (for read 1 and read 2 of each pair).");

	mates1_in = new_input_stream(argv[0]);
	if (argc > 1)
		mates2_in = new_input_stream(argv[1]);

	mkdir_p(output_dir);

	/* Open the output files.  */

	char name_buf[strlen(output_dir) + 1 + strlen(prefix) +
		      100 + strlen(out_suffix) + 1];
	char *suffix;
	suffix = name_buf + sprintf(name_buf, "%s/%s", output_dir, prefix);

	if (oparams.fmt == READ_FORMAT_TAB_DELIMITED) {
		sprintf(suffix, ".readsAndPairs.%s%s", out_filetype, out_suffix);
		out_combined = new_output_stream(out_ctype,
						 (to_stdout ? "-" : name_buf),
						 compress_prog,
						 compress_prog_args);
	} else {
		sprintf(suffix, ".extendedFrags.%s%s", out_filetype, out_suffix);
		out_combined = new_output_stream(out_ctype,
						 (to_stdout ? "-" : name_buf),
						 compress_prog,
						 compress_prog_args);

		if (!to_stdout) {
			if (interleaved_output) {
				sprintf(suffix, ".notCombined.%s%s",
					out_filetype, out_suffix);
				out_notcombined_1 = new_output_stream(out_ctype,
								      name_buf,
								      compress_prog,
								      compress_prog_args);
			} else {
				sprintf(suffix, ".notCombined_1.%s%s",
					out_filetype, out_suffix);
				out_notcombined_1 = new_output_stream(out_ctype,
								      name_buf,
								      compress_prog,
								      compress_prog_args);

				sprintf(suffix, ".notCombined_2.%s%s",
					out_filetype, out_suffix);
				out_notcombined_2 = new_output_stream(out_ctype,
								      name_buf,
								      compress_prog,
								      compress_prog_args);
			}
		}
	}

	*suffix = '\0';

	if (verbose) {
		info("Starting FLASH " VERSION_STR);
		info("Fast Length Adjustment of SHort reads");
		info(" ");
		info("Input files:");
		info("    %s", input_stream_get_name(mates1_in));
		if (mates2_in)
			info("    %s", input_stream_get_name(mates2_in));
		info(" ");
		info("Output files:");
		info("    %s", output_stream_get_name(out_combined));
		if (out_notcombined_1)
			info("    %s", output_stream_get_name(out_notcombined_1));
		if (out_notcombined_2)
			info("    %s", output_stream_get_name(out_notcombined_2));
		if (!to_stdout) {
			info("    %s.hist", name_buf);
			info("    %s.histogram", name_buf);
		}
		info(" ");
		info("Parameters:");
		info("    Min overlap:           %d",
		     alg_params.min_overlap);
		info("    Max overlap:           %d",
		     alg_params.max_overlap);
		info("    Max mismatch density:  %f",
		     alg_params.max_mismatch_density);
		info("    Cap mismatch quals:    %s",
		     alg_params.cap_mismatch_quals ? "true" : "false");
		info("    Combiner threads:      %u",
		     (unsigned)num_combiner_threads);

		char buf[256];
		info("    Input format:          %s",
		     input_format_str(buf, ARRAY_LEN(buf),
				      &iparams, interleaved_input));
		info("    Output format:         %s",
		     output_format_str(buf, ARRAY_LEN(buf),
				       &oparams, interleaved_output,
				       out_ctype,
				       compress_prog, compress_prog_args));
		info(" ");
	}

	/*
	 * We wish to do the following:
	 *
	 * "Go through each mate pair in the input files.  Determine if it can
	 * be combined, given the input parameters to the program.  If it can,
	 * write the combined read to the PREFIX.extendedFrags.fastq file.
	 * Otherwise, write the reads in the mate pair to the
	 * PREFIX.notCombined_1.fastq and PREFIX.notCombined_2.fastq files, or
	 * PREFIX.notCombined.fastq for interleaved output.  Or, if the -c /
	 * --to-stdout option is specified, write the combined reads to standard
	 * output, and ignore the uncombined reads."
	 *
	 * In the following implementation, there will be @num_combiner_threads
	 * combiner threads created that will process the reads in parallel by
	 * retrieving `struct read_set'-sized chunks of reads from the reader
	 * thread(s), and providing `struct read_set'-sized chunks of combined
	 * or uncombined reads to the writer threads.
	 */

	/* Histogram of how many combined reads have a given length.
	 *
	 * The zero index slot counts how many reads were not combined.
	 *
	 * There is a copy of the histogram for each thread, and they are
	 * combined after all the combiner threads are done. */
	struct histogram combined_read_len_hists[num_combiner_threads];
	struct histogram *combined_read_len_hist =
			&combined_read_len_hists[num_combiner_threads - 1];
	for (size_t i = 0; i < ARRAY_LEN(combined_read_len_hists); i++)
		hist_init(&combined_read_len_hists[i]);

	struct read_io_handle *iohandle =
		start_readers_and_writers(mates1_in,
					  mates2_in,
					  out_combined,
					  out_notcombined_1,
					  out_notcombined_2,
					  &iparams, &oparams,
					  num_combiner_threads, verbose);
	struct common_combiner_thread_params common = {
		.iohandle = iohandle,
		.alg_params = alg_params,
	};

	if (verbose)
		info("Starting %u combiner threads", (unsigned)num_combiner_threads);

	pthread_t other_combiner_threads[num_combiner_threads - 1];
	for (unsigned i = 0; i < num_combiner_threads; i++) {
		struct combiner_thread_params *p;

		p = xmalloc(sizeof(*p));
		p->common = &common;
		p->combined_read_len_hist = &combined_read_len_hists[i];
		if (i < num_combiner_threads - 1)
			other_combiner_threads[i] =
					create_thread(combiner_thread_proc, p);
		else
			combiner_thread_proc(p);
	}
	for (unsigned i = 0; i < num_combiner_threads - 1; i++) {
		join_thread(other_combiner_threads[i]);
		hist_combine(combined_read_len_hist,
			     &combined_read_len_hists[i]);
		hist_destroy(&combined_read_len_hists[i]);
	}
	stop_readers_and_writers(iohandle);

	if (verbose) {
		unsigned long num_combined_reads;
		unsigned long num_uncombined_reads;
		unsigned long num_total_reads;
		num_uncombined_reads = hist_total_zero(combined_read_len_hist);
		num_total_reads = hist_total(combined_read_len_hist);
		num_combined_reads = num_total_reads - num_uncombined_reads;
		info(" ");
		info("Read combination statistics:");
		info("    Total reads:      %lu", num_total_reads);
		info("    Combined reads:   %lu", num_combined_reads);
		info("    Uncombined reads: %lu", num_uncombined_reads);
		info("    Percent combined: %.2f%%", (num_total_reads) ?
			(double)num_combined_reads * 100 / num_total_reads : 0);
		info(" ");
	}

	if (!to_stdout) {
		unsigned long max_freq;
		long first_nonzero_idx;
		long last_nonzero_idx;

		if (verbose)
			info("Writing histogram files.");

		hist_stats(combined_read_len_hist,
			   &max_freq, &first_nonzero_idx, &last_nonzero_idx);

		/* Write the raw numbers of the combined read length histogram
		 * to the PREFIX.hist file. */
		strcpy(suffix, ".hist");
		write_hist_file(name_buf, combined_read_len_hist,
				first_nonzero_idx, last_nonzero_idx);
		/* Write a pretty representation of the combined read length
		 * histogram to the PREFIX.histogram file. */
		strcpy(suffix, ".histogram");
		write_histogram_file(name_buf, combined_read_len_hist,
				     first_nonzero_idx, last_nonzero_idx,
				     max_freq);
	}
	hist_destroy(combined_read_len_hist);

	if (verbose) {
		struct timeval end_time;
		gettimeofday(&end_time, NULL);
		uint64_t start_usec = (uint64_t)start_time.tv_sec * 1000000 + start_time.tv_usec;
		uint64_t end_usec = (uint64_t)end_time.tv_sec * 1000000 + end_time.tv_usec;
		info(" ");
		info("FLASH " VERSION_STR " complete!");
		info("%.3f seconds elapsed", (double)(end_usec - start_usec) / 1000000);
	}
	if (out_suffix_allocated)
		free(out_suffix);
	if (compress_prog_args_allocated)
		free(compress_prog_args);
	return 0;
}
