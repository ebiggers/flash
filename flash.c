#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <zlib.h>

#include "combine_reads.h"
#include "fastq.h"
#include "util.h"

#define VERSION_STR "v1.2.4"

static void usage()
{
	const char *usage_str =
"Usage: flash [OPTIONS] MATES_1.FASTQ MATES_2.FASTQ\n"
"\n"
"DESCRIPTION:\n"
"\n"
"FLASH (Fast Length Adjustment of SHort reads) is an accurate and fast tool\n"
"to merge paired-end reads that were generated from DNA fragments whose\n"
"lengths are shorter than twice the length of reads.  Merged read pairs result\n"
"in unpaired longer reads.  Longer reads are generally more desired in genome\n"
"assembly and genome analysis processes.\n"
"\n"
"MANDATORY INPUT:\n"
"\n"
"To run FLASH, you may provide two FASTQ files of paired-end reads\n"
"where corresponding paired reads are in the same order in both files.\n"
"Alternatively, you may provide one FASTQ file, which may be standard input,\n"
"containing interleaved paired-end reads (see the --interleaved option).\n"
"The input FASTQ files may be either plain-text or compressed with gzip.\n"
"Other compression formats for the input files are not yet supported.\n"
"\n"
"OUTPUT:\n"
"\n"
"The default output of FLASH is a FASTQ file containing the extended fragments\n"
"produced by combining read pairs, two FASTQ files containing read pairs\n"
"that were not combined, and histogram files that show the distribution of\n"
"lengths of the extended fragments.  Writing the uncombined read pairs to an\n"
"interleaved FASTQ file is also supported.  Also, writing the extended\n"
"fragments directly to standard output is supported.  Plain-text and gzip\n"
"output formats are natively supported; other output formats are supported\n"
"indirectly via the --compress-prog option.  (Note that this is all FASTQ.)\n"
"\n"
"OPTIONS:\n"
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
"                          70 bp.  If --max-overlap is specified, then the\n"
"                          specified value overrides the calculated value.\n"
"\n"
"                          If you do not know the standard deviation of the\n"
"                          fragment library, you can probably assume that the\n"
"                          standard deviation is 10% of the average fragment\n"
"                          length.\n"
"\n"
"  --interleaved-input     Instead of requiring files MATES_1.FASTQ and\n"
"                          MATES_2.FASTQ, allow a single file MATES.FASTQ that\n"
"                          has the paired-end reads interleaved.  Specify \"-\"\n"
"                          to read from standard input.\n"
"\n"
"  --interleaved-output    Write the uncombined pairs in interleaved format.\n"
"\n"
"  -I, --interleaved       Equivalent to specifying both --interleaved-input\n"
"                          and --interleaved-output.\n"
"\n"
"  -o, --output-prefix=PREFIX\n"
"                          Prefix of output files.  Default: \"out\".\n"
"\n"
"  -d, --output-directory=DIR\n"
"                          Path to directory for output files.  Default:\n"
"                          current working directory.\n"
"\n"
"  -c, --to-stdout         Write the combined reads to standard output; do not\n"
"                          write uncombined reads to anywhere.\n"
"\n"
"  -z, --compress          Compress the FASTQ output files directly with zlib.\n"
"                          Similar to specifying --compress-prog=gzip and\n"
"                          --suffix=gz, but may be slightly faster.\n"
"\n"
"  --compress-prog=PROG    Pipe the output through the compression program\n"
"                          PROG, which will be called as `PROG -c -',\n"
"                          plus any arguments specified by --compress-prog-args.\n"
"                          PROG must read uncompressed data from standard input\n"
"                          and write compressed data to standard output.\n"
"                          Examples: gzip, bzip2, xz, pigz.\n"
"\n"
"  --compress-prog-args=ARGS\n"
"                          A string of arguments that will be passed to the\n"
"                          compression program if one is specified with\n"
"                          --compress-prog.  Note: the argument -c is already\n"
"                          assumed.\n"
"\n"
"  --suffix=SUFFIX, --output-suffix=SUFFIX\n"
"                          Use SUFFIX as the suffix of the output files\n"
"                          after \".fastq\".  A dot before the suffix is assumed,\n"
"                          unless an empty suffix is provided.  Default:\n"
"                          nothing; or 'gz' if -z is specified; or PROG if\n"
"                          --compress-prog is specified.\n"
"\n"
"  -t, --threads=NTHREADS  Set the number of worker threads.  This is in\n"
"                          addition to the I/O threads.  Default: number of\n"
"                          processors.  Note: if you need FLASH's output to\n"
"                          appear deterministically or in the same order as\n"
"                          the original reads, you must specify -t 1\n"
"                          (--threads=1).\n"
"\n"
"  -q, --quiet             Do not print informational messages.  (Implied with\n"
"                          --to-stdout.)\n"
"\n"
"  -h, --help              Display this help and exit.\n"
"\n"
"  -v, --version           Display version.\n"
	;
	fputs(usage_str, stdout);
}

static void usage_short()
{
	fputs("Usage: flash [OPTIONS] MATES_1.FASTQ MATES_2.FASTQ\n"
	      "Try `flash -h' for more information.\n",
	      stderr);
}

static void version()
{
	puts("FLASH " VERSION_STR);
	puts("License:  GNU General Public License Version 3+ (http://gnu.org/licenses/gpl.html)");
	puts("Report bugs to flash.comment@gmail.com or https://sourceforge.net/p/flashpage/bugs");
}

enum {
	INTERLEAVED_INPUT_OPTION = 257,
	INTERLEAVED_OUTPUT_OPTION,
	COMPRESS_PROG_OPTION,
	COMPRESS_PROG_ARGS_OPTION,
	SUFFIX_OPTION,
};

static const char *optstring = "m:M:x:p:r:f:s:Io:d:czt:qhv";
static const struct option longopts[] = {
	{"min-overlap",          required_argument,  NULL, 'm'},
	{"max-overlap",          required_argument,  NULL, 'M'},
	{"max-mismatch-density", required_argument,  NULL, 'x'},
	{"phred-offset",         required_argument,  NULL, 'p'},
	{"read-len",             required_argument,  NULL, 'r'},
	{"fragment-len",         required_argument,  NULL, 'f'},
	{"fragment-len-stddev",  required_argument,  NULL, 's'},
	{"interleaved",          no_argument,        NULL, 'I'},
	{"interleaved-input",    no_argument,        NULL,  INTERLEAVED_INPUT_OPTION},
	{"interleaved-output",   no_argument,        NULL,  INTERLEAVED_OUTPUT_OPTION},
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


static void copy_tag(struct read *to, const struct read *from)
{
	if (to->tag_bufsz < from->tag_len + 1) {
		to->tag = xrealloc(to->tag, from->tag_len + 1);
		to->tag_bufsz = from->tag_len + 1;
	}
	to->tag_len = from->tag_len;
	memcpy(to->tag, from->tag, from->tag_len + 1);
}

/*
 * Given the FASTQ tags of two paired-end reads, find the FASTQ tag to give to
 * the combined read.
 *
 * This is done by stripping off the characters trailing the '/' (e.g. "/1" and
 * "/2"), unless there is a "barcode" beginning with the '#' character, which is
 * kept.
 */
static void get_combined_tag(const struct read *read_1,
			     const struct read *read_2,
			     struct read *combined_read)
{
	char *p;
	copy_tag(combined_read, read_1);
	for (p = &combined_read->tag[combined_read->tag_len - 1];
	     p >= combined_read->tag;
	     p--)
	{
		if (*p == '/') {
			/* Tags are different, and there's a forward slash in
			 * the first tag.  Remove everything after the forward
			 * slash, unless there's a barcode, which we keep. */
			if (*(p + 1) != '\0' && *(p + 2) == '#') {
				/* read ID has a barcode. */
				do {
					*p = *(p + 2);
				} while (*(++p + 2) != '\0');
			}
			*p = '\0';
			combined_read->tag_len = p - combined_read->tag;
			break;
		}
	}
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
	free(hist->array);
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
		if (count != 0) {
			if (fprintf(fp, "%ld\t%zu\n", i, count) < 0) {
				fatal_error_with_errno("Error writing to "
						       "the file \"%s\"",
						       hist_file);
			}
		}
	}
	xfclose(fp);
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
	xfclose(fp);
	return;
write_error:
	fatal_error_with_errno("Error writing to the file \"%s\"",
			       histogram_file);
}

struct common_combiner_thread_params {
	struct threads *threads;
	int min_overlap;
	int max_overlap;
	float max_mismatch_density;
	bool to_stdout;
	pthread_mutex_t unqueue_lock;
	pthread_mutex_t queue_lock;
};

struct combiner_thread_params {
	struct common_combiner_thread_params *common;
	struct histogram *combined_read_len_hist;
};

/* This procedure is executed in parallel by all the combiner threads. */
static void *combiner_thread_proc(void *_params)
{
	struct combiner_thread_params *params = _params;

	struct histogram *combined_read_len_hist = params->combined_read_len_hist;
	struct threads *threads = params->common->threads;
	int min_overlap = params->common->min_overlap;
	int max_overlap = params->common->max_overlap;
	float max_mismatch_density = params->common->max_mismatch_density;
	bool to_stdout = params->common->to_stdout;

	struct read_set *read_set_1;
	struct read_set *read_set_2;
	struct read_set *to_reader_1 = new_empty_read_set();
	struct read_set *to_reader_2 = new_empty_read_set();
	struct read_set *to_writer_1 = new_empty_read_set();
	struct read_set *to_writer_2 = new_empty_read_set();
	struct read_set *empty_read_set_1 = NULL;
	struct read_set *empty_read_set_2 = NULL;
	unsigned to_reader_filled = 0;
	unsigned to_writer_filled = 0;
	struct read_set *combined_read_set = read_queue_get(threads->writer_combined.free_q);
	unsigned combined_read_filled = 0;
	struct read *combined_read;
	struct read *read_1;
	struct read *read_2;
	unsigned i;
	unsigned j;
	bool combination_successful;

	while (1) { /* Process each read set, retrieved from the readers' ready queues */

		/* To make sure the read pairs are made correctly (i.e. always
		 * having read 1 correctly paired with the corresponding read 2,
		 * and not some other read 2), we need to make the retrieval of
		 * the two read sets a critical section.  Otherwise, race
		 * conditions could occur where thread A retrieves a set of
		 * read 1's, then thread B retrieves the next set of read 1's
		 * along with the read 2's that should have been received by
		 * thread A.
		 *
		 * This should not be big bottleneck because all the reads in
		 * each read set need to be processed before we get the next
		 * read set. */
		pthread_mutex_lock(&params->common->unqueue_lock);
		read_set_1 = read_queue_get(threads->reader_1.ready_q);
		read_set_2 = read_queue_get(threads->reader_2.ready_q);
		pthread_mutex_unlock(&params->common->unqueue_lock);

		/* Process each read in the read set */
		for (i = 0; i < READS_PER_READ_SET; i++) {
			if (i == READS_PER_READ_SET - 1) {
				empty_read_set_1 = read_set_1;
				empty_read_set_2 = read_set_2;
			}

			read_1 = read_set_1->reads[i];
			read_2 = read_set_2->reads[i];

			if (!read_1 || !read_2) {
				if (read_1 || read_2) {
					fatal_error("Input reads files do not "
						    "contain the same number of "
						    "reads!");
				}
				/* End-of file reached; break out of the outer
				 * loop. */
				goto no_more_reads;
			}
			reverse_complement(read_2->seq, read_2->seq_len);
			reverse(read_2->qual, read_2->seq_len);


			/* Next available combined read in the combined_read_set
			 * */
			combined_read = combined_read_set->reads[combined_read_filled];

			/* Try combining the reads. */
			combination_successful = combine_reads(read_1, read_2,
							       combined_read,
							       min_overlap,
							       max_overlap,
							       max_mismatch_density);
			if (combination_successful) {
				/* Combination was successful. */

				hist_inc(combined_read_len_hist,
					 combined_read->seq_len);

				get_combined_tag(read_1, read_2, combined_read);

				if (++combined_read_filled == READS_PER_READ_SET) {
					read_queue_put(threads->writer_combined.ready_q,
						       combined_read_set);
					combined_read_set = read_queue_get(
							threads->writer_combined.free_q);
					combined_read_filled = 0;
				}
			} else {
				hist_inc(combined_read_len_hist, 0);
			}

			if (combination_successful || to_stdout) {
				/* We do not need to write the uncombined reads,
				 * either because the reads were combined, or we
				 * are writing to stdout.  So put them back in
				 * the queues for the reader threads. */

				to_reader_1->reads[to_reader_filled] = read_1;
				to_reader_2->reads[to_reader_filled] = read_2;
				if (++to_reader_filled == READS_PER_READ_SET) {
					read_queue_put(threads->reader_1.free_q, to_reader_1);
					read_queue_put(threads->reader_2.free_q, to_reader_2);
					to_reader_filled = 0;
					to_reader_1 = empty_read_set_1;
					to_reader_2 = empty_read_set_2;
				}
			} else {
				/* We need to write the uncombined reads,
				 * because they were not combined, and we are
				 * not writing to stdout.  So enqueue them in
				 * the queues for the writer threads. */
				to_writer_1->reads[to_writer_filled] = read_1;
				to_writer_2->reads[to_writer_filled] = read_2;
				reverse_complement(read_2->seq, read_2->seq_len);
				reverse(read_2->qual, read_2->seq_len);
				if (++to_writer_filled == READS_PER_READ_SET) {
					pthread_mutex_lock(&params->common->queue_lock);
					read_queue_put(threads->writer_uncombined_1.ready_q,
						       to_writer_1);
					read_queue_put(threads->writer_uncombined_2.ready_q,
						       to_writer_2);
					pthread_mutex_unlock(&params->common->queue_lock);
					to_writer_filled = 0;
					to_writer_1 = empty_read_set_1;
					to_writer_2 = empty_read_set_2;
				}
			}
		}
	}
no_more_reads:
	/* This combiner thread has been signalled by the reader(s) that there
	 * are no more reads. */

	/* Only the first @to_reader_filled reads in the @to_reader_? read sets
	 * are valid.  Set the rest to NULL so they aren't freed. */
	assert(0 <= to_reader_filled && to_reader_filled < READS_PER_READ_SET);
	j = to_reader_filled;
	do {
		to_reader_1->reads[j] = NULL;
		to_reader_2->reads[j] = NULL;
	} while (++j != READS_PER_READ_SET);

	/* Only the first @to_writer_filled reads in the @to_writer_? read sets
	 * are valid.  Set the rest to NULL so they aren't freed, and to
	 * NULL-terminate the read sets to signal the uncombined read writer
	 * threads that there are no more reads from this combiner thread. */
	assert(0 <= to_writer_filled && to_writer_filled < READS_PER_READ_SET);
	j = to_writer_filled;
	do {
		to_writer_1->reads[j] = NULL;
		to_writer_2->reads[j] = NULL;
	} while (++j != READS_PER_READ_SET);


	/* Only the first @combined_read_filled reads in the @combined_read_set
	 * are actually ready to be written.  The rest need to be freed, but not
	 * written; and the read set needs to be NULL-terminated to signal the
	 * combined read writer that there are no more reads from this combiner
	 * thread. */
	assert(0 <= combined_read_filled && combined_read_filled < READS_PER_READ_SET);
	j = combined_read_filled;
	do {
		free_read(combined_read_set->reads[j]);
		combined_read_set->reads[j] = NULL;
	} while (++j != READS_PER_READ_SET);

	/* The first @i reads in @read_set_? are invalid because they've been
	 * moved to the @to_{reader,writer}_? sets.  Set them to NULL so they
	 * aren't freed again. */
	assert(0 <= i && i < READS_PER_READ_SET);
	for (j = 0; j < i; j++) {
		read_set_1->reads[j] = NULL;
		read_set_2->reads[j] = NULL;
	}

	/* Now, we need to dispose of the read sets properly.
	 *
	 * - @to_reader_?, @read_set_?, and @empty_read_set_? are all owned by
	 *   this thread and can be freed.
	 *
	 * - @to_writer_? and @combined_read_set must be sent to the writers
	 *   because each writer expects to receive a NULL-terminated read set
	 *   from each combiner thread; furthermore, there may be unwritten
	 *   reads in these sets.
	 */
	free_read_set(to_reader_1);
	free_read_set(to_reader_2);
	free_read_set(read_set_1);
	free_read_set(read_set_2);

	/* @empty_read_set_? are only valid if at least %READS_PER_READ_SET
	 * reads total are waiting to be sent off to the uncombined read writers
	 * and readers, total.   But although the empty read sets may be valid,
	 * the reads in them will not be, so use free() instead of
	 * free_read_set(). */
	if (to_writer_filled + to_reader_filled >= READS_PER_READ_SET) {
		free(empty_read_set_1);
		free(empty_read_set_2);
	}

	pthread_mutex_lock(&params->common->queue_lock);
	read_queue_put(threads->writer_uncombined_1.ready_q, to_writer_1);
	read_queue_put(threads->writer_uncombined_2.ready_q, to_writer_2);
	pthread_mutex_unlock(&params->common->queue_lock);

	read_queue_put(threads->writer_combined.ready_q, combined_read_set);
	free(params);
	return NULL;
}

int main(int argc, char **argv)
{
	int max_overlap            = 0;
	int min_overlap            = 10;
	float max_mismatch_density = 0.25;
	int phred_offset           = 33;
	int read_len               = 100;
	int fragment_len           = 180;
	int fragment_len_stddev    = 18;
	const char *prefix         = "out";
	const char *output_dir     = ".";
	bool to_stdout             = false;
	bool verbose               = true;
	bool interleaved_input     = false;
	bool interleaved_output    = false;
	gzFile mates1_gzf          = NULL;
	gzFile mates2_gzf          = NULL;
	void *out_combined_fp      = NULL;
	void *out_notcombined_fp_1 = NULL;
	void *out_notcombined_fp_2 = NULL;
	bool out_suffix_allocated  = false;
	char *out_suffix           = NULL;
	int num_combiner_threads   = 0;
	struct file_operations *fops = &normal_fops;
	int c;
	char *tmp;
	struct timeval start_time;
	gettimeofday(&start_time, NULL);

	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
		switch (c) {
		case 'm':
			min_overlap = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || min_overlap < 1)
				fatal_error("Minimum overlap must be a "
					    "positive integer!  Please check "
					    "option -m.");
			break;
		case 'M':
			max_overlap = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || max_overlap < 1)
				fatal_error("Maximum overlap must be "
					    "a positive integer!  Please check "
					    "option -M.");
			break;
		case 'x':
			max_mismatch_density = strtod(optarg, &tmp);
			if (tmp == optarg || *tmp || max_mismatch_density < 0.0 ||
			    max_mismatch_density > 1.0)
			{
				fatal_error("Max mismatch density must be a "
					    "number in the interval [0, 1]! "
					    "Please check option -x.");
			}
			break;
		case 'p':
			phred_offset = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp ||
			    phred_offset < 0 || phred_offset > 127)
			{
				fatal_error("Phred offset must be a in integer "
					    "in the range [0, 127]!  Please "
					    "check option -p.");
			}
			if (phred_offset != 33 && phred_offset != 64) {
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
		case 'o':
			prefix = optarg;
			break;
		case 'd':
			output_dir = optarg;
			break;
		case 'c':
			to_stdout = true;
			verbose = false;
			break;
		case 'z':
			fops = &gzip_fops;
			break;
		case COMPRESS_PROG_OPTION:
			pipe_fops.name = optarg;
			if (out_suffix == NULL) {
				out_suffix = xmalloc(strlen(optarg) + 2);
				sprintf(out_suffix, ".%s", optarg);
				out_suffix_allocated = true;
			}
			fops = &pipe_fops;
			compress_prog = optarg;
			break;
		case COMPRESS_PROG_ARGS_OPTION:
			compress_prog_args = optarg;
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
			num_combiner_threads = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || num_combiner_threads < 1) {
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
			usage();
			return 0;
		default:
			usage_short();
			return 2;
		}
	}

	if (out_suffix != NULL && out_suffix != fops->suffix)
		fops->suffix = out_suffix;

	if (max_overlap == 0)
		max_overlap = (int)(2 * read_len - fragment_len +
				    2.5 * fragment_len_stddev);

	if (max_overlap < min_overlap) {
		fatal_error(
"Maximum overlap (%d) cannot be less than the minimum overlap (%d).\n"
"Please make sure you have provided the read length and fragment length\n"
"correctly.  Or, alternatively, specify the minimum and maximum overlap\n"
"manually with the --min-overlap and --max-overlap options.",
			max_overlap, min_overlap);
	}

	if (num_combiner_threads == 0)
		num_combiner_threads = get_default_num_threads();

	argc -= optind;
	argv += optind;

	if ((interleaved_input && argc != 1) || (!interleaved_input && argc != 2)) {
		usage_short();
		return 2;
	}

	mates1_gzf = xgzopen(argv[0], "r");
	if (!interleaved_input)
		mates2_gzf = xgzopen(argv[1], "r");

	mkdir_p(output_dir);

	char name_buf[strlen(output_dir) + 1 + strlen(prefix) +
		      strlen(".notCombined_2.fastq") +
		      strlen(fops->suffix) + 1];
	char *suffix;
	suffix = name_buf + sprintf(name_buf, "%s/%s", output_dir, prefix);


	/* Open the output files. */
	if (to_stdout) {
		out_combined_fp = fops->open_file("-", "w");
	} else {
		sprintf(suffix, ".extendedFrags.fastq%s", fops->suffix);
		out_combined_fp = fops->open_file(name_buf, "w");

		if (interleaved_output) {
			sprintf(suffix, ".notCombined.fastq%s", fops->suffix);
			out_notcombined_fp_1 = fops->open_file(name_buf, "w");
		} else {
			sprintf(suffix, ".notCombined_1.fastq%s", fops->suffix);
			out_notcombined_fp_1 = fops->open_file(name_buf, "w");

			sprintf(suffix, ".notCombined_2.fastq%s", fops->suffix);
			out_notcombined_fp_2 = fops->open_file(name_buf, "w");
		}
		*suffix = '\0';
	}

	if (verbose) {
		info("Starting FLASH " VERSION_STR);
		info("Fast Length Adjustment of SHort reads");
		info(" ");
		info("Input files:");
		info("    %s", argv[0]);
		if (!interleaved_input)
			info("    %s", argv[1]);
		info(" ");
		info("Output files:");
		assert(!to_stdout);
		info("    %s.extendedFrags.fastq%s", name_buf, fops->suffix);
		if (interleaved_output) {
			info("    %s.notCombined.fastq%s", name_buf, fops->suffix);
		} else {
			info("    %s.notCombined_1.fastq%s", name_buf, fops->suffix);
			info("    %s.notCombined_2.fastq%s", name_buf, fops->suffix);
		}
		info("    %s.hist", name_buf);
		info("    %s.histogram", name_buf);
		info(" ");
		info("Parameters:");
		info("    Min overlap:          %d", min_overlap);
		info("    Max overlap:          %d", max_overlap);
		info("    Phred offset:         %d", phred_offset);
		info("    Combiner threads:     %d", num_combiner_threads);
		info("    Max mismatch density: %f", max_mismatch_density);
		info("    Output format:        %s", fops->name);
		info("    Interleaved input:    %s", interleaved_input ? "true" : "false");
		info("    Interleaved output:   %s", interleaved_output ? "true" : "false");
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

	struct threads threads;
	start_fastq_readers_and_writers(mates1_gzf, mates2_gzf, out_combined_fp,
					out_notcombined_fp_1,
					out_notcombined_fp_2, phred_offset,
					fops, &threads,
					num_combiner_threads,
					verbose);
	struct common_combiner_thread_params common = {
		.threads              = &threads,
		.min_overlap          = min_overlap,
		.max_overlap          = max_overlap,
		.max_mismatch_density = max_mismatch_density,
		.to_stdout            = to_stdout,
		.unqueue_lock         = PTHREAD_MUTEX_INITIALIZER,
		.queue_lock           = PTHREAD_MUTEX_INITIALIZER,
	};

	if (verbose)
		info("Starting %d combiner threads", num_combiner_threads);

	pthread_t other_combiner_threads[num_combiner_threads - 1];
	for (int i = 0; i < num_combiner_threads; i++) {
		struct combiner_thread_params *p;
		int ret;
		p = xmalloc(sizeof(struct combiner_thread_params));
		p->common = &common;
		p->combined_read_len_hist = &combined_read_len_hists[i];
		if (i < num_combiner_threads - 1) {
			ret = pthread_create(&other_combiner_threads[i],
					     NULL, combiner_thread_proc, p);
			if (ret != 0) {
				fatal_error_with_errno("Could not create "
						       "worker thread #%d",
						       i + 1);
			}
		} else {
			combiner_thread_proc(p);
		}
	}
	for (int i = 0; i < num_combiner_threads - 1; i++) {
		if (pthread_join(other_combiner_threads[i], NULL) != 0) {
			fatal_error_with_errno("Failed to join worker thread #%d",
					       i + 1);
		}
		hist_combine(combined_read_len_hist,
			     &combined_read_len_hists[i]);
		hist_destroy(&combined_read_len_hists[i]);
	}
	pthread_mutex_destroy(&common.unqueue_lock);
	pthread_mutex_destroy(&common.queue_lock);
	stop_fastq_readers_and_writers(&threads);

	/* The remainder is the same regardless of whether we are compiling for
	 * multiple threads or not (and we are down to one thread at this point
	 * anyway). */

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
		uint64_t start_usec = start_time.tv_sec * 1000000 + start_time.tv_usec;
		uint64_t end_usec = end_time.tv_sec * 1000000 + end_time.tv_usec;
		info(" ");
		info("FLASH " VERSION_STR " complete!");
		info("%.3f seconds elapsed", (double)(end_usec - start_usec) / 1000000);
	}
	if (out_suffix_allocated)
		free(out_suffix);
	return 0;
}
