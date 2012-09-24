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

#define VERSION_STR "v1.2.1"

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
"To run FLASH, you must provide two FASTQ files of paired-end reads.\n"
"Corresponding read pairs must be in the same order in both files.\n"
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
"  -z, --compress          Compress the FASTQ output files with zlib.\n"
"\n"
"  -t, --threads=NTHREADS  Set the number of worker threads.  This is in\n"
"                          addition to the I/O threads.  Ignored if support\n"
"                          for multiple threads was not compiled in.  Default:\n"
"                          number of processors.  Note: if you need FLASH's\n"
"                          output to appear deterministically or in the same\n"
"                          order as the original reads, you must specify\n"
"                          -t 1 (--threads=1).\n"
"\n"
"  -h, --help              Display this help and exit.\n"
"\n"
"  -v, --version           Display version.\n"
	;
	puts(usage_str);
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
}

static const char *optstring = "m:M:x:p:r:f:s:o:d:czt:qhv";
static const struct option longopts[] = {
	{"min-overlap",          required_argument,  NULL, 'm'}, 
	{"max-overlap",          required_argument,  NULL, 'M'}, 
	{"max-mismatch-density", required_argument,  NULL, 'x'}, 
	{"phred-offset",         required_argument,  NULL, 'p'}, 
	{"read_len",		 required_argument,  NULL, 'r'}, 
	{"fragment-len",         required_argument,  NULL, 'f'}, 
	{"fragment-len-stddev",  required_argument,  NULL, 's'}, 
	{"output-prefix",        required_argument,  NULL, 'o'}, 
	{"output-directory",     required_argument,  NULL, 'd'}, 
	{"to-stdout",            no_argument,        NULL, 'c'},
	{"compress",             no_argument,        NULL, 'z'},
	{"threads",              required_argument,  NULL, 't'},
	{"quiet",                no_argument,	     NULL, 'q'}, 
	{"help",                 no_argument,	     NULL, 'h'}, 
	{"version",              no_argument,	     NULL, 'v'}, 
	{NULL, 0, NULL, 0}
};


static void copy_tag(struct read *to, const struct read *from)
{
	if (to->tag_bufsz < from->tag_len) {
		to->tag = xrealloc(to->tag, from->tag_len);
		to->tag_bufsz = from->tag_len;
	}
	to->tag_len = from->tag_len;
	memcpy(to->tag, from->tag, from->tag_len);
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

#ifdef MULTITHREADED
static void hist_combine(struct histogram *hist, const struct histogram *other)
{
	for (size_t i = 0; i < other->len; i++)
		hist_add(hist, i, other->array[i]);
}
#endif

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
#ifdef MULTITHREADED

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
	struct histogram *combined_read_len_hist;
	struct common_combiner_thread_params *common;
};

/* In the multi-threaded configuration of FLASH, this procedure is executed in
 * parallel by all the combiner threads. */
static void *combiner_thread_proc(void *__params)
{
	struct combiner_thread_params *params = __params;

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
	/* This combiner thread has been signalled by the reader that there are
	 * no more reads.  We must send our pending combined and uncombined
	 * reads to the writers, terminated by a NULL read.  Actually, all the
	 * remaining reads in the uncombined and combined read sets should be
	 * set to NULL, as they have not been assigned yet, but they will be
	 * freed later.  
	 *
	 * We should also free the memory for the read_set_1 and read_set_2.  Up
	 * to but not including the read at index i, the reads have been moved
	 * to the uncombined or combined read sets, so they are taken care of
	 * and should be set to NULL.  The rest of the reads in these sets
	 * should be freed. */
	memset(&to_writer_1->reads[to_writer_filled], 0,
	       (READS_PER_READ_SET - to_writer_filled) * sizeof(struct read*));
	memset(&to_writer_2->reads[to_writer_filled], 0,
	       (READS_PER_READ_SET - to_writer_filled) * sizeof(struct read*));
	memset(&combined_read_set->reads[combined_read_filled], 0,
	       (READS_PER_READ_SET - combined_read_filled) * sizeof(struct read*));

	memset(&read_set_1->reads[0], 0, i * sizeof(struct read*));
	memset(&read_set_2->reads[0], 0, i * sizeof(struct read*));
	free_read_set(read_set_1);
	free_read_set(read_set_2);

	pthread_mutex_lock(&params->common->queue_lock);
	read_queue_put(threads->writer_uncombined_1.ready_q, to_writer_1);
	read_queue_put(threads->writer_uncombined_2.ready_q, to_writer_2);
	pthread_mutex_unlock(&params->common->queue_lock);

	read_queue_put(threads->writer_combined.ready_q, combined_read_set);
	free(params);
	return NULL;
}

#endif

int main(int argc, char **argv)
{
	int max_overlap            = 0;
	int min_overlap            = 10;
	float max_mismatch_density = 0.25;
	int phred_offset           = 33;
	int read_len               = 100;
	int fragment_len           = 180;
	int fragment_len_stddev    = 20;
	const char *prefix         = "out";
	const char *output_dir     = ".";
	bool to_stdout             = false;
	bool verbose               = true;
	gzFile mates1_gzf, mates2_gzf;
	void *out_combined_fp;
	void *out_notcombined_fp_1 = NULL;
	void *out_notcombined_fp_2 = NULL;
#ifdef MULTITHREADED
	int num_combiner_threads = 0;
#endif
	const struct file_operations *fops = &normal_fops;
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
		case 't':
			#ifdef MULTITHREADED
			num_combiner_threads = strtol(optarg, &tmp, 10);
			if (tmp == optarg || *tmp || num_combiner_threads < 1) {
				fatal_error("Number of threads must be "
					    "a positive integer!  Please "
					    "check option -t.");
			}
			#else
			warning("Multiple threads not enabled; ignoring -t "
				"option.");
			#endif
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

#ifdef MULTITHREADED
	if (num_combiner_threads == 0)
		num_combiner_threads = get_default_num_threads();
#endif

	argc -= optind;
	argv += optind;

	if (argc != 2) {
		usage_short();
		return 2;
	}

	mates1_gzf = xgzopen(argv[0], "r");
	mates2_gzf = xgzopen(argv[1], "r");

	mkdir_p(output_dir);

	char name_buf[strlen(output_dir) + 1 + strlen(prefix) +
		      sizeof(".notCombined_2.fastq") + sizeof(".gz")];
	char *suffix;
	suffix = name_buf + sprintf(name_buf, "%s/%s", output_dir, prefix);


	/* Open the output files. */
	if (to_stdout) {
		out_combined_fp = fops->open_file("-", "w");
	} else {
		sprintf(suffix, ".extendedFrags.fastq%s", fops->suffix);
		out_combined_fp = fops->open_file(name_buf, "w");

		sprintf(suffix, ".notCombined_1.fastq%s", fops->suffix);
		out_notcombined_fp_1 = fops->open_file(name_buf, "w");

		sprintf(suffix, ".notCombined_2.fastq%s", fops->suffix);
		out_notcombined_fp_2 = fops->open_file(name_buf, "w");
		*suffix = '\0';
	}

	if (verbose) {
		info("Starting FLASH " VERSION_STR);
		info("Fast Length Adjustment of SHort Reads");
		info(" ");
		info("Input files:");
		info("    %s", argv[0]);
		info("    %s", argv[1]);
		info(" ");
		info("Output files:");
		assert(!to_stdout);
		info("    %s.extendedFrags.fastq%s", name_buf, fops->suffix);
		info("    %s.notCombined_1.fastq%s", name_buf, fops->suffix);
		info("    %s.notCombined_2.fastq%s", name_buf, fops->suffix);
		info("    %s.hist", name_buf);
		info("    %s.histogram", name_buf);
		info(" ");
		info("Parameters:");
		info("    Min overlap:          %d", min_overlap);
		info("    Max overlap:          %d", max_overlap);
		info("    Phred offset:         %d", phred_offset);
	#ifdef MULTITHREADED
		info("    Combiner threads:     %d", num_combiner_threads);
	#endif
		info("    Max mismatch density: %f", max_mismatch_density);
		info("    Output format:        %s", fops->name);
		info(" ");
	}


	/* 
	 * We wish to do the following:
	 *
	 * "Go through each mate pair in the input files.  Determine if it can
	 * be combined, given the input parameters to the program.  If it can,
	 * write the combined read to the PREFIX.extendedFrags.fastq file.
	 * Otherwise, write the reads in the mate pair to the
	 * PREFIX.notCombined_1.fastq and PREFIX.notCombined_2.fastq files.  Or,
	 * if the -c / --to-stdout option is specified, write the combined reads
	 * to standard output, and ignore the uncombined reads."
	 *
	 * Two separate implementations of this follow: one for the
	 * single-threaded case, and one for the multi-threaded case.  In the
	 * multi-threaded case, there are @num_combiner_threads combiner
	 * threads created that will process the reads in parallel by retrieving
	 * `struct read_set'-sized chunks of reads from the reader threads, and
	 * providing `struct read_set'-sized chunks of combined or uncombined
	 * reads to the writer threads.
	 */

#ifdef MULTITHREADED

	/* Histogram of how many combined reads have a given length.
	 *
	 * The zero index slot counts how many reads were not combined.
	 *
	 * Here, there is a copy of the histogram for each thread, and they are
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
#else /* !MULTITHREADED */

	/* Histogram of how many combined reads have a given length.
	 *
	 * The zero index slot counts how many reads were not combined. */
	struct histogram _combined_read_len_hist;
	struct histogram *combined_read_len_hist = &_combined_read_len_hist;

	hist_init(combined_read_len_hist);

	unsigned long pair_no = 0;

	/* The single-threaded code is much simpler than the multi-threaded
	 * code; just have two `struct read's to keep using for the input reads,
	 * and another `struct read' to use for combined sequence.  */

	struct read read_1;
	struct read read_2;
	struct read combined_read;

	init_read(&read_1);
	init_read(&read_2);
	init_read(&combined_read);

	while (next_mate_pair(&read_1, &read_2,
			      mates1_gzf, mates2_gzf, phred_offset))
	{
		if (verbose && ++pair_no % 25000 == 0)
			info("Processed %lu reads", pair_no);
		if (combine_reads(&read_1, &read_2, &combined_read, min_overlap,
				  max_overlap, max_mismatch_density))
		{
			/* Combination was successful. */
			get_combined_tag(&read_1, &read_2, &combined_read);
			fops->write_read(&combined_read, out_combined_fp, 
				   	 phred_offset);
			hist_inc(combined_read_len_hist,
				 combined_read.seq_len);
		} else {
			/* Combination was unsuccessful. */
			if (!to_stdout) {
				fops->write_read(&read_1, out_notcombined_fp_1, 
						 phred_offset);
				reverse_complement(read_2.seq,
						   read_2.seq_len);
				reverse(read_2.qual, read_2.seq_len);
				fops->write_read(&read_2, out_notcombined_fp_2, 
						 phred_offset);
			}
			hist_inc(combined_read_len_hist, 0);
		}
	}
	if (verbose) {
		if (pair_no % 25000 != 0)
			info("Processed %lu reads", pair_no);
		info("Closing input and output FASTQ files");
	}
	gzclose(mates1_gzf);
	gzclose(mates2_gzf);
	fops->close_file(out_combined_fp);
	fops->close_file(out_notcombined_fp_1);
	fops->close_file(out_notcombined_fp_2);
#endif /* !MULTITHREADED */


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
	return 0;
}
