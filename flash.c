#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "combine_reads.h"
#include "fastq.h"
#include "util.h"

#define MAX_TAG_LEN 1024

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
"  -r, --read-len=LEN      Average read length.  Default: 100.\n"
"\n"
"  -f, --fragment-len=LEN  Average fragment length.  Default: 180.\n"
"\n"
"  -s, --fragment-len-stddev=LEN\n"
"                          Standard deviation of fragment lengths.  If you do\n"
"                          not know the standard deviation of the fragment\n"
"                          library, you can probably assume that the standard\n"
"                          deviation is 10% of the average fragment length.\n"
"                          Default: 20.\n"
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
"                          output to appear deterministically or in the some\n"
"                          order as the original reads, you must specify\n"
"                          -t 1 (--threads=1).\n"
"\n"
"  -h, --help              Display this help and exit.\n"
"\n"
"  -v, --version           Display version."
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
	puts("flash v1.2");
}

static const char *optstring = "m:M:x:p:r:f:s:o:d:czt:hv";
static const struct option longopts[] = {
	{"min-overlap",          required_argument,  NULL, 'm'}, 
	{"max-overlap",          required_argument,  NULL, 'M'}, 
	{"max-mismatch-density", required_argument,  NULL, 'x'}, 
	{"phred-offset",         required_argument,  NULL, 'p'}, 
	{"read-len",             required_argument,  NULL, 'r'}, 
	{"fragment-len",         required_argument,  NULL, 'f'}, 
	{"fragment-len-stddev",  required_argument,  NULL, 's'}, 
	{"output-prefix",        required_argument,  NULL, 'o'}, 
	{"output-directory",     required_argument,  NULL, 'd'}, 
	{"to-stdout",            no_argument,        NULL, 'c'},
	{"compress",             no_argument,        NULL, 'z'},
	{"threads",              required_argument,  NULL, 't'},
	{"help",                 no_argument,	     NULL, 'h'}, 
	{"version",              no_argument,	     NULL, 'v'}, 
	{NULL, 0, NULL, 0}
};


/* 
 * Given the tags of the two reads, find the tag that will be given to the
 * combined read.
 *
 * We need to strip off what trails the '/' (e.g. "/1" and "/2"), unless there
 * is a "barcode" beginning with the '#' character, which is kept.
 *
 * Instead of copying the tag, we transform @tag_1 into the combined tag.
 */
static void get_combined_tag(char *tag_1, int *tag_1_len, const char *tag_2,
			     int tag_2_len)
{
	int len = *tag_1_len;

	if (len != tag_2_len) /* Tags are the same; don't change it. */
		return;

	for (char *p = &tag_1[len - 1]; p >= tag_1; p--) {
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
			*tag_1_len = p - tag_1;
			return;
		}
	}
}

static void write_hist_file(const char *hist_file, const size_t hist[],
			    int first_nonzero_idx, int last_nonzero_idx)
{
	FILE *fp = xfopen(hist_file, "w");
	for (int i = first_nonzero_idx; i <= last_nonzero_idx; i++) {
		if (hist[i] != 0) {
			if (fprintf(fp, "%d\t%zu\n", i, hist[i]) < 0) {
				fatal_error_with_errno("Error writing to "
						       "the file \"%s\"",
						       hist_file);
			}
		}
	}
	xfclose(fp);
}

static void write_histogram_file(const char *histogram_file, 
				 const size_t hist[], int first_nonzero_idx,
				 int last_nonzero_idx, size_t max_freq)
{
	const double max_num_asterisks = 72;
	double scale = max_num_asterisks / (double)max_freq;

	FILE *fp = xfopen(histogram_file, "w");

	for (int i = first_nonzero_idx; i <= last_nonzero_idx; i++) {
		if (fprintf(fp, "%d\t", i) < 0)
			goto write_error;
		size_t num_asterisks = (size_t)(scale * (double)hist[i]);
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
	size_t *combined_read_len_hist;
	struct common_combiner_thread_params *common;
};

/* In the multi-threaded configuration of FLASH, this procedure is executed in
 * parallel by all the combiner threads. */
static void *combiner_thread_proc(void *__params)
{
	struct combiner_thread_params *params = __params;

	size_t *combined_read_len_hist = params->combined_read_len_hist;
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
	size_t combined_len;
	struct read *combined_read;
	struct read *read_1;
	struct read *read_2;
	unsigned i;
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
			combined_len = combine_reads(read_1, read_2,
						     combined_read->seq, 
						     combined_read->qual,
						     min_overlap, max_overlap,
						     max_mismatch_density);
			if (combined_len != 0) {
				/* Combination was successful. */

				/* Swap the tags of the read_1 and the
				 * combined_read so that we can set the
				 * combined_read tag by transforming read_1's
				 * tag in-place. 
				 *
				 * Note: this should only be done if the
				 * tag_bufsz's of the reads are the same. */
				char *tmp = combined_read->tag;
				combined_read->tag = read_1->tag;
				combined_read->tag_len = read_1->tag_len;
				combined_read->seq_len = combined_len;
				read_1->tag = tmp;
				get_combined_tag(combined_read->tag,
						 &combined_read->tag_len,
						 read_2->tag,
						 read_2->tag_len);
				if (++combined_read_filled == READS_PER_READ_SET) {
					read_queue_put(threads->writer_combined.ready_q,
						       combined_read_set);
					combined_read_set = read_queue_get(
							threads->writer_combined.free_q);
					combined_read_filled = 0;
				}
			} 

			if (combined_len != 0 || to_stdout) {
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
			combined_read_len_hist[combined_len]++;
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
	int max_overlap            = 70;
	bool max_overlap_specified = false;
	int min_overlap            = 10;
	float max_mismatch_density = 0.25;
	int phred_offset           = 33;
	int read_len               = 100;
	int fragment_len           = 180;
	int fragment_len_stddev    = 20;
	const char *prefix         = "out";
	const char *output_dir     = ".";
	bool to_stdout		   = false;
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

	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
		switch (c) {
		case 'm':
			min_overlap = strtol(optarg, &tmp, 10);
			if (tmp == optarg || min_overlap < 1)
				fatal_error("Minimum overlap must be a "
					    "positive integer!  Please check "
					    "option -m.");
			break;
		case 'M':
			max_overlap = strtol(optarg, &tmp, 10);
			max_overlap_specified = true;
			if (tmp == optarg || max_overlap < 1)
				fatal_error("Maximum overlap must be "
					    "a positive integer!  Please check "
					    "option -M.");
			break;
		case 'x':
			max_mismatch_density = strtod(optarg, &tmp);
			if (tmp == optarg || max_mismatch_density < 0.0 || 
			    max_mismatch_density > 1.0)
			{
				fatal_error("Max mismatch density must be a "
					    "number in the interval [0, 1]! "
					    "Please check option -x.");
			}
			break;
		case 'p':
			phred_offset = strtol(optarg, &tmp, 10);
			if (tmp == optarg || phred_offset < 0 
					  || phred_offset > 127) 
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
		case 'r':
			read_len = strtol(optarg, &tmp, 10);
			if (tmp == optarg || read_len <= 0)
				fatal_error("Read length must be a positive "
					    "integer!  Please check option -r.");
			break;
		case 'f':
			fragment_len = strtol(optarg, &tmp, 10);
			if (tmp == optarg || fragment_len <= 0)
				fatal_error("Fragment length must be a "
					    "positive integer!  Please check "
					    "option -f.");
			break;
		case 's':
			fragment_len_stddev = strtol(optarg, &tmp, 10);
			if (tmp == optarg || fragment_len_stddev <= 0)
				fatal_error("Fragment length standard deviation "
					    "must be a positive integer!  "
					    "Please check option -s.");
			break;
		case 'o':
			prefix = optarg;
			break;
		case 'd':
			output_dir = optarg;
			break;
		case 'c':
			to_stdout = true;
			break;
		case 'z':
			fops = &gzip_fops;
			break;
		case 't':
			#ifdef MULTITHREADED
			num_combiner_threads = strtol(optarg, &tmp, 10);
			if (tmp == optarg || num_combiner_threads < 1) {
				fatal_error("Number of threads must be "
					    "a positive integer!  Please "
					    "check option -t.");
			}
			#else
			warning("Multiple threads not enabled; ignoring -t "
				"option.");
			#endif
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



	if (!max_overlap_specified && ((read_len != 100) || 
				      (fragment_len != 180) || 
				      (fragment_len_stddev != 20))) {
		max_overlap = (int)(2 * read_len - fragment_len +
				    2.5 * fragment_len_stddev);
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

	char name_buf[strlen(output_dir) + 1 +
		      sizeof(".notCombined_2.fastq") + sizeof(".gz")];
	char *suffix;
	suffix = name_buf + sprintf(name_buf, "%s/%s", output_dir, prefix);


	/* Open the output files. */
	if (to_stdout) {
		out_combined_fp = fops->open_file("-", "w");
	} else {
		strcpy(suffix, ".extendedFrags.fastq");
		if (fops == &gzip_fops)
			strcat(suffix, ".gz");
		out_combined_fp = fops->open_file(name_buf, "w");

		strcpy(suffix, ".notCombined_1.fastq");
		if (fops == &gzip_fops)
			strcat(suffix, ".gz");
		out_notcombined_fp_1 = fops->open_file(name_buf, "w");

		strcpy(suffix, ".notCombined_2.fastq");
		if (fops == &gzip_fops)
			strcat(suffix, ".gz");
		out_notcombined_fp_2 = fops->open_file(name_buf, "w");
	}


	/* 
	 * We wish to do the following:
	 *
	 * "Go through each mate pair in the input file.  Determine if it can be
	 * combined, given the input parameters to the program.  If it can,
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
	size_t _combined_read_len_hist[num_combiner_threads][read_len * 2];
	ZERO_ARRAY(_combined_read_len_hist);
	size_t *combined_read_len_hist = _combined_read_len_hist[num_combiner_threads - 1];

	struct threads threads;
	start_fastq_readers_and_writers(mates1_gzf, mates2_gzf, out_combined_fp,
					out_notcombined_fp_1,
					out_notcombined_fp_2, phred_offset,
					MAX_TAG_LEN, read_len,
					fops, &threads,
					num_combiner_threads);
	struct common_combiner_thread_params common = {
		.threads              = &threads,
		.min_overlap          = min_overlap,
		.max_overlap          = max_overlap,
		.max_mismatch_density = max_mismatch_density,
		.to_stdout            = to_stdout,
		.unqueue_lock         = PTHREAD_MUTEX_INITIALIZER,
		.queue_lock           = PTHREAD_MUTEX_INITIALIZER,
	};

	pthread_t other_combiner_threads[num_combiner_threads - 1];
	for (int i = 0; i < num_combiner_threads; i++) {
		struct combiner_thread_params *p;
		int ret;
		p = xmalloc(sizeof(struct combiner_thread_params));
		p->common = &common;
		p->combined_read_len_hist = _combined_read_len_hist[i];
		if (i < num_combiner_threads - 1) {
			ret = pthread_create(&other_combiner_threads[i],
					     NULL, combiner_thread_proc, p);
			if (ret != 0) {
				fatal_error_with_errno("Could not create "
						       "worker thread #%d\n",
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
		for (unsigned j = 0; j < read_len * 2; j++)
			combined_read_len_hist[j] += _combined_read_len_hist[i][j];
	}
	pthread_mutex_destroy(&common.unqueue_lock);
	pthread_mutex_destroy(&common.queue_lock);
	stop_fastq_readers_and_writers(&threads);
#else /* !MULTITHREADED */

	/* Histogram of how many combined reads have a given length.
	 *
	 * The zero index slot counts how many reads were not combined. */
	size_t combined_read_len_hist[read_len * 2];
	ZERO_ARRAY(combined_read_len_hist);

	/* The single-threaded code is much simpler than the multi-threaded
	 * code; just have two `struct read's to keep using for the input reads,
	 * and another `struct read' to use for combined sequence.  The `struct
	 * read_queue' and `struct read_set' are not used at all. */

	struct read _read_1;
	struct read _read_2;
	struct read *read_1 = &_read_1;
	struct read *read_2 = &_read_2;
	char combined_seq[read_len * 2 + 1];
	char combined_qual[read_len * 2 + 1];
	struct read combined_read;

	/* Allocate enough memory for two reads. */
	init_read(read_1, MAX_TAG_LEN, read_len);
	init_read(read_2, MAX_TAG_LEN, read_len);

	/* Set up a `struct read' that will serve as the combined read that will
	 * be written to the output file when needed. */
	combined_read.seq = combined_seq;
	combined_read.qual = combined_qual;
	combined_read.tag = read_1->tag;
	while (next_mate_pair(read_1, read_2, mates1_gzf, mates2_gzf, phred_offset)) {
		size_t combined_len;

		combined_len = combine_reads(read_1, read_2, combined_seq, 
					     combined_qual, min_overlap, 
					     max_overlap, max_mismatch_density);
		if (combined_len != 0) {
			/* Combination was successful. */
			combined_read.tag_len = read_1->tag_len;
			combined_read.seq_len = combined_len;
			get_combined_tag(combined_read.tag, &combined_read.tag_len,
					 read_2->tag, read_2->tag_len);
			fops->write_read(&combined_read, out_combined_fp, 
				   	 phred_offset);
		} else if (!to_stdout) {
			/* Combination was unsuccessful. */
			fops->write_read(read_1, out_notcombined_fp_1, 
				   phred_offset);
			reverse_complement(read_2->seq, read_2->seq_len);
			reverse(read_2->qual, read_2->seq_len);
			fops->write_read(read_2, out_notcombined_fp_2, 
				   	 phred_offset);
		}
		combined_read_len_hist[combined_len]++;
	}
	gzclose(mates1_gzf);
	gzclose(mates2_gzf);
	fops->close_file(out_combined_fp);
	fops->close_file(out_notcombined_fp_1);
	fops->close_file(out_notcombined_fp_2);
#endif /* !MULTITHREADED */
	

	/* The remainder is the same regardless of whether we are compiling for
	 * multiple threads or not (we are down to one thread at this point
	 * anyway). */
	if (!to_stdout) {
		/* From the histogram table @combined_read_len_hist, compute:
		 * - The lowest length for which there exists a combined read of
		 *   that length.
		 * - The highest length for which there exists a combined read
		 *   of that length.
		 * - The maximum frequency of a read length in the set of
		 *   combined reads.
		 */
		size_t max_freq = 0;
		int first_nonzero_idx = -1;
		int last_nonzero_idx = 0;
		for (unsigned i = 1; i < read_len * 2; i++) {
			if (combined_read_len_hist[i] != 0) {
				if (first_nonzero_idx == -1)
					first_nonzero_idx = i;
				last_nonzero_idx = i;
			}
			if (combined_read_len_hist[i] > max_freq)
				max_freq = combined_read_len_hist[i];
		}

		if (first_nonzero_idx == -1)
			last_nonzero_idx = -2;

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
	return 0;
}
