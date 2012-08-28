#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <zlib.h>

#include "combine_reads.h"
#include "fastq.h"
#include "util.h"

#define MAX_TAG_LEN 1024

static void usage()
{
	const char *usage_str = 
"Usage: flash [OPTIONS] <mates1.fastq> <mates2.fastq>\n"
"\n"
"DESCRIPTION:\n"
"\n"
"FLASH (Fast Length Adjustment of Short reads) is an accurate and fast tool\n"
"to merge paired-end reads that were generated from DNA fragments whose\n"
"lengths are shorter than twice the length of reads. Merged read pairs result\n"
"in unpaired longer reads. Longer reads are generally more desired in genome\n"
"assembly and genome analysis processes.\n"
"\n"
"MANDATORY INPUT:\n"
"\n"
"Two fastq files of paired-end reads. Corresponding pairs must be in the same\n"
"order in both files.\n"
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
"                          standard deviation.\n"
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
"  -o, --output-prefix=PREFIX\n"
"                          Prefix of output files.  Default: \"out\".\n"
"\n"
"  -d, --output-directory=DIR\n"
"                          Path to directory for output files.  Default:\n"
"                          current working directory.\n"
"\n"
"  -r, --read-len=LEN      Average read length.  Default: 100.\n"
"\n"
"  -f, --fragment-len=LEN  Average fragment length.  Default: 180.\n"
"\n"
"  -s, --fragment-len-stddev=LEN\n"
"                          Standard deviation of fragment lengths.  If you do\n"
"                          not know standard deviation of the fragment\n"
"                          library, you can probably assume that the standard\n"
"                          deviation is 10% of the average fragment length.\n"
"                          Default: 20.\n"
"\n"
"  -h, --help              Display this help and exit.\n"
"  -v, --version           Display version.\n"
	;
	fputs(usage_str, stdout);
}

static void usage_short()
{
	puts("Usage: flash [OPTIONS] <mates1.fastq> <mates2.fastq>");
	puts("Try `flash -h' for more information.");
}

static void version()
{
	puts("flash v1.1");
}

static const char *optstring = "m:M:x:p:o:d:r:f:s:hv";
static const struct option longopts[] = {
	{"min-overlap",          required_argument,  NULL, 'm'}, 
	{"max-overlap",          required_argument,  NULL, 'M'}, 
	{"max-mismatch-density", required_argument,  NULL, 'x'}, 
	{"phred-offset",         required_argument,  NULL, 'p'}, 
	{"output-prefix",        required_argument,  NULL, 'o'}, 
	{"output-directory",     required_argument,  NULL, 'd'}, 
	{"read-len",             required_argument,  NULL, 'r'}, 
	{"fragment-len",         required_argument,  NULL, 'f'}, 
	{"fragment-len-stddev",  required_argument,  NULL, 's'}, 
	{"help",                 required_argument,  NULL, 'h'}, 
	{"version",              required_argument,  NULL, 'v'}, 
	{NULL, 0, NULL, 0}
};


static gzFile safe_gzopen(const char *filename, const char *mode)
{
	gzFile f = gzopen(filename, mode);
	if (!f) {
		fatal_error("Failed to open the file \"%s\": %s", filename,
				((errno == 0) ? 
				 	"out of memory" : strerror(errno)));
	}
	return f;
}

/* Modifies @tag_1 to contain the new tag for the combined read of a fragment.
 * */
static void get_combined_tag(char *tag_1, const char *tag_2)
{
	if (strcmp(tag_1, tag_2) != 0) {
		char *p = strrchr(tag_1, '/');
		if (p) {
			if (*(p + 1) != '\0' && *(p + 2) == '#') {
				/* read ID has a barcode. */
				do {
					*p = *(p + 2);
				} while (*(++p + 2) != '\0');
			}
			*p = '\0';
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
				fatal_error("Error writing to \"%s\": %s\n",
					    hist_file, strerror(errno));
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
	fatal_error("Error writing to file \"%s\": %s", histogram_file, 
		    strerror(errno));
}

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
	gzFile mates1_gzf, mates2_gzf;
	int c;
	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
		switch (c) {
		case 'm':
			min_overlap = atoi(optarg);
			if (min_overlap < 1)
				fatal_error("Minimum overlap must be positive!");
			break;
		case 'M':
			max_overlap = atoi(optarg);
			max_overlap_specified = true;
			if (max_overlap < 1)
				fatal_error("Maximum overlap must be positive!");
			break;
		case 'x':
			max_mismatch_density = atof(optarg);
			if (max_mismatch_density < 0.0 || 
			    max_mismatch_density > 1.0)
				fatal_error("Max mismatch density has to be in "
					    "the interval [0, 1]!");
			break;
		case 'p':
			phred_offset = atoi(optarg);
			if (phred_offset != 33 && phred_offset != 64) {
				warning("Phred offset is usually either "
				      "64 (for earlier Illumina data) or 33 "
				      "(for Sanger and later Illumina data).");
			}
			break;
		case 'o':
			prefix = optarg;
			break;
		case 'd':
			output_dir = optarg;
			break;
		case 'r':
			read_len = atoi(optarg);
			if (read_len <= 0)
				fatal_error("Read length must be positive!");
			break;
		case 'f':
			fragment_len = atoi(optarg);
			if (fragment_len <= 0)
				fatal_error("Fragment length must be positive!");
			break;
		case 's':
			fragment_len_stddev = atoi(optarg);
			if (fragment_len_stddev <= 0)
				fatal_error("Fragment length standard deviation "
					    "must be positive!");
			break;
		case 'v':
			version();
			return 0;
		case 'h':
			usage();
			return 0;
		default:
			usage_short();
			return 1;
		}
	}



	if (!max_overlap_specified && ((read_len != 100) || 
				      (fragment_len != 180) || 
				      (fragment_len_stddev != 20))) {
		max_overlap = (int)(2 * read_len - fragment_len +
				    2.5 * fragment_len_stddev);
	}

	argc -= optind - 1;
	argv += optind - 1;

	if (argc != 3) {
		usage_short();
		return 1;
	}

	mates1_gzf = safe_gzopen(argv[1], "r");
	mates2_gzf = safe_gzopen(argv[2], "r");

	mkdir_p(output_dir);

	char name_buf[strlen(output_dir) + 1 + sizeof(".notCombined_2.fastq")];
	char *suffix;
	struct mate_pair pair;
	char combined_seq[read_len * 2 + 1];
	char combined_qual[read_len * 2 + 1];
	struct read combined_read;

	suffix = name_buf + sprintf(name_buf, "%s/%s", output_dir, prefix);

	/* Allocate enough memory for two reads. */
	init_mate_pair(&pair, MAX_TAG_LEN, read_len);

	/* Set up a `struct read' that will serve as the combined read that will
	 * be written to the output file when needed. */
	combined_read.seq = combined_seq;
	combined_read.qual = combined_qual;
	combined_read.tag = pair.read_1.tag;

	/* Open the output files. */
	strcpy(suffix, ".extendedFrags.fastq");
	FILE *out_combined_fp = xfopen(name_buf, "w");

	strcpy(suffix, ".notCombined_1.fastq");
	FILE *out_notcombined_fp_1 = xfopen(name_buf, "w");

	strcpy(suffix, ".notCombined_2.fastq");
	FILE *out_notcombined_fp_2 = xfopen(name_buf, "w");

	/* Histogram of how many combined reads have a given length.
	 *
	 * The zero index slot counts how many reads were not combined. */
	size_t combined_read_len_hist[read_len * 2];
	ZERO_ARRAY(combined_read_len_hist);

	/* Go through each mate pair in the input file.  Determine if it can be
	 * combined with the provided paramaters.  If it can, write the combined
	 * read to the PREFIX.extendedFrags.fastq file.  Otherwise, write the
	 * reads in the mate pair to the PREFIX.notCombined_1.fastq and
	 * PREFIX.notCombined_2.fastq files. */
	while (next_mate_pair(&pair, mates1_gzf, mates2_gzf, phred_offset)) {
		size_t combined_len;

		combined_len = combine_reads(&pair, combined_seq, 
					     combined_qual, min_overlap, 
					     max_overlap, max_mismatch_density);
		if (combined_len != 0) {
			/* Combination was successful. */
			get_combined_tag(pair.read_1.tag, pair.read_2.tag);
			combined_read.seq_len = combined_len;
			write_read(&combined_read, out_combined_fp, 
				   phred_offset);
		} else {
			/* Combination was unsuccessful. */
			write_read(&pair.read_1, out_notcombined_fp_1, 
				   phred_offset);
			reverse_complement(pair.read_2.seq, 
					   pair.read_2.seq_len);
			reverse(pair.read_2.qual, pair.read_2.seq_len);
			write_read(&pair.read_2, out_notcombined_fp_2, 
				   phred_offset);
		}
		combined_read_len_hist[combined_len]++;
	}

	gzclose(mates1_gzf);
	gzclose(mates2_gzf);
	xfclose(out_combined_fp);
	xfclose(out_notcombined_fp_1);
	xfclose(out_notcombined_fp_2);

	/* From the histogram table @combined_read_len_hist, compute:
	 * - The lowest length for which there exists a combined read of that
	 *   length.
	 * - The highest length for which there exists a combined read of that
	 *   length.
	 * - The maximum frequency of a read length in the set of combined
	 *   reads.
	 */
	size_t max_freq = 0;
	int first_nonzero_idx = 0;
	int last_nonzero_idx = 0;
	for (int i = 1; i < ARRAY_LEN(combined_read_len_hist); i++) {
		if (combined_read_len_hist[i] != 0) {
			if (first_nonzero_idx == 0)
				first_nonzero_idx = i;
			last_nonzero_idx = i;
		}
		if (combined_read_len_hist[i] > max_freq)
			max_freq = combined_read_len_hist[i];
	}
	
	/* Write the raw numbers of the combined read length histogram to the
	 * PREFIX.hist file. */
	strcpy(suffix, ".hist");
	write_hist_file(name_buf, combined_read_len_hist, first_nonzero_idx,
			last_nonzero_idx);

	/* Write a pretty representation of the combined read length histogram
	 * to the PREFIX.histogram file. */
	strcpy(suffix, ".histogram");
	write_histogram_file(name_buf, combined_read_len_hist, 
			     first_nonzero_idx, last_nonzero_idx, 
			     max_freq);
	return 0;
}
