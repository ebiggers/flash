#include "fastq.h"
#include "util.h"
#include <errno.h>

static void init_read(struct read *read, int max_tag_len, 
		      int max_read_len)
{
	max_tag_len += 10;
	max_read_len += 10;
	read->tag_bufsz  = max_tag_len;
	read->read_bufsz = max_read_len;
	read->tag        = xmalloc(max_tag_len);
	read->seq        = xmalloc(max_read_len);
	read->separator  = xmalloc(max_read_len);
	read->qual       = xmalloc(max_read_len);

	read->tag[max_tag_len - 2] = '\0';
	read->seq[max_read_len - 2] = '\0';
	read->qual[max_read_len - 2] = '\0';
}

void init_mate_pair(struct mate_pair *pair, int max_tag_len, 
		    int max_read_len)
{
	init_read(&pair->read_1, max_tag_len, max_read_len);
	init_read(&pair->read_2, max_tag_len, max_read_len);
}

/* 
 * Reads the next read from the FASTQ file @mates_gzf into the @read structure.
 *
 * Whitespace is stripped from the end of the sequence, tag, and quality scores.
 *
 * The sequence is translated into only the characters A, C, G, T, and N, while
 * the quality values are re-scaled to start at 0.
 *
 * Returns true on success, false on EOF.  Aborts on read error.
 */
static bool next_read(struct read *read, gzFile mates_gzf, int phred_offset)
{
	int errnum;

	if (!gzgets(mates_gzf, read->tag, read->tag_bufsz))
		goto out_eofok;
	if (read->tag[read->tag_bufsz - 2] != '\0')
		goto line_too_long;

	if (!gzgets(mates_gzf, read->seq, read->read_bufsz))
		goto out_eofnotok;
	if (read->seq[read->read_bufsz - 2] != '\0')
		goto line_too_long;

	if (!gzgets(mates_gzf, read->separator, read->read_bufsz))
		goto out_eofnotok;
	if (read->separator[0] != '+')
		fatal_error("Expected '+' character!");

	if (!gzgets(mates_gzf, read->qual, read->read_bufsz))
		goto out_eofnotok;
	if (read->qual[read->read_bufsz - 2] != '\0')
		goto line_too_long;

	read->seq_len = trim(read->seq);
	trim(read->tag);

	if (trim(read->qual) != read->seq_len)
		fatal_error("Qual string length not the same as sequence "
			    "length!");

	for (int i = 0; i < read->seq_len; i++)
		read->seq[i] = canonical_ascii_char(read->seq[i]);

	for (int i = 0; i < read->seq_len; i++) {
		if (read->qual[i] < phred_offset) {
			fatal_error("Qual string contains character under "
				    "phred_offset = %d!", phred_offset);
		}
		read->qual[i] -= phred_offset;
	}
	return true;

out_eofnotok:
	if (gzeof(mates_gzf))
		fatal_error("Unexpected EOF reading input file!");
out_eofok:
	gzerror(mates_gzf, &errnum);
	if (errnum != Z_OK)
		fatal_error("Error reading input file!");
	return false;
line_too_long:
	fatal_error("Line too long!  Make sure you have specified the "
		    "read length correctly.");
}

/* Reads the next mate pair from the FASTQ files @mates1_gzf and @mates2_gzf.  
 *
 * Whitespace is stripped from the ends of the sequences, tags, and quality
 * scores.
 *
 * The sequences are translated into only the characters A, C, G, T, and N,
 * while the quality values are re-scaled to start at 0.
 *
 * The sequence of the second read is reverse complemented, so that the reads go
 * in the same direction rather than being oriented in opposite directions on
 * opposite strands.
 *
 * Returns true on success, false on EOF.  Aborts on read error.
 */
bool next_mate_pair(struct mate_pair *pair, gzFile mates1_gzf, 
		    gzFile mates2_gzf, int phred_offset)
{
	if (!next_read(&pair->read_1, mates1_gzf, phred_offset))
		return false;
	if (!next_read(&pair->read_2, mates2_gzf, phred_offset))
		return false;
	reverse_complement(pair->read_2.seq, pair->read_2.seq_len);
	reverse(pair->read_2.qual, pair->read_2.seq_len);
	return true;
}

/* 
 * Writes a read to an output uncompressed FASTQ file.
 *
 * @read:
 * 	The read to write (sequence, tag, and quality values), where the quality
 * 	values are scaled to start at 0.
 * @fp:
 * 	FASTQ file opened for writing.
 * @phred_offset:
 * 	Offset for quality values in the output.  
 *
 * Aborts on write error.
 */
void write_read(struct read *read, FILE *fp, int phred_offset)
{
	int i;
	int ret;
	for (i = 0; i < read->seq_len; i++)
		read->qual[i] += phred_offset;

	ret = fprintf(fp, "%s\n%s\n+\n%s\n", read->tag, read->seq, read->qual);
	if (ret < 0)
		fatal_error("Error writing to output file: %s", 
			    strerror(errno));
}
