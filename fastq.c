#include <stdio.h>
#include <stdlib.h>

#include "fastq.h"
#include "util.h"

static ssize_t gzip_getline(gzFile gz_fp, char **lineptr, size_t *n)
{
	size_t offset = 0;
	ssize_t ret;
	if (!*lineptr) {
		*n = 128;
		*lineptr = (char*)xmalloc(*n);
	}
	while (1) {
		char *line = *lineptr;
		if (gzgets(gz_fp, line + offset, *n - offset)) {
			ret = strlen(line);
			if (line[ret - 1] == '\n')
				return ret;
		} else {
			const char *error_str;
			int errnum;

			if (gzeof(gz_fp))
				return -1;
			error_str = gzerror(gz_fp, &errnum);
			if (errnum == Z_ERRNO)
				fatal_error_with_errno("Error reading reads file");
			else
				fatal_error("zlib error while reading reads "
					    "file: %s", error_str);
		}
		offset = *n - 1;
		*n *= 2;
		*lineptr = xrealloc(line, *n);
	}
}

/* 
 * Reads the next read from the FASTQ file @mates_gzf into the @read structure.
 *
 * Whitespace is stripped from the end of the sequence, tag, and quality scores.
 *
 * The sequence is translated into only the characters A, C, G, T, and N, while
 * the quality values are re-scaled to start at 0.
 *
 * Returns true on success, false on EOF.  Aborts on read error, or if the data
 * is invalid.
 */
static bool next_read(struct read *read, gzFile mates_gzf, int phred_offset)
{
	ssize_t tag_len, seq_len, qual_len;

	tag_len = gzip_getline(mates_gzf, &read->tag, &read->tag_bufsz);
	if (tag_len == -1)
		return false;

	seq_len = gzip_getline(mates_gzf, &read->seq, &read->seq_bufsz);
	if (seq_len == -1)
		fatal_error("Unexpected EOF reading input file!");

	qual_len = gzip_getline(mates_gzf, &read->qual, &read->qual_bufsz);
	if (qual_len == -1)
		fatal_error("Unexpected EOF reading input file!");

	if (read->qual[0] != '+')
		fatal_error("Expected '+' character in FASTQ separator!");

	qual_len = gzip_getline(mates_gzf, &read->qual, &read->qual_bufsz);
	if (qual_len == -1)
		fatal_error("Unexpected EOF reading input file!");

	seq_len = trim(read->seq, seq_len);
	tag_len = trim(read->tag, tag_len);
	qual_len = trim(read->qual, qual_len);

	if (qual_len != seq_len)
		fatal_error("Qual string length (%zu) not the same as sequence "
			    "length (%zu)!", qual_len, seq_len);

	for (size_t i = 0; i < seq_len; i++)
		read->seq[i] = canonical_ascii_char(read->seq[i]);

	for (size_t i = 0; i < seq_len; i++) {
		if (read->qual[i] < phred_offset) {
			fatal_error("Qual string contains character under "
				    "phred_offset = %d!", phred_offset);
		}
		read->qual[i] -= phred_offset;
	}
	read->seq_len = seq_len;
	read->tag_len = tag_len;
	return true;
}


static void sprint_read(char buf[], const struct read *read, int phred_offset)
{
	char *p = buf;
	memcpy(p, read->tag, read->tag_len);
	p += read->tag_len;
	*p++ = '\n';
	memcpy(p, read->seq, read->seq_len);
	p += read->seq_len;
	*p++ = '\n';
	*p++ = '+';
	*p++ = '\n';
	for (int i = 0; i < read->seq_len; i++)
		p[i] = read->qual[i] + phred_offset;
	p += read->seq_len;
	*p++ = '\n';
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
void write_read_uncompressed(struct read *read, void *_fp, int phred_offset)
{
	FILE *fp = (FILE*)_fp;
	size_t size = read->tag_len + 1 + read->seq_len + 3 + read->seq_len + 1;
	char buf[size];
	sprint_read(buf, read, phred_offset);

	size_t ret = fwrite(buf, 1, size, fp);
	if (ret != size)
		fatal_error_with_errno("Error writing to output file");
}



/* 
 * Writes a read to an output zlib-compressed FASTQ file.
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
void write_read_compressed(struct read *read, void *_fp, int phred_offset)
{
	gzFile fp = (gzFile)_fp;

	size_t size = read->tag_len + 1 + read->seq_len + 3 + read->seq_len + 1;
	char buf[size];
	sprint_read(buf, read, phred_offset);

	int ret = gzwrite(fp, buf, size);
	if (ret != size) {
		int errnum;
		const char *err_str = gzerror(fp, &errnum);
		if (errnum == Z_ERRNO)
			fatal_error_with_errno("Error writing to output file");
		else
			fatal_error("zlib error writing to output file: %s",
				    err_str);
	}
}



const struct file_operations gzip_fops = {
	.name       = "gzip",
	.suffix     = ".gz",
	.write_read = write_read_compressed,
	.open_file  = xgzopen,
	.close_file = xgzclose,
};

const struct file_operations normal_fops = {
	.name       = "text",
	.suffix     = "",
	.write_read = write_read_uncompressed,
	.open_file  = xfopen,
	.close_file = xfclose,
};

#ifdef MULTITHREADED

static struct read *new_read()
{
	struct read *r = xmalloc(sizeof (struct read));
	init_read(r);
	return r;
}

static void free_read(struct read *r)
{
	if (r) {
		destroy_read(r);
		free(r);
	}
}

static struct read_set *new_read_set()
{
	struct read_set *p = new_empty_read_set();
	for (size_t i = 0; i < READS_PER_READ_SET; i++)
		p->reads[i] = new_read();
	return p;
}


void free_read_set(struct read_set *p)
{
	if (p) {
		for (size_t i = 0; i < READS_PER_READ_SET; i++)
			free_read(p->reads[i]);
		free(p);
	}
}

#include <semaphore.h>

/* 
 * Producer-consumer queue; it holds pointers to `struct read_sets', which can
 * be added or removed from the queue in a thread-safe manner using
 * read_queue_put() and read_queue_get(), respectively.
 */
struct read_queue {
	sem_t		  q_filled_slots;
	sem_t		  q_empty_slots;
	sem_t		  q_lock;
	unsigned	  q_front;
	unsigned	  q_back;
	struct read_set **q_reads;
	size_t	          q_size;
};

static struct read_queue *new_read_queue(size_t size, bool full)
{
	struct read_queue *q = xmalloc(sizeof(struct read_queue));
	q->q_reads = xmalloc(size * sizeof(struct read_set*));
	if (full)
		for (size_t i = 0; i < size; i++)
			q->q_reads[i] = new_read_set();
	else
		for (size_t i = 0; i < size; i++)
			q->q_reads[i] = NULL;

	sem_init(&q->q_filled_slots, 0, (full) ? size : 0);
	sem_init(&q->q_empty_slots, 0, (full) ? 0 : size);
	sem_init(&q->q_lock, 0, 1);
	q->q_front = 0;
	q->q_back = size - 1;
	q->q_size = size;
	return q;
}

static void free_read_queue(struct read_queue *q)
{
	if (q) {
		int filled_slots;
		sem_getvalue(&q->q_filled_slots, &filled_slots);
		sem_destroy(&q->q_filled_slots);
		sem_destroy(&q->q_empty_slots);
		sem_destroy(&q->q_lock);
		int i = q->q_front;
		while (filled_slots-- > 0) {
			free_read_set(q->q_reads[i]);
			i = (i + 1) % q->q_size;
		}
		free(q->q_reads);
		free(q);
	}
}

/* Retrieve the next available read set from the queue, blocking until one is
 * available. */
struct read_set *read_queue_get(struct read_queue *q)
{
	sem_wait(&q->q_filled_slots);
	sem_wait(&q->q_lock);

	struct read_set *r = q->q_reads[q->q_front];
	q->q_front = (q->q_front + 1) % q->q_size;

	sem_post(&q->q_empty_slots);
	sem_post(&q->q_lock);
	return r;
}

/* Put a read set into the queue, blocking until there is an empty space
 * available. */
void read_queue_put(struct read_queue *q, struct read_set *r)
{
	sem_wait(&q->q_empty_slots);
	sem_wait(&q->q_lock);

	q->q_back = (q->q_back + 1) % q->q_size;
	q->q_reads[q->q_back] = r;

	sem_post(&q->q_filled_slots);
	sem_post(&q->q_lock);
}


/* Parameters structure to pass to the reader and writer threads */
struct reader_writer_params {
	struct read_queue *free_q;
	struct read_queue *ready_q;
	void *fp;
	const struct file_operations *fops;
	int phred_offset;
	int num_combiner_threads;
	bool verbose;
};

static struct reader_writer_params *
new_reader_writer_params(size_t q_size, void *fp,
			 const struct file_operations *fops, int phred_offset,
			 struct read_queue *free_q, int num_combiner_threads,
			 bool verbose)
{
	struct reader_writer_params *p;
	
	p = xmalloc(sizeof(struct reader_writer_params));
	if (free_q)
		p->free_q = free_q;
	else
		p->free_q = new_read_queue(q_size, true);
	p->ready_q = new_read_queue(q_size, false);
	p->fp = fp;
	p->fops = fops;
	p->phred_offset = phred_offset;
	p->num_combiner_threads = num_combiner_threads;
	p->verbose = verbose;
	return p;
}

/* 
 * This is the procedure that is executed by the reader threads.
 *
 * Each one of these threads repeatedly reads the next read from a FASTQ file.
 * It retrieves a `struct read_set' from the "free" queue, then fills it up with
 * these reads.  When the read set it is full, it places it the on "ready" queue
 * for the combiner thread to retrieve.
 *
 * In the event of EOF, we need to signal each combiner thread by providing a
 * `struct read_set' that contains a NULL pointer as one of these reads.
 */
static void *fastq_reader_thread_proc(void *p)
{
	struct reader_writer_params *params = p;

	struct read_queue *free_q = params->free_q;
	struct read_queue *ready_q = params->ready_q;
	gzFile fp = (gzFile)params->fp;
	int phred_offset = params->phred_offset;
	struct read_set *r;
	unsigned i;
	unsigned long pair_no = 0;

	while (1) {
		r = read_queue_get(free_q);
		for (i = 0; i < READS_PER_READ_SET; i++) {
			if (!next_read(r->reads[i], fp, phred_offset))
				goto out;
			if (params->verbose && ++pair_no % 25000 == 0)
				info("Processed %lu reads", pair_no);
		}
		read_queue_put(ready_q, r);
	}
out:
	if (params->verbose && pair_no % 25000 != 0)
		info("Processed %lu reads", pair_no);
	free_read(r->reads[i]);
	r->reads[i] = NULL;
	read_queue_put(ready_q, r);
	for (int j = 1; j < params->num_combiner_threads; j++) {
		r = read_queue_get(free_q);
		free_read(r->reads[0]);
		r->reads[0] = NULL;
		read_queue_put(ready_q, r);
	}
	gzclose(fp);
	free(params);
	return NULL;
}

/* 
 * This is the procedure that is executed by the writer threads.
 *
 * Each of these threads repeatedly retrieves a read set from the "ready" queue
 * and writes the reads contained in it to a FASTQ file.
 *
 * After writing the reads contained in a read set, the read set is returned to
 * the "free" queue so that the memory is recycled.
 *
 * If a read is NULL, this is the signal for the last read; when this has
 * happened one time for each combiner thread, the output file is closed and the
 * writer thread exits.
 */
static void *fastq_writer_thread_proc(void *p)
{
	struct reader_writer_params *params = p;

	struct read_queue *free_q = params->free_q;
	struct read_queue *ready_q = params->ready_q;
	void *fp = params->fp;
	int phred_offset = params->phred_offset;

	write_read_t write_read = params->fops->write_read;

	int num_poisons = 0;
	struct read_set *r;
	while (1) {
		r = read_queue_get(ready_q);
		for (unsigned i = 0; i < READS_PER_READ_SET; i++) {
			if (!r->reads[i]) {
				if (++num_poisons == params->num_combiner_threads)
					goto out;
				else
					break;
			}
			write_read(r->reads[i], fp, phred_offset);
		}
		read_queue_put(free_q, r);
	}
out:
	read_queue_put(free_q, r);
	params->fops->close_file(fp);
	free(params);
	return NULL;
}

/* Launches a FASTQ reader thread. */
static void start_fastq_reader(void *fp, int phred_offset,
			       struct thread *thread,
			       int num_combiner_threads,
			       size_t queue_size,
			       bool verbose)
{
	struct reader_writer_params *p;
	int ret;
	
	p = new_reader_writer_params(queue_size, fp, NULL, phred_offset,
				     NULL, num_combiner_threads, verbose);
	thread->free_q = p->free_q;
	thread->ready_q = p->ready_q;
	ret = pthread_create(&thread->pthread, NULL,
			     fastq_reader_thread_proc, p);
	if (ret != 0)
		fatal_error_with_errno("Could not create FASTQ reader thread");
}

/* Launches a FASTQ writer thread. */
static void start_fastq_writer(void *fp, int phred_offset,
			       const struct file_operations *fops,
			       struct thread *thread,
			       struct read_queue *free_q, size_t queue_size,
			       int num_combiner_threads)
{
	struct reader_writer_params *p;
	int ret;
	
	p = new_reader_writer_params(queue_size, fp, fops, phred_offset,
				     free_q, num_combiner_threads, false);
	thread->free_q = p->free_q;
	thread->ready_q = p->ready_q;
	ret = pthread_create(&thread->pthread, NULL,
			     fastq_writer_thread_proc, p);
	if (ret != 0)
		fatal_error_with_errno("Could not create FASTQ writer thread");
}

/* Starts the FASTQ readers and writers needed for the FLASH program to run. */
void start_fastq_readers_and_writers(void *mates1_gzf, void *mates2_gzf,
				     void *out_combined_fp,
				     void *out_notcombined_fp_1,
				     void *out_notcombined_fp_2,
				     int phred_offset,
				     const struct file_operations *fops,
				     struct threads *threads,
				     int num_combiner_threads,
				     bool verbose)
{
	size_t queue_size = num_combiner_threads * QUEUE_SIZE_PER_THREAD;

	if (verbose)
		info("Starting FASTQ readers and writer threads");

	start_fastq_reader(mates1_gzf, phred_offset, &threads->reader_1,
			   num_combiner_threads, queue_size, verbose);

	start_fastq_reader(mates2_gzf, phred_offset, &threads->reader_2,
			   num_combiner_threads, queue_size, false);

	start_fastq_writer(out_combined_fp, phred_offset, fops,
			   &threads->writer_combined, NULL, queue_size,
			   num_combiner_threads);

	start_fastq_writer(out_notcombined_fp_1, phred_offset, fops,
			   &threads->writer_uncombined_1,
			   threads->reader_1.free_q, queue_size,
			   num_combiner_threads);

	start_fastq_writer(out_notcombined_fp_2, phred_offset, fops,
			   &threads->writer_uncombined_2,
			   threads->reader_2.free_q, queue_size,
			   num_combiner_threads);
}

/* Waits for all the reader and writer threads to exit, then free the read
 * queues.  */
void stop_fastq_readers_and_writers(const struct threads *threads)
{
	if (pthread_join(threads->reader_1.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join first reader thread");
	if (pthread_join(threads->reader_2.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join second reader thread");
	if (pthread_join(threads->writer_combined.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join extended fragments "
				       "writer thread");
	if (pthread_join(threads->writer_uncombined_1.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join uncombined fragments "
				       "writer thread #1");
	if (pthread_join(threads->writer_uncombined_2.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join uncombined fragments "
				       "writer thread #2");
	free_read_queue(threads->reader_1.free_q);
	free_read_queue(threads->reader_2.free_q);
	free_read_queue(threads->reader_1.ready_q);
	free_read_queue(threads->reader_2.ready_q);
	free_read_queue(threads->writer_uncombined_1.ready_q);
	free_read_queue(threads->writer_uncombined_2.ready_q);
	free_read_queue(threads->writer_combined.free_q);
	free_read_queue(threads->writer_combined.ready_q);
}

#else /* MULTITHREADED */

/* 
 * Reads the next mate pair from the FASTQ files @mates1_gzf and @mates2_gzf.  
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
bool next_mate_pair(struct read *read_1, struct read *read_2, gzFile mates1_gzf, 
		    gzFile mates2_gzf, int phred_offset)
{
	if (!next_read(read_1, mates1_gzf, phred_offset))
		return false;
	if (!next_read(read_2, mates2_gzf, phred_offset))
		return false;
	reverse_complement(read_2->seq, read_2->seq_len);
	reverse(read_2->qual, read_2->seq_len);
	return true;
}

#endif /* !MULTITHREADED */
