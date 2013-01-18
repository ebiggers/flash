#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "fastq.h"
#include "util.h"

/* Similar to getline(), but read from a gzFile, and abort on read error. */
static ssize_t gzip_getline(gzFile gz_fp, char **lineptr, size_t *n)
{
	size_t offset = 0;
	if (!*lineptr) {
		*n = 128;
		*lineptr = (char*)xmalloc(*n);
	}
	while (1) {
		char *line = *lineptr;
		assert(*n - offset > 0);
		if (gzgets(gz_fp, line + offset, *n - offset)) {
			/* gzgets() reads until '\n' or until *n - offset - 1
			 * characters have been read.  In both cases the buffer
			 * is terminated with a '\0' character. */
			char *nl_p = memchr(line + offset, '\n', *n - offset - 1);
			if (nl_p) {
				/* Read the whole line.  Return its length,
				 * including the '\n'. */
				assert(*(nl_p + 1) == '\0');
				return nl_p - line + 1;
			} else {
				/* Didn't read the whole line.  Increase the
				 * buffer size. */
				assert(line[*n - 1] == '\0');
				offset = *n - 1;
				*n *= 2;
				*lineptr = xrealloc(line, *n);
				continue;
			}
		} else {
			/* gzgets() returned NULL, so EOF or error occurred. */
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



struct file_operations gzip_fops = {
	.name       = "gzip",
	.suffix     = ".gz",
	.write_read = write_read_compressed,
	.open_file  = xgzopen,
	.close_file = xgzclose,
};

struct file_operations normal_fops = {
	.name       = "text",
	.suffix     = "",
	.write_read = write_read_uncompressed,
	.open_file  = xfopen,
	.close_file = xfclose,
};

struct file_operations pipe_fops = {
	.name       = NULL,
	.suffix     = NULL,
	.write_read = write_read_uncompressed,
	.open_file  = xpopen,
	.close_file = xpclose,
};

char *compress_prog = NULL;
char *compress_prog_args = "";

void destroy_read(struct read *r)
{
	free(r->tag);
	free(r->seq);
	free(r->qual);
}

void free_read(struct read *r)
{
	if (r) {
		destroy_read(r);
		free(r);
	}
}

static struct read *new_read()
{
	struct read *r = xmalloc(sizeof (struct read));
	init_read(r);
	return r;
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

/*
 * Producer-consumer queue; it holds pointers to `struct read_sets', which can
 * be added or removed from the queue in a thread-safe manner using
 * read_queue_put() and read_queue_get(), respectively.
 */
struct read_queue {
	size_t size;
	size_t front;
	size_t back;
	size_t filled_slots;
	struct read_set **read_sets;
	pthread_mutex_t lock;
	pthread_cond_t read_set_avail_cond;
	pthread_cond_t space_avail_cond;
};

static struct read_queue *new_read_queue(size_t size, bool full)
{
	struct read_queue *q = xmalloc(sizeof(struct read_queue));
	q->read_sets = xmalloc(size * sizeof(struct read_set*));
	if (full) {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = new_read_set();
		q->filled_slots = size;
	} else {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = NULL;
		q->filled_slots = 0;
	}
	q->front = 0;
	q->back = size - 1;
	q->size = size;
	if (pthread_mutex_init(&q->lock, NULL))
		fatal_error_with_errno("Failed to initialize mutex");
	if (pthread_cond_init(&q->read_set_avail_cond, NULL))
		fatal_error_with_errno("Failed to initialize condition variable");
	if (pthread_cond_init(&q->space_avail_cond, NULL))
		fatal_error_with_errno("Failed to initialize condition variable");
	return q;
}

static void free_read_queue(struct read_queue *q)
{
	if (q) {
		size_t filled_slots = q->filled_slots;
		size_t i = q->front;

		while (filled_slots--) {
			free_read_set(q->read_sets[i]);
			i = (i + 1) % q->size;
		}

		pthread_mutex_destroy(&q->lock);
		pthread_cond_destroy(&q->read_set_avail_cond);
		pthread_cond_destroy(&q->space_avail_cond);
		free(q->read_sets);
	}
}

/* Retrieve the next available read set from the queue, blocking until one is
 * available. */
struct read_set *read_queue_get(struct read_queue *q)
{
	struct read_set *r;

	pthread_mutex_lock(&q->lock);
	while (q->filled_slots == 0)
		pthread_cond_wait(&q->read_set_avail_cond, &q->lock);

	r = q->read_sets[q->front];
	q->front = (q->front + 1) % q->size;
	q->filled_slots--;

	pthread_cond_broadcast(&q->space_avail_cond);
	pthread_mutex_unlock(&q->lock);
	return r;
}

/* Put a read set into the queue, blocking until there is an empty space
 * available. */
void read_queue_put(struct read_queue *q, struct read_set *r)
{
	pthread_mutex_lock(&q->lock);
	while (q->filled_slots == q->size)
		pthread_cond_wait(&q->space_avail_cond, &q->lock);

	q->back = (q->back + 1) % q->size;
	q->read_sets[q->back] = r;
	q->filled_slots++;

	pthread_cond_broadcast(&q->read_set_avail_cond);
	pthread_mutex_unlock(&q->lock);
}

struct common_reader_params {
	int phred_offset;
	int num_combiner_threads;
	bool verbose;
	gzFile fp;
};

struct reader_params {
	struct common_reader_params common;
	struct read_queue *free_q;
	struct read_queue *ready_q;
};

struct interleaved_reader_params {
	struct common_reader_params common;
	struct read_queue *free_q_1;
	struct read_queue *ready_q_1;
	struct read_queue *free_q_2;
	struct read_queue *ready_q_2;
};

struct common_writer_params {
	int phred_offset;
	int num_combiner_threads;
	void *fp;
	const struct file_operations *fops;
};

struct writer_params {
	struct common_writer_params common;
	struct read_queue *free_q;
	struct read_queue *ready_q;
};

struct interleaved_writer_params {
	struct common_writer_params common;
	struct read_queue *free_q_1;
	struct read_queue *ready_q_1;
	struct read_queue *free_q_2;
	struct read_queue *ready_q_2;
};


static void init_common_reader_params(struct common_reader_params *params,
				      gzFile fp, int phred_offset,
				      int num_combiner_threads, bool verbose)
{
	params->fp = fp;
	params->phred_offset = phred_offset;
	params->num_combiner_threads = num_combiner_threads;
	params->verbose = verbose;
}

static void init_common_writer_params(struct common_writer_params *params,
				      void *fp, const struct file_operations *fops,
				      int phred_offset, int num_combiner_threads)
{
	params->fp = fp;
	params->fops = fops;
	params->phred_offset = phred_offset;
	params->num_combiner_threads = num_combiner_threads;
}

static void *fastq_reader_thread_proc(void *_params)
{
	struct reader_params *params = _params;
	void *fp = params->common.fp;
	int phred_offset = params->common.phred_offset;
	bool verbose = params->common.verbose;
	unsigned i;
	unsigned long pair_no = 0;
	struct read_set *r;

	while (1) {
		r = read_queue_get(params->free_q);
		for (i = 0; i < READS_PER_READ_SET; i++) {
			if (!next_read(r->reads[i], fp, phred_offset))
				goto out;
			if (verbose && ++pair_no % 25000 == 0)
				info("Processed %lu read pairs", pair_no);
		}
		read_queue_put(params->ready_q, r);
	}
out:
	if (verbose && pair_no % 25000 != 0)
		info("Processed %lu read pairs", pair_no);

	free_read(r->reads[i]);
	r->reads[i] = NULL;
	read_queue_put(params->ready_q, r);
	for (int j = 1; j < params->common.num_combiner_threads; j++) {
		r = read_queue_get(params->free_q);
		free_read(r->reads[0]);
		r->reads[0] = NULL;
		read_queue_put(params->ready_q, r);
	}
	gzclose(fp);
	free(params);
	return NULL;
}

static void *fastq_interleaved_reader_thread_proc(void *_params)
{
	struct interleaved_reader_params *params = _params;
	void *fp = params->common.fp;
	int phred_offset = params->common.phred_offset;
	bool verbose = params->common.verbose;
	unsigned i;
	unsigned long pair_no = 0;
	struct read_set *r1, *r2;

	while (1) {
		r1 = read_queue_get(params->free_q_1);
		r2 = read_queue_get(params->free_q_2);
		for (i = 0; i < READS_PER_READ_SET; i++) {
			if (!next_read(r1->reads[i], fp, phred_offset))
				goto out;
			if (!next_read(r2->reads[i], fp, phred_offset)) {
				fatal_error("Found odd number of reads "
					    "in interleaved reads file");
			}
			if (verbose && ++pair_no % 25000 == 0)
				info("Processed %lu read pairs", pair_no);
		}
		read_queue_put(params->ready_q_1, r1);
		read_queue_put(params->ready_q_2, r2);
	}

out:
	if (verbose && pair_no % 25000 != 0)
		info("Processed %lu read pairs", pair_no);

	free_read(r1->reads[i]);
	free_read(r2->reads[i]);
	r1->reads[i] = NULL;
	r2->reads[i] = NULL;
	read_queue_put(params->ready_q_1, r1);
	read_queue_put(params->ready_q_2, r2);
	for (int j = 1; j < params->common.num_combiner_threads; j++) {
		r1 = read_queue_get(params->free_q_1);
		free_read(r1->reads[0]);
		r1->reads[0] = NULL;
		read_queue_put(params->ready_q_1, r1);

		r2 = read_queue_get(params->free_q_2);
		free_read(r2->reads[0]);
		r2->reads[0] = NULL;
		read_queue_put(params->ready_q_2, r2);
	}
	gzclose(fp);
	free(params);
	return NULL;
}

static void *fastq_writer_thread_proc(void *_params)
{
	struct writer_params *params = _params;
	write_read_t write_read = params->common.fops->write_read;
	int phred_offset = params->common.phred_offset;
	void *fp = params->common.fp;
	int num_poisons = 0;
	struct read_set *r;

	while (1) {
		r = read_queue_get(params->ready_q);
		for (unsigned i = 0; i < READS_PER_READ_SET; i++) {
			if (!r->reads[i]) {
				if (++num_poisons == params->common.num_combiner_threads)
					goto out;
				else
					break;
			}
			write_read(r->reads[i], fp, phred_offset);
		}
		read_queue_put(params->free_q, r);
	}
out:
	read_queue_put(params->free_q, r);
	params->common.fops->close_file(fp);
	free(params);
	return NULL;
}

static void *fastq_interleaved_writer_thread_proc(void *_params)
{
	struct interleaved_writer_params *params = _params;
	write_read_t write_read = params->common.fops->write_read;
	int phred_offset = params->common.phred_offset;
	void *fp = params->common.fp;
	int num_poisons = 0;
	struct read_set *r1, *r2;

	while (1) {
		r1 = read_queue_get(params->ready_q_1);
		r2 = read_queue_get(params->ready_q_2);
		for (unsigned i = 0; i < READS_PER_READ_SET; i++) {
			if (r1->reads[i] == NULL) {
				assert(r2->reads[i] == NULL);
				if (++num_poisons == params->common.num_combiner_threads)
					goto out;
				else
					break;
			}
			assert(r2->reads[i] != NULL);
			write_read(r1->reads[i], fp, phred_offset);
			write_read(r2->reads[i], fp, phred_offset);
		}
		read_queue_put(params->free_q_1, r1);
		read_queue_put(params->free_q_2, r2);
	}
out:
	read_queue_put(params->free_q_1, r1);
	read_queue_put(params->free_q_2, r2);
	params->common.fops->close_file(fp);
	free(params);
	return NULL;
}

static void start_noninterleaved_fastq_reader(gzFile fp,
					      int phred_offset,
					      int num_combiner_threads,
					      bool verbose,
					      struct thread *thread,
					      struct read_queue *free_q,
					      struct read_queue *ready_q)
{
	struct reader_params *params = xmalloc(sizeof *params);
	init_common_reader_params(&params->common, fp, phred_offset,
				  num_combiner_threads, verbose);
	params->free_q = free_q;
	params->ready_q = ready_q;
	if (pthread_create(&thread->pthread, NULL,
			   fastq_reader_thread_proc, params))
	{
		fatal_error_with_errno("Could not create FASTQ reader thread");
	}
}

static void start_interleaved_fastq_reader(gzFile fp,
					   int phred_offset,
					   int num_combiner_threads,
					   bool verbose,
					   struct thread *thread,
					   struct read_queue *free_q_1,
					   struct read_queue *ready_q_1,
					   struct read_queue *free_q_2,
					   struct read_queue *ready_q_2)
{
	struct interleaved_reader_params *params = xmalloc(sizeof *params);
	init_common_reader_params(&params->common, fp, phred_offset,
				  num_combiner_threads, verbose);
	params->free_q_1 = free_q_1;
	params->ready_q_1 = ready_q_1;
	params->free_q_2 = free_q_2;
	params->ready_q_2 = ready_q_2;
	if (pthread_create(&thread->pthread, NULL,
			   fastq_interleaved_reader_thread_proc, params))
	{
		fatal_error_with_errno("Could not create FASTQ reader thread");
	}
}


static void start_noninterleaved_fastq_writer(void *fp,
					      const struct file_operations *fops,
					      int phred_offset,
					      int num_combiner_threads,
					      struct thread *thread,
					      struct read_queue *free_q,
					      struct read_queue *ready_q)
{
	struct writer_params *params = xmalloc(sizeof *params);
	init_common_writer_params(&params->common,
				  fp, fops, phred_offset, num_combiner_threads);
	params->free_q = free_q;
	params->ready_q = ready_q;
	if (pthread_create(&thread->pthread, NULL,
			   fastq_writer_thread_proc, params))
	{
		fatal_error_with_errno("Could not create FASTQ writer thread");
	}
}

static void start_interleaved_fastq_writer(void *fp,
					   const struct file_operations *fops,
					   int phred_offset,
					   int num_combiner_threads,
					   struct thread *thread,
					   struct read_queue *free_q_1,
					   struct read_queue *ready_q_1,
					   struct read_queue *free_q_2,
					   struct read_queue *ready_q_2)
{
	struct interleaved_writer_params *params = xmalloc(sizeof *params);
	init_common_writer_params(&params->common,
				  fp, fops, phred_offset, num_combiner_threads);
	params->free_q_1 = free_q_1;
	params->ready_q_1 = ready_q_1;
	params->free_q_2 = free_q_2;
	params->ready_q_2 = ready_q_2;

	if (pthread_create(&thread->pthread, NULL,
			   fastq_interleaved_writer_thread_proc, params))
	{
		fatal_error_with_errno("Could not create FASTQ writer thread");
	}
}


/* Starts the FASTQ readers and writers needed for the FLASH program to run. */
void start_fastq_readers_and_writers(gzFile mates1_gzf,
				     gzFile mates2_gzf,
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

	threads->reader_1.free_q             = new_read_queue(queue_size, true);
	threads->reader_1.ready_q            = new_read_queue(queue_size, false);
	threads->reader_2.free_q             = new_read_queue(queue_size, true);
	threads->reader_2.ready_q            = new_read_queue(queue_size, false);
	threads->writer_combined.free_q      = new_read_queue(queue_size, true);
	threads->writer_combined.ready_q     = new_read_queue(queue_size, false);
	threads->writer_uncombined_1.free_q  = threads->reader_1.free_q;
	threads->writer_uncombined_1.ready_q = new_read_queue(queue_size, false);
	threads->writer_uncombined_2.free_q  = threads->reader_2.free_q;
	threads->writer_uncombined_2.ready_q = new_read_queue(queue_size, false);

	if (mates2_gzf) {
		start_noninterleaved_fastq_reader(mates1_gzf, phred_offset,
						  num_combiner_threads, verbose,
						  &threads->reader_1,
						  threads->reader_1.free_q,
						  threads->reader_1.ready_q);

		start_noninterleaved_fastq_reader(mates2_gzf, phred_offset,
						  num_combiner_threads, false,
						  &threads->reader_2,
						  threads->reader_2.free_q,
						  threads->reader_2.ready_q);
	} else {
		start_interleaved_fastq_reader(mates1_gzf, phred_offset,
					       num_combiner_threads, verbose,
					       &threads->reader_1,
					       threads->reader_1.free_q,
					       threads->reader_1.ready_q,
					       threads->reader_2.free_q,
					       threads->reader_2.ready_q);
		threads->reader_2.pthread = (pthread_t)-1;
	}

	/* combined reads are not paired and therefore are non-interleaved */
	start_noninterleaved_fastq_writer(out_combined_fp, fops, phred_offset,
					  num_combiner_threads,
					  &threads->writer_combined,
					  threads->writer_combined.free_q,
					  threads->writer_combined.ready_q);

	if (out_notcombined_fp_2) {
		/* non-interleaved output for uncombined reads */
		start_noninterleaved_fastq_writer(out_notcombined_fp_1, fops,
						  phred_offset, num_combiner_threads,
						  &threads->writer_uncombined_1,
						  threads->writer_uncombined_1.free_q,
						  threads->writer_uncombined_1.ready_q);

		start_noninterleaved_fastq_writer(out_notcombined_fp_2, fops,
						  phred_offset, num_combiner_threads,
						  &threads->writer_uncombined_2,
						  threads->writer_uncombined_2.free_q,
						  threads->writer_uncombined_2.ready_q);
	} else {
		/* interleaved output for uncombined reads */
		start_interleaved_fastq_writer(out_notcombined_fp_1, fops,
					       phred_offset, num_combiner_threads,
					       &threads->writer_uncombined_1,
					       threads->writer_uncombined_1.free_q,
					       threads->writer_uncombined_1.ready_q,
					       threads->writer_uncombined_2.free_q,
					       threads->writer_uncombined_2.ready_q);
		threads->writer_uncombined_2.pthread = (pthread_t)-1;
	}
}

/* Waits for all the reader and writer threads to exit, then free the read
 * queues.  */
void stop_fastq_readers_and_writers(const struct threads *threads)
{
	if (pthread_join(threads->reader_1.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join first reader thread");

	if (threads->reader_2.pthread != (pthread_t)-1) {
		if (pthread_join(threads->reader_2.pthread, NULL) != 0)
			fatal_error_with_errno("Failed to join second reader thread");
	} else {
		/* Interleaved reads--- only had one reader thread */
	}

	if (pthread_join(threads->writer_combined.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join extended fragments "
				       "writer thread");
	if (pthread_join(threads->writer_uncombined_1.pthread, NULL) != 0)
		fatal_error_with_errno("Failed to join uncombined fragments "
				       "writer thread #1");

	if (threads->writer_uncombined_2.pthread != (pthread_t)-1) {
		if (pthread_join(threads->writer_uncombined_2.pthread, NULL) != 0)
			fatal_error_with_errno("Failed to join uncombined fragments "
					       "writer thread #2");
	} else {
		/* Interleaved reads--- only had one uncombined read writer
		 * thread */
	}
	free_read_queue(threads->reader_1.free_q);
	free_read_queue(threads->reader_2.free_q);
	free_read_queue(threads->reader_1.ready_q);
	free_read_queue(threads->reader_2.ready_q);
	free_read_queue(threads->writer_uncombined_1.ready_q);
	free_read_queue(threads->writer_uncombined_2.ready_q);
	free_read_queue(threads->writer_combined.free_q);
	free_read_queue(threads->writer_combined.ready_q);
}
