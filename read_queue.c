/*
 * read_queue.c: Code to set up reader/writer threads and shared queues to pass
 * reads between threads in memory.
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

#include "iostream.h"
#include "read.h"
#include "read_io.h"
#include "read_queue.h"
#include "util.h"

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

static struct read *
new_read(void)
{
	return xzalloc(sizeof(struct read));
}

void
free_read(struct read *r)
{
	if (r) {
		xfree(r->tag, r->tag_bufsz);
		xfree(r->seq, r->seq_bufsz);
		xfree(r->qual, r->qual_bufsz);
		xfree(r, sizeof(*r));
	}
}

struct read_set *
new_empty_read_set(void)
{
	return xzalloc(sizeof(struct read_set));
}

static struct read_set *
new_full_read_set(void)
{
	struct read_set *s = xmalloc(sizeof(*s));
	for (size_t i = 0; i < READS_PER_READ_SET; i++)
		s->reads[i] = new_read();
	return s;
}

void
free_read_set(struct read_set *s)
{
	if (s) {
		for (size_t i = 0; i < READS_PER_READ_SET; i++)
			free_read(s->reads[i]);
		xfree(s, sizeof(*s));
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
	size_t filled;
	struct read_set **read_sets;
	pthread_mutex_t lock;
	pthread_cond_t read_set_avail_cond;
	pthread_cond_t space_avail_cond;
};

static struct read_queue *
new_read_queue(size_t size, bool full)
{
	struct read_queue *q = xmalloc(sizeof(*q));

	q->read_sets = xmalloc(size * sizeof(q->read_sets[0]));
	q->front = 0;
	q->size = size;
	if (full) {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = new_full_read_set();
		q->filled = size;
	} else {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = NULL;
		q->filled = 0;
	}
	if (pthread_mutex_init(&q->lock, NULL))
		fatal_error_with_errno("Failed to initialize mutex");
	if (pthread_cond_init(&q->space_avail_cond, NULL))
		fatal_error_with_errno("Failed to initialize condition variable");
	if (pthread_cond_init(&q->read_set_avail_cond, NULL))
		fatal_error_with_errno("Failed to initialize condition variable");
	return q;
}

static void
free_read_queue(struct read_queue *q)
{
	if (q) {
		size_t filled = q->filled;
		size_t i = q->front;

		while (filled--) {
			free_read_set(q->read_sets[i]);
			i = (i + 1) % q->size;
		}

		pthread_mutex_destroy(&q->lock);
		pthread_cond_destroy(&q->read_set_avail_cond);
		pthread_cond_destroy(&q->space_avail_cond);

		xfree(q->read_sets, q->size * sizeof(q->read_sets[0]));
		xfree(q, sizeof(*q));
	}
}

/* Retrieve the next available read set from the queue, blocking until one is
 * available. */
static struct read_set *
read_queue_get(struct read_queue *q)
{
	struct read_set *s;

	pthread_mutex_lock(&q->lock);
	while (q->filled == 0)
		pthread_cond_wait(&q->read_set_avail_cond, &q->lock);

	s = q->read_sets[q->front];
	q->front = (q->front + 1) % q->size;
	q->filled--;

	pthread_cond_signal(&q->space_avail_cond);
	pthread_mutex_unlock(&q->lock);
	return s;
}

/* Put a read set into the queue, blocking until there is an empty space
 * available. */
static void
read_queue_put(struct read_queue *q, struct read_set *s)
{
	pthread_mutex_lock(&q->lock);
	while (q->filled == q->size)
		pthread_cond_wait(&q->space_avail_cond, &q->lock);

	q->read_sets[(q->front + q->filled) % q->size] = s;
	q->filled++;

	pthread_cond_signal(&q->read_set_avail_cond);
	pthread_mutex_unlock(&q->lock);
}

struct reader_params {
	struct input_stream *in;
	const struct read_format_params *iparams;
	unsigned num_combiner_threads;
	bool verbose;
	struct read_queue *avail_read_1_q;
	struct read_queue *avail_read_2_q;
	struct read_queue *unprocessed_read_1_q;
	struct read_queue *unprocessed_read_2_q;
};

struct writer_params {
	struct output_stream *out;
	const struct read_format_params *oparams;
	unsigned poison_count;
	struct read_queue *to_write_queue_1;
	struct read_queue *to_write_queue_2;
	struct read_queue *avail_queue_1;
	struct read_queue *avail_queue_2;
};

static void *
reader_proc(void *_params)
{
	struct reader_params *params = _params;
	struct read_set *s1, *s2;
	uint64_t pair_no = 0;
	uint64_t line_no = 1;
	size_t i;

	for (;;) {
		s1 = read_queue_get(params->avail_read_1_q);
		if (params->avail_read_2_q)
			s2 = read_queue_get(params->avail_read_2_q);
		else
			s2 = NULL;

		for (i = 0; i < READS_PER_READ_SET; i++) {
			if (s2) {
				if (!load_read_pair(params->in,
						    params->iparams,
						    s1->reads[i],
						    s2->reads[i],
						    &line_no))
					goto eof_reached;
			} else {
				if (!load_read(params->in, params->iparams,
					       s1->reads[i],
					       &line_no))
					goto eof_reached;
			}

			if (params->verbose && ++pair_no % 25000 == 0)
				info("Processed %"PRIu64" read pairs", pair_no);
		}
		read_queue_put(params->unprocessed_read_1_q, s1);
		if (s2)
			read_queue_put(params->unprocessed_read_2_q, s2);
	}

eof_reached:
	if (params->verbose && pair_no % 25000 != 0)
		info("Processed %"PRIu64" read pairs", pair_no);

	free_read(s1->reads[i]);
	s1->reads[i] = NULL;
	read_queue_put(params->unprocessed_read_1_q, s1);
	if (s2) {
		free_read(s2->reads[i]);
		s2->reads[i] = NULL;
		read_queue_put(params->unprocessed_read_2_q, s2);
	}

	for (unsigned k = 1; k < params->num_combiner_threads; k++) {
		struct read_set *s;

		s = read_queue_get(params->avail_read_1_q);
		free_read(s->reads[0]);
		s->reads[0] = NULL;
		read_queue_put(params->unprocessed_read_1_q, s);

		if (params->avail_read_2_q) {
			s = read_queue_get(params->avail_read_2_q);
			free_read(s->reads[0]);
			s->reads[0] = NULL;
			read_queue_put(params->unprocessed_read_2_q, s);
		}
	}
	free_input_stream(params->in);
	xfree(params, sizeof(*params));
	return NULL;
}

static void *
writer_proc(void *_params)
{
	struct writer_params *params = _params;
	unsigned poisons_remaining = params->poison_count;
	struct read_set *s1, *s2;

	while (poisons_remaining) {
		s1 = read_queue_get(params->to_write_queue_1);
		if (s1->paired && params->to_write_queue_2)
			s2 = read_queue_get(params->to_write_queue_2);
		else
			s2 = NULL;

		for (size_t i = 0; i < READS_PER_READ_SET; i++) {
			if (!s1->reads[i]) {
				assert(!s2 || !s2->reads[i]);
				--poisons_remaining;
				break;
			}
			if (s2)
				write_read_pair(params->out, params->oparams,
						s1->reads[i], s2->reads[i]);
			else
				write_read(params->out, params->oparams,
					   s1->reads[i]);
		}
		read_queue_put(params->avail_queue_1, s1);
		if (s2)
			read_queue_put(params->avail_queue_2, s2);
	}
	free_output_stream(params->out);
	xfree(params, sizeof(*params));
	return NULL;
}

static pthread_t
start_reader2(struct input_stream *in,
	      const struct read_format_params *iparams,
	      unsigned num_combiner_threads,
	      bool verbose,
	      struct read_queue *avail_read_1_q,
	      struct read_queue *avail_read_2_q,
	      struct read_queue *unprocessed_read_1_q,
	      struct read_queue *unprocessed_read_2_q)
{
	struct reader_params *params = xmalloc(sizeof(*params));

	params->in = in;
	params->iparams = iparams;
	params->num_combiner_threads = num_combiner_threads;
	params->verbose = verbose;
	params->avail_read_1_q = avail_read_1_q;
	params->avail_read_2_q = avail_read_2_q;
	params->unprocessed_read_1_q = unprocessed_read_1_q;
	params->unprocessed_read_2_q = unprocessed_read_2_q;

	return create_thread(reader_proc, params);
}

static pthread_t
start_reader1(struct input_stream *in,
	      const struct read_format_params *iparams,
	      unsigned num_combiner_threads,
	      bool verbose,
	      struct read_queue *avail_read_q,
	      struct read_queue *unprocessed_read_q)
{
	return start_reader2(in, iparams, num_combiner_threads, verbose,
			     avail_read_q, NULL, unprocessed_read_q, NULL);
}

static pthread_t
start_writer2(struct output_stream *out,
	      const struct read_format_params *oparams,
	      unsigned poison_count,
	      struct read_queue *to_write_queue_1,
	      struct read_queue *to_write_queue_2,
	      struct read_queue *avail_queue_1,
	      struct read_queue *avail_queue_2)
{
	struct writer_params *params = xmalloc(sizeof(*params));

	params->out = out;
	params->oparams = oparams;
	params->poison_count = poison_count;
	params->to_write_queue_1 = to_write_queue_1;
	params->to_write_queue_2 = to_write_queue_2;
	params->avail_queue_1 = avail_queue_1;
	params->avail_queue_2 = avail_queue_2;

	return create_thread(writer_proc, params);
}

static pthread_t
start_writer1(struct output_stream *out,
	      const struct read_format_params *oparams,
	      unsigned poison_count,
	      struct read_queue *to_write_queue,
	      struct read_queue *avail_queue)
{
	return start_writer2(out, oparams, poison_count,
			     to_write_queue, NULL, avail_queue, NULL);
}


struct read_io_handle {

	pthread_t reader_1;
	pthread_t reader_2;
	pthread_t writer_1;
	pthread_t writer_2;
	pthread_t writer_3;
	bool reader_1_started;
	bool reader_2_started;
	bool writer_1_started;
	bool writer_2_started;
	bool writer_3_started;

	struct read_queue *avail_read_q;
	struct read_queue *unprocessed_read_1_q;
	struct read_queue *unprocessed_read_2_q;

	struct read_queue *combined_read_q;
	struct read_queue *uncombined_read_1_q;
	struct read_queue *uncombined_read_2_q;

	pthread_mutex_t get_unprocessed_pair_mutex;
	pthread_mutex_t put_uncombined_pair_mutex;
};

/* Retrieve a set of unprocessed read pairs from the I/O layer.  */
void
get_unprocessed_read_pairs(struct read_io_handle *h, struct read_set **s1_p,
			   struct read_set **s2_p)
{
	/* get_unprocessed_pair_mutex ensures the reads are paired up correctly.
	 */

	pthread_mutex_lock(&h->get_unprocessed_pair_mutex);

	*s1_p = read_queue_get(h->unprocessed_read_1_q);
	*s2_p = read_queue_get(h->unprocessed_read_2_q);

	pthread_mutex_unlock(&h->get_unprocessed_pair_mutex);
}

/* Submits a set of combined reads to the I/O layer to be written.  */
void
put_combined_reads(struct read_io_handle *h, struct read_set *s)
{
	s->paired = false;

	read_queue_put(h->combined_read_q, s);
}

/* Submits a set of uncombined read pairs to the I/O layer to be written.  */
void
put_uncombined_read_pairs(struct read_io_handle *h,
			  struct read_set *s1, struct read_set *s2)
{
	s1->paired = true;
	s2->paired = true;

	/* put_unprocessed_pair_mutex ensures the reads are paired up correctly.
	 */

	pthread_mutex_lock(&h->put_uncombined_pair_mutex);

	read_queue_put(h->uncombined_read_1_q, s1);
	read_queue_put(h->uncombined_read_2_q, s2);

	pthread_mutex_unlock(&h->put_uncombined_pair_mutex);
}


/* Retrieve a read set that is ready to be reused.  */
struct read_set *
get_empty_read_set(struct read_io_handle *h)
{
	return read_queue_get(h->avail_read_q);
}

/* Return a set of read pairs to the pool for reuse.  */
void
put_empty_read_pairs(struct read_io_handle *h,
		     struct read_set *s1, struct read_set *s2)
{
	read_queue_put(h->avail_read_q, s1);
	read_queue_put(h->avail_read_q, s2);
}

/* Starts the FLASH I/O layer, which is responsible for input/output of reads.
 *
 * If @in_2 is not NULL, then @in_1 and @in_2 are the input files for read 1 and
 * read 2 of the pairs, respectively.  Otherwise @in_1 contains both read 1 and
 * read 2 of the pairs interleaved.
 *
 * Either 1, 2, or 3 output files may be specified --- see below for more
 * details.  */
struct read_io_handle *
start_readers_and_writers(struct input_stream *in_1,
			  struct input_stream *in_2,
			  struct output_stream *out_combined,
			  struct output_stream *out_uncombined_1,
			  struct output_stream *out_uncombined_2,
			  const struct read_format_params *iparams,
			  const struct read_format_params *oparams,
			  unsigned num_combiner_threads,
			  bool verbose)
{
	assert(in_1 != NULL);
	assert(out_combined != NULL &&
	       (out_uncombined_1 != NULL || out_uncombined_2 == NULL));
	assert(iparams != NULL);
	assert(oparams != NULL);
	assert(num_combiner_threads > 0);

	if (verbose)
		info("Starting reader and writer threads");

	struct read_io_handle *h = xzalloc(sizeof(*h));

	size_t queue_size = num_combiner_threads * QUEUE_SIZE_PER_THREAD;

	h->avail_read_q = new_read_queue(queue_size * 3, true);
	h->unprocessed_read_1_q = new_read_queue(queue_size, false);
	h->unprocessed_read_2_q = new_read_queue(queue_size, false);
	h->combined_read_q = new_read_queue(queue_size, false);

	if (pthread_mutex_init(&h->get_unprocessed_pair_mutex, NULL))
		fatal_error_with_errno("Failed to initialize mutex");
	if (pthread_mutex_init(&h->put_uncombined_pair_mutex, NULL))
		fatal_error_with_errno("Failed to initialize mutex");

	if (in_2) {
		/* Two input files:  read 1 in each pair comes from the first
		 * file, and read 2 in each pair comes from the second file.
		 *
		 * Only set @verbose for one.  */
		h->reader_1 = start_reader1(in_1,
					    iparams,
					    num_combiner_threads,
					    verbose,
					    h->avail_read_q,
					    h->unprocessed_read_1_q);
		h->reader_1_started = true;

		h->reader_2 = start_reader1(in_2,
					    iparams,
					    num_combiner_threads,
					    false,
					    h->avail_read_q,
					    h->unprocessed_read_2_q);
		h->reader_2_started = true;
	} else {
		/* One input file:  both reads in each pair come from the same
		 * file.  */
		h->reader_1 = start_reader2(in_1,
					    iparams,
					    num_combiner_threads,
					    verbose,
					    h->avail_read_q,
					    h->avail_read_q,
					    h->unprocessed_read_1_q,
					    h->unprocessed_read_2_q);
		h->reader_1_started = true;
	}

	if (out_uncombined_2) {
		/* All 3 output files specified: one for combined reads, one for
		 * read 1 of uncombined pairs, and one for read 2 of uncombined
		 * pairs.  */

		h->uncombined_read_1_q = new_read_queue(queue_size, false);
		h->uncombined_read_2_q = new_read_queue(queue_size, false);

		h->writer_1 = start_writer1(out_combined, oparams,
					    num_combiner_threads,
					    h->combined_read_q,
					    h->avail_read_q);
		h->writer_1_started = true;

		h->writer_2 = start_writer1(out_uncombined_1, oparams,
					    num_combiner_threads,
					    h->uncombined_read_1_q,
					    h->avail_read_q);
		h->writer_2_started = true;

		h->writer_3 = start_writer1(out_uncombined_2, oparams,
					    num_combiner_threads,
					    h->uncombined_read_2_q,
					    h->avail_read_q);
		h->writer_3_started = true;
	} else if (out_uncombined_1) {
		/* 2 output files specified: one for combined reads and one for
		 * uncombined pairs.  */

		h->uncombined_read_1_q = new_read_queue(queue_size, false);
		h->uncombined_read_2_q = new_read_queue(queue_size, false);

		h->writer_1 = start_writer1(out_combined, oparams,
					    num_combiner_threads,
					    h->combined_read_q,
					    h->avail_read_q);
		h->writer_1_started = true;

		h->writer_2 = start_writer2(out_uncombined_1, oparams,
					    num_combiner_threads,
					    h->uncombined_read_1_q,
					    h->uncombined_read_2_q,
					    h->avail_read_q,
					    h->avail_read_q);
		h->writer_2_started = true;
	} else {
		/* 1 output file specified: combined reads, plus optionally
		 * uncombined pairs if supported by the format.  */

		if (read_format_supports_mixed_reads(oparams)) {
			h->uncombined_read_1_q = h->combined_read_q;
			h->uncombined_read_2_q = new_read_queue(queue_size, false);

			h->writer_1 = start_writer2(out_combined, oparams,
						    num_combiner_threads * 2,
						    h->combined_read_q,
						    h->uncombined_read_2_q,
						    h->avail_read_q,
						    h->avail_read_q);
			h->writer_1_started = true;
		} else {
			/* Can only output combined reads.
			 * Reroute uncombined reads back to the queue of
			 * available (for reuse) reads.  */
			h->uncombined_read_1_q = h->avail_read_q;
			h->uncombined_read_2_q = h->avail_read_q;

			h->writer_1 = start_writer1(out_combined, oparams,
						    num_combiner_threads,
						    h->combined_read_q,
						    h->avail_read_q);
			h->writer_1_started = true;
		}
	}

	return h;
}

/* Terminates the FLASH I/O layer, which is responsible for input/output of
 * reads.
 */
void
stop_readers_and_writers(struct read_io_handle *h)
{
	if (h->reader_1_started)
		join_thread(h->reader_1);
	if (h->reader_2_started)
		join_thread(h->reader_2);
	if (h->writer_1_started)
		join_thread(h->writer_1);
	if (h->writer_2_started)
		join_thread(h->writer_2);
	if (h->writer_3_started)
		join_thread(h->writer_3);

	free_read_queue(h->avail_read_q);
	free_read_queue(h->unprocessed_read_1_q);
	free_read_queue(h->unprocessed_read_2_q);
	free_read_queue(h->combined_read_q);

	if (h->uncombined_read_1_q != h->avail_read_q &&
	    h->uncombined_read_1_q != h->combined_read_q)
		free_read_queue(h->uncombined_read_1_q);

	if (h->uncombined_read_2_q != h->avail_read_q)
		free_read_queue(h->uncombined_read_2_q);

	xfree(h, sizeof(*h));
}
