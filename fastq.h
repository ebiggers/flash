#ifndef _FASTQ_H
#define _FASTQ_H

#include <zlib.h>
#include <stdbool.h>
#include "util.h"
#include <stdlib.h>

struct file_operations;

/* Read from a FASTQ file. */
struct read {
	char *tag;
	char *seq;
	char *qual;
	size_t seq_bufsz;
	size_t tag_bufsz;
	size_t qual_bufsz;
	int seq_len;
	int tag_len;
};

#include <pthread.h>
#define READS_PER_READ_SET 32
#define QUEUE_SIZE_PER_THREAD 8

struct read_queue;

struct thread {
	pthread_t pthread;
	struct read_queue *free_q;
	struct read_queue *ready_q;
};
struct threads {
	struct thread reader_1;
	struct thread reader_2;
	struct thread writer_combined;
	struct thread writer_uncombined_1;
	struct thread writer_uncombined_2;
};

struct read_set {
	struct read *reads[READS_PER_READ_SET];
};


extern void
start_fastq_readers_and_writers(gzFile mates1_gzf,
				gzFile mates2_gzf,
				void *out_combined_fp,
				void *out_notcombined_fp_1,
				void *out_notcombined_fp_2,
				int phred_offset,
				const struct file_operations *fops,
				struct threads *threads,
				int num_worker_threads,
				bool verbose);

extern void stop_fastq_readers_and_writers(const struct threads *threads);

extern struct read_set *read_queue_get(struct read_queue *q);
extern void read_queue_put(struct read_queue *q, struct read_set *r);

static inline struct read_set *new_empty_read_set()
{
	return xmalloc(sizeof(struct read_set));
}

extern void free_read_set(struct read_set *p);

static inline void init_read(struct read *read)
{
	memset(read, 0, sizeof(*read));
}

extern void destroy_read(struct read *r);

extern void free_read(struct read *r);

extern void write_read_uncompressed(struct read *read, void *fp,
				    int phred_offset);
extern void write_read_compressed(struct read *read, void *fp,
				  int phred_offset);

#endif /* _FASTQ_H */
