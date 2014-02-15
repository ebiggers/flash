#ifndef _FLASH_READ_QUEUE_H_
#define _FLASH_READ_QUEUE_H_

#include <stdbool.h>

struct read;
struct read_format_params;
struct input_stream;
struct output_stream;

#define READS_PER_READ_SET 32
#define QUEUE_SIZE_PER_THREAD 8

struct read_io_handle;

struct read_set {
	struct read *reads[READS_PER_READ_SET];
	size_t filled;
	enum {
		READS_UNCOMBINED,
		READS_COMBINED,
		READS_UNPAIRED,
	} type;
};

extern struct read_io_handle *
start_readers_and_writers(struct input_stream *in_1,
			  struct input_stream *in_2,
			  struct output_stream *out_combined,
			  struct output_stream *out_uncombined_1,
			  struct output_stream *out_uncombined_2,
			  const struct read_format_params *iparams,
			  const struct read_format_params *oparams,
			  unsigned num_combiner_threads,
			  bool verbose);

extern struct read_set *
get_avail_read_set(struct read_io_handle *handle);

extern bool
get_unprocessed_read_pairs(struct read_io_handle *handle,
			   struct read_set **s1_ret, struct read_set **s2_ret);

extern void
put_combined_reads(struct read_io_handle *handle, struct read_set *s);

extern void
put_uncombined_read_pairs(struct read_io_handle *handle,
			  struct read_set *s1, struct read_set *s2);

extern void
put_avail_read_pairs(struct read_io_handle *handle,
		     struct read_set *s1, struct read_set *s2);

extern void
notify_combiner_terminated(struct read_io_handle *h);

extern struct read_set *
new_empty_read_set(void);

extern void
free_read_set(struct read_set *s);

extern void
stop_readers_and_writers(struct read_io_handle *h);

#endif /* _FLASH_READ_QUEUE_H_ */
