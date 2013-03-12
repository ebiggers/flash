CC       := cc
CFLAGS   := -O2 -Wall -std=c99 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS := -pthread
LDFLAGS  := -pthread
LDLIBS   := -lz
OBJ      := combine_reads.o fastq.o flash.o util.o
EXE      := flash

$(EXE):$(OBJ)

clean:
	rm -f $(OBJ) $(EXE)

.PHONY: clean
