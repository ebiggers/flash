CC      := cc
CFLAGS  := -O2 -Wall -std=c99
LDFLAGS := -lz
OBJ     := combine_reads.o fastq.o flash.o util.o
EXE     := flash

$(EXE):$(OBJ)

clean:
	rm -f $(OBJ) $(EXE)

.PHONY: clean
