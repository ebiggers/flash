MULTITHREADED := yes

CC      := cc
CFLAGS  := -O2 -Wall -std=c99
LDFLAGS := -lz
OBJ     := combine_reads.o fastq.o flash.o util.o
EXE     := flash

ifeq ($(MULTITHREADED),yes)
	CFLAGS += -DMULTITHREADED
	LDFLAGS += -lpthread
endif

$(EXE):$(OBJ)

clean:
	rm -f $(OBJ) $(EXE)

.PHONY: clean
