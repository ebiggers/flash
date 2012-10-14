MULTITHREADED := yes

CC       := cc
CFLAGS   := -O2 -Wall -std=c99
LDLIBS   := -lz
OBJ      := combine_reads.o fastq.o flash.o util.o
EXE      := flash

ifeq ($(MULTITHREADED),yes)
	CFLAGS   += -DMULTITHREADED
	CPPFLAGS += -pthread
	LDFLAGS  += -pthread
endif

$(EXE):$(OBJ)

clean:
	rm -f $(OBJ) $(EXE)

.PHONY: clean
