CC ?= gcc
MPICC ?= mpicc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
MPI_CFLAGS ?= $(CFLAGS)
INCLUDES := -Isrc

SERIAL_SRCS := src/tsp.c src/ga_common.c src/ga_serial.c src/cli.c src/serial_main.c
PARALLEL_SRCS := src/tsp.c src/ga_common.c src/parallel_ga.c src/cli.c src/main.c

SERIAL_OBJS := $(patsubst src/%.c, build/serial/%.o, $(SERIAL_SRCS))
PARALLEL_OBJS := $(patsubst src/%.c, build/parallel/%.o, $(PARALLEL_SRCS))

.PHONY: all serial parallel clean

all: serial

serial: build/serial_ga

parallel: build/parallel_ga

build/serial_ga: $(SERIAL_OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(SERIAL_OBJS) -lm -o $@

build/parallel_ga: $(PARALLEL_OBJS)
	@mkdir -p $(dir $@)
	$(MPICC) $(MPI_CFLAGS) $(PARALLEL_OBJS) -lm -o $@

build/serial/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

build/parallel/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(MPICC) $(MPI_CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
