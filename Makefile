CC ?= mpicc
CFLAGS ?= -std=c11 -Wall -Wextra -O2

SRC_DIR := src
SRCS := $(SRC_DIR)/main.c \
        $(SRC_DIR)/ga.c \
        $(SRC_DIR)/parallel_ga.c \
        $(SRC_DIR)/tsp_utils.c \
        $(SRC_DIR)/random_utils.c

OBJS := $(SRCS:$(SRC_DIR)/%.c=build/%.o)

.PHONY: all clean

all: tsp_ga

tsp_ga: $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(OBJS) -lm -o $@

build/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build tsp_ga
