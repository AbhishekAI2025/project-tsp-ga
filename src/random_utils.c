#include "random_utils.h"

#include <stdint.h>
#include <stdlib.h>

static uint32_t rng_state = 1u;

static uint32_t lcg_next(void) {
    rng_state = 1664525u * rng_state + 1013904223u;
    return rng_state;
}

void seed_random(unsigned int seed) {
    if (seed == 0u) {
        seed = 1u;
    }
    rng_state = seed;
}

double rand_double(void) {
    return (lcg_next() + 1.0) / (UINT32_MAX + 1.0);
}

int rand_int(int min_inclusive, int max_exclusive) {
    if (max_exclusive <= min_inclusive) {
        return min_inclusive;
    }
    const uint32_t span = (uint32_t)(max_exclusive - min_inclusive);
    return min_inclusive + (int)(lcg_next() % span);
}

void shuffle_array(int *array, int length) {
    if (!array || length <= 1) {
        return;
    }
    for (int i = length - 1; i > 0; --i) {
        const int j = rand_int(0, i + 1);
        const int tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }
}
