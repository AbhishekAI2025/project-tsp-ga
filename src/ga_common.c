#include "ga_common.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint32_t rng_state = 1u;

static uint32_t rng_next(void) {
    rng_state = 1664525u * rng_state + 1013904223u;
    return rng_state;
}

void ga_seed_rng(unsigned int seed) {
    if (seed == 0) {
        seed = 1u;
    }
    rng_state = seed;
}

double ga_random_unit(void) {
    return (rng_next() + 1.0) / (UINT32_MAX + 1.0);
}

int ga_random_int(int min_inclusive, int max_exclusive) {
    const uint32_t span = (uint32_t)(max_exclusive - min_inclusive);
    return min_inclusive + (int)(rng_next() % span);
}

void ga_shuffle(int *array, int n) {
    for (int i = n - 1; i > 0; --i) {
        const int j = ga_random_int(0, i + 1);
        const int tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }
}

void ga_initialize_individual(Individual *ind, const TSPInstance *instance) {
    const int n = instance->dimension;
    for (int i = 0; i < n; ++i) {
        ind->tour[i] = i;
    }
    ga_shuffle(ind->tour, n);
    ga_evaluate_individual(ind, instance);
}

void ga_evaluate_individual(Individual *ind, const TSPInstance *instance) {
    ind->length = tsp_tour_length(instance, ind->tour);
    ind->fitness = 1.0 / (ind->length + 1e-9);
}

void ga_copy_individual(const Individual *src, Individual *dst, int dimension) {
    memcpy(dst->tour, src->tour, (size_t)dimension * sizeof(int));
    dst->length = src->length;
    dst->fitness = src->fitness;
}

void ga_pmx_crossover(const int *parent1, const int *parent2, int *child, int dimension) {
    if (dimension < 2) {
        if (dimension == 1) {
            child[0] = parent1[0];
        }
        return;
    }
    int cut1 = ga_random_int(0, dimension - 1);
    int cut2 = ga_random_int(cut1 + 1, dimension);

    for (int i = 0; i < dimension; ++i) {
        child[i] = -1;
    }

    int *index_in_parent1 = malloc((size_t)dimension * sizeof(int));
    int *index_in_parent2 = malloc((size_t)dimension * sizeof(int));
    int *in_child = calloc((size_t)dimension, sizeof(int));
    if (!index_in_parent1 || !index_in_parent2 || !in_child) {
        /* fall back to copying parent1 to avoid undefined behaviour */
        for (int i = 0; i < dimension; ++i) {
            child[i] = parent1[i];
        }
        free(index_in_parent1);
        free(index_in_parent2);
        free(in_child);
        return;
    }

    for (int i = 0; i < dimension; ++i) {
        index_in_parent1[parent1[i]] = i;
        index_in_parent2[parent2[i]] = i;
    }

    for (int i = cut1; i <= cut2; ++i) {
        const int gene = parent1[i];
        child[i] = gene;
        in_child[gene] = 1;
    }

    for (int i = cut1; i <= cut2; ++i) {
        int gene = parent2[i];
        if (in_child[gene]) {
            continue;
        }
        int position = i;
        while (child[position] != -1) {
            const int mapped_gene = parent1[position];
            position = index_in_parent2[mapped_gene];
        }
        child[position] = gene;
        in_child[gene] = 1;
    }

    for (int i = 0; i < dimension; ++i) {
        if (child[i] == -1) {
            const int gene = parent2[i];
            child[i] = gene;
            in_child[gene] = 1;
        }
    }

    free(index_in_parent1);
    free(index_in_parent2);
    free(in_child);
}

void ga_inversion_mutation(int *tour, int dimension) {
    const int i = ga_random_int(0, dimension);
    const int j = ga_random_int(i, dimension);
    int left = i;
    int right = j;
    while (left < right) {
        const int tmp = tour[left];
        tour[left] = tour[right];
        tour[right] = tmp;
        ++left;
        --right;
    }
}

int ga_tournament_select(const Individual *population, int population_size, int tournament_size) {
    int best_idx = -1;
    double best_fitness = -DBL_MAX;
    for (int i = 0; i < tournament_size; ++i) {
        const int idx = ga_random_int(0, population_size);
        if (population[idx].fitness > best_fitness) {
            best_idx = idx;
            best_fitness = population[idx].fitness;
        }
    }
    return best_idx;
}

int ga_two_opt_local_search(int *tour, const TSPInstance *instance, int max_iterations) {
    const int n = instance->dimension;
    int improved = 0;

    for (int iter = 0; iter < max_iterations; ++iter) {
        int delta_found = 0;
        for (int i = 0; i < n - 1; ++i) {
            const int a = tour[i];
            const int b = tour[i + 1];
            for (int j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1) {
                    continue; /* would break the loop */
                }
                const int c = tour[j];
                const int d = tour[(j + 1) % n];

                const double current = tsp_distance(instance, a, b) + tsp_distance(instance, c, d);
                const double proposed = tsp_distance(instance, a, c) + tsp_distance(instance, b, d);
                if (proposed + 1e-9 < current) {
                    int left = i + 1;
                    int right = j;
                    while (left < right) {
                        const int tmp = tour[left];
                        tour[left] = tour[right];
                        tour[right] = tmp;
                        ++left;
                        --right;
                    }
                    improved = 1;
                    delta_found = 1;
                    break;
                }
            }
            if (delta_found) {
                break;
            }
        }
        if (!delta_found) {
            break;
        }
    }

    return improved;
}
