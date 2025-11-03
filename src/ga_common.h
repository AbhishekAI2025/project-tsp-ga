#ifndef GA_COMMON_H
#define GA_COMMON_H

#include "tsp.h"

typedef struct {
    int population_size;
    int generations;
    double crossover_rate;
    double mutation_rate;
    int tournament_size;
    int two_opt_iterations;
    unsigned int seed;
} GAParams;

typedef struct {
    int *tour;
    double length;
    double fitness;
} Individual;

void ga_seed_rng(unsigned int seed);
double ga_random_unit(void);
int ga_random_int(int min_inclusive, int max_exclusive);
void ga_shuffle(int *array, int n);
void ga_initialize_individual(Individual *ind, const TSPInstance *instance);
void ga_evaluate_individual(Individual *ind, const TSPInstance *instance);
void ga_copy_individual(const Individual *src, Individual *dst, int dimension);
void ga_pmx_crossover(const int *parent1, const int *parent2, int *child, int dimension);
void ga_inversion_mutation(int *tour, int dimension);
int ga_tournament_select(const Individual *population, int population_size, int tournament_size);
int ga_two_opt_local_search(int *tour, const TSPInstance *instance, int max_iterations);

#endif
