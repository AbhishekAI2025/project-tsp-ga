#ifndef GA_H
#define GA_H

#include "tsp_utils.h"

typedef struct {
    int population_size;
    int generations;
    double crossover_rate;
    double mutation_rate;
    int tournament_size;
    int two_opt_max_swaps;
    unsigned int seed;
    int sync_interval;  // MPI best-tour synchronization interval (generations)
} GAParams;

double evaluate_tour(int *tour, double **dist_matrix, int num_cities);
void tournament_selection(int **population, double *fitness, int population_size, int k, int *selected_index);
void pmx_crossover(const int *parent1, const int *parent2, int *child, int num_cities);
void inversion_mutation(int *tour, int num_cities, double mutation_rate);
void two_opt(int *tour, double **dist_matrix, int num_cities);
void set_two_opt_limit(int limit);
int run_serial_ga(const City *cities, int num_cities, double **dist_matrix, const GAParams *params, int *best_tour,
                  double *best_length);

#endif
