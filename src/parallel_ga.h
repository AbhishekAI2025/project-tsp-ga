#ifndef PARALLEL_GA_H
#define PARALLEL_GA_H

#include "ga.h"

#include <mpi.h>

int run_parallel_ga(const City *cities, int num_cities, double **dist_matrix, const GAParams *params, int *best_tour,
                    double *best_length, MPI_Comm comm);
void migrate_elites(int **population, double *lengths, double *fitness, int population_size, int num_cities,
                    int rank, int size, MPI_Comm comm);

#endif
