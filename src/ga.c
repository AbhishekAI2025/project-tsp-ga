#include "ga.h"

#include "random_utils.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_two_opt_limit = 0;

void set_two_opt_limit(int limit) {
    g_two_opt_limit = limit;
}

static int **allocate_population(int population_size, int num_cities) {
    int **population = calloc((size_t)population_size, sizeof(int *));
    int *storage = malloc((size_t)population_size * (size_t)num_cities * sizeof(int));
    if (!population || !storage) {
        free(population);
        free(storage);
        return NULL;
    }
    for (int i = 0; i < population_size; ++i) {
        population[i] = storage + (size_t)i * (size_t)num_cities;
    }
    return population;
}

static void free_population(int **population) {
    if (!population) {
        return;
    }
    free(population[0]);
    free(population);
}

double evaluate_tour(int *tour, double **dist_matrix, int num_cities) {
    double length = 0.0;
    for (int i = 0; i < num_cities - 1; ++i) {
        length += dist_matrix[tour[i]][tour[i + 1]];
    }
    length += dist_matrix[tour[num_cities - 1]][tour[0]];
    return length;
}

void tournament_selection(int **population, double *fitness, int population_size, int k, int *selected_index) {
    if (!selected_index) {
        return;
    }
    int best_idx = rand_int(0, population_size);
    double best_fit = fitness[best_idx];
    const int iterations = k < population_size ? k : population_size;
    for (int i = 1; i < iterations; ++i) {
        const int idx = rand_int(0, population_size);
        if (fitness[idx] > best_fit) {
            best_fit = fitness[idx];
            best_idx = idx;
        }
    }
    *selected_index = best_idx;
}

void pmx_crossover(const int *parent1, const int *parent2, int *child, int num_cities) {
    if (num_cities <= 1) {
        if (num_cities == 1) {
            child[0] = parent1[0];
        }
        return;
    }

    int cut1 = rand_int(0, num_cities - 1);
    int cut2 = rand_int(cut1 + 1, num_cities);

    int *index_in_parent1 = malloc((size_t)num_cities * sizeof(int));
    int *index_in_parent2 = malloc((size_t)num_cities * sizeof(int));
    int *in_child = calloc((size_t)num_cities, sizeof(int));
    if (!index_in_parent1 || !index_in_parent2 || !in_child) {
        memcpy(child, parent1, (size_t)num_cities * sizeof(int));
        free(index_in_parent1);
        free(index_in_parent2);
        free(in_child);
        return;
    }

    for (int i = 0; i < num_cities; ++i) {
        child[i] = -1;
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

    for (int i = 0; i < num_cities; ++i) {
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

void inversion_mutation(int *tour, int num_cities, double mutation_rate) {
    if (rand_double() >= mutation_rate) {
        return;
    }
    const int start = rand_int(0, num_cities);
    const int end = rand_int(start, num_cities);
    int left = start;
    int right = end;
    while (left < right) {
        const int tmp = tour[left];
        tour[left] = tour[right];
        tour[right] = tmp;
        ++left;
        --right;
    }
}

void two_opt(int *tour, double **dist_matrix, int num_cities) {
    if (num_cities < 4) {
        return;
    }
    int improved = 1;
    int iteration = 0;
    const int max_iterations = g_two_opt_limit > 0 ? g_two_opt_limit : num_cities * num_cities;
    while (improved && iteration < max_iterations) {
        improved = 0;
        for (int i = 0; i < num_cities - 1; ++i) {
            int a = tour[i];
            int b = tour[i + 1];
            for (int j = i + 2; j < num_cities; ++j) {
                if (i == 0 && j == num_cities - 1) {
                    continue;
                }
                int c = tour[j];
                int d = tour[(j + 1) % num_cities];
                double current = dist_matrix[a][b] + dist_matrix[c][d];
                double swapped = dist_matrix[a][c] + dist_matrix[b][d];
                if (swapped + 1e-9 < current) {
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
                    break;
                }
            }
            if (improved) {
                break;
            }
        }
        ++iteration;
    }
}

static void initialize_individual(int *tour, int num_cities) {
    for (int i = 0; i < num_cities; ++i) {
        tour[i] = i;
    }
    shuffle_array(tour, num_cities);
}

int run_serial_ga(const City *cities, int num_cities, double **dist_matrix, const GAParams *params, int *best_tour,
                  double *best_length) {
    if (!cities || num_cities <= 0 || !dist_matrix || !params) {
        return -1;
    }
    if (params->population_size < 2 || params->generations <= 0) {
        return -1;
    }

    (void)cities;
    seed_random(params->seed);
    set_two_opt_limit(params->two_opt_max_swaps);

    int **population = allocate_population(params->population_size, num_cities);
    int **next_population = allocate_population(params->population_size, num_cities);
    double *fitness = calloc((size_t)params->population_size, sizeof(double));
    double *lengths = calloc((size_t)params->population_size, sizeof(double));
    double *next_fitness = calloc((size_t)params->population_size, sizeof(double));
    double *next_lengths = calloc((size_t)params->population_size, sizeof(double));
    if (!population || !next_population || !fitness || !lengths || !next_fitness || !next_lengths) {
        free_population(population);
        free_population(next_population);
        free(fitness);
        free(lengths);
        free(next_fitness);
        free(next_lengths);
        return -1;
    }

    double global_best_length = DBL_MAX;
    int *global_best_tour = malloc((size_t)num_cities * sizeof(int));
    if (!global_best_tour) {
        free_population(population);
        free_population(next_population);
        free(fitness);
        free(lengths);
        free(next_fitness);
        free(next_lengths);
        return -1;
    }

    for (int i = 0; i < params->population_size; ++i) {
        initialize_individual(population[i], num_cities);
        two_opt(population[i], dist_matrix, num_cities);
        lengths[i] = evaluate_tour(population[i], dist_matrix, num_cities);
        fitness[i] = 1.0 / (lengths[i] + 1e-9);
        if (lengths[i] < global_best_length) {
            global_best_length = lengths[i];
            memcpy(global_best_tour, population[i], (size_t)num_cities * sizeof(int));
        }
    }

    for (int generation = 0; generation < params->generations; ++generation) {
        int best_idx = 0;
        double best_fit = fitness[0];
        for (int i = 1; i < params->population_size; ++i) {
            if (fitness[i] > best_fit) {
                best_fit = fitness[i];
                best_idx = i;
            }
        }

        memcpy(next_population[0], population[best_idx], (size_t)num_cities * sizeof(int));
        next_lengths[0] = lengths[best_idx];
        next_fitness[0] = fitness[best_idx];

        for (int i = 1; i < params->population_size; ++i) {
            int parent_a_idx = 0;
            int parent_b_idx = 0;
            tournament_selection(population, fitness, params->population_size, params->tournament_size, &parent_a_idx);
            tournament_selection(population, fitness, params->population_size, params->tournament_size, &parent_b_idx);

            if (rand_double() < params->crossover_rate) {
                pmx_crossover(population[parent_a_idx], population[parent_b_idx], next_population[i], num_cities);
            } else {
                memcpy(next_population[i], population[parent_a_idx], (size_t)num_cities * sizeof(int));
            }

            inversion_mutation(next_population[i], num_cities, params->mutation_rate);
            two_opt(next_population[i], dist_matrix, num_cities);
            next_lengths[i] = evaluate_tour(next_population[i], dist_matrix, num_cities);
            next_fitness[i] = 1.0 / (next_lengths[i] + 1e-9);

            if (next_lengths[i] < global_best_length) {
                global_best_length = next_lengths[i];
                memcpy(global_best_tour, next_population[i], (size_t)num_cities * sizeof(int));
            }
        }

        int **tmp_pop = population;
        population = next_population;
        next_population = tmp_pop;

        double *tmp_fit = fitness;
        fitness = next_fitness;
        next_fitness = tmp_fit;

        double *tmp_len = lengths;
        lengths = next_lengths;
        next_lengths = tmp_len;
    }

    if (best_tour) {
        memcpy(best_tour, global_best_tour, (size_t)num_cities * sizeof(int));
    }
    if (best_length) {
        *best_length = global_best_length;
    }

    free(global_best_tour);
    free_population(population);
    free_population(next_population);
    free(fitness);
    free(lengths);
    free(next_fitness);
    free(next_lengths);
    return 0;
}
