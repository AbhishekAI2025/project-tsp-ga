#include "ga_common.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    Individual *population;
    Individual *next_population;
    int *population_storage;
    int *next_storage;
} PopulationBuffers;

static int allocate_population_buffers(PopulationBuffers *buffers, int population_size, int dimension) {
    buffers->population = calloc((size_t)population_size, sizeof(Individual));
    buffers->next_population = calloc((size_t)population_size, sizeof(Individual));
    if (!buffers->population || !buffers->next_population) {
        free(buffers->population);
        free(buffers->next_population);
        buffers->population = NULL;
        buffers->next_population = NULL;
        return -1;
    }

    buffers->population_storage = malloc((size_t)population_size * (size_t)dimension * sizeof(int));
    buffers->next_storage = malloc((size_t)population_size * (size_t)dimension * sizeof(int));
    if (!buffers->population_storage || !buffers->next_storage) {
        free(buffers->population_storage);
        free(buffers->next_storage);
        free(buffers->population);
        free(buffers->next_population);
        buffers->population = NULL;
        buffers->next_population = NULL;
        buffers->population_storage = NULL;
        buffers->next_storage = NULL;
        return -1;
    }

    for (int i = 0; i < population_size; ++i) {
        buffers->population[i].tour = buffers->population_storage + (size_t)i * (size_t)dimension;
        buffers->next_population[i].tour = buffers->next_storage + (size_t)i * (size_t)dimension;
    }

    return 0;
}

static void free_population_buffers(PopulationBuffers *buffers) {
    free(buffers->population);
    free(buffers->next_population);
    free(buffers->population_storage);
    free(buffers->next_storage);
    buffers->population = NULL;
    buffers->next_population = NULL;
    buffers->population_storage = NULL;
    buffers->next_storage = NULL;
}

int ga_run_serial(const TSPInstance *instance, const GAParams *params, Individual *best_out, int log_interval,
                  FILE *log_stream) {
    const int dimension = instance->dimension;
    const int population_size = params->population_size;

    if (population_size < 2) {
        fprintf(stderr, "Population size must be at least 2\n");
        return -1;
    }

    PopulationBuffers buffers = {0};
    if (allocate_population_buffers(&buffers, population_size, dimension) != 0) {
        fprintf(stderr, "Failed to allocate population buffers\n");
        free_population_buffers(&buffers);
        return -1;
    }

    ga_seed_rng(params->seed);

    double best_length = DBL_MAX;
    int best_index = -1;

    for (int i = 0; i < population_size; ++i) {
        ga_initialize_individual(&buffers.population[i], instance);
        if (buffers.population[i].length < best_length) {
            best_length = buffers.population[i].length;
            best_index = i;
        }
    }

    Individual global_best = {0};
    global_best.tour = malloc((size_t)dimension * sizeof(int));
    if (!global_best.tour) {
        fprintf(stderr, "Out of memory allocating best tour buffer\n");
        free_population_buffers(&buffers);
        return -1;
    }
    ga_copy_individual(&buffers.population[best_index], &global_best, dimension);

    for (int generation = 0; generation < params->generations; ++generation) {
        ga_copy_individual(&global_best, &buffers.next_population[0], dimension);

        int fill_index = 1;
        while (fill_index < population_size) {
            const int parent_a_idx = ga_tournament_select(buffers.population, population_size, params->tournament_size);
            const int parent_b_idx = ga_tournament_select(buffers.population, population_size, params->tournament_size);
            const Individual *parent_a = &buffers.population[parent_a_idx];
            const Individual *parent_b = &buffers.population[parent_b_idx];

            Individual *child1 = &buffers.next_population[fill_index];
            Individual *child2 = NULL;
            if (fill_index + 1 < population_size) {
                child2 = &buffers.next_population[fill_index + 1];
            }

            const double rand_val = ga_random_unit();
            if (rand_val < params->crossover_rate) {
                ga_pmx_crossover(parent_a->tour, parent_b->tour, child1->tour, dimension);
                if (child2) {
                    ga_pmx_crossover(parent_b->tour, parent_a->tour, child2->tour, dimension);
                }
            } else {
                ga_copy_individual(parent_a, child1, dimension);
                if (child2) {
                    ga_copy_individual(parent_b, child2, dimension);
                }
            }

            if (ga_random_unit() < params->mutation_rate) {
                ga_inversion_mutation(child1->tour, dimension);
            }
            ga_two_opt_local_search(child1->tour, instance, params->two_opt_iterations);
            ga_evaluate_individual(child1, instance);

            if (child2) {
                if (ga_random_unit() < params->mutation_rate) {
                    ga_inversion_mutation(child2->tour, dimension);
                }
                ga_two_opt_local_search(child2->tour, instance, params->two_opt_iterations);
                ga_evaluate_individual(child2, instance);
            }

            fill_index += (child2 ? 2 : 1);
        }

        Individual *tmp = buffers.population;
        buffers.population = buffers.next_population;
        buffers.next_population = tmp;

        double best_generation_length = global_best.length;
        const Individual *best_generation_individual = NULL;
        for (int i = 0; i < population_size; ++i) {
            if (buffers.population[i].length < best_generation_length) {
                best_generation_length = buffers.population[i].length;
                best_generation_individual = &buffers.population[i];
            }
        }
        if (best_generation_individual) {
            ga_copy_individual(best_generation_individual, &global_best, dimension);
        }

        if (log_stream && (generation % log_interval == 0 || generation == params->generations - 1)) {
            fprintf(log_stream, "Generation %d: best length %.3f\n", generation, global_best.length);
            fflush(log_stream);
        }
    }

    if (best_out) {
        ga_copy_individual(&global_best, best_out, dimension);
    }

    free(global_best.tour);
    free_population_buffers(&buffers);
    return 0;
}
