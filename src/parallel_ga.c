#include "parallel_ga.h"

#include "random_utils.h"

#include <float.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIGRATION_INTERVAL 50

typedef struct {
    int size;
    int **population;
    int **next_population;
    int *population_storage;
    int *next_storage;
    double *lengths;
    double *next_lengths;
    double *fitness;
    double *next_fitness;
} PopulationBuffers;

static double **g_dist_matrix = NULL;
static int g_num_cities = 0;
static int g_total_population = 0;

static int **allocate_population_array(int population_size, int num_cities, int **storage_out) {
    int **array = calloc((size_t)population_size, sizeof(int *));
    int *storage = malloc((size_t)population_size * (size_t)num_cities * sizeof(int));
    if (!array || !storage) {
        free(array);
        free(storage);
        return NULL;
    }
    for (int i = 0; i < population_size; ++i) {
        array[i] = storage + (size_t)i * (size_t)num_cities;
    }
    *storage_out = storage;
    return array;
}

static int allocate_population_buffers(PopulationBuffers *buffers, int population_size, int num_cities) {
    memset(buffers, 0, sizeof(*buffers));
    buffers->size = population_size;
    buffers->population = allocate_population_array(population_size, num_cities, &buffers->population_storage);
    buffers->next_population = allocate_population_array(population_size, num_cities, &buffers->next_storage);
    buffers->lengths = calloc((size_t)population_size, sizeof(double));
    buffers->next_lengths = calloc((size_t)population_size, sizeof(double));
    buffers->fitness = calloc((size_t)population_size, sizeof(double));
    buffers->next_fitness = calloc((size_t)population_size, sizeof(double));
    if (!buffers->population || !buffers->next_population || !buffers->lengths || !buffers->next_lengths ||
        !buffers->fitness || !buffers->next_fitness) {
        free(buffers->population);
        free(buffers->next_population);
        free(buffers->population_storage);
        free(buffers->next_storage);
        free(buffers->lengths);
        free(buffers->next_lengths);
        free(buffers->fitness);
        free(buffers->next_fitness);
        memset(buffers, 0, sizeof(*buffers));
        return -1;
    }
    return 0;
}

static void free_population_buffers(PopulationBuffers *buffers) {
    if (!buffers) {
        return;
    }
    free(buffers->population);
    free(buffers->next_population);
    free(buffers->population_storage);
    free(buffers->next_storage);
    free(buffers->lengths);
    free(buffers->next_lengths);
    free(buffers->fitness);
    free(buffers->next_fitness);
    memset(buffers, 0, sizeof(*buffers));
}

static void initialize_individual(int *tour, int num_cities) {
    for (int i = 0; i < num_cities; ++i) {
        tour[i] = i;
    }
    shuffle_array(tour, num_cities);
}

static void evaluate_individual(const int *tour, double *length_out, double *fitness_out) {
    double length = evaluate_tour((int *)tour, g_dist_matrix, g_num_cities);
    if (length_out) {
        *length_out = length;
    }
    if (fitness_out) {
        *fitness_out = 1.0 / (length + 1e-9);
    }
}

static void select_best_indices(const double *lengths, int population_size, int count, int *indices) {
    if (count <= 0) {
        return;
    }
    char *selected = calloc((size_t)population_size, sizeof(char));
    if (!selected) {
        return;
    }
    for (int i = 0; i < count; ++i) {
        double best_length = DBL_MAX;
        int best_idx = 0;
        for (int j = 0; j < population_size; ++j) {
            if (selected[j]) {
                continue;
            }
            if (lengths[j] < best_length) {
                best_length = lengths[j];
                best_idx = j;
            }
        }
        indices[i] = best_idx;
        selected[best_idx] = 1;
    }
    free(selected);
}

static void select_worst_indices(const double *lengths, int population_size, int count, int *indices) {
    if (count <= 0) {
        return;
    }
    char *selected = calloc((size_t)population_size, sizeof(char));
    if (!selected) {
        return;
    }
    for (int i = 0; i < count; ++i) {
        double worst_length = -DBL_MAX;
        int worst_idx = 0;
        for (int j = 0; j < population_size; ++j) {
            if (selected[j]) {
                continue;
            }
            if (lengths[j] > worst_length) {
                worst_length = lengths[j];
                worst_idx = j;
            }
        }
        indices[i] = worst_idx;
        selected[worst_idx] = 1;
    }
    free(selected);
}

static int resize_population_buffers(PopulationBuffers *buffers, int new_size, int num_cities) {
    if (new_size == buffers->size) {
        return 0;
    }

    PopulationBuffers new_buffers = {0};
    if (allocate_population_buffers(&new_buffers, new_size, num_cities) != 0) {
        return -1;
    }

    int keep = buffers->size < new_size ? buffers->size : new_size;
    if (keep > 0 && buffers->population) {
        int *best_indices = malloc((size_t)keep * sizeof(int));
        if (!best_indices) {
            free_population_buffers(&new_buffers);
            return -1;
        }
        select_best_indices(buffers->lengths, buffers->size, keep, best_indices);
        for (int i = 0; i < keep; ++i) {
            int src_idx = best_indices[i];
            memcpy(new_buffers.population[i], buffers->population[src_idx],
                   (size_t)num_cities * sizeof(int));
            new_buffers.lengths[i] = buffers->lengths[src_idx];
            new_buffers.fitness[i] = buffers->fitness[src_idx];
        }
        free(best_indices);
    }

    for (int i = keep; i < new_size; ++i) {
        initialize_individual(new_buffers.population[i], num_cities);
        two_opt(new_buffers.population[i], g_dist_matrix, num_cities);
        evaluate_individual(new_buffers.population[i], &new_buffers.lengths[i], &new_buffers.fitness[i]);
    }

    free_population_buffers(buffers);
    *buffers = new_buffers;
    return 0;
}

void migrate_elites(int **population, double *lengths, double *fitness, int population_size, int num_cities,
                    int rank, int size, MPI_Comm comm) {
    int local_elite_count = (int)(0.05 * population_size);
    if (local_elite_count < 1) {
        local_elite_count = 1;
    }

    int *elite_indices = malloc((size_t)local_elite_count * sizeof(int));
    if (!elite_indices) {
        return;
    }
    select_best_indices(lengths, population_size, local_elite_count, elite_indices);

    int *send_tours = malloc((size_t)local_elite_count * (size_t)num_cities * sizeof(int));
    double *send_lengths = malloc((size_t)local_elite_count * sizeof(double));
    if (!send_tours || !send_lengths) {
        free(elite_indices);
        free(send_tours);
        free(send_lengths);
        return;
    }
    for (int i = 0; i < local_elite_count; ++i) {
        int idx = elite_indices[i];
        memcpy(send_tours + (size_t)i * (size_t)num_cities, population[idx], (size_t)num_cities * sizeof(int));
        send_lengths[i] = lengths[idx];
    }
    free(elite_indices);

    int *counts = NULL;
    int *tour_counts = NULL;
    int *tour_displs = NULL;
    int *length_displs = NULL;
    int total_elites = 0;
    int total_tour_entries = 0;
    if (rank == 0) {
        counts = calloc((size_t)size, sizeof(int));
        tour_counts = calloc((size_t)size, sizeof(int));
        tour_displs = calloc((size_t)size, sizeof(int));
        length_displs = calloc((size_t)size, sizeof(int));
    }

    int send_count = local_elite_count;
    MPI_Gather(&send_count, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);

    if (rank == 0) {
        int offset_tour = 0;
        int offset_length = 0;
        for (int i = 0; i < size; ++i) {
            tour_counts[i] = counts[i] * num_cities;
            tour_displs[i] = offset_tour;
            length_displs[i] = offset_length;
            offset_tour += tour_counts[i];
            offset_length += counts[i];
        }
        total_tour_entries = offset_tour;
        total_elites = offset_length;
    }

    int *all_tours = NULL;
    double *all_lengths = NULL;
    if (rank == 0 && total_tour_entries > 0) {
        all_tours = malloc((size_t)total_tour_entries * sizeof(int));
    }
    if (rank == 0 && total_elites > 0) {
        all_lengths = malloc((size_t)total_elites * sizeof(double));
    }

    MPI_Gatherv(send_tours, local_elite_count * num_cities, MPI_INT, all_tours, tour_counts, tour_displs, MPI_INT, 0,
                comm);
    MPI_Gatherv(send_lengths, local_elite_count, MPI_DOUBLE, all_lengths, counts, length_displs, MPI_DOUBLE, 0, comm);

    free(send_tours);
    free(send_lengths);

    int global_elite_count = (int)(0.05 * g_total_population);
    if (global_elite_count < 1) {
        global_elite_count = 1;
    }
    if (rank == 0 && total_elites < global_elite_count) {
        global_elite_count = total_elites;
    }

    int *broadcast_tours = NULL;
    double *broadcast_lengths = NULL;
    if (rank == 0 && global_elite_count > 0) {
        broadcast_tours = malloc((size_t)global_elite_count * (size_t)num_cities * sizeof(int));
        broadcast_lengths = malloc((size_t)global_elite_count * sizeof(double));
        if (!broadcast_tours || !broadcast_lengths) {
            free(broadcast_tours);
            free(broadcast_lengths);
            global_elite_count = 0;
        }
    }

    if (rank == 0 && global_elite_count > 0) {
        int *selected_indices = calloc((size_t)global_elite_count, sizeof(int));
        if (!selected_indices) {
            global_elite_count = 0;
        } else {
            char *used = calloc((size_t)total_elites, sizeof(char));
            if (!used) {
                global_elite_count = 0;
            } else {
                for (int i = 0; i < global_elite_count; ++i) {
                    double best_length = DBL_MAX;
                    int best_idx = 0;
                    for (int j = 0; j < total_elites; ++j) {
                        if (used[j]) {
                            continue;
                        }
                        if (all_lengths[j] < best_length) {
                            best_length = all_lengths[j];
                            best_idx = j;
                        }
                    }
                    used[best_idx] = 1;
                    selected_indices[i] = best_idx;
                    memcpy(broadcast_tours + (size_t)i * (size_t)num_cities,
                           all_tours + (size_t)best_idx * (size_t)num_cities,
                           (size_t)num_cities * sizeof(int));
                    broadcast_lengths[i] = all_lengths[best_idx];
                }
                free(used);
            }
            free(selected_indices);
        }
    }

    MPI_Bcast(&global_elite_count, 1, MPI_INT, 0, comm);

    if (global_elite_count > 0) {
        if (rank != 0) {
            broadcast_tours = malloc((size_t)global_elite_count * (size_t)num_cities * sizeof(int));
            broadcast_lengths = malloc((size_t)global_elite_count * sizeof(double));
        }
        MPI_Bcast(broadcast_tours, global_elite_count * num_cities, MPI_INT, 0, comm);
        MPI_Bcast(broadcast_lengths, global_elite_count, MPI_DOUBLE, 0, comm);

        int replace_count = global_elite_count < population_size ? global_elite_count : population_size;
        int *worst_indices = malloc((size_t)replace_count * sizeof(int));
        if (worst_indices) {
            select_worst_indices(lengths, population_size, replace_count, worst_indices);
            for (int i = 0; i < replace_count; ++i) {
                int dest = worst_indices[i];
                memcpy(population[dest], broadcast_tours + (size_t)i * (size_t)num_cities,
                       (size_t)num_cities * sizeof(int));
                evaluate_individual(population[dest], &lengths[dest], &fitness[dest]);
            }
            free(worst_indices);
        }
        free(broadcast_tours);
        free(broadcast_lengths);
    }

    if (rank == 0) {
        free(all_tours);
        free(all_lengths);
        free(counts);
        free(tour_counts);
        free(tour_displs);
        free(length_displs);
    }
}

int run_parallel_ga(const City *cities, int num_cities, double **dist_matrix, const GAParams *params, int *best_tour,
                    double *best_length, MPI_Comm comm) {
    if (!cities || num_cities <= 0 || !dist_matrix || !params) {
        return -1;
    }

    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank(comm, &world_rank);
    MPI_Comm_size(comm, &world_size);

    g_dist_matrix = dist_matrix;
    g_num_cities = num_cities;
    g_total_population = params->population_size;
    const int sync_interval = params->sync_interval > 0 ? params->sync_interval : 1;

    int base = params->population_size / world_size;
    int remainder = params->population_size % world_size;
    int local_size = base + (world_rank < remainder ? 1 : 0);
    if (local_size < 2) {
        local_size = 2;
    }

    seed_random(params->seed + (unsigned int)world_rank * 97u);
    set_two_opt_limit(params->two_opt_max_swaps);

    PopulationBuffers buffers = {0};
    if (allocate_population_buffers(&buffers, local_size, num_cities) != 0) {
        return -1;
    }

    for (int i = 0; i < buffers.size; ++i) {
        initialize_individual(buffers.population[i], num_cities);
        two_opt(buffers.population[i], dist_matrix, num_cities);
        evaluate_individual(buffers.population[i], &buffers.lengths[i], &buffers.fitness[i]);
    }

    double local_best_length = DBL_MAX;
    int *local_best_tour = malloc((size_t)num_cities * sizeof(int));
    int *global_best_tour = malloc((size_t)num_cities * sizeof(int));
    if (!local_best_tour || !global_best_tour) {
        free(local_best_tour);
        free(global_best_tour);
        free_population_buffers(&buffers);
        return -1;
    }

    for (int i = 0; i < buffers.size; ++i) {
        if (buffers.lengths[i] < local_best_length) {
            local_best_length = buffers.lengths[i];
            memcpy(local_best_tour, buffers.population[i], (size_t)num_cities * sizeof(int));
        }
    }

    struct {
        double value;
        int rank;
    } local_info = {local_best_length, world_rank}, global_info = {0.0, 0};

    MPI_Allreduce(&local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    if (world_rank == global_info.rank) {
        memcpy(global_best_tour, local_best_tour, (size_t)num_cities * sizeof(int));
    }
    MPI_Bcast(global_best_tour, num_cities, MPI_INT, global_info.rank, comm);
    double global_best_length = 0.0;
    if (world_rank == global_info.rank) {
        global_best_length = local_best_length;
    }
    MPI_Bcast(&global_best_length, 1, MPI_DOUBLE, global_info.rank, comm);

    double block_start = MPI_Wtime();

    for (int generation = 0; generation < params->generations; ++generation) {
        memcpy(buffers.next_population[0], global_best_tour, (size_t)num_cities * sizeof(int));
        buffers.next_lengths[0] = global_best_length;
        buffers.next_fitness[0] = 1.0 / (global_best_length + 1e-9);

        for (int i = 1; i < buffers.size; ++i) {
            int parent_a_idx = 0;
            int parent_b_idx = 0;
            tournament_selection(buffers.population, buffers.fitness, buffers.size, params->tournament_size,
                                 &parent_a_idx);
            tournament_selection(buffers.population, buffers.fitness, buffers.size, params->tournament_size,
                                 &parent_b_idx);

            if (rand_double() < params->crossover_rate) {
                pmx_crossover(buffers.population[parent_a_idx], buffers.population[parent_b_idx],
                              buffers.next_population[i], num_cities);
            } else {
                memcpy(buffers.next_population[i], buffers.population[parent_a_idx], (size_t)num_cities * sizeof(int));
            }

            inversion_mutation(buffers.next_population[i], num_cities, params->mutation_rate);
            two_opt(buffers.next_population[i], dist_matrix, num_cities);
            evaluate_individual(buffers.next_population[i], &buffers.next_lengths[i], &buffers.next_fitness[i]);
        }

        int **tmp_pop = buffers.population;
        buffers.population = buffers.next_population;
        buffers.next_population = tmp_pop;

        double *tmp_lengths = buffers.lengths;
        buffers.lengths = buffers.next_lengths;
        buffers.next_lengths = tmp_lengths;

        double *tmp_fitness = buffers.fitness;
        buffers.fitness = buffers.next_fitness;
        buffers.next_fitness = tmp_fitness;

        local_best_length = DBL_MAX;
        for (int i = 0; i < buffers.size; ++i) {
            if (buffers.lengths[i] < local_best_length) {
                local_best_length = buffers.lengths[i];
                memcpy(local_best_tour, buffers.population[i], (size_t)num_cities * sizeof(int));
            }
        }

        int sync_now = ((generation + 1) % sync_interval == 0) || (generation == params->generations - 1);
        if (sync_now) {
            local_info.value = local_best_length;
            local_info.rank = world_rank;
            MPI_Allreduce(&local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

            if (world_rank == global_info.rank) {
                memcpy(global_best_tour, local_best_tour, (size_t)num_cities * sizeof(int));
                global_best_length = local_best_length;
            }
            MPI_Bcast(global_best_tour, num_cities, MPI_INT, global_info.rank, comm);
            MPI_Bcast(&global_best_length, 1, MPI_DOUBLE, global_info.rank, comm);
        } else if (local_best_length < global_best_length) {
            global_best_length = local_best_length;
            memcpy(global_best_tour, local_best_tour, (size_t)num_cities * sizeof(int));
        }

        if ((generation + 1) % MIGRATION_INTERVAL == 0) {
            migrate_elites(buffers.population, buffers.lengths, buffers.fitness, buffers.size, num_cities, world_rank,
                           world_size, comm);

            double block_end = MPI_Wtime();
            double local_block_time = block_end - block_start;
            double *all_block_times = malloc((size_t)world_size * sizeof(double));
            double *speeds = malloc((size_t)world_size * sizeof(double));
            int *new_sizes = malloc((size_t)world_size * sizeof(int));
            if (all_block_times && speeds && new_sizes) {
                MPI_Allgather(&local_block_time, 1, MPI_DOUBLE, all_block_times, 1, MPI_DOUBLE, comm);
                double total_speed = 0.0;
                for (int i = 0; i < world_size; ++i) {
                    speeds[i] = all_block_times[i] > 0.0 ? 1.0 / all_block_times[i] : 1.0;
                    total_speed += speeds[i];
                }
                int assigned = 0;
                for (int i = 0; i < world_size; ++i) {
                    double fraction = speeds[i] / total_speed;
                    new_sizes[i] = (int)(g_total_population * fraction);
                    if (new_sizes[i] < 2) {
                        new_sizes[i] = 2;
                    }
                    assigned += new_sizes[i];
                }
                while (assigned < g_total_population) {
                    int best_idx = 0;
                    double best_speed = speeds[0];
                    for (int i = 1; i < world_size; ++i) {
                        if (speeds[i] > best_speed) {
                            best_speed = speeds[i];
                            best_idx = i;
                        }
                    }
                    new_sizes[best_idx] += 1;
                    assigned += 1;
                }
                while (assigned > g_total_population) {
                    int worst_idx = 0;
                    double worst_speed = speeds[0];
                    for (int i = 1; i < world_size; ++i) {
                        if (speeds[i] < worst_speed && new_sizes[i] > 2) {
                            worst_speed = speeds[i];
                            worst_idx = i;
                        }
                    }
                    if (new_sizes[worst_idx] > 2) {
                        new_sizes[worst_idx] -= 1;
                        assigned -= 1;
                    } else {
                        break;
                    }
                }

                MPI_Bcast(new_sizes, world_size, MPI_INT, 0, comm);
                int desired_size = new_sizes[world_rank];
                if (desired_size != buffers.size) {
                    if (resize_population_buffers(&buffers, desired_size, num_cities) != 0) {
                        free(all_block_times);
                        free(speeds);
                        free(new_sizes);
                        free(local_best_tour);
                        free(global_best_tour);
                        free_population_buffers(&buffers);
                        return -1;
                    }
                }
            }
            free(all_block_times);
            free(speeds);
            free(new_sizes);
            block_start = MPI_Wtime();
        }
    }

    local_best_length = DBL_MAX;
    for (int i = 0; i < buffers.size; ++i) {
        if (buffers.lengths[i] < local_best_length) {
            local_best_length = buffers.lengths[i];
            memcpy(local_best_tour, buffers.population[i], (size_t)num_cities * sizeof(int));
        }
    }
    local_info.value = local_best_length;
    local_info.rank = world_rank;
    MPI_Allreduce(&local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    if (world_rank == global_info.rank) {
        memcpy(global_best_tour, local_best_tour, (size_t)num_cities * sizeof(int));
        global_best_length = local_best_length;
    }
    MPI_Bcast(global_best_tour, num_cities, MPI_INT, global_info.rank, comm);
    MPI_Bcast(&global_best_length, 1, MPI_DOUBLE, global_info.rank, comm);

    if (best_tour && world_rank == 0) {
        memcpy(best_tour, global_best_tour, (size_t)num_cities * sizeof(int));
    }
    if (best_length && world_rank == 0) {
        *best_length = global_best_length;
    }

    free(local_best_tour);
    free(global_best_tour);
    free_population_buffers(&buffers);
    return 0;
}
