#include "parallel_ga.h"

#include <float.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define MIGRATION_INTERVAL 50

typedef struct {
    Individual *population;
    Individual *next_population;
    int *population_storage;
    int *next_storage;
    int size;
} PopulationBuffers;

static void select_best_indices(const Individual *population, int population_size, int count, int *indices) {
    for (int i = 0; i < count; ++i) {
        indices[i] = -1;
    }
    double *lengths = malloc((size_t)count * sizeof(double));
    if (!lengths) {
        return;
    }
    for (int i = 0; i < count; ++i) {
        lengths[i] = DBL_MAX;
    }
    for (int i = 0; i < population_size; ++i) {
        const double length = population[i].length;
        for (int pos = 0; pos < count; ++pos) {
            if (length < lengths[pos]) {
                for (int shift = count - 1; shift > pos; --shift) {
                    lengths[shift] = lengths[shift - 1];
                    indices[shift] = indices[shift - 1];
                }
                lengths[pos] = length;
                indices[pos] = i;
                break;
            }
        }
    }
    free(lengths);
}

static void select_worst_indices(const Individual *population, int population_size, int count, int *indices) {
    for (int i = 0; i < count; ++i) {
        indices[i] = -1;
    }
    double *lengths = malloc((size_t)count * sizeof(double));
    if (!lengths) {
        return;
    }
    for (int i = 0; i < count; ++i) {
        lengths[i] = -DBL_MAX;
    }
    for (int i = 0; i < population_size; ++i) {
        const double length = population[i].length;
        for (int pos = 0; pos < count; ++pos) {
            if (length > lengths[pos]) {
                for (int shift = count - 1; shift > pos; --shift) {
                    lengths[shift] = lengths[shift - 1];
                    indices[shift] = indices[shift - 1];
                }
                lengths[pos] = length;
                indices[pos] = i;
                break;
            }
        }
    }
    free(lengths);
}

static int allocate_population_buffers(PopulationBuffers *buffers, int population_size, int dimension) {
    buffers->size = population_size;
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
    buffers->size = 0;
}

static int resize_population_buffers(PopulationBuffers *buffers, int new_size, int dimension,
                                     const TSPInstance *instance) {
    if (new_size == buffers->size) {
        return 0;
    }

    Individual *old_population = buffers->population;
    Individual *old_next = buffers->next_population;
    int *old_storage = buffers->population_storage;
    int *old_next_storage = buffers->next_storage;
    int old_size = buffers->size;

    PopulationBuffers new_buffers = {0};
    if (allocate_population_buffers(&new_buffers, new_size, dimension) != 0) {
        return -1;
    }

    const int keep = old_population ? (old_size < new_size ? old_size : new_size) : 0;
    if (old_population && keep > 0) {
        int *indices = malloc((size_t)keep * sizeof(int));
        int *selected = calloc((size_t)old_size, sizeof(int));
        if (!indices || !selected) {
            free(indices);
            free(selected);
            free_population_buffers(&new_buffers);
            return -1;
        }
        for (int i = 0; i < keep; ++i) {
            double best_length = DBL_MAX;
            int best_idx = -1;
            for (int j = 0; j < old_size; ++j) {
                if (selected[j]) {
                    continue;
                }
                if (old_population[j].length < best_length) {
                    best_length = old_population[j].length;
                    best_idx = j;
                }
            }
            if (best_idx >= 0) {
                selected[best_idx] = 1;
                indices[i] = best_idx;
            } else {
                indices[i] = 0;
            }
        }
        free(selected);
        for (int i = 0; i < keep; ++i) {
            ga_copy_individual(&old_population[indices[i]], &new_buffers.population[i], dimension);
        }
        free(indices);
    }

    for (int i = keep; i < new_size; ++i) {
        ga_initialize_individual(&new_buffers.population[i], instance);
    }

    for (int i = 0; i < new_size; ++i) {
        /* ensure derived metrics are up to date */
        ga_evaluate_individual(&new_buffers.population[i], instance);
    }

    buffers->population = new_buffers.population;
    buffers->next_population = new_buffers.next_population;
    buffers->population_storage = new_buffers.population_storage;
    buffers->next_storage = new_buffers.next_storage;
    buffers->size = new_buffers.size;

    free(old_population);
    free(old_next);
    free(old_storage);
    free(old_next_storage);
    return 0;
}

static void integrate_elites(PopulationBuffers *buffers, const TSPInstance *instance, const int *elite_tours,
                             const double *elite_lengths, int elite_count, int dimension) {
    if (elite_count <= 0 || buffers->size == 0) {
        return;
    }
    const int replace_count = elite_count < buffers->size ? elite_count : buffers->size;
    int *worst_indices = malloc((size_t)replace_count * sizeof(int));
    if (!worst_indices) {
        return;
    }
    select_worst_indices(buffers->population, buffers->size, replace_count, worst_indices);
    for (int i = 0; i < replace_count; ++i) {
        const int dest_idx = worst_indices[i] >= 0 ? worst_indices[i] : i;
        memcpy(buffers->population[dest_idx].tour, elite_tours + (size_t)i * (size_t)dimension,
               (size_t)dimension * sizeof(int));
        buffers->population[dest_idx].length = elite_lengths[i];
        buffers->population[dest_idx].fitness = 1.0 / (elite_lengths[i] + 1e-9);
    }
    free(worst_indices);
}

int ga_run_parallel(const TSPInstance *instance, const GAParams *params, Individual *best_out, int log_interval,
                    MPI_Comm comm, FILE *log_stream) {
    const int dimension = instance->dimension;
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    int base = params->population_size / world_size;
    int remainder = params->population_size % world_size;
    int local_size = base + (world_rank < remainder ? 1 : 0);

    if (local_size <= 0) {
        local_size = 1;
    }

    ga_seed_rng(params->seed + (unsigned int)world_rank * 17u);

    PopulationBuffers buffers = {0};
    if (allocate_population_buffers(&buffers, local_size, dimension) != 0) {
        if (log_stream && world_rank == 0) {
            fprintf(log_stream, "Failed to allocate population buffers for rank %d\n", world_rank);
        }
        return -1;
    }

    for (int i = 0; i < local_size; ++i) {
        ga_initialize_individual(&buffers.population[i], instance);
    }

    Individual global_best = {0};
    global_best.tour = malloc((size_t)dimension * sizeof(int));
    if (!global_best.tour) {
        free_population_buffers(&buffers);
        return -1;
    }
    global_best.length = DBL_MAX;
    global_best.fitness = 0.0;

    Individual local_best = {0};
    local_best.tour = malloc((size_t)dimension * sizeof(int));
    if (!local_best.tour) {
        free(global_best.tour);
        free_population_buffers(&buffers);
        return -1;
    }
    local_best.length = DBL_MAX;
    local_best.fitness = 0.0;

    for (int i = 0; i < buffers.size; ++i) {
        if (buffers.population[i].length < local_best.length) {
            ga_copy_individual(&buffers.population[i], &local_best, dimension);
        }
    }

    struct {
        double length;
        int rank;
    } local_info = {local_best.length, world_rank};
    struct {
        double length;
        int rank;
    } global_info = {0.0, 0};
    MPI_Allreduce(&local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    if (world_rank == global_info.rank) {
        ga_copy_individual(&local_best, &global_best, dimension);
    }
    MPI_Bcast(global_best.tour, dimension, MPI_INT, global_info.rank, comm);
    MPI_Bcast(&global_best.length, 1, MPI_DOUBLE, global_info.rank, comm);
    global_best.fitness = 1.0 / (global_best.length + 1e-9);

    double block_start = MPI_Wtime();

    for (int generation = 0; generation < params->generations; ++generation) {
        ga_copy_individual(&global_best, &buffers.next_population[0], dimension);
        int fill_index = 1;
        while (fill_index < buffers.size) {
            const int parent_a_idx = ga_tournament_select(buffers.population, buffers.size, params->tournament_size);
            const int parent_b_idx = ga_tournament_select(buffers.population, buffers.size, params->tournament_size);
            const Individual *parent_a = &buffers.population[parent_a_idx];
            const Individual *parent_b = &buffers.population[parent_b_idx];

            Individual *child1 = &buffers.next_population[fill_index];
            Individual *child2 = NULL;
            if (fill_index + 1 < buffers.size) {
                child2 = &buffers.next_population[fill_index + 1];
            }

            if (ga_random_unit() < params->crossover_rate) {
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

        Individual *tmp_pop = buffers.population;
        buffers.population = buffers.next_population;
        buffers.next_population = tmp_pop;

        local_best.length = DBL_MAX;
        for (int i = 0; i < buffers.size; ++i) {
            if (buffers.population[i].length < local_best.length) {
                ga_copy_individual(&buffers.population[i], &local_best, dimension);
            }
        }

        local_info.length = local_best.length;
        local_info.rank = world_rank;
        global_info.length = 0.0;
        global_info.rank = 0;
        MPI_Allreduce(&local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

        if (world_rank == global_info.rank) {
            ga_copy_individual(&local_best, &global_best, dimension);
        }
        MPI_Bcast(global_best.tour, dimension, MPI_INT, global_info.rank, comm);
        MPI_Bcast(&global_best.length, 1, MPI_DOUBLE, global_info.rank, comm);
        global_best.fitness = 1.0 / (global_best.length + 1e-9);

        if ((generation + 1) % MIGRATION_INTERVAL == 0) {
            int local_elite_count = (int)(0.05 * buffers.size);
            if (local_elite_count < 1) {
                local_elite_count = 1;
            }
            int *elite_indices = malloc((size_t)local_elite_count * sizeof(int));
            if (!elite_indices) {
                free(global_best.tour);
                free(local_best.tour);
                free_population_buffers(&buffers);
                return -1;
            }
            select_best_indices(buffers.population, buffers.size, local_elite_count, elite_indices);

            int send_count = local_elite_count * dimension;
            int *send_buffer = malloc((size_t)send_count * sizeof(int));
            double *send_lengths = malloc((size_t)local_elite_count * sizeof(double));
            if (!send_buffer || !send_lengths) {
                free(elite_indices);
                free(send_buffer);
                free(send_lengths);
                free(global_best.tour);
                free(local_best.tour);
                free_population_buffers(&buffers);
                return -1;
            }
            for (int i = 0; i < local_elite_count; ++i) {
                const int idx = elite_indices[i] >= 0 ? elite_indices[i] : 0;
                memcpy(send_buffer + (size_t)i * (size_t)dimension, buffers.population[idx].tour,
                       (size_t)dimension * sizeof(int));
                send_lengths[i] = buffers.population[idx].length;
            }
            free(elite_indices);

            int *elite_counts = malloc((size_t)world_size * sizeof(int));
            int *recv_counts_int = malloc((size_t)world_size * sizeof(int));
            int *displs_int = malloc((size_t)world_size * sizeof(int));
            int *length_displs = malloc((size_t)world_size * sizeof(int));
            int *all_tours = NULL;
            double *all_lengths = NULL;
            int total_elites = 0;

            if (!elite_counts || !recv_counts_int || !displs_int || !length_displs) {
                free(elite_counts);
                free(recv_counts_int);
                free(displs_int);
                free(length_displs);
                free(send_buffer);
                free(send_lengths);
                free(global_best.tour);
                free(local_best.tour);
                free_population_buffers(&buffers);
                return -1;
            }

            for (int i = 0; i < world_size; ++i) {
                elite_counts[i] = 0;
                recv_counts_int[i] = 0;
                displs_int[i] = 0;
                length_displs[i] = 0;
            }

            int send_elites_metric = local_elite_count;
            MPI_Gather(&send_elites_metric, 1, MPI_INT, elite_counts, 1, MPI_INT, 0, comm);

            if (world_rank == 0) {
                int offset_int = 0;
                int offset_len = 0;
                for (int i = 0; i < world_size; ++i) {
                    const int count = elite_counts[i];
                    recv_counts_int[i] = count * dimension;
                    displs_int[i] = offset_int;
                    length_displs[i] = offset_len;
                    offset_int += recv_counts_int[i];
                    offset_len += count;
                    total_elites += count;
                }
                if (offset_int > 0) {
                    all_tours = malloc((size_t)offset_int * sizeof(int));
                }
                if (total_elites > 0) {
                    all_lengths = malloc((size_t)total_elites * sizeof(double));
                }
            }

            MPI_Gatherv(send_buffer, send_elites_metric * dimension, MPI_INT, all_tours, recv_counts_int,
                        displs_int, MPI_INT, 0, comm);
            MPI_Gatherv(send_lengths, send_elites_metric, MPI_DOUBLE, all_lengths, elite_counts,
                        length_displs, MPI_DOUBLE, 0, comm);

            free(send_buffer);
            free(send_lengths);

            int global_elite_count = (int)(0.05 * params->population_size);
            if (global_elite_count < 1) {
                global_elite_count = 1;
            }
            int *global_elite_tours = NULL;
            double *global_elite_lengths = NULL;

            if (world_rank == 0) {
                if (total_elites < global_elite_count) {
                    global_elite_count = total_elites;
                }
                if (global_elite_count > 0) {
                    global_elite_tours =
                        malloc((size_t)global_elite_count * (size_t)dimension * sizeof(int));
                    global_elite_lengths = malloc((size_t)global_elite_count * sizeof(double));
                }
                int *best_indices = global_elite_count > 0
                                        ? malloc((size_t)global_elite_count * sizeof(int))
                                        : NULL;
                if ((global_elite_count > 0) &&
                    (!global_elite_tours || !global_elite_lengths || !best_indices)) {
                    free(global_elite_tours);
                    free(global_elite_lengths);
                    free(best_indices);
                    global_elite_tours = NULL;
                    global_elite_lengths = NULL;
                    global_elite_count = 0;
                } else if (global_elite_count > 0) {
                    for (int i = 0; i < global_elite_count; ++i) {
                        double best_length = DBL_MAX;
                        int best_idx = -1;
                        for (int j = 0; j < total_elites; ++j) {
                            int already_selected = 0;
                            for (int k = 0; k < i; ++k) {
                                if (best_indices[k] == j) {
                                    already_selected = 1;
                                    break;
                                }
                            }
                            if (already_selected) {
                                continue;
                            }
                            if (all_lengths[j] < best_length) {
                                best_length = all_lengths[j];
                                best_idx = j;
                            }
                        }
                        if (best_idx >= 0) {
                            best_indices[i] = best_idx;
                            memcpy(global_elite_tours + (size_t)i * (size_t)dimension,
                                   all_tours + (size_t)best_idx * (size_t)dimension, (size_t)dimension * sizeof(int));
                            global_elite_lengths[i] = all_lengths[best_idx];
                        }
                    }
                    free(best_indices);
                }
            }

            MPI_Bcast(&global_elite_count, 1, MPI_INT, 0, comm);
            if (global_elite_count > 0) {
                if (world_rank != 0) {
                    global_elite_tours = malloc((size_t)global_elite_count * (size_t)dimension * sizeof(int));
                    global_elite_lengths = malloc((size_t)global_elite_count * sizeof(double));
                }
                MPI_Bcast(global_elite_tours, global_elite_count * dimension, MPI_INT, 0, comm);
                MPI_Bcast(global_elite_lengths, global_elite_count, MPI_DOUBLE, 0, comm);
                integrate_elites(&buffers, instance, global_elite_tours, global_elite_lengths, global_elite_count,
                                  dimension);
                free(global_elite_tours);
                free(global_elite_lengths);
            }

            if (world_rank == 0) {
                free(all_tours);
                free(all_lengths);
            }
            free(elite_counts);
            free(recv_counts_int);
            free(displs_int);
            free(length_displs);

            double block_end = MPI_Wtime();
            double local_block_time = block_end - block_start;
            double *all_block_times = malloc((size_t)world_size * sizeof(double));
            MPI_Allgather(&local_block_time, 1, MPI_DOUBLE, all_block_times, 1, MPI_DOUBLE, comm);

            double total_speed = 0.0;
            double *speed = malloc((size_t)world_size * sizeof(double));
            for (int i = 0; i < world_size; ++i) {
                speed[i] = all_block_times[i] > 0.0 ? (1.0 / all_block_times[i]) : 1.0;
                total_speed += speed[i];
            }

            int *new_sizes = malloc((size_t)world_size * sizeof(int));
            int assigned = 0;
            for (int i = 0; i < world_size; ++i) {
                double fraction = speed[i] / total_speed;
                new_sizes[i] = (int)(params->population_size * fraction);
                if (new_sizes[i] < 1) {
                    new_sizes[i] = 1;
                }
                assigned += new_sizes[i];
            }
            while (assigned < params->population_size) {
                int best_idx = 0;
                double best_speed = speed[0];
                for (int i = 1; i < world_size; ++i) {
                    if (speed[i] > best_speed) {
                        best_speed = speed[i];
                        best_idx = i;
                    }
                }
                new_sizes[best_idx] += 1;
                assigned += 1;
            }
            while (assigned > params->population_size) {
                int worst_idx = 0;
                double worst_speed = speed[0];
                for (int i = 1; i < world_size; ++i) {
                    if (speed[i] < worst_speed && new_sizes[i] > 1) {
                        worst_speed = speed[i];
                        worst_idx = i;
                    }
                }
                if (new_sizes[worst_idx] > 1) {
                    new_sizes[worst_idx] -= 1;
                    assigned -= 1;
                } else {
                    break;
                }
            }

            int new_size = local_size;
            if (new_sizes) {
                new_size = new_sizes[world_rank];
            }

            block_start = MPI_Wtime();

            MPI_Bcast(new_sizes, world_size, MPI_INT, 0, comm);
            new_size = new_sizes[world_rank];

            free(all_block_times);
            free(speed);
            free(new_sizes);

            if (resize_population_buffers(&buffers, new_size, dimension, instance) != 0) {
                free(global_best.tour);
                free(local_best.tour);
                free_population_buffers(&buffers);
                return -1;
            }
            local_size = new_size;
        }

        if (log_stream && world_rank == 0 &&
            (generation % log_interval == 0 || generation == params->generations - 1)) {
            fprintf(log_stream, "[MPI] Generation %d global best %.3f\n", generation, global_best.length);
        }
    }

    if (best_out && best_out->tour) {
        ga_copy_individual(&global_best, best_out, dimension);
    }

    free(global_best.tour);
    free(local_best.tour);
    free_population_buffers(&buffers);
    return 0;
}
