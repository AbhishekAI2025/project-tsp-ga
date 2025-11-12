#include "ga.h"
#include "parallel_ga.h"
#include "tsp_utils.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(void) {
    fprintf(stderr,
            "Usage: tsp_ga <instance.tsp> [options]\n"
            "Options:\n"
            "  --population <int>    Population size (default 256)\n"
            "  --generations <int>   Number of generations (default 500)\n"
            "  --crossover <double>  Crossover rate (default 0.9)\n"
            "  --mutation <double>   Mutation rate (default 0.05)\n"
            "  --tournament <int>    Tournament size (default 4)\n"
            "  --two-opt <int>       Two-opt max swaps per tour (default 2000)\n"
            "  --sync-interval <int> Generations between MPI best-tour syncs (default 10)\n"
            "  --seed <uint>         RNG seed (default 42)\n"
            "  --help                Show this message\n");
}

static int parse_arguments(int argc, char **argv, GAParams *params, char **instance_path, int rank) {
    if (argc < 2) {
        if (rank == 0) {
            print_usage();
        }
        return -1;
    }

    *instance_path = argv[1];
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            if (rank == 0) {
                print_usage();
            }
            return 1;
        } else if (strcmp(argv[i], "--population") == 0 && i + 1 < argc) {
            params->population_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--generations") == 0 && i + 1 < argc) {
            params->generations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--crossover") == 0 && i + 1 < argc) {
            params->crossover_rate = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mutation") == 0 && i + 1 < argc) {
            params->mutation_rate = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tournament") == 0 && i + 1 < argc) {
            params->tournament_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--two-opt") == 0 && i + 1 < argc) {
            params->two_opt_max_swaps = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--sync-interval") == 0 && i + 1 < argc) {
            params->sync_interval = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            params->seed = (unsigned int)strtoul(argv[++i], NULL, 10);
        } else {
            if (rank == 0) {
                fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
                print_usage();
            }
            return -1;
        }
    }

    if (params->population_size < 2 || params->generations <= 0 || params->tournament_size < 2 ||
        params->sync_interval <= 0) {
        if (rank == 0) {
            fprintf(stderr, "Invalid GA parameters.\n");
        }
        return -1;
    }
    return 0;
}

static void print_tour(const int *tour, int num_cities) {
    for (int i = 0; i < num_cities; ++i) {
        printf("%d ", tour[i] + 1);
    }
    printf("%d\n", tour[0] + 1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    GAParams params = {
        .population_size = 256,
        .generations = 500,
        .crossover_rate = 0.9,
        .mutation_rate = 0.05,
        .tournament_size = 4,
        .two_opt_max_swaps = 2000,
        .seed = 42u,
        .sync_interval = 10,
    };

    char *instance_path = NULL;
    int parse_result = 0;
    if (rank == 0) {
        parse_result = parse_arguments(argc, argv, &params, &instance_path, rank);
    }
    MPI_Bcast(&parse_result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (parse_result != 0) {
        MPI_Finalize();
        return parse_result > 0 ? 0 : 1;
    }

    int path_length = 0;
    if (rank == 0) {
        path_length = (int)strlen(instance_path) + 1;
    }
    MPI_Bcast(&path_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        instance_path = malloc((size_t)path_length);
        if (!instance_path) {
            fprintf(stderr, "Rank %d failed to allocate instance path buffer\n", rank);
            MPI_Finalize();
            return 1;
        }
    }
    MPI_Bcast(instance_path, path_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params, sizeof(GAParams), MPI_BYTE, 0, MPI_COMM_WORLD);

    City *cities = NULL;
    int num_cities = 0;
    int read_status = 0;
    if (rank == 0) {
        cities = read_tsp_file(instance_path, &num_cities);
        if (!cities) {
            read_status = -1;
        }
    }
    MPI_Bcast(&read_status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (read_status != 0) {
        if (rank == 0) {
            fprintf(stderr, "Failed to load TSP instance %s\n", instance_path);
        }
        if (rank != 0) {
            free(instance_path);
        }
        MPI_Finalize();
        return 1;
    }

    MPI_Bcast(&num_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        cities = calloc((size_t)num_cities, sizeof(City));
        if (!cities) {
            fprintf(stderr, "Rank %d failed to allocate cities buffer\n", rank);
            if (rank != 0) {
                free(instance_path);
            }
            MPI_Finalize();
            return 1;
        }
    }

    MPI_Bcast(cities, num_cities * (int)sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

    double **dist_matrix = compute_distance_matrix(cities, num_cities);
    if (!dist_matrix) {
        if (rank == 0) {
            fprintf(stderr, "Failed to compute distance matrix for %s\n", instance_path);
        }
        free(cities);
        if (rank != 0 && instance_path) {
            free(instance_path);
        }
        MPI_Finalize();
        return 1;
    }

    int *best_tour = NULL;
    double best_length = 0.0;
    if (rank == 0) {
        best_tour = malloc((size_t)num_cities * sizeof(int));
        if (!best_tour) {
            fprintf(stderr, "Failed to allocate best tour buffer\n");
            free_distance_matrix(dist_matrix, num_cities);
            free(cities);
            if (rank != 0) {
                free(instance_path);
            }
            MPI_Finalize();
            return 1;
        }
    }

    int status = 0;
    if (size == 1) {
        status = run_serial_ga(cities, num_cities, dist_matrix, &params, best_tour, &best_length);
    } else {
        status = run_parallel_ga(cities, num_cities, dist_matrix, &params, rank == 0 ? best_tour : NULL,
                                 rank == 0 ? &best_length : NULL, MPI_COMM_WORLD);
    }

    if (status != 0) {
        if (rank == 0) {
            fprintf(stderr, "Genetic algorithm execution failed.\n");
        }
        free(best_tour);
        free_distance_matrix(dist_matrix, num_cities);
        free(cities);
        if (rank != 0 && instance_path) {
            free(instance_path);
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        printf("Instance: %s\n", instance_path);
        printf("Cities: %d\n", num_cities);
        printf("Generations: %d | Population: %d | Mutation: %.3f | Crossover: %.3f\n", params.generations,
               params.population_size, params.mutation_rate, params.crossover_rate);
        printf("Best tour length: %.3f\n", best_length);
        printf("Tour (1-based): ");
        print_tour(best_tour, num_cities);
    }

    free(best_tour);
    free_distance_matrix(dist_matrix, num_cities);
    free(cities);
    if (rank != 0 && instance_path) {
        free(instance_path);
    }

    MPI_Finalize();
    return 0;
}
