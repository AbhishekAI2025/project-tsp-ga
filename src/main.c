#include "cli.h"
#include "parallel_ga.h"
#include "tsp.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

static void print_tour(const int *tour, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        printf("%d ", tour[i] + 1);
    }
    printf("%d\n", tour[0] + 1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CLIOptions options;
    int parse_status = 0;
    if (world_rank == 0) {
        parse_status = cli_parse_args(argc, argv, &options, stdout);
    }
    MPI_Bcast(&parse_status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (parse_status != 0) {
        MPI_Finalize();
        return parse_status > 0 ? 0 : 1;
    }
    if (world_rank != 0) {
        cli_set_defaults(&options);
    }
    MPI_Bcast(&options, sizeof(CLIOptions), MPI_BYTE, 0, MPI_COMM_WORLD);

    TSPInstance instance = {0};
    int local_load_ok = (tsp_load(options.instance_path, &instance) == 0) ? 1 : 0;
    int global_load_ok = 0;
    MPI_Allreduce(&local_load_ok, &global_load_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_load_ok) {
        if (world_rank == 0) {
            fprintf(stderr, "Failed to load TSP instance %s\n", options.instance_path);
        }
        if (local_load_ok) {
            tsp_free(&instance);
        }
        MPI_Finalize();
        return 1;
    }

    Individual best = {0};
    best.tour = malloc((size_t)instance.dimension * sizeof(int));
    if (!best.tour) {
        if (world_rank == 0) {
            fprintf(stderr, "Failed to allocate buffer for best tour\n");
        }
        tsp_free(&instance);
        MPI_Finalize();
        return 1;
    }

    if (world_rank == 0) {
        printf("Running parallel GA with %d MPI ranks on %s (%d cities)\n", world_size, instance.name,
               instance.dimension);
    }

    int status = ga_run_parallel(&instance, &options.params, &best, options.log_interval, MPI_COMM_WORLD,
                                 world_rank == 0 ? stdout : NULL);
    if (status != 0) {
        if (world_rank == 0) {
            fprintf(stderr, "Parallel GA execution failed\n");
        }
        free(best.tour);
        tsp_free(&instance);
        MPI_Finalize();
        return 1;
    }

    if (world_rank == 0) {
        printf("Best tour length: %.3f\n", best.length);
        printf("Tour sequence (1-based): ");
        print_tour(best.tour, instance.dimension);
    }

    free(best.tour);
    tsp_free(&instance);
    MPI_Finalize();
    return 0;
}
