#include "cli.h"
#include "ga_serial.h"
#include "tsp.h"

#include <stdio.h>
#include <stdlib.h>

static void print_tour(const int *tour, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        printf("%d ", tour[i] + 1);
    }
    printf("%d\n", tour[0] + 1);
}

int main(int argc, char **argv) {
    CLIOptions options;
    const int parse_status = cli_parse_args(argc, argv, &options, stdout);
    if (parse_status != 0) {
        return parse_status > 0 ? 0 : 1;
    }

    TSPInstance instance = {0};
    if (tsp_load(options.instance_path, &instance) != 0) {
        return 1;
    }

    Individual best = {0};
    best.tour = malloc((size_t)instance.dimension * sizeof(int));
    if (!best.tour) {
        fprintf(stderr, "Failed to allocate buffer for best tour\n");
        tsp_free(&instance);
        return 1;
    }

    printf("Running serial GA on %s (%d cities)\n", instance.name, instance.dimension);

    const int status =
        ga_run_serial(&instance, &options.params, &best, options.log_interval, stdout);
    if (status != 0) {
        fprintf(stderr, "Serial GA execution failed\n");
        free(best.tour);
        tsp_free(&instance);
        return 1;
    }

    printf("Best tour length: %.3f\n", best.length);
    printf("Tour sequence (1-based): ");
    print_tour(best.tour, instance.dimension);

    free(best.tour);
    tsp_free(&instance);
    return 0;
}
