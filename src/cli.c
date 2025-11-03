#include "cli.h"

#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

static double parse_double(const char *value, double fallback) {
    if (!value) {
        return fallback;
    }
    char *endptr = NULL;
    errno = 0;
    const double parsed = strtod(value, &endptr);
    if (errno != 0 || endptr == value) {
        return fallback;
    }
    return parsed;
}

static int parse_int(const char *value, int fallback) {
    if (!value) {
        return fallback;
    }
    char *endptr = NULL;
    errno = 0;
    const long parsed = strtol(value, &endptr, 10);
    if (errno != 0 || endptr == value) {
        return fallback;
    }
    return (int)parsed;
}

void cli_set_defaults(CLIOptions *options) {
    memset(options, 0, sizeof(*options));
    options->params.population_size = 200;
    options->params.generations = 1000;
    options->params.crossover_rate = 0.8;
    options->params.mutation_rate = 0.05;
    options->params.tournament_size = 4;
    options->params.two_opt_iterations = 20;
    options->params.seed = 42;
    options->log_interval = 10;
}

void cli_print_usage(const char *program, FILE *output) {
    fprintf(output,
            "Usage: %s --instance <path> [options]\n"
            "Options:\n"
            "  -i, --instance PATH       Path to TSPLIB .tsp instance (required)\n"
            "  -p, --population N        Population size (default 200)\n"
            "  -g, --generations N       Number of generations (default 1000)\n"
            "  -c, --crossover RATE      Crossover rate in [0,1] (default 0.8)\n"
            "  -m, --mutation RATE       Mutation rate in [0,1] (default 0.05)\n"
            "  -t, --tournament K        Tournament size (default 4)\n"
            "      --two-opt N           Max 2-opt iterations per child (default 20)\n"
            "  -s, --seed N              Random seed (default 42)\n"
            "      --log-interval N      Generations between log messages (default 10)\n"
            "  -h, --help                Show this help message\n",
            program);
}

int cli_parse_args(int argc, char **argv, CLIOptions *options, FILE *output) {
    cli_set_defaults(options);

    static struct option long_options[] = {
        {"instance", required_argument, 0, 'i'},
        {"population", required_argument, 0, 'p'},
        {"generations", required_argument, 0, 'g'},
        {"crossover", required_argument, 0, 'c'},
        {"mutation", required_argument, 0, 'm'},
        {"tournament", required_argument, 0, 't'},
        {"two-opt", required_argument, 0, 1000},
        {"seed", required_argument, 0, 's'},
        {"log-interval", required_argument, 0, 1001},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};

    int option_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, "i:p:g:c:m:t:s:h", long_options, &option_index)) != -1) {
        switch (opt) {
        case 'i':
            strncpy(options->instance_path, optarg, sizeof(options->instance_path) - 1);
            break;
        case 'p':
            options->params.population_size = parse_int(optarg, options->params.population_size);
            break;
        case 'g':
            options->params.generations = parse_int(optarg, options->params.generations);
            break;
        case 'c':
            options->params.crossover_rate = parse_double(optarg, options->params.crossover_rate);
            break;
        case 'm':
            options->params.mutation_rate = parse_double(optarg, options->params.mutation_rate);
            break;
        case 't':
            options->params.tournament_size = parse_int(optarg, options->params.tournament_size);
            break;
        case 's':
            options->params.seed = (unsigned int)parse_int(optarg, (int)options->params.seed);
            break;
        case 1000:
            options->params.two_opt_iterations = parse_int(optarg, options->params.two_opt_iterations);
            break;
        case 1001:
            options->log_interval = parse_int(optarg, options->log_interval);
            break;
        case 'h':
            cli_print_usage(argv[0], output);
            return 1;
        default:
            cli_print_usage(argv[0], output);
            return -1;
        }
    }

    if (options->instance_path[0] == '\0') {
        fprintf(output, "Missing required --instance argument\n");
        cli_print_usage(argv[0], output);
        return -1;
    }

    if (options->params.population_size < 2) {
        options->params.population_size = 2;
    }
    if (options->params.generations < 1) {
        options->params.generations = 1;
    }
    if (options->params.tournament_size < 2) {
        options->params.tournament_size = 2;
    }
    if (options->params.two_opt_iterations < 0) {
        options->params.two_opt_iterations = 0;
    }
    if (options->params.crossover_rate < 0.0) {
        options->params.crossover_rate = 0.0;
    } else if (options->params.crossover_rate > 1.0) {
        options->params.crossover_rate = 1.0;
    }
    if (options->params.mutation_rate < 0.0) {
        options->params.mutation_rate = 0.0;
    } else if (options->params.mutation_rate > 1.0) {
        options->params.mutation_rate = 1.0;
    }
    if (options->log_interval < 1) {
        options->log_interval = 1;
    }

    return 0;
}
