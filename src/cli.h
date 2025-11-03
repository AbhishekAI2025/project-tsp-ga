#ifndef CLI_H
#define CLI_H

#include "ga_common.h"

#include <stdio.h>

typedef struct {
    char instance_path[512];
    GAParams params;
    int log_interval;
} CLIOptions;

void cli_set_defaults(CLIOptions *options);
int cli_parse_args(int argc, char **argv, CLIOptions *options, FILE *output);
void cli_print_usage(const char *program, FILE *output);

#endif
