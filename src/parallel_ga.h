#ifndef PARALLEL_GA_H
#define PARALLEL_GA_H

#include "ga_common.h"

#include <mpi.h>
#include <stdio.h>

int ga_run_parallel(const TSPInstance *instance, const GAParams *params, Individual *best_out, int log_interval,
                    MPI_Comm comm, FILE *log_stream);

#endif
