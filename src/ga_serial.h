#ifndef GA_SERIAL_H
#define GA_SERIAL_H

#include "ga_common.h"
#include <stdio.h>

int ga_run_serial(const TSPInstance *instance, const GAParams *params, Individual *best_out, int log_interval,
                  FILE *log_stream);

#endif
