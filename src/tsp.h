#ifndef TSP_H
#define TSP_H

#include <stddef.h>

typedef struct {
    char name[128];
    int dimension;
    double (*coords)[2];
    double *distance_matrix;
} TSPInstance;

int tsp_load(const char *path, TSPInstance *instance);
void tsp_free(TSPInstance *instance);
double tsp_distance(const TSPInstance *instance, int city_a, int city_b);
double tsp_tour_length(const TSPInstance *instance, const int *tour);

#endif
