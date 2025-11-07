#ifndef TSP_UTILS_H
#define TSP_UTILS_H

typedef struct {
    int id;
    double x;
    double y;
} City;

City *read_tsp_file(const char *filename, int *num_cities);
double **compute_distance_matrix(const City *cities, int num_cities);
void free_distance_matrix(double **matrix, int num_cities);

#endif
