#include "tsp_utils.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void trim(char *line) {
    char *end = line + strlen(line);
    while (end > line && isspace((unsigned char)end[-1])) {
        --end;
    }
    *end = '\0';
}

static double euclidean(double x1, double y1, double x2, double y2) {
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy);
}

City *read_tsp_file(const char *filename, int *num_cities) {
    if (!filename || !num_cities) {
        return NULL;
    }

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open TSP file %s: %s\n", filename, strerror(errno));
        return NULL;
    }

    char line[256];
    int dimension = -1;
    int found_section = 0;

    while (fgets(line, sizeof(line), fp)) {
        trim(line);
        if (line[0] == '\0') {
            continue;
        }
        if (strncmp(line, "DIMENSION", 9) == 0) {
            const char *sep = strchr(line, ':');
            if (!sep) {
                fprintf(stderr, "Invalid DIMENSION line in %s\n", filename);
                fclose(fp);
                return NULL;
            }
            dimension = atoi(sep + 1);
        } else if (strcmp(line, "NODE_COORD_SECTION") == 0) {
            found_section = 1;
            break;
        }
    }

    if (dimension <= 0 || !found_section) {
        fprintf(stderr, "Invalid TSP file %s: missing DIMENSION or NODE_COORD_SECTION\n", filename);
        fclose(fp);
        return NULL;
    }

    City *cities = calloc((size_t)dimension, sizeof(City));
    if (!cities) {
        fprintf(stderr, "Out of memory allocating %d cities\n", dimension);
        fclose(fp);
        return NULL;
    }

    for (int i = 0; i < dimension; ++i) {
        int id = 0;
        double x = 0.0;
        double y = 0.0;
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "Unexpected EOF while reading coordinates in %s\n", filename);
            free(cities);
            fclose(fp);
            return NULL;
        }
        if (sscanf(line, "%d %lf %lf", &id, &x, &y) != 3) {
            fprintf(stderr, "Invalid coordinate line: %s\n", line);
            free(cities);
            fclose(fp);
            return NULL;
        }
        if (id <= 0 || id > dimension) {
            fprintf(stderr, "City index out of range (%d) in %s\n", id, filename);
            free(cities);
            fclose(fp);
            return NULL;
        }
        cities[id - 1].id = id;
        cities[id - 1].x = x;
        cities[id - 1].y = y;
    }

    fclose(fp);
    *num_cities = dimension;
    return cities;
}

double **compute_distance_matrix(const City *cities, int num_cities) {
    if (!cities || num_cities <= 0) {
        return NULL;
    }

    double **matrix = calloc((size_t)num_cities, sizeof(double *));
    double *storage = calloc((size_t)num_cities * (size_t)num_cities, sizeof(double));
    if (!matrix || !storage) {
        fprintf(stderr, "Failed to allocate distance matrix for %d cities\n", num_cities);
        free(matrix);
        free(storage);
        return NULL;
    }

    for (int i = 0; i < num_cities; ++i) {
        matrix[i] = storage + (size_t)i * (size_t)num_cities;
        for (int j = 0; j < num_cities; ++j) {
            matrix[i][j] = euclidean(cities[i].x, cities[i].y, cities[j].x, cities[j].y);
        }
    }

    return matrix;
}

void free_distance_matrix(double **matrix, int num_cities) {
    if (!matrix) {
        return;
    }
    if (num_cities > 0 && matrix[0]) {
        free(matrix[0]);
    }
    free(matrix);
}
