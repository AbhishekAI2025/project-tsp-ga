#include "tsp.h"

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

static double compute_distance(double x1, double y1, double x2, double y2) {
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy);
}

int tsp_load(const char *path, TSPInstance *instance) {
    memset(instance, 0, sizeof(*instance));

    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Failed to open TSP file %s: %s\n", path, strerror(errno));
        return -1;
    }

    char line[256];
    int dimension = -1;
    int in_coord_section = 0;
    char name_buf[128] = {0};

    while (fgets(line, sizeof(line), f)) {
        trim(line);
        if (line[0] == '\0') {
            continue;
        }
        if (strncmp(line, "NAME", 4) == 0) {
            const char *sep = strchr(line, ':');
            if (sep) {
                while (*sep == ':' || isspace((unsigned char)*sep)) {
                    ++sep;
                }
                strncpy(name_buf, sep, sizeof(name_buf) - 1);
            }
        } else if (strncmp(line, "DIMENSION", 9) == 0) {
            const char *sep = strchr(line, ':');
            if (!sep) {
                fprintf(stderr, "Invalid DIMENSION line in %s\n", path);
                fclose(f);
                return -1;
            }
            dimension = atoi(sep + 1);
        } else if (strcmp(line, "NODE_COORD_SECTION") == 0) {
            in_coord_section = 1;
            break;
        }
    }

    if (dimension <= 0 || !in_coord_section) {
        fprintf(stderr, "Invalid TSP file %s: missing DIMENSION or NODE_COORD_SECTION\n", path);
        fclose(f);
        return -1;
    }

    double (*coords)[2] = calloc((size_t)dimension, sizeof(*coords));
    if (!coords) {
        fprintf(stderr, "Out of memory allocating coordinates\n");
        fclose(f);
        return -1;
    }

    for (int i = 0; i < dimension; ++i) {
        int idx = 0;
        double x = 0.0, y = 0.0;
        if (!fgets(line, sizeof(line), f)) {
            fprintf(stderr, "Unexpected EOF reading coordinates in %s\n", path);
            free(coords);
            fclose(f);
            return -1;
        }
        if (sscanf(line, "%d %lf %lf", &idx, &x, &y) != 3) {
            fprintf(stderr, "Invalid coordinate line: %s\n", line);
            free(coords);
            fclose(f);
            return -1;
        }
        if (idx < 1 || idx > dimension) {
            fprintf(stderr, "Node index out of range: %d\n", idx);
            free(coords);
            fclose(f);
            return -1;
        }
        coords[idx - 1][0] = x;
        coords[idx - 1][1] = y;
    }

    fclose(f);

    double *distance_matrix = calloc((size_t)dimension * (size_t)dimension, sizeof(double));
    if (!distance_matrix) {
        fprintf(stderr, "Out of memory allocating distance matrix\n");
        free(coords);
        return -1;
    }

    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            distance_matrix[(size_t)i * (size_t)dimension + (size_t)j] =
                compute_distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
        }
    }

    strncpy(instance->name, name_buf[0] ? name_buf : "unknown", sizeof(instance->name) - 1);
    instance->dimension = dimension;
    instance->coords = coords;
    instance->distance_matrix = distance_matrix;
    return 0;
}

void tsp_free(TSPInstance *instance) {
    if (!instance) {
        return;
    }
    free(instance->coords);
    instance->coords = NULL;
    free(instance->distance_matrix);
    instance->distance_matrix = NULL;
    instance->dimension = 0;
    instance->name[0] = '\0';
}

double tsp_distance(const TSPInstance *instance, int city_a, int city_b) {
    const int n = instance->dimension;
    return instance->distance_matrix[(size_t)city_a * (size_t)n + (size_t)city_b];
}

double tsp_tour_length(const TSPInstance *instance, const int *tour) {
    const int n = instance->dimension;
    double length = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        length += tsp_distance(instance, tour[i], tour[i + 1]);
    }
    length += tsp_distance(instance, tour[n - 1], tour[0]);
    return length;
}
