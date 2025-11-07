#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

void seed_random(unsigned int seed);
double rand_double(void);
int rand_int(int min_inclusive, int max_exclusive);
void shuffle_array(int *array, int length);

#endif
