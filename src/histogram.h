#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "saxs_config.h"
#include "saxs_common.h"

#include "utils.h"

typedef enum histogram_type_t
{
  INT_HIST, FLOAT_HIST
} histogram_type;

typedef struct histogram_t
{
  int num_bins;
  histogram_type type;
  union
  {
    int *int_bin;
    double *float_bin;
  } data;
  void *delta_x;
} histogram;

histogram *histogram_(int nbins, void *delta_x, histogram_type ht);
histogram *create_histogram(int size, histogram_type ht);
void destroy_histogram(histogram *hist);

double delta_x_d(histogram *hist);
int delta_x_i(histogram *hist);

void print_histogram(histogram *hist);
void fprint_histogram(FILE *fp, histogram *hist);

#endif /* _HISTOGRAM_H_ */
