#ifndef _RADIAL_DISTRIBUTION_H_
#define _RADIAL_DISTRIBUTION_H_

#include "saxs_profile.h"
#include "saxs_common.h"

typedef struct radial_distribution_t
{
  double bin_size;

  int nbins;
  double *values;

  // cached values
  double one_over_bin_size;
} radial_distribution;


void radial_distribution_create(radial_distribution *dist,
                                double bin_size, double max_dist);
void radial_distribution_reset(radial_distribution *dist);
void radial_distribution_destroy(radial_distribution *dist);

void radial_distributions_to_partials(
  saxs_profile *profile,
  int ndists,
  radial_distribution *dists);

static inline void add2distribution(radial_distribution *dist,
                                    double distance, double val)
{
  int bin = (int) (dist->one_over_bin_size * distance);
  dist->values[bin] += val;
}

#endif /* _RADIAL_DISTRIBUTION_H_ */
