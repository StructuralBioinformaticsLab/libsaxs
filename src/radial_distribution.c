#include "radial_distribution.h"

#include "saxs_utils.h"

void radial_distribution_create(
  radial_distribution *dist, double bin_size, double max_dist)
{
  dist->bin_size = bin_size;
  dist->one_over_bin_size = 1 / bin_size;
  dist->nbins = (int) (max_dist * dist->one_over_bin_size) + 1;
  dist->values = calloc(dist->nbins, sizeof(double));
}

void radial_distribution_reset(radial_distribution *dist)
{
  memset(dist->values, 0, dist->nbins * sizeof(double));
}

void radial_distribution_destroy(radial_distribution *dist)
{
  if (dist != NULL && dist->values != NULL)
    free(dist->values);
}

void radial_distributions_to_partials(
  saxs_profile *profile,
  int ndists,
  radial_distribution *r_distributions)
{
  int nbins = r_distributions[0].nbins;
  double delta_x = r_distributions[0].bin_size;
  double q, qd, x;
  double scaling_factor;

  for (int iq = 0; iq < profile->nsamples; ++iq)
  {
    q = profile->q[iq];
    for (int r = 0; r < nbins; ++r)
    {
      qd = r * delta_x * q; // r * delta_x = dist
      x = sinc(qd);
      if (r_distributions[0].values[r] > 0.0)
      {
        profile->vac_vac[iq] += r_distributions[0].values[r] * x;
        profile->dum_dum[iq] += r_distributions[1].values[r] * x;
        profile->vac_dum[iq] += r_distributions[2].values[r] * x;

        if (ndists == 6)
        {
          profile->h2o_h2o[iq]+= r_distributions[3].values[r] * x;
          profile->vac_h2o[iq] += r_distributions[4].values[r] * x;
          profile->dum_h2o[iq] += r_distributions[5].values[r] * x;
        }
      }
    }

    scaling_factor = exp(-modulation_function_parameter * profile->q[iq] * profile->q[iq]);
    profile->vac_vac[iq] *= scaling_factor;
    profile->vac_dum[iq] *= scaling_factor;
    profile->dum_dum[iq] *= scaling_factor;

    if (ndists == 6)
    {
      profile->vac_h2o[iq] *= scaling_factor;
      profile->dum_h2o[iq] *= scaling_factor;
      profile->h2o_h2o[iq] *= scaling_factor;
    }
  }
}
