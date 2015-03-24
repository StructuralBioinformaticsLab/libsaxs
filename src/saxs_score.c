#include "saxs_score.h"

#include "saxs_utils.h"

#include <assert.h>


double chi_score(saxs_profile *expt_profile, saxs_profile *comp_profile)
{
  double chi = 0.0;
  double err;

  /* Find c that minimizes chi by finding where derivative is zero*/
  double numer = 0.0, denom = 0.0;
  for (int i = 0; i < expt_profile->nsamples; ++i)
  {
    if (expt_profile->err != NULL)
    {
      err = expt_profile->err[i];
    }
    else
    {
      err = 0.05 * expt_profile->In[i];
    }

    if (!( err == 0 || isnan(err) || isinf(err) ))
    {

      numer += (expt_profile->In[i] * comp_profile->In[i])/(err * err);
      denom += (comp_profile->In[i] * comp_profile->In[i])/(err * err);
      assert(!isinf(numer));
      assert(!isinf(denom));
    }
  }
  double c = numer / denom;

  /* Compute chi */
  for (int i = 0; i < expt_profile->nsamples; ++i)
  {
    if (expt_profile->err != NULL)
    {
      err = expt_profile->err[i];
    }
    else
    {
      err = 0.05 * expt_profile->In[i];
    }

    if (! (err == 0 || isnan(err) || isinf(err)))
      chi += pow((expt_profile->In[i] - c * comp_profile->In[i])/err, 2);
  }

  chi /= expt_profile->nsamples;
  chi = sqrt(chi);

  return chi;
}

static double chi_fit_score_helper(
  saxs_profile *expt_profile, saxs_profile *comp_profile,
  double *c1, double *c2, double initial_chi,
  double min_c1, double max_c1,
  double min_c2, double max_c2,
  int depth)
{
  static int max_depth = 3;
  if (depth >= max_depth)
    return initial_chi;

  int nc1cells = 10, nc2cells = 10;
  double delta_c1 = (max_c1 - min_c1) / nc1cells;
  double delta_c2 = (max_c2 - min_c2) / nc2cells;

  double cur_c1, cur_c2;
  double best_c1 = *c1, best_c2 = *c2;
  double best_chi = initial_chi, cur_chi;

  for (int i = 0; i < nc1cells; ++i)
  {
    cur_c1 = min_c1 + i * delta_c1;

    for (int j = 0; j < nc2cells; ++j)
    {
      cur_c2 = min_c2 + j * delta_c2;

      saxs_profile_sum_partials(comp_profile, cur_c1, cur_c2);
      cur_chi = chi_score(expt_profile, comp_profile);

      if (cur_chi < best_chi)
      {
        best_c1 = cur_c1;
        best_c2 = cur_c2;
        best_chi = cur_chi;
      }

      *c1 = best_c1;
      *c2 = best_c2;
    }
  }

  if ( (initial_chi - best_chi) > 0.0001)
  {
    min_c1 = max(best_c1 - delta_c1, min_c1);
    max_c1 = min(best_c1 + delta_c1, max_c1);
    min_c2 = max(best_c2 - delta_c2, min_c2);
    max_c2 = min(best_c2 + delta_c2, max_c2);

    return chi_fit_score_helper(expt_profile, comp_profile,
                                c1, c2, best_chi,
                                min_c1, max_c1,
                                min_c2, max_c2,
                                depth + 1);
  }

  saxs_profile_sum_partials(comp_profile, best_c1, best_c2);
  return best_chi;
}

double chi_fit_score_atom_group(saxs_profile *expt_profile, mol_atom_group *ag,
                                form_factor_table *ff_table, struct prm *prm,
                                double *c1, double *c2,
                                double min_c1, double max_c1,
                                double min_c2, double max_c2)
{
  saxs_profile comp_profile;
  saxs_profile_clone(&comp_profile, expt_profile);
  saxs_profile_create_partials(&comp_profile, true);
  comp_profile.average_radius = mol_atom_group_average_radius(ag, prm);

  float *saxs_sa;
  double *vacuum_ff;
  double *dummy_ff;
  double max_d = max_dist(ag) + 2;

  faccs(saxs_sa, ag, prm, 1.4); // Assumes msur_k = 1.0
  get_dummy_ff(dummy_ff, ag, ff_table, prm);
  get_vacuum_ff(vacuum_ff, ag, ff_table, prm);

  compute_profile_partials(&comp_profile,
                           ag,
                           saxs_sa, vacuum_ff, dummy_ff,
                           max_d, ff_table);

  return chi_fit_score(expt_profile, &comp_profile,
                       c1, c2,
                       min_c1, max_c1,
                       min_c2, max_c2);
}

double chi_fit_score(
  saxs_profile *expt_profile, saxs_profile *comp_profile,
  double *c1, double *c2,
  double min_c1, double max_c1,
  double min_c2, double max_c2)
{
  saxs_profile_sum_partials(comp_profile, *c1, *c2);
  double initial_chi = chi_score(expt_profile, comp_profile);

  return chi_fit_score_helper(expt_profile, comp_profile,
                              c1, c2, initial_chi,
                              min_c1, max_c1,
                              min_c2, max_c2,
                              0);
}

