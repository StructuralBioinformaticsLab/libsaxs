/** @file
 * Contains all functions related to computing SAXS profiles.
 */
#ifndef _SAXS_PROFILE_H_
#define _SAXS_PROFILE_H_

#include "form_factor_table.h"
#include "saxs_atom.h"

#include "saxs_common.h"
#include "saxs_config.h"

/**
 */
typedef struct saxs_profile_t
{
  int nsamples;
  double *q;
  double *In;
  double *err;

  double rgyration;
  double average_radius;

  // Partials
  int npartials;

  double *vac_vac;
  double *vac_dum;
  double *dum_dum;
  double *vac_h2o;
  double *dum_h2o;
  double *h2o_h2o;
} saxs_profile;


// Initialization

/**
 * Initializes a new profile. Profile will be left untouched if `q_max`
 * is less than `q_min`.
 *
 * @param profile Profile to be initialized.
 */
void saxs_profile_create(saxs_profile *profile,
                         double q_min, double q_max, double q_step);

/**
 * Initializes a profile from a given filename. Will fail if path of
 * filename is not readable. The file is expected to contain 2 or 3
 * columns of numerical data, with the columns being `q`, `In`, and
 * optionally `err`. Blank lines and lines starting with '#' are
 * ignored.
 */
void saxs_profile_from_path(saxs_profile *profile, char *filename);

/**
 * Initializes a profile from a file pointer. Will fail if `file`
 * is in a bad state. Format should be same as described in
 * saxs_profile_from_path().
 */
void saxs_profile_from_file(saxs_profile *profile, FILE *file);

/**
 * Clones a profile by copying the `nsamples` and `q`. Mainly
 * used to initialize a computed profile from a experimental
 * profile.
 */
void saxs_profile_clone(saxs_profile *dst, saxs_profile *src);

/**
 *
 */
void saxs_profile_create_partials(saxs_profile *profile, bool fit_water);

// Output

/**
 *
 */
void saxs_profile_fprint(FILE *fp, saxs_profile *profile);

/**
 *
 */
void saxs_profile_print(saxs_profile *profile);


// Cleanup

/**
 *
 */
void saxs_profile_destroy(saxs_profile *profile);

/**
 *
 */
void saxs_profile_reset_In(saxs_profile *profile);

/**
 *
 */
void saxs_profile_reset_partials(saxs_profile *profile);


// Actions

/**
 *
 */
void saxs_profile_scale(saxs_profile *profile, double scale);

/**
 * Sums the intensities from `p1` and `p2` and places the result in `dst`.
 * Assumes the q-values all match up, and weights the intensities by
 * `w1` and `w2` respectively.
 */
void saxs_profile_add(saxs_profile *dst,
                      saxs_profile *p1, saxs_profile *p2,
                      double w1, double w2);

/**
 * Computes the profile using the partials method, ignoring the hydration layer.
 * This is equivalent to calling compute_profile_partials() and then
 * saxs_profile_sum_partials() with `c1 = 1.0` and `c2 = 0.0`.
 */
void compute_profile(
  saxs_profile *profile,
  mol_atom_group *structure,
  double *vacuum_ff, double *dummy_ff,
  double max_d, form_factor_table *ff_table);

/**
 * Computes the partial profiles and stores them in the profile.
 * @param[in] max_d The expected maximum distance to be encountered between
 *                  atoms in `structure`
 */
void compute_profile_partials(
  saxs_profile *profile,
  mol_atom_group *structure,
  float *saxs_sa, double *vacuum_ff, double *dummy_ff,
  double max_d, form_factor_table *ff_table);

/**
 * Adds the partials for a profile up and updates the profile.
 */
void saxs_profile_sum_partials(saxs_profile *profile, double c1, double c2);

/**
 * Computes only the parts of the profile needed, assuming the internal
 * profile for `pa` and `pb` have been computed.
 * @deprecated Needs to be updated
 */
void compute_profile_reclig(
  saxs_profile *profile,
  mol_atom_group *pa, mol_atom_group *pb,
  saxs_profile *pa_profile, saxs_profile *pb_profile,
  double max_dist, form_factor_table *ff_table, struct prm* prms);

/**
 * Compute profile using the Debye formula directly. Slower but more likely
 * to give correct results.
 */
void compute_profile_debye(
  saxs_profile *profile,
  mol_atom_group *structure,
  float *saxs_sa, double *vacuum_ff, double *dummy_ff,
  double c1, double c2, struct prm* prm, double water_ff);

// Internal

void saxs_partials_create(saxs_profile *profile,
                          bool fit_water);

#endif /* _SAXS_PROFILE_H_ */
