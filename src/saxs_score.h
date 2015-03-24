/** @file
 * Contains functions related to scoring profiles.
 */
#ifndef _SAXS_SCORE_H_
#define _SAXS_SCORE_H_

#include "saxs_profile.h"

#include "saxs_common.h"

// default c1/c2 limits
#define _DEFAULT_MIN_C1  0.95
#define _DEFAULT_MAX_C1  1.12
#define _DEFAULT_MIN_C2 -2.0
#define _DEFAULT_MAX_C2  4.0

// Simple interface

/**
 * Returns chi score between an experimental profile and a computed profile.
 * The experimental profile
 */
double chi_score(saxs_profile *expt_profile, saxs_profile *comp_profile);

double chi_fit_score_atom_group(saxs_profile *expt_profile, mol_atom_group *ag,
                                form_factor_table *ff_table, struct prm *prm,
                                double *c1, double *c2,
                                double min_c1, double max_c1,
                                double min_c2, double max_c2);


// Advanced interface

/**
 * Returns the chi for the fitted profile. The value of `c1` and
 * `c2` are used as the initial conditions.
 * @param[in,out] c1
 * @param[in,out] c2
 */
double chi_fit_score(saxs_profile *expt_profile, saxs_profile *comp_profile,
                     double *c1, double *c2,
                     double min_c1, double max_c1,
                     double min_c2, double max_c2);

#endif /* _SAXS_SCORE_H_ */
