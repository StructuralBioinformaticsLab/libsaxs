/** @mainpage
    @section examples Examples
    The following is a simple example of reading a profile from a file, and
    then fitting the computed profile from a `mol_atom_group` onto the
    experimental profile.

        char *expt_profile_path;
        saxs_profile expt_profile;
        saxs_profile_from_path(&expt_profile, expt_profile_path);
        mol_atom_group *ag = read_pdb(ag_path, prm);
        double c1 = 1.0, c2 = 0.0;
        double chi_score;
        chi_score = chi_fit_score_atom_group(expt_profile, ag,
                                             &c1, &c2,
                                             _DEFAULT_MIN_C1, _DEFAULT_MAX_C1,
                                             _DEFAULT_MIN_C2, _DEFAULT_MAX_C2);


    This is a slightly more involved example, which you would want
    if you are going to be computing scores for many profiles.

        char *expt_profile_path;
        saxs_profile expt_profile;
        saxs_profile_from_path(&expt_profile, expt_profile_path);
        mol_atom_group *ag = read_pdb(ag_path, prm);

    In all cases, including this header file should be sufficient to have
    visibility of all SAXS related functions and data structures.
 */
#ifndef _SAXS_LIB_H_
#define _SAXS_LIB_H_

// Users will call functions from these files
#include "saxs_profile.h"
#include "saxs_score.h"
#include "form_factor_table.h"
#include "rgyration.h"

// Very rarely need to peek inside these
#include "saxs_atom.h"
#include "radial_distribution.h"
#include "saxs_utils.h"

#endif /* _SAXS_LIB_H_ */
