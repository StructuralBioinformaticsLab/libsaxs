#include "saxs_lib.h"

#include "mol.0.0.6/pdb.h"

#include <getopt.h>

void print_usage(char *exe_name)
{
  fprintf(stderr,
          "Usage: %s MAPPING_PRM ATOMPRM SAXS_PROFILE REC [LIG]\n",
          basename(exe_name));
}

int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    print_usage(argv[0]);
    exit(1);
  }

  char *mapping_file_path = argv[1];
  char *prms_file_path    = argv[2];
  char *profile_path      = argv[3];
  char *rec_file_path     = argv[4];
  char *lig_file_path = NULL;
  if (argc > 5)
    lig_file_path         = argv[5];

  form_factor_table *ff_table = default_ff_table(mapping_file_path);
  struct prm *prm = read_prm(prms_file_path, "0.0.6");

  mol_atom_group *pa = read_pdb(rec_file_path, prm);
  pa = remove_hydrogens(pa, prm);
  if (argc > 5) {
    mol_atom_group *lig = read_pdb(lig_file_path, prm);
    lig = remove_hydrogens(lig, prm);
    pa = join_2atomgrps(pa, lig);
  }

  saxs_profile *expt_profile = calloc(1, sizeof(saxs_profile));
  saxs_profile_from_path(expt_profile, profile_path);

  saxs_profile *comp_profile = calloc(1, sizeof(saxs_profile));
  saxs_profile_clone(comp_profile, expt_profile);
  saxs_profile_create_partials(comp_profile, true);
  comp_profile->average_radius = mol_atom_group_average_radius(pa, prm);

  double c1 = 1.0;
  double c2 = 0.0;

  float *saxs_sa = malloc(pa->natoms * sizeof(float));
  double *vacuum_ff = malloc(pa->natoms * sizeof(double));
  double *dummy_ff = malloc(pa->natoms * sizeof(double));
  get_vacuum_ff(vacuum_ff, pa, ff_table, prm);
  get_dummy_ff(dummy_ff, pa, ff_table, prm);
  faccs(saxs_sa, pa, prm, 1.4);

  compute_profile_partials(comp_profile, pa, saxs_sa,
                           vacuum_ff, dummy_ff,
                           max_dist(pa), ff_table);
  saxs_profile_sum_partials(comp_profile, 1.0, 0.0);

  double default_chi = chi_score(expt_profile, comp_profile);

  double fitted_chi = chi_fit_score(expt_profile, comp_profile,
                                    &c1, &c2,
                                    _DEFAULT_MIN_C1, _DEFAULT_MAX_C1,
                                    _DEFAULT_MIN_C2, _DEFAULT_MAX_C2);

  printf("%f %f %f %f\n", fitted_chi, c1, c2, default_chi);

  return 0;
}
