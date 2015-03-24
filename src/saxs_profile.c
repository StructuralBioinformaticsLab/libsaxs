#include "saxs_profile.h"

#include "radial_distribution.h"
#include "saxs_utils.h"

#include "mol.0.0.6/myhelpers.h"

// Initialization
void saxs_profile_create(saxs_profile *profile,
                         double q_min, double q_max, double q_step)
{
  if (q_min > q_max)
    return;

  profile->nsamples = (q_max - q_min) / q_step;
  profile->q = malloc(profile->nsamples * sizeof(double));
  profile->In = malloc(profile->nsamples * sizeof(double));
  profile->err = NULL;

  profile->q[0] = q_min;
  for (int i = 1; i < profile->nsamples; ++i)
  {
    profile->q[i] = profile->q[i-1] + q_step;
  }

  profile->rgyration = 0.0;
  profile->npartials = 0;
  profile->vac_vac = NULL;
  profile->vac_dum = NULL;
  profile->dum_dum = NULL;
  profile->vac_h2o = NULL;
  profile->dum_h2o = NULL;
  profile->h2o_h2o = NULL;
}

void saxs_profile_from_path(saxs_profile *profile, char *filename)
{
  FILE *fp = fopen(filename, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "no such file %s\n", filename);
    return;
  }

  saxs_profile_from_file(profile, fp);
  fclose(fp);
}

void saxs_profile_from_file(saxs_profile *profile, FILE *fp)
{
  int arr_size = 16;

  profile->nsamples = 0;
  profile->q = malloc(arr_size * sizeof(double));
  profile->In = malloc(arr_size * sizeof(double));
  profile->err = malloc(arr_size * sizeof(double));
  profile->npartials = 0;
  profile->vac_vac = NULL;
  profile->vac_dum = NULL;
  profile->dum_dum = NULL;
  profile->vac_h2o = NULL;
  profile->dum_h2o = NULL;
  profile->h2o_h2o = NULL;

  char *line = NULL;
  size_t len;
  double q, val, err;
  int tmp;
  int matched_fields;

  while (getline(&line, &len, fp) != -1)
  {
    if (line[0] == '#' || iswhiteline(line))
      continue;

    if (profile->nsamples >= arr_size)
    {
      arr_size *= 2;  // Double size of arrays
      profile->q = realloc(profile->q, arr_size * sizeof(double));
      profile->In = realloc(profile->In, arr_size * sizeof(double));
      profile->err = realloc(profile->err, arr_size * sizeof(double));
    }

    matched_fields = sscanf(line, "%lf %lf %lf", &q, &val, &err);
    if (matched_fields < 2)
    {
      fprintf(stderr, "Error parsing saxs profile file\n");
      exit(1); // error parsing
    }
    else
    {
      profile->q[profile->nsamples] = q;
      profile->In[profile->nsamples] = val;
    }
    if (matched_fields == 3)
      profile->err[profile->nsamples] = err;
    else
      // Fill in error
      profile->err[profile->nsamples] = val * 0.05;

    ++(profile->nsamples);
  }

  profile->q = realloc(profile->q, profile->nsamples * sizeof(double));
  profile->In = realloc(profile->In, profile->nsamples * sizeof(double));
  profile->err = realloc(profile->err, profile->nsamples * sizeof(double));

  profile->rgyration = 0.0;
}

void saxs_profile_clone(saxs_profile *dst, saxs_profile *src)
{
  dst->nsamples = src->nsamples;
  dst->q = (double *) malloc(src->nsamples * sizeof(double));
  dst->In = (double *) malloc(src->nsamples * sizeof(double));

  memcpy(dst->q, src->q, src->nsamples * sizeof(double));
  if (src->err != NULL)
  {
    dst->err = (double *) malloc(src->nsamples * sizeof(double));
    memcpy(dst->err, src->err, src->nsamples * sizeof(double));
  }

  memcpy(dst->In, src->In, src->nsamples * sizeof(double));

  dst->npartials = 0;
  dst->vac_vac = NULL;
  dst->vac_dum = NULL;
  dst->dum_dum = NULL;
  dst->vac_h2o = NULL;
  dst->dum_h2o = NULL;
  dst->h2o_h2o = NULL;
}

void saxs_profile_create_partials(saxs_profile *profile,
                                  bool fit_water)
{
  int nsamples = profile->nsamples;

  profile->vac_vac = calloc(nsamples, sizeof(double));
  profile->vac_dum = calloc(nsamples, sizeof(double));
  profile->dum_dum = calloc(nsamples, sizeof(double));

  if (fit_water)
  {
    profile->vac_h2o = calloc(nsamples, sizeof(double));
    profile->dum_h2o = calloc(nsamples, sizeof(double));
    profile->h2o_h2o = calloc(nsamples, sizeof(double));

    profile->npartials = 6;
  }
  else
  {
    profile->vac_h2o = profile->dum_h2o = profile->h2o_h2o = NULL;

    profile->npartials = 3;
  }
}


// Output
void saxs_profile_fprint(FILE *fp, saxs_profile *profile)
{
  for (int i = 0; i < profile->nsamples; ++i)
  {
    double q = profile->q[i];
    double In = profile->In[i];

    if (profile->err != NULL)
      fprintf(fp, "%f %f %f\n", q, In, profile->err[i]);
    else
      fprintf(fp, "%f %f\n", q, In);
  }
}

void saxs_profile_print(saxs_profile *profile)
{
  saxs_profile_fprint(stdout, profile);
}


// Cleanup
void saxs_profile_destroy(saxs_profile *profile)
{
  if (profile != NULL)
  {
    if (profile->q != NULL)
      free(profile->q);
    if (profile->In != NULL)
      free(profile->In);
    if (profile->err != NULL)
      free(profile->err);

    if (profile->vac_vac != NULL)
      free(profile->vac_vac);
    if (profile->vac_dum != NULL)
      free(profile->vac_dum);
    if (profile->dum_dum != NULL)
      free(profile->dum_dum);
    if (profile->vac_h2o != NULL)
      free(profile->vac_h2o);
    if (profile->dum_h2o != NULL)
      free(profile->dum_h2o);
    if (profile->h2o_h2o != NULL)
      free(profile->h2o_h2o);

    free(profile);
  }
}

void saxs_profile_reset_In(saxs_profile *profile)
{
  memset(profile->In, 0, profile->nsamples * sizeof(double));
}


// Actions
void saxs_profile_scale(saxs_profile *profile, double scale)
{
  for (int iq = 0; iq < profile->nsamples; ++iq)
  {
    profile->In[iq] *= scale;
  }
}

void saxs_profile_add(saxs_profile *dst,
                      saxs_profile *p1, saxs_profile *p2,
                      double w1, double w2)
{
  if (p1->nsamples != p2->nsamples)
    return;

  for (int iq = 0; iq < p1->nsamples; ++iq)
  {
    dst->In[iq] = p1->In[iq] * w1 + p2->In[iq] * w2;
  }
}

void compute_profile(
  saxs_profile *profile,
  mol_atom_group *structure,
  double *vacuum_ff, double *dummy_ff,
  double max_dist, form_factor_table *ff_table)
{
  if(profile->npartials == 0)
    saxs_profile_create_partials(profile, false);

  compute_profile_partials(
    profile,
    structure,
    NULL, vacuum_ff, dummy_ff,
    max_dist, ff_table);
  saxs_profile_sum_partials(profile, 1.0, 0.0); // c1 = 1.0, c2 = 0.0
}

void compute_profile_partials(
  saxs_profile *profile,
  mol_atom_group *structure,
  float *saxs_sa, double *vacuum_ff, double *dummy_ff,
  double max_d, form_factor_table *ff_table)
{
  int ndists = saxs_sa == NULL ? 3 : 6;

  assert(ndists <= profile->npartials);
  if (ndists > profile->npartials)
    return;

  radial_distribution r_distributions[ndists];
  for (int i = 0; i < ndists; ++i)
    radial_distribution_create(&r_distributions[i], 0.5, max_d);

  memset(profile->vac_vac, 0, profile->nsamples * sizeof(double));
  memset(profile->vac_dum, 0, profile->nsamples * sizeof(double));
  memset(profile->dum_dum, 0, profile->nsamples * sizeof(double));
  if (ndists == 6)
  {
    memset(profile->vac_h2o, 0, profile->nsamples * sizeof(double));
    memset(profile->dum_h2o, 0, profile->nsamples * sizeof(double));
    memset(profile->h2o_h2o, 0, profile->nsamples * sizeof(double));
  }

  int natoms = structure->natoms;
  double vac_ff_i, dum_ff_i, h2o_ff_i;
  double vac_ff_j, dum_ff_j, h2o_ff_j;
  double h2o_ff = water_ff(ff_table);

  double dist;
  struct atom *ai, *aj;
  {
    for (int i = 0; i < natoms; ++i)
    {
      ai = &(structure->atoms[i]);
      vac_ff_i = vacuum_ff[i];
      dum_ff_i = dummy_ff[i];
      if (saxs_sa != NULL)
        h2o_ff_i = h2o_ff * saxs_sa[i];

      for (int j = i + 1; j < natoms; ++j)
      {
        aj = &(structure->atoms[j]);
        vac_ff_j = vacuum_ff[j];
        dum_ff_j = dummy_ff[j];
        dist = distance_3d(ai, aj);

        add2distribution(&(r_distributions[0]), dist,
                            2 * vac_ff_i * vac_ff_j); // constant
        add2distribution(&(r_distributions[1]), dist,
                            2 * dum_ff_i * dum_ff_j); // c1^2
        add2distribution(&(r_distributions[2]), dist,
                            2 * (vac_ff_i * dum_ff_j +
                                 vac_ff_j * dum_ff_i)); // -c1
        if (saxs_sa != NULL)
        {
          h2o_ff_j = h2o_ff * saxs_sa[j];

          add2distribution(&(r_distributions[3]), dist,
                              2 * h2o_ff_i * h2o_ff_j); // c2^2
          add2distribution(&(r_distributions[4]), dist,
                              2 * (vac_ff_i * h2o_ff_j +
                                   vac_ff_j * h2o_ff_i)); // c2
          add2distribution(&(r_distributions[5]), dist,
                              2 * (h2o_ff_i * dum_ff_j +
                                   h2o_ff_j * dum_ff_i)); // -c1*c2
        }
      }

      // Autocorrelation
      add2distribution(&(r_distributions[0]), 0,
                          vac_ff_i * vac_ff_i); // constant
      add2distribution(&(r_distributions[1]), 0,
                          dum_ff_i * dum_ff_i); // c1^2
      add2distribution(&(r_distributions[2]), 0,
                          2 * vac_ff_i * dum_ff_i); // -c
      if (saxs_sa != NULL)
      {
        add2distribution(&(r_distributions[3]), 0,
                            h2o_ff_i * h2o_ff_i); // c2^2
        add2distribution(&(r_distributions[4]), 0,
                            2 * vac_ff_i * h2o_ff_i);// c2
        add2distribution(&(r_distributions[5]), 0,
                            2 * h2o_ff_i * dum_ff_i); // -c1*c2
      }
    }
  }

  radial_distributions_to_partials(profile, ndists, r_distributions);

  for (int i = 0; i < ndists; ++i)
  {
    radial_distribution_destroy(&(r_distributions[i]));
  }
}

void saxs_profile_sum_partials(saxs_profile *profile, double c1, double c2)
{
  double rm = profile->average_radius;
  double coeff = -pow(4.0 * M_PI / 3.0, 3.0/2.0) * (c1 * c1 - 1.0) / (16*M_PI);
  coeff *= (rm * rm);
  int npartials = profile->npartials;

  if (npartials < 3)
    return;

  // Initialize profile
  saxs_profile_reset_In(profile);

  // Add profiles
  double q, G_q;
  for (int iq = 0; iq < profile->nsamples; ++iq)
  {
    q = profile->q[iq];
    G_q = (c1*c1*c1)*exp(coeff*q*q);

    profile->In[iq] +=
      profile->vac_vac[iq] +
      profile->dum_dum[iq] * (G_q * G_q) +
      profile->vac_dum[iq] * (-G_q);

    if (npartials == 6)
    {
      profile->In[iq] +=
        profile->h2o_h2o[iq] * (c2 * c2) +
        profile->vac_h2o[iq] * (c2) +
        profile->dum_h2o[iq] * (-G_q * c2);
    }
  }
}

void compute_docked_profile(
  saxs_profile *profile,
  mol_atom_group *pa, mol_atom_group *pb,
  saxs_profile *pa_profile, saxs_profile *pb_profile,
  double max_dist, form_factor_table *ff_table, struct prm* prm)
{
  mol_atom_group *mol_pa = pa;
  mol_atom_group *mol_pb = pb;

  int ndists = 3;
  double delta_x = 0.5; // TODO
  radial_distribution *r_distributions =
    (radial_distribution *) malloc(ndists * sizeof(radial_distribution));
  for (int i = 0; i < ndists; ++i)
  {
    radial_distribution_create(&(r_distributions[i]), delta_x, max_dist);
  }

  double dist;
  struct atom *ai, *aj;
  double vac_ff_ai, vac_ff_bj;
  double dum_ff_ai, dum_ff_bj;
  for (int i = 0; i < mol_pa->natoms; ++i)
  {
    ai = &(mol_pa->atoms[i]);
    for (int j = 0; j < mol_pb->natoms; ++j)
    {
      aj = &(mol_pb->atoms[j]);
      dist = distance_3d(ai, aj);

      add2distribution(&(r_distributions[0]), dist,
                       2 * vac_ff_ai * vac_ff_bj); // constant
      add2distribution(&(r_distributions[1]), dist,
                       2 * dum_ff_ai * dum_ff_bj); // c1^2
      add2distribution(&(r_distributions[2]), dist,
                       2 * (vac_ff_ai * dum_ff_bj +
                            vac_ff_bj * dum_ff_ai)); // -c1
    }
  }

  // Create partial profiles
  // radial_distributions2partials(r_distributions, ndists, profile);

  // sum_partials(profile, 1.0, 0.0);

  saxs_profile_add(profile, profile, pa_profile, 1.0, 1.0);
  saxs_profile_add(profile, profile, pb_profile, 1.0, 1.0);
}

void compute_profile_debye(
  saxs_profile *profile,
  mol_atom_group *structure,
  float *saxs_sa, double *vacuum_ff, double *dummy_ff,
  double c1, double c2, struct prm* prm, double water_ff)
{
  double q;
  mol_atom *ai, *aj;
  double fi, fj;
  double d;

  for (int i = 0; i < structure->natoms; ++i)
  {
    ai = &(structure->atoms[i]);
    fj = vacuum_ff[i] - c1*dummy_ff[i] + c2*saxs_sa[i]*water_ff;

    for (int j = i + 1; j < structure->natoms; ++j)
    {
      aj = &(structure->atoms[j]);
      fj = vacuum_ff[j] - c1*dummy_ff[j] + c2*saxs_sa[j]*water_ff;
      d = distance_3d(ai, aj);

      for (int iq = 0; iq < profile->nsamples; ++iq)
      {
        q = profile->q[iq];

        profile->In[iq] += fi * fj * sinc(q*d);
      }
    }
  }
}
