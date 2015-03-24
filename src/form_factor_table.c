#include "form_factor_table.h"

#include "saxs_utils.h"

static double default_zero_form_factors[] = {
  -0.720147, -0.720228,
   //   H       He
  1.591, 2.591, 3.591, 0.50824, 6.16294, 4.94998, 7.591, 6.993,
  // Li     Be      B     C       N        O       F      Ne
  7.9864, 8.9805, 9.984, 10.984, 13.0855, 9.36656, 13.984, 16.591,
  //  Na      Mg     Al     Si      P        S       Cl     Ar
  15.984, 14.9965, 20.984, 21.984, 20.9946, 23.984,
  // K     Ca2+     Cr      Mn      Fe2+      Co
  24.984, 25.984, 24.9936, 30.9825, 31.984, 49.16,
  // Ni   Cu      Zn2+       Se      Br      I
  70.35676, 71.35676, 72.324,  73.35676,
  // Ir         Pt      Au      Hg
  -0.211907, -0.932054, -1.6522, 5.44279, 4.72265,4.0025,4.22983,3.50968,8.64641
  //  CH        CH2        CH3     NH       NH2       NH3     OH       OH2   SH
};

static double default_vacuum_zero_form_factors[] = {
  0.999953, 0.999872, 2.99, 3.99, 4.99, 5.9992, 6.9946, 7.9994, 8.99, 9.999,
  //   H       He      Li     Be     B    C       N       O      F     Ne
  10.9924, 11.9865, 12.99, 13.99, 14.9993, 15.9998, 16.99, 17.99,
  //  Na      Mg     Al     Si      P        S       Cl     Ar
  18.99, 18.0025,  23.99, 24.99,  24.0006, 26.99,
  // K     Ca2+     Cr     Mn      Fe2+      Co
  27.99, 28.99, 27.9996, 33.99, 34.99, 52.99, 76.99, 77.99, 78.9572, 79.99,
  // Ni   Cu      Zn2+    Se     Br     I       Ir     Pt      Au     Hg
  6.99915, 7.99911, 8.99906, 7.99455, 8.99451, 9.99446, 8.99935, 9.9993, 16.9998
  //  CH      CH2     CH3     NH       NH2       NH3     OH      OH2      SH
};

static double default_dummy_zero_form_factors[] = {
  1.7201, 1.7201, 1.399, 1.399, 1.399 , 5.49096, 0.83166, 3.04942, 1.399, 3.006,
  //  H     He     Li?    Be?    B?       C        N        O      F?     Ne
  3.006, 3.006, 3.006, 3.006, 1.91382, 6.63324, 3.006, 1.399,
  // Na     Mg    Al?    Si?      P        S      Cl?    Ar?
  3.006, 3.006, 3.006, 3.006, 3.006, 3.006,
  // K?   Ca2+    Cr?    Mn?   Fe2+   Co?
  3.006, 3.006, 3.006, 3.006, 3.006, 3.83, 6.63324, 6.63324, 6.63324, 6.63324,
  // Ni?   Cu?   Zn2+    Se     Br?     I?   Ir?      Pt?       Au      Hg
  7.21106, 8.93116, 10.6513, 2.55176, 4.27186, 5.99196, 4.76952, 6.48962,8.35334
  //  CH       CH2      CH3     NH       NH2       NH3     OH       OH2   SH
};

static int compatype_map (void* a1, void *a2)
{
  atype_map* atype1 = (atype_map *) a1;
  atype_map* atype2 = (atype_map *) a2;
  return (strcmp (atype1->key, atype2->key));
}

static ff_type string_to_ff_type(char *name)
{
  static char cmp[8];
  strcpy(cmp, name);
  upcase(cmp);

  if (!strcmp(name, "H"))     return H;
  if (!strcmp(name, "HE"))    return He;
  if (!strcmp(name, "C"))     return C;
  if (!strcmp(name, "N"))     return N;
  if (!strcmp(name, "O"))     return O;
  if (!strcmp(name, "NE"))    return Ne;
  if (!strcmp(name, "SOD+"))  return Na;
  if (!strcmp(name, "MG2+"))  return Mg;
  if (!strcmp(name, "P"))     return P;
  if (!strcmp(name, "S"))     return S;
  if (!strcmp(name, "K"))     return K;
  if (!strcmp(name, "CAL2+")) return Ca;
  if (!strcmp(name, "FE2+"))  return Fe;
  if (!strcmp(name, "ZN2+"))  return Zn;
  if (!strcmp(name, "SE"))    return Se;
  if (!strcmp(name, "AU"))    return Au;
  if (!strcmp(name, "CH"))    return CH;
  if (!strcmp(name, "CH2"))   return CH2;
  if (!strcmp(name, "CH3"))   return CH3;
  if (!strcmp(name, "NH"))    return NH;
  if (!strcmp(name, "NH2"))   return NH2;
  if (!strcmp(name, "NH3"))   return NH3;
  if (!strcmp(name, "OH"))    return OH;
  if (!strcmp(name, "SH"))    return SH;
  return s_UNK;
}

static const char *ff_type_to_sring(ff_type type)
{
  if (type == H) return "H";
  if (type == He) return "HE";
  if (type == C) return "C";
  if (type == N) return "N";
  if (type == O) return "O";
  if (type == Ne) return "NE";
  if (type == Na) return "SOD+";
  if (type == Mg) return "MG2+";
  if (type == P) return "P";
  if (type == S) return "S";
  if (type == K) return "K";
  if (type == Ca) return "CAL2+";
  if (type == Fe) return "FE2+";
  if (type == Zn) return "ZN2+";
  if (type == Se) return "SE";
  if (type == Au) return "AU";
  if (type == CH) return "CH";
  if (type == CH2) return "CH2";
  if (type == CH3) return "CH3";
  if (type == NH) return "NH";
  if (type == NH2) return "NH2";
  if (type == NH3) return "NH3";
  if (type == OH) return "OH";
  if (type == SH) return "SH";
  return NULL;
}

form_factor_table *default_ff_table(char *type_mapping_file)
{
  static form_factor_table *table = NULL;
  static char space[] = " ";

  if (table != NULL)
    return table;

  FILE *tm_fp = fopen(type_mapping_file, "r");
  if (tm_fp == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", type_mapping_file);
    exit(1);
  }

  table = (form_factor_table *) malloc(sizeof(form_factor_table));

  table->natom_types = 0;
  table->factors = (form_factor *) malloc(HEAVY_ATOM_SIZE * sizeof(form_factor));
  table->coeffs = (form_factor_coeffs *) malloc(HEAVY_ATOM_SIZE * sizeof(form_factor_coeffs));
  for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
  {
    table->factors[i].group_id = i;
    table->coeffs[i].group_id = i;
    table->factors[i].values = NULL;
    table->factors[i].zero_ff = default_zero_form_factors[i];
    table->factors[i].vacuum_ff = default_vacuum_zero_form_factors[i];
    table->factors[i].dummy_ff = default_dummy_zero_form_factors[i];
  }

  int arr_size = 128;
  table->map_size = 0;
  table->type_map = (atype_map *) malloc(arr_size * sizeof(atype_map));

  char *line = (char *) malloc(256 * sizeof(char));
  size_t len;
  //                      res name atom
  char *tm_line_format = "%s  %s   %s";
  char atypemaj[16], atypemin[16], ff_type[16];

  // Create search key (so we don't have to malloc everytime)
  table->search_key = (atype_map *) malloc(sizeof(atype_map));
  table->search_key->key = (char *) malloc(64 * sizeof(char)); // Allocate 64 bytes for key

  while (getline(&line, &len, tm_fp) != -1)
  {
    if (line[0] == '\0' || line[0] == '#') // Skip empty lines
    {
      continue;
    }

    if (table->map_size >= arr_size)
    {
      arr_size *= 2;
      table->type_map = (atype_map *) realloc(table->type_map, arr_size * sizeof(atype_map));
    }

    sscanf(line, tm_line_format, atypemaj, atypemin, ff_type);

    atype_map *atkey = &(table->type_map[table->map_size]);
    atkey->key = (char*) malloc (64 * sizeof (char));
    atkey->key = strcpy (atkey->key, atypemaj);
    atkey->key = strcat (atkey->key, space);
    atkey->key = strcat (atkey->key, atypemin);
    atkey->val = string_to_ff_type(ff_type);

    // Make factor valid
    if (atkey->val != s_UNK)
      table->factors[atkey->val].group_id = atkey->val;

    ++(table->map_size);
  }
  fclose(tm_fp);

  // Repack tables
  table->type_map = (atype_map *) realloc(table->type_map, table->map_size * sizeof(atype_map));
  free(line);

  return table;
}

form_factor_table *read_ff_table(char *ff_file, char *type_mapping_file)
{
  FILE *ff_fp = fopen(ff_file, "r");
  if (ff_fp == NULL)
  {
    return NULL;
  }

  form_factor_table *table = (form_factor_table *) malloc(sizeof(form_factor_table));

  table->natom_types = 0;
  table->factors = (form_factor *) malloc(HEAVY_ATOM_SIZE * sizeof(form_factor));
  table->coeffs = (form_factor_coeffs *) malloc(HEAVY_ATOM_SIZE * sizeof(form_factor_coeffs));
  for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
  {
    table->factors[i].group_id = s_UNK;
    table->coeffs[i].group_id = s_UNK;
  }

  char *line = (char *) malloc(256 * sizeof(char));
  size_t len;
  char group_name[16];
  double a1, a2, a3, a4, a5;
  double b1, b2, b3, b4, b5;
  double c, excl_vol;
  //                   gn a1  a2  a3  a4  a5  c   b1  b2  b3  b4  b4  vol
  char *line_format = "%15s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  form_factor_coeffs *cc;

  while (getline(&line, &len, ff_fp) != -1)
  {
    if (line[0] == '#')
    {
      continue;
    }

    sscanf(line, line_format, group_name,
           &a1, &a2, &a3, &a3, &a4,
           &c,
           &b1, &b2, &b3, &b3, &b4,
           &excl_vol);
    ff_type type = string_to_ff_type(group_name);
    if (type != s_UNK)
    {
      // Fill in coefficient data
      cc = &(table->coeffs[type]);
      cc->a = (double *) malloc(5 * sizeof(double));
      cc->b = (double *) malloc(5 * sizeof(double));
      cc->a[0] = a1; cc->a[1] = a2; cc->a[2] = a3; cc->a[3] = a4; cc->a[4] = a5;
      cc->b[0] = b1; cc->b[1] = b2; cc->b[2] = b3; cc->b[3] = b4; cc->b[4] = b5;
      cc->group_name = (char *) malloc((strlen(group_name) + 1) * sizeof(char));
      cc->excl_vol = excl_vol;
      strcpy(cc->group_name, group_name);
      cc->group_id = type;
      ++(table->natom_types);

      // Make the form_factor valid
      table->factors[type].group_id = type;
    }
  }
  fclose(ff_fp);

  // Read type mappings
  FILE *tm_fp = fopen(type_mapping_file, "r");

  int arr_size = 128;
  table->map_size = 0;
  table->type_map = (atype_map *) malloc(arr_size * sizeof(atype_map));
  char space[] = " ";
  //                      res name atom
  char *tm_line_format = "%s  %s   %s";
  char atypemaj[16], atypemin[16], ff_type[16];

  // Create search key (so we don't have to malloc everytime)
  table->search_key = (atype_map *) malloc(sizeof(atype_map));
  table->search_key->key = (char *) malloc(64 * sizeof(char)); // Allocate 64 bytes for key

  while (getline(&line, &len, tm_fp) != -1)
  {
    if (line[0] == '\0' || line[0] == '#') // Skip empty lines
    {
      continue;
    }

    if (table->map_size >= arr_size)
    {
      arr_size *= 2;
      table->type_map = (atype_map *) realloc(table->type_map, arr_size * sizeof(atype_map));
    }

    sscanf(line, tm_line_format, atypemaj, atypemin, ff_type);

    atype_map *atkey = &(table->type_map[table->map_size]);
    atkey->key = (char*) malloc (64 * sizeof (char));
    atkey->key = strcpy (atkey->key, atypemaj);
    atkey->key = strcat (atkey->key, space);
    atkey->key = strcat (atkey->key, atypemin);
    atkey->val = string_to_ff_type(ff_type);

    // Make factor valid
    table->factors[atkey->val].group_id = atkey->val;

    ++(table->map_size);
  }
  fclose(tm_fp);

  // Repack tables
  table->type_map = realloc(table->type_map,
                            table->map_size * sizeof(atype_map));

  return table;
}

void ff_table_destroy(form_factor_table *table)
{
  if (table != NULL)
  {
    if (table->coeffs != NULL)
    {
      for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
      {
        form_factor_coeffs *coeffs = &(table->coeffs[i]);
        if (coeffs->group_name !=NULL)
        {
          free(coeffs->group_name);
        }
        if (coeffs->a !=NULL)
        {
          free(coeffs->a);
        }
        if (coeffs->b !=NULL)
        {
          free(coeffs->a);
        }
      }
      free(table->coeffs);
    }

    if (table->factors != NULL)
    {
      for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
      {

      }
      free(table->factors);
    }

    free(table);
  }
}

form_factor *get_ff(struct atom *a, form_factor_table *table, struct prm *prm)
{
  static char space[] = " ";

  char *atypemaj = prm->atoms[a->atom_typen].typemaj;
  char *atypemin = prm->atoms[a->atom_typen].typemin;

  if (strlen(atypemaj) == 4) {
    atypemaj[3] = '\0';
  }

  char key[16];
  atype_map search_key;
  search_key.key = key;

  atype_map *atkey = &search_key;
  atype_map *res;

  atkey->key = strcpy (atkey->key, atypemaj);
  atkey->key = strcat (atkey->key, space);
  atkey->key = strcat (atkey->key, atypemin);

  res = bsearch (atkey, table->type_map, table->map_size,
                 sizeof(struct atype_map_t), (void *) compatype_map);

  if (res == NULL)
  {
    fprintf(stderr, "mapping for %s not found\n", atkey->key);
    return NULL;
  }
  else
  {
    if (table->factors[res->val].group_id != s_UNK)
    {
      return &(table->factors[res->val]);
    }
    else
    {
      printf("%s maps to UNK\n", atkey->key);
      return NULL;
    }
  }
}

double compute_form_factor(form_factor_table *table, ff_type type, double q)
{
  form_factor_coeffs *coeff = &(table->coeffs[type]);

  double f = coeff->c;
  // f = c + (SUM_i a_i*e^(-b_i*(q^2)))
  for (int i = 0; i < 5; ++i)
  {
    f += coeff->a[i] * exp(-coeff->b[i] * (q*q));
  }

  return 0.0;
}

void build_form_factor_tables(form_factor_table *table, double *q_vals,
                              size_t n_qvals)
{
  static double two_thirds = 2.0/3.0;
  static double one_over_four_pi = 1.0/(4.0*M_PI);
  static double rho = 0.334; // TODO: replace with electron density of solvent

  // Precompute qq and ss;
  double *qq = (double *) malloc(n_qvals * sizeof(double));
  for (size_t i = 0; i < n_qvals; ++i)
  {
    qq[i] = q_vals[i] * q_vals[i];
  }
  double *ss = (double *) malloc(n_qvals * sizeof(double));
  for (size_t i = 0; i < n_qvals; ++i)
  {
    ss[i] = qq[i] * one_over_four_pi * one_over_four_pi;
  }

  // Create tables for single atoms
  for (int ia = 0; ia < SINGLE_ATOM_SIZE; ++ia)
  {
    form_factor_coeffs *coeff = &(table->coeffs[ia]);
    if (coeff->group_id == s_UNK)
      continue;
    form_factor *factor = &(table->factors[ia]);
    factor->values = (double *) malloc(n_qvals * sizeof(double));

    for (int iq = 0; iq < n_qvals; ++iq)
    {
      factor->values[iq] = coeff->c;
      // f = c + (SUM_i a_i*e^(-b_i*(q^2)))
      for (int j = 0; j < 5; ++j)
      {
        factor->values[iq] += coeff->a[j] * exp(-coeff->b[j] * (ss[iq]));
      }

      double vol_c = - pow(coeff->excl_vol, two_thirds) * one_over_four_pi;
      factor->values[iq] -= rho * coeff->excl_vol * exp(vol_c * qq[iq]);
    }
    factor->zero_ff = coeff->c;
    for (int j = 0; j < 5; ++j)
    {
      factor->zero_ff += coeff->a[j];
    }
    factor->vacuum_ff = factor->zero_ff;
    factor->dummy_ff = rho * coeff->excl_vol;
    factor->zero_ff -= rho * coeff->excl_vol;
  }

  // Create tables for compound groups
  form_factor *CH_factor = &(table->factors[CH]);
  form_factor *CH2_factor = &(table->factors[CH2]);
  form_factor *CH3_factor = &(table->factors[CH3]);
  form_factor *NH_factor = &(table->factors[NH]);
  form_factor *NH2_factor = &(table->factors[NH2]);
  form_factor *NH3_factor = &(table->factors[NH3]);
  form_factor *OH_factor = &(table->factors[OH]);
  form_factor *OH2_factor = &(table->factors[s_OH2]);
  form_factor *SH_factor = &(table->factors[SH]);
  form_factor *PO4_factor = &(table->factors[PO4]);

  CH_factor->values = (double *) malloc(n_qvals * sizeof(double));
  CH2_factor->values = (double *) malloc(n_qvals * sizeof(double));
  CH3_factor->values = (double *) malloc(n_qvals * sizeof(double));
  NH_factor->values = (double *) malloc(n_qvals * sizeof(double));
  NH2_factor->values = (double *) malloc(n_qvals * sizeof(double));
  NH3_factor->values = (double *) malloc(n_qvals * sizeof(double));
  OH_factor->values = (double *) malloc(n_qvals * sizeof(double));
  OH2_factor->values = (double *) malloc(n_qvals * sizeof(double));
  SH_factor->values = (double *) malloc(n_qvals * sizeof(double));
  PO4_factor->values = (double *) malloc(n_qvals * sizeof(double));

  // Fill in factor table by summing component atoms
  for (int iq = 0; iq < n_qvals; ++iq)
  {
    double q = q_vals[iq];

    CH_factor->values[iq] =
      table->factors[C].values[iq] + table->factors[H].values[iq];
    CH2_factor->values[iq] =
      table->factors[C].values[iq] + 2 * table->factors[H].values[iq];
    CH3_factor->values[iq] =
      table->factors[C].values[iq] + 3 * table->factors[H].values[iq];
    NH_factor->values[iq] =
      table->factors[N].values[iq] + table->factors[H].values[iq];
    NH2_factor->values[iq] =
      table->factors[N].values[iq] + 2 * table->factors[H].values[iq];
    NH3_factor->values[iq] =
      table->factors[N].values[iq] + 3 * table->factors[H].values[iq];
    OH_factor->values[iq] =
      table->factors[O].values[iq] + table->factors[H].values[iq];
    OH2_factor->values[iq] =
      table->factors[O].values[iq] + 2 * table->factors[H].values[iq];
    SH_factor->values[iq] =
      table->factors[S].values[iq] + table->factors[H].values[iq];
    PO4_factor->values[iq] =
      table->factors[P].values[iq] + 4 * table->factors[O].values[iq];
  }

  // Compute zero, dummy, and vacuum by summing
  CH_factor->zero_ff = table->factors[C].zero_ff + table->factors[H].zero_ff;
  CH2_factor->zero_ff = table->factors[C].zero_ff + 2 * table->factors[H].zero_ff;
  CH3_factor->zero_ff = table->factors[C].zero_ff + 3 * table->factors[H].zero_ff;
  NH_factor->zero_ff = table->factors[N].zero_ff + table->factors[H].zero_ff;
  NH2_factor->zero_ff = table->factors[N].zero_ff + 2 * table->factors[H].zero_ff;
  NH3_factor->zero_ff = table->factors[N].zero_ff + 3 * table->factors[H].zero_ff;
  OH_factor->zero_ff = table->factors[O].zero_ff + table->factors[H].zero_ff;
  OH2_factor->zero_ff = table->factors[O].zero_ff + 2 * table->factors[H].zero_ff;
  SH_factor->zero_ff = table->factors[S].zero_ff + table->factors[H].zero_ff;
  PO4_factor->zero_ff = table->factors[P].zero_ff + 4 * table->factors[O].zero_ff;

  CH_factor->dummy_ff = table->factors[C].dummy_ff + table->factors[H].dummy_ff;
  CH2_factor->dummy_ff = table->factors[C].dummy_ff + 2 * table->factors[H].dummy_ff;
  CH3_factor->dummy_ff = table->factors[C].dummy_ff + 3 * table->factors[H].dummy_ff;
  NH_factor->dummy_ff = table->factors[N].dummy_ff + table->factors[H].dummy_ff;
  NH2_factor->dummy_ff = table->factors[N].dummy_ff + 2 * table->factors[H].dummy_ff;
  NH3_factor->dummy_ff = table->factors[N].dummy_ff + 3 * table->factors[H].dummy_ff;
  OH_factor->dummy_ff = table->factors[O].dummy_ff + table->factors[H].dummy_ff;
  OH2_factor->dummy_ff = table->factors[O].dummy_ff + 2 * table->factors[H].dummy_ff;
  SH_factor->dummy_ff = table->factors[S].dummy_ff + table->factors[H].dummy_ff;
  PO4_factor->dummy_ff = table->factors[P].dummy_ff + 4 * table->factors[O].dummy_ff;

  CH_factor->vacuum_ff =
    table->factors[C].vacuum_ff + table->factors[H].vacuum_ff;
  CH2_factor->vacuum_ff =
    table->factors[C].vacuum_ff + 2 * table->factors[H].vacuum_ff;
  CH3_factor->vacuum_ff =
    table->factors[C].vacuum_ff + 3 * table->factors[H].vacuum_ff;
  NH_factor->vacuum_ff =
    table->factors[N].vacuum_ff + table->factors[H].vacuum_ff;
  NH2_factor->vacuum_ff =
    table->factors[N].vacuum_ff + 2 * table->factors[H].vacuum_ff;
  NH3_factor->vacuum_ff =
    table->factors[N].vacuum_ff + 3 * table->factors[H].vacuum_ff;
  OH_factor->vacuum_ff =
    table->factors[O].vacuum_ff + table->factors[H].vacuum_ff;
  OH2_factor->vacuum_ff =
    table->factors[O].vacuum_ff + 2 * table->factors[H].vacuum_ff;
  SH_factor->vacuum_ff =
    table->factors[S].vacuum_ff + table->factors[H].vacuum_ff;
  PO4_factor->vacuum_ff =
    table->factors[P].vacuum_ff + 4 * table->factors[O].vacuum_ff;
}

void dump_mappings(form_factor_table *table)
{
  for (int i = 0; i < table->map_size; ++i)
  {
    atype_map *mapping = &(table->type_map[i]);
    printf("%s %d\n", mapping->key, mapping->val);
  }
}

void print_ff_table(form_factor_table *table)
{
  for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
  {
    fprintf(stderr, "FFTYPE %2d zero_ff: %f vacuum_ff: %f dummy_ff: %f\n", i,
            table->factors[i].zero_ff,
            table->factors[i].vacuum_ff,
            table->factors[i].dummy_ff);
  }
}

ff_type atom2ff_type(struct atom* a, form_factor_table *ff_table,
                     struct prm *prm)
{
  static char space[] = " ";

  char *atypemaj = prm->atoms[a->atom_typen].typemaj;
  char *atypemin = prm->atoms[a->atom_typen].typemin;

  atype_map *atkey;
  atype_map *res;

#pragma omp critical(ff_update_search_key)
  {
    atkey = ff_table->search_key;
    atkey->key = strcpy (atkey->key, atypemaj);
    atkey->key = strcat (atkey->key, space);
    atkey->key = strcat (atkey->key, atypemin);

    res = bsearch (atkey, ff_table->type_map, ff_table->map_size,
                   sizeof(struct atype_map_t), (void*) compatype_map);
  }

  if (res == NULL)
  {
    fprintf(stderr, "mapping for %s not found\n", atkey->key);
    return s_UNK;
  }

  return res->val;
}

double dummy_ff(struct atom *a, form_factor_table *table, struct prm *prm)
{
  return get_ff(a, table, prm)->dummy_ff;
}

double vacuum_ff(struct atom *a, form_factor_table *table, struct prm *prm)
{
  return get_ff(a, table, prm)->vacuum_ff;
}

void get_dummy_ff(double *dummy_ffs,
                  mol_atom_group *ag, form_factor_table *table, struct prm *prm)
{
  for (int i = 0; i < ag->natoms; ++i)
  {
    dummy_ffs[i] = dummy_ff(&(ag->atoms[i]), table, prm);
  }
}

void get_vacuum_ff(double *vacuum_ffs,
                   mol_atom_group *ag, form_factor_table *table, struct prm *prm)
{
  for (int i = 0; i < ag->natoms; ++i)
  {
   vacuum_ffs[i] = vacuum_ff(&(ag->atoms[i]), table, prm);
  }
}
