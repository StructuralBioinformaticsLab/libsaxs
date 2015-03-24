#ifndef _FORM_FACTOR_H_
#define _FORM_FACTOR_H_

#include "mol.0.0.6/atom.h"
#include "mol.0.0.6/atom_group.h"
#include "mol.0.0.6/prms.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


typedef enum ff_type_t
{
  H, He,
  Li, Be, B, C, N, O, F, Ne, // periodic table, lines 1-2 (10)
  Na, Mg, Al, Si, P, S, Cl, Ar, // line 3 (8)
  K, Ca, Cr, Mn, Fe, Co,
  Ni, Cu, Zn, Se, Br, // line 4 (11)
  Io, Ir, Pt, Au, Hg, SINGLE_ATOM_SIZE = 34,
  CH=34, CH2=35, CH3=36, NH=37, NH2=38, NH3=39, OH=40, s_OH2=41, SH=42, PO4=43,
  HEAVY_ATOM_SIZE=44, s_UNK=99
} ff_type;

typedef struct atype_map_t
{
  char* key; /**< "atypemaj atypemin" */
  ff_type val; /** ff_type */
} atype_map;

// f[q] = c + (sum_i a_i*e^(-b_i*q(q^2)))
typedef struct form_factor_coeffs_t
{
  ff_type group_id;
  char *group_name;
  double *a;
  double *b;
  double c;
  double excl_vol;
} form_factor_coeffs;

typedef struct form_factor_t
{
  ff_type group_id;
  double *values;
  // Need vacuum_ff, dummy_ff
  double zero_ff;
  double vacuum_ff;
  double dummy_ff;
} form_factor;

typedef struct form_factor_table_t
{
  int natom_types;
  form_factor *factors;
  form_factor_coeffs *coeffs;

  int map_size;
  atype_map *type_map;

  // Temp variables
  atype_map *search_key; // Used for searching
  form_factor *water_ff;
} form_factor_table;


form_factor_table *default_ff_table(char *type_mapping_file);
form_factor_table *read_ff_table(char *ff_file, char *type_mapping_file);
void ff_table_destroy(form_factor_table *table);

form_factor *get_ff(struct atom *a, form_factor_table *table, struct prm *prms);
double dummy_ff(struct atom *a, form_factor_table *table, struct prm *prms);
double vacuum_ff(struct atom *a, form_factor_table *table, struct prm *prms);
double compute_form_factor(form_factor_table *table, ff_type type, double q);
void build_form_factor_tables(form_factor_table *table,
                              double *q_vals, size_t n_qvals);
ff_type atom2ff_type(struct atom *a, form_factor_table *ff_table,
                       struct prm *prm);

void get_dummy_ff(double *dummy_ffs,
                  mol_atom_group *ag,
                  form_factor_table *table, struct prm *prm);
void get_vacuum_ff(double *vacuum_ffs,
                   mol_atom_group *ag,
                   form_factor_table *table, struct prm *prm);
static inline double water_ff(form_factor_table *table)
{
  return table->factors[s_OH2].zero_ff;
}


// debugging functions
void dump_mappings(form_factor_table *table);
void print_ff_table(form_factor_table *table);

#endif /* _FORM_FACTOR_H_ */
