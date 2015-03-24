#include "saxs_atom.h"
#include "saxs_common.h"

#include "saxs_utils.h"

void faccs(float *sa,
           mol_atom_group *ag, struct prm *prm, float r_solv)
{
  static const short cont_acc = 1;
  static const short rpr = 1;
  static const double fourpi = 4 * M_PI;

  accs(ag, prm, r_solv, cont_acc, rpr, sa);

  for (int i = 0; i < ag->natoms; ++i)
  {
    // sa expected by saxs is the fraction of accessible surface area
    double ri = prm->atoms[ag->atoms[i].atom_typen].r;
    if (ri == 0 || isnan(sa[i]))
      sa[i] = 0.0;
    else
      sa[i] = sa[i] / (fourpi * ri*ri);
  }
}

mol_atom_group *remove_hydrogens(struct atomgrp *src, struct prm *prm)
{
  int atomn;
  int natoms;

  mol_atom_group* exag = calloc (1, sizeof (mol_atom_group));
  exag->atoms = malloc (sizeof (struct atom) * src->natoms);

  natoms = 0;
  for (atomn = 0; atomn < src->natoms; atomn++)
  {
    if ( prm->atoms[src->atoms[atomn].atom_typen].typemin[0] != 'H' )
    {
      exag->atoms[natoms].atom_typen = src->atoms[atomn].atom_typen;
      exag->atoms[natoms].sa = src->atoms[atomn].sa;
      exag->atoms[natoms].X = src->atoms[atomn].X;
      exag->atoms[natoms].Y = src->atoms[atomn].Y;
      exag->atoms[natoms].Z = src->atoms[atomn].Z;
      natoms++;
    }
  }

  if (natoms == 0)
    return NULL;

  exag->atoms = realloc (exag->atoms, sizeof (struct atom) * natoms);
  exag->natoms = natoms;

  return exag;
}

double mol_atom_group_average_radius(mol_atom_group *ag, struct prm *prm)
{
  double av_r = 0.0;

  for (int i = 0; i < ag->natoms; ++i)
  {
    av_r += prm->atoms[ag->atoms[i].atom_typen].r;
  }

  return av_r / ag->natoms;
}
