#ifndef _SAXS_ATOM_H_
#define _SAXS_ATOM_H_

#include "saxs_common.h"
#include "form_factor_table.h"

#include "mol.0.0.6/atom.h"
#include "mol.0.0.6/atom_group.h"
#include "mol.0.0.6/matrix.h"
#include "mol.0.0.6/vector.h"

void faccs(float *sa,
           mol_atom_group *ag, struct prm *prm, float r_solv);

mol_atom_group *remove_hydrogens(struct atomgrp *src, struct prm *prm);

double mol_atom_group_average_radius(mol_atom_group *ag, struct prm *prm);

#endif /* _SAXS_ATOM_H_ */
