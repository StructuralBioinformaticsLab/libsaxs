#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "rgyration.h"

#include "saxs_common.h"

double rgyration_structure(struct atomgrp *structure)
{
  double centroid[] = {0.0, 0.0, 0.0};
  double rgyr = 0.0;

  // Compute centroid
  for (int i = 0; i < structure->natoms; ++i)
  {
    centroid[0] += structure->atoms[i].X;
    centroid[1] += structure->atoms[i].Y;
    centroid[2] += structure->atoms[i].Z;
  }
  centroid[0] /= structure->natoms;
  centroid[1] /= structure->natoms;
  centroid[2] /= structure->natoms;

  // Compute rgyr = sqrt( (sum_k (r_k - r_c)^2) / N )
  double dx, dy, dz;
  for (int i = 0; i < structure->natoms; ++i)
  {
    dx = centroid[0] - structure->atoms[i].X;
    dy = centroid[1] - structure->atoms[i].Y;
    dz = centroid[2] - structure->atoms[i].Z;
    rgyr += dx*dx + dy*dy + dz+dz;
  }

  rgyr = sqrt(rgyr/structure->natoms);

  return rgyr;
}

double rgyration_complex(struct atomgrp *rec, struct atomgrp *lig)
{
  struct atomgrp *cmplx = join_2atomgrps(rec, lig);
  double rgyr = rgyration_structure(cmplx);
  free_atomgrp(cmplx);
  return rgyr;
}
