#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "saxs_utils.h"

#include "saxs_common.h"

double max_dist(struct atomgrp *ag)
{
  double max_d = 0.0;

  for (int i = 0; i < ag->natoms; ++i)
  {
    struct atom *ai = &(ag->atoms[i]);
    for (int j = i + 1; j < ag->natoms; ++j)
    {
      struct atom *aj = &(ag->atoms[j]);

      double d = distance_3d_squared(ai, aj);

      max_d = d > max_d ? d : max_d; // one line max
    }
  }

  return sqrt(max_d);
}

void upcase(char *s)
{
  int len = strlen(s);

  for (int i = 0; i < len; ++i)
  {
    s[i] = toupper(s[i]);
  }
}
