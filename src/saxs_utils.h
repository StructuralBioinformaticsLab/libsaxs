#ifndef _SAXS_UTILS_H_
#define _SAXS_UTILS_H_

#include "saxs_common.h"

#include "mol.0.0.6/atom_group.h"
#include "mol.0.0.6/sasa.h"

#ifndef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

#ifndef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif

typedef struct dec_ftresult_t
{
  struct ftresult *res;
  double rgyr;
  double chi;
  int rank;
} dec_ftresult;

double max_dist(struct atomgrp *ag);
void upcase(char *s);

static inline double distance_3d_squared(struct atom *ai, struct atom *aj)
{
  return
    (ai->X - aj->X) * (ai->X - aj->X) +
    (ai->Y - aj->Y) * (ai->Y - aj->Y) +
    (ai->Z - aj->Z) * (ai->Z - aj->Z);
}

static inline double distance_3d(struct atom *ai, struct atom *aj)
{
  return sqrt(
    (ai->X - aj->X) * (ai->X - aj->X) +
    (ai->Y - aj->Y) * (ai->Y - aj->Y) +
    (ai->Z - aj->Z) * (ai->Z - aj->Z));
}

static inline double sinc(double x)
{
  return x == 0.0 ? 1 : sin(x)/x;
}

static inline double square(double x)
{
  return x*x;
}

#endif /* _SAXS_UTILS_H_ */
