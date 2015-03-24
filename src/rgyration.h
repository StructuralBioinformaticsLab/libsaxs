#ifndef _RGYRATION_H_
#define _RGYRATION_H_

#include "saxs_profile.h"

double rgyration_structure(struct atomgrp *structure);
double rgyration_complex(struct atomgrp *rec, struct atomgrp *lig);

#endif /* _RGYRATION_H_ */
