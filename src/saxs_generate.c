#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "fftw3.h"
#include "mol.0.0.4.h"

#include "arg.h"
#include "form_factor_table.h"
#include "profile.h"
#include "test.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  static struct option long_options[] =
    {
      {"help",     no_argument,       0,    'h'},
      {"factors",  optional_argument, 0,    'f'},
      {"ft",       optional_argument, 0,      0},
      {"mappings", required_argument, 0,    'm'},
      {"prms",     required_argument, 0,    'p'},
      {"rots",     required_argument, 0,    'r'},
      {"nrots",    required_argument, 0,    'R'},
      {"pa",       required_argument, 0,      0},
      {"pb",       optional_argument, 0,      0},
      {"test",     no_argument,       0,      0},
      {"verbose",  no_argument,       0,    'v'},
      {"version",  no_argument,       0,    'V'},
      {"cluster",  required_argument, 0,      0},
      {"profile",  optional_argument, 0,    'P'},
      {"command",  required_argument, 0,    'C'},
      {0, 0, 0, 0}
    };

  return 0;
}
