libsaxs
=======

libsaxs contains functions for computing SAXS profiles, chi-squared
scores, and other things related to SAXS.

Installation
------------

Installation requires libmol v0.0.6, cmake >= 2.8, and gcc. The
following should work if you have these:

    mkdir -p build && cd build/
    cmake ..
    make

To get libmol and install libmol, follow the instructions at https://github.com/StructuralBioinformaticsLab/libmol.

Usage
-----

There are example files for computing the chi score from a pdb structure
and an experimental profile in the `example` directory. After compiling, you can
run this example as follows from the `build` directory:

    ./tools/pdb_chi_score ../example/pdb_formfactor_mapping.prm ../example/atoms.0.0.6.prm ../example/iofq_orig.dat ../example/rec.pdb
