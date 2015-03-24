# Try to find libmol automatically

if(LIBMOL_LIBRARIES AND LIBMOL_INCLUDE_DIRS)
  set(LIBMOL_FOUND TRUE)
else(LIBMOL_LIBRARIES AND LIBMOL_INCLUDE_DIRS)
  find_path(LIBMOL_INCLUDE_DIR
    NAMES
      mol.0.0.6.h
    PATHS
      /usr/include
      /usr/local/include
      /opt/local/include
      /sw/include
      $ENV{HOME}/include
      $ENV{HOME}/usr/include)
  find_library(LIBMOL_LIBRARY
    NAMES
      mol.0.0.6
    PATHS
      /usr/lib
      /usr/local/lib
      /opt/local/lib
      /sw/lib
      $ENV{HOME}/lib
      $ENV{HOME}/usr/lib)

  set(LIBMOL_INCLUDE_DIRS ${LIBMOL_INCLUDE_DIR})

  if(LIBMOL_LIBRARY)
    set(LIBMOL_LIBRARIES
      ${LIBMOL_LIBRARIES}
      ${LIBMOL_LIBRARY})
  endif(LIBMOL_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(libmol DEFAULT_MSG
    LIBMOL_LIBRARIES LIBMOL_INCLUDE_DIRS)
endif(LIBMOL_LIBRARIES AND LIBMOL_INCLUDE_DIRS)

mark_as_advanced( LIBMOL_INCLUDE_DIRS LIBMOL_LIBRARIES )