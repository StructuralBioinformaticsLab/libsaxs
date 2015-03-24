include_directories(../src)

# Executables
add_executable(gen_saxs_profile gen_saxs_profile.c)
target_link_libraries(gen_saxs_profile saxs)

add_executable(pdb_chi_score pdb_chi_score.c)
target_link_libraries(pdb_chi_score saxs)
