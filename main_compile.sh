gfortran -Wall -fcheck=all  types.f95 parameter_shared.f95 atom_module.f95 wannier_center_module.f95 count_time.f95 traj.f95 read_input.f95 1a_sample.f95 1b_label.f95 1c7_coord_translation.f95 2a7_calculate_dipole_moment.f95 main.f95 -o main 

rm *.mod

# -fbacktrace
## to run:
#./main < input_sample
