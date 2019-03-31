gfortran -Wall -fcheck=all  types.f95 parameter_shared.f95 atom_module.f95 wannier_center_module.f95 count_time.f95 traj.f95 read_input.f95 1a_sample.f95 1b_label.f95 1c6_coord_translation.f95 2a_calculate_dipole_moment.f95 main1_6.f95 -o main1_6 

rm *.mod

# -fbacktrace
## to run:
#./main < input_sample
