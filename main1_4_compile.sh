gfortran -Wall -fcheck=all  types.f95 parameter_shared.f95 atom_module.f95 wannier_center_module.f95 count_time.f95 traj.f95 read_input.f95 1a_sample.f95 1b_label.f95 1c4_coord_translation.f95  main1_4.f95 -o main1_4 

rm *.mod

# -fbacktrace
## to run:
#./main < input_sample
