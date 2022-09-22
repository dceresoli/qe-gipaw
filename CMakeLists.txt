set(src_gipaw
  src/atomic.f90
  src/biot_savart.f90
  src/cg_psi.f90
  src/cgsolve_all.f90
  src/ch_psi_all.f90
  src/compute_u_kq.f90
  src/core_relax.f90
  src/efg.f90
  src/epr_routines.f90
  src/gipaw_main.f90
  src/gipaw_module.f90
  src/gipaw_results.f90
  src/gipaw_routines_bands.f90
  src/gipaw_routines.f90
  src/gipaw_setup.f90
  src/gipaw_version.f90
  src/greenfunction.f90
  src/h_psiq.f90
  src/hyperfine.f90
  src/init_gipaw_1.f90
  src/init_gipaw_2.f90
  src/init_us_2_no_phase.f90
  src/j_para.f90
  src/knight_shift.f90
  src/mossbauer.f90
  src/nmr_routines.f90
  src/orthoatwfc1.f90
  src/output_magres.f90
  src/paw_gipaw.f90
  src/stop_code.f90
  src/suscept_crystal.f90
  src/symmetrize_field.f90
  src/symmetrize_rho.f90
  src/util.f90
  src/velocity.f90
  src/write_tensor_field.f90
  src/xml_routines.f90)


###########################################################
# gipaw.x
###########################################################
qe_add_executable(qe_gipaw_exe ${src_gipaw})
qe_install_targets(qe_gipaw_exe)
set_target_properties(qe_gipaw_exe
    PROPERTIES
        OUTPUT_NAME gipaw.x)
target_link_libraries(qe_gipaw_exe
    PRIVATE
        qe_hdf5_fortran
        qe_openmp_fortran
        qe_openacc_fortran
        qe_mpi_fortran
        qe_pw
        qe_modules
        qe_upflib
        qe_fftx
        qe_fox
        qe_lax
        qe_upflib
        qe_xclib)

