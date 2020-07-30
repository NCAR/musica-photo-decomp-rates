################################################################################
# musica-photo-decomp-rates library

set(PHOTO_SRC oxygen_xsect/la_srb_mod.F90
              oxygen_xsect/molec_ox_xsect.F90
              photo_rates/module_prates_tuv.F90
              photo_rates/module_rxn.F90
              photo_rates/tuv_photolysis.F90
              rad_xfer/rad_abs_xsect.F90
              rad_xfer/rad_trans.F90
              rad_xfer/tuv_radiation_transfer.F90
              rad_xfer/tuv_subs.F90
              share/module_xsections.F90
              share/numer.F90
              share/params_mod.F90
              share/phot_kind_mod.F90
              share/phot_util_mod.F90
              share/wavelength_grid.F90)

add_library(musica-photolib ${PHOTO_SRC})

target_link_libraries(musica-photolib ${NETCDF_LIBS})

set_target_properties(musica-photolib PROPERTIES OUTPUT_NAME musica-photo)

################################################################################
