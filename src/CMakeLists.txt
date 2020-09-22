################################################################################
# musica-photo-decomp-rates library

set(PHOTO_SRC
    photo_rates/module_prates_tuv.F90
    photo_rates/module_rxn.F90
    photo_rates/tuv_photolysis.F90
   )

add_library(musica-photolib ${PHOTO_SRC})

target_link_libraries(musica-photolib ${NETCDF_LIBS})

set_target_properties(musica-photolib PROPERTIES OUTPUT_NAME musica-photo)

################################################################################
