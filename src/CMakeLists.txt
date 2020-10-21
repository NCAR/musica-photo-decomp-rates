################################################################################
# musica-photo-decomp-rates library

set(PHOTO_SRC
    phot_specs_mod.F90
    photo_utils.F90
    etf_mod.F90
    dissociation_class.F90
    o3_dissociation_class.F90
    jrates_mod.F90
   )

add_library(musica-photolib ${PHOTO_SRC})

target_link_libraries(musica-photolib ${NETCDF_LIBS})

set_target_properties(musica-photolib PROPERTIES OUTPUT_NAME musica-photo)

################################################################################
