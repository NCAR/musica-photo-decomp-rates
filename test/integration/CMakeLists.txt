################################################################################
# Integration smoke test

set(INTEGRATION_SRC ccpp_kinds.F90
                    driver.F90)

add_executable(integration_test ${INTEGRATION_SRC})

target_include_directories(integration_test PUBLIC ${CMAKE_BINARY_DIR}/src)

target_link_libraries(integration_test musica-photolib musica-rad-xfer-lib musica-photo-utils-lib)

add_test(NAME integration COMMAND integration_test)

add_test(NAME compare_out COMMAND compare_output.sh)

################################################################################
