################################################################################
# Copy needed data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# macro to link the photo rate library and add the test to the suite

macro(add_unit_test test_name)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} musica-photolib musica-rad-xfer-lib musica-photo-utils-lib)
  add_test(NAME ${test_name} COMMAND ${test_name})
endmacro(add_unit_test)

################################################################################
# Add all the tests to the test suite

add_subdirectory(integration)

################################################################################
# shared module tests

add_executable(share_numer share/numer_test.F90)
add_unit_test(share_numer)

add_executable(empty_atmos empty_atmos/test.F90)
add_unit_test(empty_atmos)

################################################################################
