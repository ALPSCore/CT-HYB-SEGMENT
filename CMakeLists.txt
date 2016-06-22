cmake_minimum_required (VERSION 2.8.12)
# Not strictly necessary, but will create executables that are
# aware of ALPSCore location
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(alps_cthyb)

# Disable in-source builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Request the ALPSCore package:
find_package(ALPSCore COMPONENTS accumulators params mc hdf5)
find_package(LAPACK)
find_package(MPI)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    ${MPI_CXX_INCLUDE_PATH}
)

set (ALPS_CTHYB_SOURCE
   hybsim.cpp
   hyblocal.cpp
   hybint.cpp
   hybfun.cpp
   hybretintfun.cpp
   hybmatrix.cpp
   hybmatrix_ft.cpp
   hybconfig.cpp
   hybupdates.cpp
   hybevaluate.cpp
   hybmeasurements.cpp)
add_library(alps_hyb_exp ${ALPS_CTHYB_SOURCE})
target_link_libraries(alps_hyb_exp ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES})
set_target_properties(alps_hyb_exp PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(EXEC_NAME alps_cthyb)
add_executable(${EXEC_NAME} hybmain.cpp)
target_link_libraries(${EXEC_NAME} alps_hyb_exp ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES})
install (TARGETS ${EXEC_NAME} DESTINATION bin )

### Additional installation files
# configuration file for module managers, such as LMOD
configure_file("${CMAKE_SOURCE_DIR}/misc/alps_cthyb.lmod.in" "${CMAKE_BINARY_DIR}/alps_cthyb.lmod")
