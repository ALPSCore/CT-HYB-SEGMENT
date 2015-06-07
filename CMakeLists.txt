# Not strictly necessary, but will create executables that are
# aware of ALPSCore location
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

cmake_minimum_required (VERSION 2.8.12)

project(alps-hybridization)

# Disable in-source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Request the ALPSCore package:
#find_package(ALPSCore)
# (optionally, components can be also listed, like this:) 
find_package(ALPSCore COMPONENTS accumulators params)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
)

set (QMC_HYBRIDIZATION_SOURCE
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
add_library(alps_hyb_exp STATIC ${QMC_HYBRIDIZATION_SOURCE})
target_link_libraries(alps_hyb_exp ${ALPSCore_LIBRARIES})

# fixme
#add_executable(hybridization hybmain.cpp ${QMC_HYBRIDIZATION_SOURCE})

# Use ALPSCore_LIBRARIES variable to link to ALPSCore 
#target_link_libraries(hybridization ${ALPSCore_LIBRARIES})
