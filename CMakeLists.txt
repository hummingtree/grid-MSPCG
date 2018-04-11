cmake_minimum_required(VERSION 2.8.0)
project (MeasurementPackage)


SET( CMAKE_C_COMPILER cc )
SET( CMAKE_CXX_COMPILER CC )

FILE ( GLOB SRC *.cc )

set ( GRID_VERSION "develop")

set ( GRID_DIR "/global/homes/j/jiquntu/grid-build" )

set ( QLAT_DIR "/global/homes/j/jiquntu/qlat-build/1.0" )

set ( FLAGS "-xmic-avx512 -O3 -mkl -std=c++11 -fopenmp -axMIC-AVX512,CORE-AVX2 -g" )

SET( CMAKE_CXX_FLAGS ${FLAGS} )
SET( CMAKE_C_FLAGS ${FLAGS} )

set ( LIBS libGrid.a )
link_directories( "${GRID_DIR}/${GRID_VERSION}/build/lib" )

include_directories( "${GRID_DIR}/${GRID_VERSION}/build/include" "${QLAT_DIR}/include" "${QLAT_DIR}/include/utils" "${QLAT_DIR}/include/eigen3")

add_executable ( knl.x ${SRC} )
target_link_libraries(knl.x libGrid.a )
