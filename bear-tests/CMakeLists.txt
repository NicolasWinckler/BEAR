

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/bear-utils/src
${CMAKE_SOURCE_DIR}/bear-equations/src
${CMAKE_SOURCE_DIR}/bear-matrix-operations/src
${CMAKE_SOURCE_DIR}/bear-tests/policy-impl
)

Set(SYSTEM_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${Boost_LIBRARY_DIRS})

link_directories(${LINK_DIRECTORIES})


Set(EXE_NAME runSolveSysAtEqui)
Set(SRCS run/runSolveSystemAtEquilibrium.cxx)
Set(DEPENDENCIES bear_utils blas lapack gfortran)
GENERATE_EXECUTABLE()