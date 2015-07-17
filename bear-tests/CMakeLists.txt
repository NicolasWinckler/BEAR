
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/solveSteadyEq.sh.in ${CMAKE_BINARY_DIR}/bin/solveSteadyEq.sh )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/config/bear-config-steady.cfg.in ${CMAKE_SOURCE_DIR}/data/config/bear-config-steady.cfg )
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/solveSteadyEq.sh")

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/bear-utils/src
${CMAKE_SOURCE_DIR}/bear-equations/src
${CMAKE_SOURCE_DIR}/bear-matrix-operations/src
${CMAKE_SOURCE_DIR}/bear-tests/policy-impl
)

Set(SYSTEM_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})

if(LAPACK_FOUND AND BNB_FOUND)
  set(INCLUDE_DIRECTORIES 
    ${BNB_INCLUDE_DIR}
    ${INCLUDE_DIRECTORIES}
  )

endif(LAPACK_FOUND AND BNB_FOUND)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${Boost_LIBRARY_DIRS})

link_directories(${LINK_DIRECTORIES})


Set(EXE_NAME runSolveSteadyEq)
Set(SRCS run/runSolveSystemAtEquilibrium.cxx)
Set(DEPENDENCIES bear_utils blas lapack gfortran)
GENERATE_EXECUTABLE()

if(LAPACK_FOUND AND BNB_FOUND)
  Set(EXE_NAME runSolveSteadyEqLapack)
  Set(SRCS 
    run/runSolveSystemAtEquilibriumLapack.cxx
    ${SRC}
  )
  #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  #set(LIBS_MAP_LAPACK "LAPACK")
  Set(DEPENDENCIES bear_utils blas lapack gfortran ${LAPACK_LIBRARIES})
  GENERATE_EXECUTABLE()

endif(LAPACK_FOUND AND BNB_FOUND)

