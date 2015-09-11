

set(SCRIPT_CMD1 "\${filename##*.}")
set(SCRIPT_CMD2 "\${filename%.*}")

configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/solveSteadyEq.sh.in ${CMAKE_BINARY_DIR}/bin/solveSteadyEq.sh )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/solveDynamicEq.sh.in ${CMAKE_BINARY_DIR}/bin/solveDynamicEq.sh )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/config/bear-config-steady.cfg.in ${CMAKE_SOURCE_DIR}/data/config/bear-config-steady.cfg )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/solveDynamicEqRootAll.sh.in ${CMAKE_BINARY_DIR}/bin/solveDynamicEqRootAll.sh)


configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/bear-eq.sh.in ${CMAKE_BINARY_DIR}/bin/bear-eq.sh )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/macro/bear-plot.sh.in ${CMAKE_BINARY_DIR}/bin/bear-plot.sh )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/config/bear-config-save-to-file.cfg.in ${CMAKE_SOURCE_DIR}/data/config/bear-config-save-to-file.cfg )
configure_file( ${CMAKE_SOURCE_DIR}/bear-tests/config/bear-config-no-save-to-file.cfg.in ${CMAKE_SOURCE_DIR}/data/config/bear-config-no-save-to-file.cfg )


EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/bear-eq.sh")
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/bear-plot.sh")

EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/solveSteadyEq.sh")
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/solveDynamicEq.sh")
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/solveDynamicEqRootAll.sh")

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
  if(ROOT_FOUND)
        set(INCLUDE_DIRECTORIES 
            ${ROOT_INCLUDE_DIR}
            ${INCLUDE_DIRECTORIES}
           )
  endif(ROOT_FOUND)
endif(LAPACK_FOUND AND BNB_FOUND)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${Boost_LIBRARY_DIRS})
set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${LINK_DIRECTORIES}) 

link_directories(${LINK_DIRECTORIES})


Set(EXE_NAME runSolveSteadyEq)
Set(SRCS run/runSolveSystemAtEquilibrium.cxx)
Set(DEPENDENCIES bear_utils blas lapack gfortran)
GENERATE_EXECUTABLE()

if(LAPACK_FOUND AND BNB_FOUND)
  Set(EXE_NAME runSolveSteadyEqLapack)
  Set(SRCS 
    run/runSolveSystemAtEquilibriumLapack.cxx
  )
  #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  #set(LIBS_MAP_LAPACK "LAPACK")
  Set(DEPENDENCIES bear_utils blas lapack gfortran ${LAPACK_LIBRARIES})
  GENERATE_EXECUTABLE()

  Set(EXE_NAME rundiagoTest)
  Set(SRCS 
    run/test_diagonalization.cxx
  )
  #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  #set(LIBS_MAP_LAPACK "LAPACK")
  Set(DEPENDENCIES bear_utils blas lapack gfortran ${LAPACK_LIBRARIES})
  GENERATE_EXECUTABLE()


  ## ROOT GUI
  if(ROOT_FOUND)
  
      
      Set(EXE_NAME runSolveDynEqRoot)
      Set(SRCS 
        run/runSolveDynEqRoot.cxx
      )
      #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
      #set(LIBS_MAP_LAPACK "LAPACK")
      Set(DEPENDENCIES bear_utils blas lapack gfortran ${LAPACK_LIBRARIES}
         Hist Graf Gpad RIO Cint Core
         )
      GENERATE_EXECUTABLE()



      Set(EXE_NAME runSolveDynEqRootRK)
      Set(SRCS 
        run/runSolveDynEqRoot_RK.cxx
      )
      #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
      #set(LIBS_MAP_LAPACK "LAPACK")
      Set(DEPENDENCIES bear_utils blas lapack gfortran ${LAPACK_LIBRARIES}
         Hist Graf Gpad RIO Cint Core
         )
      GENERATE_EXECUTABLE()


  endif(ROOT_FOUND)

endif(LAPACK_FOUND AND BNB_FOUND)

