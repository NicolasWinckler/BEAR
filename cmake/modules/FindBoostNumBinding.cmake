# - Try to find the BNB library.
#
# The following are set after configuration is done: 
#  BNB_FOUND
#  BNB_INCLUDE_DIRS
#  BNB_LIBRARY_DIRS
#  BNB_LIBRARIES

message(STATUS "Looking for Boost numeric bindings...")
set(BNB_H syevd.hpp)
find_path(BNB_INCLUDE_DIR NAMES ${BNB_H}
  PATHS ${SIMPATH}/include/boost/numeric/bindings/lapack
  NO_DEFAULT_PATH
  DOC   "Path to Boost numeric bindings include header files."
)
#find_library(BNB_LIBRARY NAMES libbnb bnb)

#message("BNB include dir = ${BNB_INCLUDE_DIR}")
#message("BNB lib = ${BNB_LIBRARY}")

set(BNB_LIBRARIES ${BNB_LIBRARY})
set(BNB_INCLUDE_DIRS ${BNB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set the BNB_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BNB DEFAULT_MSG
                                  #BNB_LIBRARY 
                                  BNB_INCLUDE_DIR)



mark_as_advanced(BNB_INCLUDE_DIR BNB_LIBRARY)

if(BNB_FOUND)
    message(STATUS "Looking for boost numeric bindings... - found")
endif(BNB_FOUND)
