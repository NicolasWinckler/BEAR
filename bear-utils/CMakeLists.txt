 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# Create a library called "libbear_utils" 

set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/bear-utils/src)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${SYSTEM_INCLUDE_DIRECTORIES}
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${Boost_LIBRARY_DIRS}
) 

link_directories(${LINK_DIRECTORIES})

Set(SRCS
      src/options_manager.cxx
      src/logger.cxx
    )

Set(LIBRARY_NAME bear_utils)
Set(DEPENDENCIES boost_system boost_filesystem boost_program_options)

GENERATE_LIBRARY()
