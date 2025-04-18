######################################################################
## CMakeLists.txt --- gsKLShell
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Hugo Verhelst
## Copyright (C) 2023
######################################################################

## gsKLShell module
project(gsKLShell)

# Apply G+Smo config
include(gsConfig)

## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules" )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" )

# add filedata folder
add_definitions(-D${PROJECT_NAME}_DATA="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

# add example files
if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(examples)
  add_subdirectory(tutorials)
else()
  add_subdirectory(examples EXCLUDE_FROM_ALL)
  add_subdirectory(tutorials EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)