include(ExternalProject)

# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

set(proj GTest)

set(GTEST_TARGET_VERSION 1.7.0)
set(GTEST_DOWNLOAD_SOURCE_HASH "2d6ec8ccdf5c46b05ba54a9fd1d130d7")

# follow the standard EP_PREFIX locations
set(GTEST_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-build)
set(GTEST_source_dir ${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj})
set(GTEST_install_dir ${CMAKE_CURRENT_BINARY_DIR}/${proj})

set(${proj}_ARCHIVE_OUTPUT_DIRECTORY "<BINARY_DIR>/lib")
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  set(${proj}_ARCHIVE_OUTPUT_DIRECTORY "<BINARY_DIR>/lib/$<CONFIG>")
endif()

#
# Function which converts a list a cmake cache variable into a list of
# "-Dvar:type=value;" suitable for command line initialization.
#
function( VariableListToArgs var_list args )
  foreach( var IN LISTS ${var_list} )
    if( DEFINED ${var} AND NOT ${var} STREQUAL "" ) # if variable has been set
      get_property( type CACHE ${var} PROPERTY TYPE )
      if (NOT "${type}" STREQUAL "")
        set(type ":${type}")
      else()
        set(type ":UNINITIALIZED")
      endif()
      set(value ${${var}})
      STRING( REPLACE ";" "$<SEMICOLON>" value "${value}" )
      list( APPEND _args "-D${var}${type}=${value}" )
    endif()
  endforeach()
  set( ${args} "${_args}" PARENT_SCOPE)
endfunction( )

list( APPEND ep_common_list
  CMAKE_BUILD_TYPE
  CMAKE_MAKE_PROGRAM

  CMAKE_C_COMPILER
  CMAKE_C_COMPILER_ARG1

  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEBINFO

  CMAKE_CXX_COMPILER
  CMAKE_CXX_COMPILER_ARG1

  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_RELWITHDEBINFO

  CMAKE_LINKER

  CMAKE_EXE_LINKER_FLAGS
  CMAKE_EXE_LINKER_FLAGS_DEBUG
  CMAKE_EXE_LINKER_FLAGS_MINSIZEREL
  CMAKE_EXE_LINKER_FLAGS_RELEASE
  CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_MODULE_LINKER_FLAGS
  CMAKE_MODULE_LINKER_FLAGS_DEBUG
  CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL
  CMAKE_MODULE_LINKER_FLAGS_RELEASE
  CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_DEBUG
  CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO

  CMAKE_PREFIX_PATH
  CMAKE_FRAMEWORK_PATH
  CMAKE_SYSTEM_PREFIX_PATH
  CMAKE_SYSTEM_INCLUDE_PATH
  CMAKE_SYSTEM_LIBRARY_PATH
  CMAKE_SYSTEM_PROGRAM_PATH
  CMAKE_SYSTEM_IGNORE_PATH

  CMAKE_GENERATOR
  CMAKE_EXTRA_GENERATOR
 )

if( APPLE )
  list( APPEND ep_common_list
    CMAKE_OSX_SYSROOT
    CMAKE_OSX_DEPLOYMENT_TARGET
    CMAKE_OSX_ARCHITECTURES )
endif()


VariableListToArgs( ep_common_list ep_common_args )

set(ep_extra_args)
if(MSVC_VERSION EQUAL 1700)
  # Tuples are limited by _VARIADIC_MAX in VS11. The variadic
  # templates are not deep enough by default. We are not currently
  # using the GTest features which require tuple, so just disable them
  # and hope that upstream premanetly addresses the problem, with out
  # required more CMake core for compiler issues.
  set(ep_extra_args ${ep_extra_args} -DCMAKE_CXX_FLAGS=-DGTEST_HAS_TR1_TUPLE=0 ${CMAKE_CXX_FLAGS})
endif()

if(MSVC)
  set(ep_extra_args ${ep_extra_args} -Dgtest_force_shared_crt:BOOL=ON)
endif()


ExternalProject_Add(${proj}
  URL "https://midas3.kitware.com/midas/download/item/318120/gtest-1.7.0.zip"
  URL_MD5 ${GTEST_DOWNLOAD_SOURCE_HASH}
  INSTALL_DIR ${GTEST_install_dir}
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    --no-warn-unused-cli
  CMAKE_CACHE_ARGS
    ${ep_common_args}
    ${ep_extra_args}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=<BINARY_DIR>/lib
  INSTALL_COMMAND
      ${CMAKE_COMMAND} -E copy_directory ${${proj}_ARCHIVE_OUTPUT_DIRECTORY} <INSTALL_DIR>/lib
    COMMAND
      ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include
)


set(GTEST_ROOT ${GTEST_install_dir})
set(GTEST_INCLUDE_DIRS "${GTEST_install_dir}/include")
set(GTEST_LIBRARIES_DIR "${GTEST_install_dir}/lib/")
set(GTEST_LIBRARIES "gtest")
set(GTEST_MAIN_LIBRARIES "gtest_main")
set(GTEST_BOTH_LIBRARIES ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})

