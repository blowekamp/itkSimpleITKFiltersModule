#
# Example Dashboard script for CI services
#
# Used environment variable:
#  ITK_MODULE_NAME - name of the ITK external/remote module
#  ITK_SRC - ITK source directory ( already upto date )
#  ITK_REPOSITORY - local or remote repository
#  ITK_TAG - name of ITK tag or branch
set(itk_module "$ENV{ITK_MODULE_NAME}")

set(CTEST_BUILD_NAME "${itk_module}")
if("$ENV{TRAVIS}" STREQUAL "true" )
  set(CTEST_SITE "travisci-$ENV{TRAVIS_BUILD_ID}")
  string(SUBSTRING "$ENV{TRAVIS_COMMIT}" 0 6 short_sha1)
  set(git_moniker "-")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-$ENV{TRAVIS_OS_NAME}-$ENV{TRAVIS_JOB_NUMBER}-${short_sha1}")
endif()


set(CTEST_BUILD_CONFIGURATION "Debug")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(dashboard_no_update 0)
if(DEFINED ENV{ITK_SRC})
  set(CTEST_SOURCE_DIRECTORY "$ENV{ITK_SRC}")
  set(dashboard_no_update 1)
elseif(DEFINED ENV{ITK_REPOSITORY})
  set(dashboard_git_url "$ENV{ITK_REPOSITORY}")
  set(dashboard_git_branch "$ENV{ITK_TAG}")
  set(CTEST_DASHBOARD_ROOT "$ENV{HOME}/dash")
endif()

set(CTEST_BUILD_TARGET "${itk_module}-all")
set(CTEST_TEST_ARGS INCLUDE_LABEL ${itk_module})

set(dashboard_model "Experimental")
set(dashboard_track "Remote")

include( ProcessorCount )
ProcessorCount( PROCESSOR_COUNT )
if(PROCESSOR_COUNT)
  set( CTEST_BUILD_FLAGS -j${PROCESSOR_COUNT})
  set( CTEST_TEST_ARGS ${CTEST_TEST_ARGS} PARALLEL_LEVEL ${PROCESSOR_COUNT} )
endif()


# this is the initial cache to use for the binary tree.
SET (dashboard_cache "

    BUILD_DOCUMENTATION:BOOL=OFF
    BUILD_EXAMPLES:BOOL=OFF
    BUILD_SHARED_LIBS:BOOL=OFF
    BUILD_TESTING:BOOL=ON
    ITK_USE_KWSTYLE:BOOL=OFF

    ITK_BUILD_DEFAULT_MODULES:BOOL=OFF
    Module_${itk_module}:BOOL=ON
" )

include(${CTEST_SCRIPT_DIRECTORY}/../dashboard/itk_common.cmake)

