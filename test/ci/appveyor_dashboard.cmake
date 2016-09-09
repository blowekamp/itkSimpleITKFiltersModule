# Client maintainer: blowekamp@mail.nih.gov
set(CTEST_SITE "appveyor")
set(CTEST_DASHBOARD_ROOT $ENV{APPVEYOR_BUILD_FOLDER}/..)
string(SUBSTRING $ENV{APPVEYOR_REPO_COMMIT} 0 7 commit)

# Extract major/minor/patch  versions

set(what "$ENV{APPVEYOR_PULL_REQUEST_TITLE}_#$ENV{APPVEYOR_PULL_REQUEST_NUMBER}")
if("$ENV{APPVEYOR_PULL_REQUEST_NUMBER}" STREQUAL "")
  set(what "$ENV{APPVEYOR_REPO_BRANCH}")
endif()
set(CTEST_CONFIGURATION_TYPE $ENV{CONFIGURATION})
set(CTEST_CMAKE_GENERATOR "$ENV{GENERATOR}")
if("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  set(CTEST_CMAKE_GENERATOR "Visual Studio 9 2008")
endif()

#set(CTEST_BUILD_NAME "VS-${platform}-$ENV{CONFIGURATION}_${what}_${commit}")
set(CTEST_BUILD_FLAGS "")


set(dashboard_binary_name "python-cmake-buildsystem/build")
set(dashboard_model Experimental)
set(dashboard_track AppVeyor-CI)
