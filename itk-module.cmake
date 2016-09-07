get_filename_component( MY_CURENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file( READ "${MY_CURENT_DIR}/README.rst" DOCUMENTATION )


# ITK version 4.5 changed it from EXCLUDE_FROM_ALL to EXCLUDE_FROM_DEFAULT
set( _EXCLUDE "EXCLUDE_FROM_ALL" )
if (NOT "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_MINOR_PATCH}" VERSION_LESS "4.5")
  set( _EXCLUDE "EXCLUDE_FROM_DEFAULT" )
endif()

# itk_module() defines the module dependencies in SimpleITKFiltersModule
# SimpleITKFiltersModule depends on ITKCommon
# The testing module in SimpleITKFiltersModule depends on ITKTestKernel
# and ITKMetaIO(besides SimpleITKFiltersModule and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(SimpleITKFiltersModule
  DEPENDS
    ITKImageFeature
  TEST_DEPENDS
    ITKTestKernel
    ITKImageSources
  DESCRIPTION
    "${DOCUMENTATION}"
  ${_EXCLUDE}
)
