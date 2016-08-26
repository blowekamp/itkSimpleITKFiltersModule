itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}Tests
  itkObjectnessMeasureImageFilterTest.cxx
)


set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}Tests}")

itk_add_test(NAME ObjectnessMeasureImageFilterTest1
  COMMAND ${itk-module}TestDriver itkObjectnessMeasureImageFilterTest
    DATA{Input/itkHessianToObjectnessMeasureImageFilterTest.mha}
    ${TEMP}/ObjectnessMeasureImageFilterTest1.nii
    )
