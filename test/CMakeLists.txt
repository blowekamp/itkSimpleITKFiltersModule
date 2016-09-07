itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}Tests
  itkObjectnessMeasureImageFilterTest.cxx
  itkHessianImageFilterTest.cxx
)


set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}Tests}")

itk_add_test(NAME itkObjectnessMeasureImageFilterTest1
  COMMAND ${itk-module}TestDriver
    --compareIntensityTolerance .001
    --compare
      DATA{Baseline/ObjectnessMeasureImageFilterTest1.mha}
      ${TEMP}/ObjectnessMeasureImageFilterTestOutput1.nii
    itkObjectnessMeasureImageFilterTest
      DATA{Input/DSA.png}
      ${TEMP}/ObjectnessMeasureImageFilterTestOutput1.nii
      1 0
    )

itk_add_test(NAME itkObjectnessMeasureImageFilterTest2
  COMMAND ${itk-module}TestDriver
    --compareIntensityTolerance .001
    --compare
      DATA{Baseline/ObjectnessMeasureImageFilterTest2.mha}
      ${TEMP}/ObjectnessMeasureImageFilterTestOutput2.nii
    itkObjectnessMeasureImageFilterTest
      DATA{Input/DSA.png}
      ${TEMP}/ObjectnessMeasureImageFilterTestOutput2.nii
      0 0
    )

add_test(NAME itkHessianImageFilterTest
      COMMAND ${itk-module}TestDriver itkHessianImageFilterTest )
