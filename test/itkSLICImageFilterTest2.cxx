/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkSLICImageFilter.h"
#include "itkVectorImage.h"
#include "itkRandomImageSource.h"


namespace
{
}

int itkSLICImageFilterTest2(int, char *[])
{
  const unsigned int gridSize = 50;
  const float proximityWeight = 10.0;

  const unsigned int VDimension = 2;
  typedef itk::VectorImage<float, VDimension>  InputImageType;
  typedef itk::Image<unsigned int, VDimension> OutputImageType;

  InputImageType::Pointer input = InputImageType::New();

  InputImageType::RegionType region;
  InputImageType::SizeType size = {{gridSize+1, gridSize+1}};
  region.SetSize( size );
  input->SetRegions(region);
  input->SetVectorLength(3);
  input->Allocate();

  typedef itk::SLICImageFilter< InputImageType, OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetSpatialProximityWeight(proximityWeight);
  filter->SetInput(input);


  // check case where grid size is bigger than image
  filter->SetSuperGridSize(gridSize);
  filter->Update();

  // check case where more threads than slices
  input->Modified();
  filter->SetNumberOfThreads(125);
  filter->Update();

  // check 1x1 image
  InputImageType::SizeType size2 = {{1,1}};
  region.SetSize( size2 );
  input->SetRegions(region);
  input->Allocate();
  filter->Update();

  filter->GetOutput()->Print(std::cout);

  return 0;
}
