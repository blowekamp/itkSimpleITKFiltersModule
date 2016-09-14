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

#include "itkSliceImageFilter.h"

#include "itkBitwiseNotFunctor.h"
#include "itkDivideFloorFunctor.h"
#include "itkDivideRealFunctor.h"

#include "gtest/gtest.h"


TEST(BitwiseNotFunctorTest, Test1)
{


  itk::Functor::BitwiseNot<unsigned char, unsigned char> f;
  itk::Functor::BitwiseNot<unsigned char, unsigned char> f2;


  EXPECT_FALSE(f != f2);
  EXPECT_TRUE(f == f2);
  EXPECT_TRUE(f == f);

  EXPECT_EQ(254, f(1));
  EXPECT_EQ(1, f(f(1)));


  itk::Functor::BitwiseNot<signed char, signed char> f3;

  EXPECT_EQ(0,f3(-1));
  EXPECT_EQ(-1,f3(f3(-1)));
}


TEST(DivideFloorFunctorTest, Test1)
{
  itk::Functor::DivFloor<int,int,int> f;

  EXPECT_TRUE(f==f);
  EXPECT_FALSE(f!=f);

  EXPECT_EQ(2, f(5,2));
  EXPECT_EQ(-2, f(8,-7));
}


TEST(DivideRealFunctorTest, Test1)
{

  itk::Functor::DivReal<int,int,double> f;

  EXPECT_TRUE(f==f);
  EXPECT_FALSE(f!=f);

  EXPECT_FLOAT_EQ(2.5, f(5,2));
}
