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
#ifndef itkSLICImageFilter_hxx
#define itkSLICImageFilter_hxx

#include "itkSLICImageFilter.h"


#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include <numeric>
#include <functional>


namespace itk
{

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SLICImageFilter()
  : m_MaximumNumberOfIterations( (ImageDimension > 2) ? 5 : 10),
    m_SpatialProximityWeight( 10.0 ),
    m_LabelConnectivityEnforce(true),
    m_LabelConnectivityMinimumSize(0.25),
    m_LabelConnectivityRelabelSequential(false),
    m_NumberOfThreadsUsed(1),
    m_Barrier(Barrier::New())
{
  m_SuperGridSize.Fill(50);
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::~SLICImageFilter()
{
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SetSuperGridSize(unsigned int factor)
{
  unsigned int i;
  for (i = 0; i < ImageDimension; ++i)
    {
    if (factor != m_SuperGridSize[i])
      {
      break;
      }
    }
  if ( i < ImageDimension )
    {
    this->Modified();
    m_SuperGridSize.Fill(factor);
    }
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SetSuperGridSize(unsigned int i, unsigned int factor)
{
  if (m_SuperGridSize[i] == factor)
    {
    return;
    }

  this->Modified();
  m_SuperGridSize[i] = factor;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "SuperGridSize: " << m_SuperGridSize << std::endl;
  os << indent << "MaximumNumberOfIterations: " << m_MaximumNumberOfIterations << std::endl;
  os << indent << "SpatialProximityWeight: " << m_SpatialProximityWeight << std::endl;
  os << indent << "LabelConnectivityEnforce: " << m_LabelConnectivityEnforce << std::endl;
  os << indent << "LabelConnectivityMinimumSize: " << m_LabelConnectivityMinimumSize << std::endl;
  os << indent << "LabelConnectivityRelabelSequential: " << m_LabelConnectivityRelabelSequential << std::endl;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::VerifyInputInformation ()
{
  Superclass::VerifyInputInformation();

  const InputImageType *inputImage = this->GetInput();

  typename InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();

  size_t numberOfClusters = 1u;

  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    numberOfClusters *= Math::Ceil<size_t>( double(size[i])/m_SuperGridSize[i] );
    }

  if (numberOfClusters >= static_cast<size_t>(itk::NumericTraits<LabelPixelType>::max()))
    {
    itkExceptionMacro( "Too many clusters for output pixel type!" );
    }

}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::BeforeThreadedGenerateData()
{
  itkDebugMacro("Starting BeforeThreadedGenerateData");


  // Compute actual number of threads used
  {
  ThreadIdType numberOfThreads = this->GetNumberOfThreads();

  if ( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 )
    {
    numberOfThreads = std::min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
    }

  typename TOutputImage::RegionType splitRegion;  // dummy region - just to call the following method

  m_NumberOfThreadsUsed = this->SplitRequestedRegion(0, numberOfThreads, splitRegion);
  }

  m_Barrier->Initialize(m_NumberOfThreadsUsed);

  const InputImageType *inputImage = this->GetInput();


  itkDebugMacro("Initializing Clusters");


  typename InputImageType::SizeType  strips, size, totalErr, accErr;
  typename InputImageType::IndexType startIdx, idx;

  typename InputImageType::RegionType region = inputImage->GetLargestPossibleRegion();

  size = region.GetSize();

  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    // number of super pixels
    strips[i] = size[i]/m_SuperGridSize[i];

    // the remainder of the pixels
    totalErr[i] = size[i]%m_SuperGridSize[i];

    // the starting superpixel index
    startIdx[i] = region.GetIndex()[i]+m_SuperGridSize[i]/2 + totalErr[i]/(strips[i]*2);
    idx[i] = startIdx[i];

    // with integer math keep track of the remaining odd pixel.
    // accErr/strips is the fractional pixels missing per superpixel
    // from even division.
    accErr[i] = totalErr[i]%(strips[i]*2);
    }


  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;
  const size_t numberOfClusters =  std::accumulate( &strips.m_Size[0], &strips.m_Size[0]+ImageDimension, size_t(1), std::multiplies<size_t>() );

   itkDebugMacro("numberOfClusters: " << numberOfClusters );

  // allocate array of scalars
  std::vector<ClusterComponentType>(numberOfClusters*numberOfClusterComponents).swap(m_Clusters);
  std::vector<ClusterComponentType>(numberOfClusters*numberOfClusterComponents).swap(m_OldClusters);


  size_t cnt = 0;
  while( idx[ImageDimension-1] < region.GetUpperIndex()[ImageDimension-1]  )
    {


    for( SizeValueType i = 0; i < strips[0] - 1; ++i )
      {
      RefClusterType cluster( numberOfClusterComponents, &m_Clusters[cnt*numberOfClusterComponents] );
      ++cnt;

      CreateClusterPoint(inputImage->GetPixel(idx),
                         cluster,
                         numberOfComponents,
                         inputImage,
                         idx );
      itkDebugMacro("Initial cluster " << cnt-1 << " : " << cluster << " idx: " << idx );
      accErr[0] += totalErr[0];
      idx[0] += m_SuperGridSize[0] + accErr[0]/strips[0];
      accErr[0] %= strips[0];
      }

      RefClusterType cluster( numberOfClusterComponents, &m_Clusters[cnt*numberOfClusterComponents] );
      ++cnt;

      CreateClusterPoint(inputImage->GetPixel(idx),
                         cluster,
                         numberOfComponents,
                         inputImage,
                         idx );
      itkDebugMacro("Initial cluster " << cnt-1<< " : " << cluster );

    // increment the startIdx to next line on sample grid
    idx[0] = startIdx[0];
    accErr[0] = totalErr[0]%(strips[0]*2);
    for ( unsigned int i = 1; i < ImageDimension; ++i )
      {

      accErr[i] += totalErr[i];
      idx[i] += m_SuperGridSize[0] + accErr[i]/strips[i];
      accErr[i] %= strips[i];

      if (idx[i] < region.GetUpperIndex()[i]
          || i == ImageDimension-1)
        {
        break;
        }
      idx[i] = startIdx[i];
      accErr[i] = totalErr[i]%(strips[i]*2);
      }
    }

  itkDebugMacro("Initial Clustering Completed");


  m_DistanceImage = DistanceImageType::New();
  m_DistanceImage->CopyInformation(inputImage);
  m_DistanceImage->SetBufferedRegion( region );
  m_DistanceImage->Allocate();

  m_MarkerImage = ITK_NULLPTR;

  const typename InputImageType::SpacingType spacing = inputImage->GetSpacing();

  for (unsigned int i = 0; i < ImageDimension; ++i)
    {
    const double physicalGridSize = m_SuperGridSize[i]*spacing[i];
    m_DistanceScales[i] = 1.0/physicalGridSize;
    }


  m_UpdateClusterPerThread.resize(m_NumberOfThreadsUsed);
  std::vector<std::list<LabelPixelType> >(m_NumberOfThreadsUsed).swap(m_MissedLabelsPerThread);

  this->Superclass::BeforeThreadedGenerateData();
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedUpdateDistanceAndLabel(const OutputImageRegionType & outputRegionForThread, ThreadIdType  itkNotUsed(threadId))
{
  // This method modifies the OutputImage and the DistanceImage only
  // in the outputRegionForThread. It searches for any cluster, whose
  // search radius is within the output region for the thread. Then it
  // updates DistnaceImage with the minimum distance and the
  // corresponding label id in the output image.
  //
  typedef ImageScanlineConstIterator< InputImageType > InputConstIteratorType;
  typedef ImageScanlineIterator< DistanceImageType >   DistanceIteratorType;

  const InputImageType *inputImage = this->GetInput();
  OutputImageType *outputImage = this->GetOutput();
  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;

  typename InputImageType::SizeType searchRadius;
  for (unsigned int i = 0; i < ImageDimension; ++i)
    {
    searchRadius[i] = m_SuperGridSize[i];
    }

  for (size_t i = 0; i*numberOfClusterComponents < m_Clusters.size(); ++i)
    {
    const RefClusterType cluster(numberOfClusterComponents, &m_Clusters[i*numberOfClusterComponents]);
    typename InputImageType::RegionType localRegion;
    typename InputImageType::PointType pt;
    IndexType idx;

    for (unsigned int d = 0; d < ImageDimension; ++d)
      {
      pt[d] = cluster[numberOfComponents+d];
      }

    inputImage->TransformPhysicalPointToIndex(pt, idx);

    localRegion.SetIndex(idx);
    localRegion.GetModifiableSize().Fill(1u);
    localRegion.PadByRadius(searchRadius);

    // Check cluster is in the output region for this thread.
    if (!localRegion.Crop(outputRegionForThread))
      {
      continue;
      }


    const size_t         ln =  localRegion.GetSize(0);

    InputConstIteratorType inputIter(inputImage, localRegion);
    DistanceIteratorType   distanceIter(m_DistanceImage, localRegion);


    while ( !inputIter.IsAtEnd() )
      {
      for( size_t x = 0; x < ln; ++x )
        {
        const IndexType &currentIdx = inputIter.GetIndex();

        inputImage->TransformIndexToPhysicalPoint(currentIdx, pt);
        const double distance = this->Distance(cluster,
                                               inputIter.Get(),
                                               pt );
        if (distance < distanceIter.Get() )
          {
          distanceIter.Set(distance);
          outputImage->SetPixel(currentIdx, static_cast<LabelPixelType>(i));
          }

        ++distanceIter;
        ++inputIter;
        }
      inputIter.NextLine();
      distanceIter.NextLine();
      }

    // for neighborhood iterator size S
    }

}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedUpdateClusters(const OutputImageRegionType & updateRegionForThread, ThreadIdType threadId)
{
  const InputImageType *inputImage = this->GetInput();
  OutputImageType *outputImage = this->GetOutput();

  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;

  typedef ImageScanlineConstIterator< InputImageType > InputConstIteratorType;
  typedef ImageScanlineIterator< OutputImageType >     OutputIteratorType;

  UpdateClusterMap &clusterMap = m_UpdateClusterPerThread[threadId];
  clusterMap.clear();

  itkDebugMacro("Estimating Centers");
  // calculate new centers
  OutputIteratorType itOut = OutputIteratorType(outputImage, updateRegionForThread);
  InputConstIteratorType itIn = InputConstIteratorType(inputImage, updateRegionForThread);
  ClusterType incr_cluster(numberOfClusterComponents);

  while(!itOut.IsAtEnd() )
    {
    const size_t ln =  updateRegionForThread.GetSize(0);
    for (size_t x = 0; x < ln; ++x)
      {
      const LabelPixelType l = itOut.Get();

      std::pair<typename UpdateClusterMap::iterator, bool> r =  clusterMap.insert(std::make_pair(l,UpdateCluster()));
      ClusterType &cluster = r.first->second.cluster;
      if (r.second)
        {
        cluster.set_size(numberOfClusterComponents);
        cluster.fill(0.0);
        r.first->second.count = 0;
        }
      ++r.first->second.count;

      // create cluster point
      CreateClusterPoint(itIn.Get(),
                         incr_cluster,
                         numberOfComponents,
                         inputImage,
                         itOut.GetIndex() );

      cluster += incr_cluster;

      ++itIn;
      ++itOut;
      }
    itIn.NextLine();
    itOut.NextLine();
    }
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedPerturbClusters(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId )
{
  // Update the m_Clusters array by spiting the threads over the
  // cluster indexes, moving cluster center to the
  // lowest gradient position in a 1-radius neighborhood.

  const InputImageType *inputImage = this->GetInput();

  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;
  const size_t numberOfClusters = m_Clusters.size()/numberOfClusterComponents;

  itk::Size<ImageDimension> radius;
  radius.Fill( 1 );
  unsigned long center;
  unsigned long stride[ImageDimension];


  typename InputImageType::SizeType searchRadius;
  searchRadius.Fill(1);


  typedef ConstNeighborhoodIterator< TInputImage > NeighborhoodType;

  // get center and dimension strides for iterator neighborhoods
  NeighborhoodType it( radius, inputImage, outputRegionForThread);
  center = it.Size()/2;
  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    stride[i] = it.GetStride(i);
    }

  const typename InputImageType::SpacingType spacing = inputImage->GetSpacing();

  // ceiling of number of clusters divided by actual number of threads
  const size_t strideCluster = 1 + ((numberOfClusters - 1) / m_NumberOfThreadsUsed);
  size_t clusterIndex = strideCluster*threadId;
  const size_t stopCluster = std::min(numberOfClusters, clusterIndex+strideCluster);

  for (; clusterIndex < stopCluster; ++clusterIndex)
    {
    // cluster is a reference to array
    RefClusterType cluster(numberOfClusterComponents, &m_Clusters[clusterIndex*numberOfClusterComponents]);
    typename InputImageType::RegionType localRegion;
    typename InputImageType::PointType pt;
    IndexType idx;

    for (unsigned int d = 0; d < ImageDimension; ++d)
      {
      pt[d] = cluster[numberOfComponents+d];
      }
    inputImage->TransformPhysicalPointToIndex(pt, idx);

    localRegion.SetIndex(idx);
    localRegion.GetModifiableSize().Fill(1u);
    localRegion.PadByRadius(searchRadius);

    it.SetRegion( localRegion );

    double minG = NumericTraits<double>::max();

    IndexType minIdx = idx;

    double gNorm;
    vnl_vector<ClusterComponentType> A(numberOfComponents), B(numberOfComponents);

    while ( !it.IsAtEnd() )
      {

      // Use the 1-norm of the Jacobian for the "Gradient Magnitude"
      // of single or multi-component images.
      gNorm = 0;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        NumericTraits<InputPixelType>::AssignToArray(it.GetPixel(center + stride[i]), A);
        NumericTraits<InputPixelType>::AssignToArray(it.GetPixel(center - stride[i]), B);
        A -= B;
        A /= spacing[i]; // omitting constant 2
        gNorm += A.one_norm();
        }

      if ( gNorm < minG)
        {
        minG = gNorm;
        minIdx = it.GetIndex();
        }
      ++it;
      }

    // create cluster point
     CreateClusterPoint(inputImage->GetPixel(minIdx),
                        cluster,
                        numberOfComponents,
                        inputImage,
                        minIdx );

    }

}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
{
  const InputImageType *inputImage = this->GetInput();

  const typename InputImageType::RegionType region = inputImage->GetBufferedRegion();
  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;

  itkDebugMacro("Perturb cluster centers");
  ThreadedPerturbClusters(outputRegionForThread,threadId);

  itkDebugMacro("Entering Main Loop");
  for(unsigned int loopCnt = 0;  loopCnt<m_MaximumNumberOfIterations; ++loopCnt)
    {

    if (threadId == 0)
      {
      itkDebugMacro("Iteration :" << loopCnt);
      m_DistanceImage->FillBuffer(NumericTraits<typename DistanceImageType::PixelType>::max());
      }
    m_Barrier->Wait();

    ThreadedUpdateDistanceAndLabel(outputRegionForThread,threadId);

    m_Barrier->Wait();


    ThreadedUpdateClusters(outputRegionForThread, threadId);

    m_Barrier->Wait();

    if (threadId==0)
      {
      // prepare to update clusters
      swap(m_Clusters, m_OldClusters);
      std::fill(m_Clusters.begin(), m_Clusters.end(), 0.0);
      std::vector<size_t> clusterCount(m_Clusters.size()/numberOfClusterComponents, 0);

      // reduce the produce cluster maps per-thread into m_Cluster array
      for(size_t i = 0; i < m_UpdateClusterPerThread.size(); ++i)
        {
        const UpdateClusterMap &clusterMap = m_UpdateClusterPerThread[i];
        for(typename UpdateClusterMap::const_iterator clusterIter = clusterMap.begin(); clusterIter != clusterMap.end(); ++clusterIter)
          {
          const size_t clusterIdx = clusterIter->first;
          clusterCount[clusterIdx] += clusterIter->second.count;

          RefClusterType cluster(numberOfClusterComponents, &m_Clusters[clusterIdx*numberOfClusterComponents]);
          cluster += clusterIter->second.cluster;
          }
        }

      // average
      for (size_t i = 0; i*numberOfClusterComponents < m_Clusters.size(); ++i)
        {
        if (clusterCount[i] != 0)
          {
          RefClusterType cluster(numberOfClusterComponents,&m_Clusters[i*numberOfClusterComponents]);
          cluster /= clusterCount[i];
          }
        }

     // residual
#if !defined NDEBUG
      double l1Residual = 0.0;
      for (size_t i = 0; i*numberOfClusterComponents < m_Clusters.size(); ++i)
        {

        RefClusterType cluster(numberOfClusterComponents,&m_Clusters[i*numberOfClusterComponents]);
        RefClusterType oldCluster(numberOfClusterComponents, &m_OldClusters[i*numberOfClusterComponents]);
        l1Residual += Distance(cluster,oldCluster);
        }
      itkDebugMacro( << "L1 residual: " << std::sqrt(l1Residual) );
#endif
      }
    // while error <= threshold
    }

  if(m_LabelConnectivityEnforce && !m_LabelConnectivityRelabelSequential)
    {

    m_Barrier->Wait();
    if (threadId == 0)
      {
      m_DistanceImage = ITK_NULLPTR;

      m_MarkerImage = MarkerImageType::New();
      m_MarkerImage->CopyInformation(inputImage);
      m_MarkerImage->SetBufferedRegion( region );
      m_MarkerImage->Allocate();
      m_MarkerImage->FillBuffer(-2.0);
      }
    m_Barrier->Wait();

    const size_t superGirdArea =
      std::accumulate( &m_SuperGridSize[0],
                       &m_SuperGridSize[0]+ImageDimension, size_t(1),
                       std::multiplies<size_t>() );

    // Split the clusters up over threads, updating m_Marker map
    // with visited regions, centered and connected with the cluster

    // ceiling of number of clusters divided by actual number of threads
    const size_t numberOfClusters = m_Clusters.size() / numberOfClusterComponents;
    const size_t stride = 1 + ((numberOfClusters - 1) / m_NumberOfThreadsUsed );


    std::list<LabelPixelType> &missedLabels = this->m_MissedLabelsPerThread[threadId];

    size_t i = stride*threadId;
    const size_t stopCluster = std::min(numberOfClusters, i+stride);

    for (; i < stopCluster; ++i)
      {
      RefClusterType cluster(numberOfClusterComponents,&m_Clusters[i*numberOfClusterComponents]);
      typename InputImageType::PointType pt;
      IndexType idx;

      for (unsigned int d = 0; d < ImageDimension; ++d)
        {
        pt[d] = cluster[numberOfComponents+d];
        }

      m_MarkerImage->TransformPhysicalPointToIndex(pt,idx);

      // count
      const size_t count = this->FillCluster(idx, i, -1);

      if ( count > m_LabelConnectivityMinimumSize*superGirdArea )
        {
        this->FillCluster(idx, i, 1, i);
        }
      else
        {
        missedLabels.push_back(i);
        }
      }

    }


}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
size_t
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::FillCluster(const IndexType &idx,
              size_t label,
              int relabel,
              LabelPixelType outLabel )
{

    OutputImageType *outputImage = this->GetOutput();

    std::vector<IndexType> idxStack;
    // todo reserve

    const typename OutputImageType::RegionType region = outputImage->GetBufferedRegion();

    size_t count = 0;
    IndexType tempIdx = idx;
    if ( m_MarkerImage->GetPixel(tempIdx) < relabel
         && outputImage->GetPixel(tempIdx) == label)
      {
      ++count;
      m_MarkerImage->SetPixel(tempIdx, relabel);
      idxStack.push_back(tempIdx);
      if (relabel > 0)
        {
        outputImage->SetPixel(tempIdx, outLabel);
        }
      }

    while( !idxStack.empty() )
      {
      tempIdx = idxStack.back();
      idxStack.pop_back();

      for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          for ( int j = -1; j <= 1; j += 2 )
            {
            tempIdx[i] += j;
            if (region.IsInside(tempIdx))
              {
               LabelPixelType &currLabel = outputImage->GetPixel(tempIdx);
              if (currLabel == label
                   && m_MarkerImage->GetPixel(tempIdx) < relabel )
                {
                ++count;
                m_MarkerImage->SetPixel(tempIdx, relabel);
                idxStack.push_back(tempIdx);
                if (relabel > 0)
                  {
                  currLabel = outLabel;
                  }
                }

              }
            tempIdx[i] -= j;
            }
          }
      }
    return count;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::AfterThreadedGenerateData()
{

  itkDebugMacro("Starting AfterThreadedGenerateData");

  m_DistanceImage = ITK_NULLPTR;

  const InputImageType *inputImage = this->GetInput();
  OutputImageType      *outputImage = this->GetOutput();

  const typename OutputImageType::RegionType region = outputImage->GetBufferedRegion();

  if (m_LabelConnectivityEnforce)
    {

    if (!m_MarkerImage)
      {
      m_MarkerImage = MarkerImageType::New();
      m_MarkerImage->CopyInformation(inputImage);
      m_MarkerImage->SetBufferedRegion( region );
      m_MarkerImage->Allocate();
      }

    const size_t superGirdArea =
      std::accumulate( &m_SuperGridSize[0],
                       &m_SuperGridSize[0]+ImageDimension, size_t(1),
                       std::multiplies<size_t>() );


    std::list<LabelPixelType> missedLabels;
    if (m_LabelConnectivityRelabelSequential)
      {
      missedLabels.push_back(0);
      m_MarkerImage->FillBuffer(-1.0);
      }
    else
      {
      // merge missing labels, from threaded initial labeling
      for (typename std::vector<std::list<LabelPixelType> >::iterator i =  m_MissedLabelsPerThread.begin();
           i != m_MissedLabelsPerThread.end();
           ++i)
        {
        missedLabels.splice(missedLabels.end(), *i);
        }
      const unsigned int   numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
      const unsigned int   numberOfClusterComponents = numberOfComponents+ImageDimension;
      const LabelPixelType numberOfClusters = static_cast<LabelPixelType>(m_Clusters.size() / numberOfClusterComponents);
      missedLabels.push_back(numberOfClusters);
      }

    typedef ImageScanlineIterator< OutputImageType >      OutputIteratorType;
    typedef ImageScanlineConstIterator< MarkerImageType > MarkerIteratorType;


    OutputIteratorType outputIter(outputImage, region);
    MarkerIteratorType markerIter(m_MarkerImage, region);
    const size_t       ln =  region.GetSize(0);

    // the next label to use for relabeling
    LabelPixelType nextLabel = missedLabels.front();
    missedLabels.pop_front();

    while ( !outputIter.IsAtEnd() )
      {
      for( size_t x = 0; x < ln; ++x )
        {
        if ( markerIter.Get() == 0 )
          {
          std::cout << "unexpected 0 mark @ " << markerIter.GetIndex() << " label: " << outputIter.Get() << std::endl;
          }
        if ( markerIter.Get() < 0)
          {
          const LabelPixelType currLabel = outputIter.Get();

          const size_t countLabel = this->FillCluster(outputIter.GetIndex(), currLabel, 0);

          if (countLabel < m_LabelConnectivityMinimumSize*superGirdArea)
            {

            LabelPixelType replaceLabel =  0;
            IndexType tempIdx = outputIter.GetIndex();
            for ( unsigned int i = 0; i < ImageDimension; i++ )
              {
              for ( int j = -1; j <= 1; j += 2 )
                {
                tempIdx[i] += j;
                 if (region.IsInside(tempIdx))
                   {
                   const LabelPixelType l = outputImage->GetPixel(tempIdx);
                   if (m_MarkerImage->GetPixel(tempIdx) > 0)
                     {
                     replaceLabel = l;
                     }
                   }
                 tempIdx[i] -= j;
                }
              }

            itkDebugMacro("Replacing island with neighbor: " << outputIter.Get() << "->" << replaceLabel<<  " " << outputIter.GetIndex());
            this->FillCluster(outputIter.GetIndex(), currLabel, 1, replaceLabel);
            }
          else
            {

            itkDebugMacro("Relabling big island of: " << m[outputIter.Get()] <<  " with new label: " << outputIter.Get() << "->" << nextLabel);
            this->FillCluster(outputIter.GetIndex(), currLabel, 1, nextLabel);

            if ( nextLabel != NumericTraits<LabelPixelType>::max() )
              {
              if (!missedLabels.empty())
                {
                nextLabel = missedLabels.front();
                missedLabels.pop_front();
                }
              else
                {
                ++nextLabel;
                }
              }
            }

          }
        ++outputIter;
        ++markerIter;
        }
      outputIter.NextLine();
      markerIter.NextLine();
      }
  }

  // Clean up all algorithm variables
  m_MarkerImage = ITK_NULLPTR;

  // cleanup
  std::vector<ClusterComponentType>().swap(m_Clusters);
  std::vector<ClusterComponentType>().swap(m_OldClusters);
  std::vector<std::list<LabelPixelType> >().swap(m_MissedLabelsPerThread);
  for(unsigned int i = 0; i < m_UpdateClusterPerThread.size(); ++i)
    {
    UpdateClusterMap().swap(m_UpdateClusterPerThread[i]);
    }
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
typename SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>::DistanceType
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::Distance(const ClusterType &cluster1, const ClusterType &cluster2)
{
  const unsigned int s = cluster1.size();
  DistanceType d1 = 0.0;
  DistanceType d2 = 0.0;
  unsigned int i = 0;
  for (; i<s-ImageDimension; ++i)
    {
    const DistanceType d = (cluster1[i] - cluster2[i]);
    d1 += d*d;
    }
  //d1 = std::sqrt(d1);

  for (unsigned int j = 0; j < ImageDimension; ++j, ++i)
    {
    const DistanceType d = (cluster1[i] - cluster2[i]) * m_DistanceScales[j];
    d2 += d*d;
    }
  d2 *= m_SpatialProximityWeight * m_SpatialProximityWeight;
  //d2 = std::sqrt(d2);
  return d1+d2;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
typename SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>::DistanceType
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::Distance(const ClusterType &cluster, const InputPixelType &v, const PointType &pt)
{
  const unsigned int s = cluster.size();
  DistanceType d1 = 0.0;
  DistanceType d2 = 0.0;
  unsigned int i = 0;
  for (; i<s-ImageDimension; ++i)
    {
    const DistanceType d = (cluster[i] - v[i]);
    d1 += d*d;
    }
  //d1 = std::sqrt(d1);

  for (unsigned int j = 0; j < ImageDimension; ++j, ++i)
    {
    const DistanceType d = (cluster[i] - pt[j])  * m_DistanceScales[j];
    d2 += d*d;
    }
  d2 *= m_SpatialProximityWeight * m_SpatialProximityWeight;
  //d2 = std::sqrt(d2);
  return d1+d2;
}

} // end namespace itk

#endif // itkSLICImageFilter_hxx
