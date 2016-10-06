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
#ifndef itkSLICImageFilter_h
#define itkSLICImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkIsSame.h"

#include "itkBarrier.h"


#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageScanlineIterator.h"

#include "itkMath.h"

namespace itk
{

/** \class SLICImageFilter
 * \brief Simple Linear Iterative Clustering (SLIC) Superpixel algorithm
 *
 * The Simple Linear Iterative Clustering (SLIC) Superpixel performs
 * joint domain ( image intensity and physical location ) clustering
 * of the input image to form a superpixel labeled output image.
 *
 * This implementation is multi-threaded, works in n-dimensions with
 * m-component images. The filter works with VectorImage, scalar
 * Images, and Images of FixedArrays.
 *
 * R. Achanta, A. Shaji, K. Smith, and A. Lucchi. Slic superpixels. Technical report, 2010.
 *
 * \ingroup SimpleITKFiltersModule
 */
template< typename TInputImage, typename TOutputImage, typename TDistancePixel = float>
class SLICImageFilter:
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef SLICImageFilter                                 Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SLICImageFilter, ImageToImageFilter);

  /** Image type information. */
  typedef TInputImage                         InputImageType;
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType LabelPixelType;
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);
  typedef TDistancePixel                      DistanceType;
  typedef Image<DistanceType, ImageDimension> DistanceImageType;
  typedef Image<signed char, ImageDimension>  MarkerImageType;

  typedef typename InputImageType::IndexType IndexType;
  typedef typename InputImageType::PointType PointType;
  // assume variable length vector right now
  typedef double                               ClusterComponentType;
  typedef vnl_vector<ClusterComponentType>     ClusterType;
  typedef vnl_vector_ref<ClusterComponentType> RefClusterType;

  typedef typename OutputImageType::RegionType   OutputImageRegionType;

  typedef FixedArray< unsigned int, ImageDimension > SuperGridSizeType;

  itkSetMacro( SpatialProximityWeight, double );
  itkGetConstMacro( SpatialProximityWeight, double );

  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

  itkSetMacro(SuperGridSize, SuperGridSizeType);
  void SetSuperGridSize(unsigned int factor);
  void SetSuperGridSize(unsigned int i, unsigned int factor);

  itkSetMacro(LabelConnectivityEnforce, bool);
  itkGetMacro(LabelConnectivityEnforce, bool);
  itkBooleanMacro(LabelConnectivityEnforce);

  itkSetMacro(LabelConnectivityMinimumSize, float);
  itkGetMacro(LabelConnectivityMinimumSize, float);

  itkSetMacro(LabelConnectivityRelabelSequential, bool);
  itkGetMacro(LabelConnectivityRelabelSequential, bool);
  itkBooleanMacro(LabelConnectivityRelabelSequential);

protected:
  SLICImageFilter();
  ~SLICImageFilter();

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void VerifyInputInformation ()  ITK_OVERRIDE;

  /** Generate full output and require full input */
  void EnlargeOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;

  void ThreadedUpdateDistanceAndLabel(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedUpdateClusters(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedPerturbClusters(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) ITK_OVERRIDE;

  void AfterThreadedGenerateData() ITK_OVERRIDE;

  size_t FillCluster(const IndexType &idx, size_t label, int fill=0, LabelPixelType outLabel=0);

  DistanceType Distance(const ClusterType &cluster1,
                        const ClusterType &cluster2);

  DistanceType Distance(const ClusterType &cluster,
                        const InputPixelType &v,
                        const IndexType &idx)
    {
      return Self::DistanceDispatched(cluster, v, idx);
    }

  inline static void CreateClusterPoint( const InputPixelType &v,
                                         ClusterType &outCluster,
                                         const unsigned int numberOfComponents,
                                         const IndexType &idx )
    {
      NumericTraits<InputPixelType>::AssignToArray(v, outCluster);
      for(unsigned int i = 0; i < ImageDimension; ++i)
        {
        outCluster[numberOfComponents+i] = idx[i];
        }
    }

private:
  SLICImageFilter(const Self &);    //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  template<typename TPixelType>
  inline DistanceType DistanceDispatched(const ClusterType &cluster,
                                         const TPixelType &v,
                                         const IndexType &idx,
                                         mpl::FalseType isScalar = typename IsSame<TPixelType, typename itk::NumericTraits<InputPixelType>::ValueType>::Type() )
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

      for (unsigned int j = 0; j < ImageDimension; ++j, ++i)
        {
        const DistanceType d = (cluster[i] - idx[j])  * m_DistanceScales[j];
        d2 += d*d;
        }
      d2 *= m_SpatialProximityWeight * m_SpatialProximityWeight;

      return d1+d2;
    }


  inline DistanceType DistanceDispatched(const ClusterType &cluster,
                                         const typename itk::NumericTraits<InputPixelType>::ValueType &v,
                                         const IndexType &idx )
    {
      const unsigned int s = cluster.size();
      DistanceType d1 = 0.0;
      DistanceType d2 = 0.0;
      unsigned int i = 0;

      {
      const DistanceType d = (cluster[i] - v);
      d1 += d*d;
      ++i;
      }

      for (unsigned int j = 0; j < ImageDimension; ++j, ++i)
        {
        const DistanceType d = (cluster[i] - idx[j])  * m_DistanceScales[j];
        d2 += d*d;
        }
      d2 *= m_SpatialProximityWeight * m_SpatialProximityWeight;

      return d1+d2;
    }


  SuperGridSizeType m_SuperGridSize;
  unsigned int      m_MaximumNumberOfIterations;
  double            m_SpatialProximityWeight;
  bool              m_LabelConnectivityEnforce;
  float             m_LabelConnectivityMinimumSize;
  bool              m_LabelConnectivityRelabelSequential;

  FixedArray<double,ImageDimension> m_DistanceScales;
  std::vector<ClusterComponentType> m_Clusters;
  std::vector<ClusterComponentType> m_OldClusters;

  struct UpdateCluster
  {
    size_t count;
    vnl_vector<ClusterComponentType> cluster;
  };

  typedef std::map<LabelPixelType, UpdateCluster> UpdateClusterMap;

  std::vector<UpdateClusterMap> m_UpdateClusterPerThread;

  std::vector<std::list<LabelPixelType> > m_MissedLabelsPerThread;

  ThreadIdType m_NumberOfThreadsUsed;

  typename Barrier::Pointer           m_Barrier;
  typename DistanceImageType::Pointer m_DistanceImage;
  typename MarkerImageType::Pointer   m_MarkerImage;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSLICImageFilter.hxx"
#endif

#endif //itkSLICImageFilter_h
