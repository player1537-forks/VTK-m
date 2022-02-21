//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_particleadvection_ParticleAdvectionAlgorithm_h
#define vtk_m_filter_particleadvection_ParticleAdvectionAlgorithm_h

#include <vtkm/filter/particleadvection/AdvectorBaseAlgorithm.h>
#include <vtkm/filter/particleadvection/AdvectorBaseThreadedAlgorithm.h>
#include <vtkm/filter/particleadvection/DataSetIntegrator.h>

namespace vtkm
{
namespace filter
{
namespace particleadvection
{

using DSIType = vtkm::filter::particleadvection::DataSetIntegrator;
using TDSIType = vtkm::filter::particleadvection::TemporalDataSetIntegrator;

//-------------------------------------------------------------------------------------------
//Steady state advection algorithms

//Particle advection
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT ParticleAdvectionAlgorithm
  : public AdvectorBaseAlgorithm<DSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>
{
public:
  ParticleAdvectionAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                             const std::vector<DSIType>& blocks)
    : AdvectorBaseAlgorithm<DSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>(bm,
                                                                                           blocks)
  {
  }
};

//Threaded particle advection
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT ParticleAdvectionThreadedAlgorithm
  : public AdvectorBaseThreadedAlgorithm<DSIType,
                                         vtkm::worklet::ParticleAdvectionResult,
                                         ParticleType>
{
public:
  ParticleAdvectionThreadedAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                                     const std::vector<DSIType>& blocks)
    : AdvectorBaseThreadedAlgorithm<DSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>(
        bm,
        blocks)
  {
  }
};

//Streamline
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT StreamlineAlgorithm
  : public AdvectorBaseAlgorithm<DSIType, vtkm::worklet::StreamlineResult, ParticleType>
{
public:
  StreamlineAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                      const std::vector<DSIType>& blocks)
    : AdvectorBaseAlgorithm<DSIType, vtkm::worklet::StreamlineResult, ParticleType>(bm, blocks)
  {
  }
};

//Threaded streamline
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT StreamlineThreadedAlgorithm
  : public AdvectorBaseThreadedAlgorithm<DSIType, vtkm::worklet::StreamlineResult, ParticleType>
{
public:
  StreamlineThreadedAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                              const std::vector<DSIType>& blocks)
    : AdvectorBaseThreadedAlgorithm<DSIType, vtkm::worklet::StreamlineResult, ParticleType>(bm,
                                                                                            blocks)
  {
  }
};

//-------------------------------------------------------------------------------------------
//Unsteady state advection algorithms

//PathParticle
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT PathParticleAlgorithm
  : public AdvectorBaseAlgorithm<TDSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>
{
public:
  PathParticleAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                        const std::vector<TDSIType>& blocks)
    : AdvectorBaseAlgorithm<TDSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>(bm,
                                                                                            blocks)
  {
  }
};


//Threaded PathParticle
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT PathParticleThreadedAlgorithm
  : public AdvectorBaseThreadedAlgorithm<TDSIType,
                                         vtkm::worklet::ParticleAdvectionResult,
                                         ParticleType>
{
public:
  PathParticleThreadedAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                                const std::vector<TDSIType>& blocks)
    : AdvectorBaseThreadedAlgorithm<TDSIType, vtkm::worklet::ParticleAdvectionResult, ParticleType>(
        bm,
        blocks)
  {
  }
};

//Pathline
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT PathlineAlgorithm
  : public AdvectorBaseAlgorithm<TDSIType, vtkm::worklet::StreamlineResult, ParticleType>
{
public:
  PathlineAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                    const std::vector<TDSIType>& blocks)
    : AdvectorBaseAlgorithm<TDSIType, vtkm::worklet::StreamlineResult, ParticleType>(bm, blocks)
  {
  }
};

//Threaded pathline
template <typename ParticleType>
class VTKM_ALWAYS_EXPORT PathlineThreadedAlgorithm
  : public AdvectorBaseThreadedAlgorithm<TDSIType, vtkm::worklet::StreamlineResult, ParticleType>
{
public:
  PathlineThreadedAlgorithm(const vtkm::filter::particleadvection::BoundsMap& bm,
                            const std::vector<TDSIType>& blocks)
    : AdvectorBaseThreadedAlgorithm<TDSIType, vtkm::worklet::StreamlineResult, ParticleType>(bm,
                                                                                             blocks)
  {
  }
};

}
}
} // namespace vtkm::filter::particleadvection

#endif //vtk_m_filter_particleadvection_ParticleAdvectionAlgorithm_h
