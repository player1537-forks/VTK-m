//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_DataSetIntegrator_hxx
#define vtk_m_filter_DataSetIntegrator_hxx

namespace vtkm
{
namespace filter
{
namespace particleadvection
{

using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<
  vtkm::worklet::particleadvection::VelocityField<vtkm::cont::ArrayHandle<vtkm::Vec3f>>>;
using TemporalGridEvalType = vtkm::worklet::particleadvection::TemporalGridEvaluator<
  vtkm::worklet::particleadvection::VelocityField<vtkm::cont::ArrayHandle<vtkm::Vec3f>>>;
using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;
using TemporalRK4Type = vtkm::worklet::particleadvection::RK4Integrator<TemporalGridEvalType>;
using TemporalStepper =
  vtkm::worklet::particleadvection::Stepper<TemporalRK4Type, TemporalGridEvalType>;

#define DEFINE_DOADVECT(GridEvalType, ParticleType, ResultType, WorkletType) \
  template <>                                                                \
  template <>                                                                \
  inline void DataSetIntegratorBase<GridEvalType>::DoAdvect(                 \
    vtkm::cont::ArrayHandle<ParticleType>& seeds,                            \
    const Stepper& stepper,                                                  \
    vtkm::Id maxSteps,                                                       \
    ResultType<ParticleType>& result) const                                  \
  {                                                                          \
    WorkletType Worklet;                                                     \
    result = Worklet.Run(stepper, seeds, maxSteps);                          \
  }

//-----
// Specialization for different worklets.

// ParticleAdvection worklet
DEFINE_DOADVECT(GridEvalType,
                vtkm::Particle,
                vtkm::worklet::ParticleAdvectionResult,
                vtkm::worklet::ParticleAdvection)
DEFINE_DOADVECT(GridEvalType,
                vtkm::ChargedParticle,
                vtkm::worklet::ParticleAdvectionResult,
                vtkm::worklet::ParticleAdvection)

// Streamline worklet
DEFINE_DOADVECT(GridEvalType,
                vtkm::Particle,
                vtkm::worklet::StreamlineResult,
                vtkm::worklet::Streamline)
DEFINE_DOADVECT(GridEvalType,
                vtkm::ChargedParticle,
                vtkm::worklet::StreamlineResult,
                vtkm::worklet::Streamline)

// PathParticle worklet
DEFINE_DOADVECT(TemporalGridEvalType,
                vtkm::Particle,
                vtkm::worklet::ParticleAdvectionResult,
                vtkm::worklet::ParticleAdvection)
DEFINE_DOADVECT(TemporalGridEvalType,
                vtkm::ChargedParticle,
                vtkm::worklet::ParticleAdvectionResult,
                vtkm::worklet::ParticleAdvection)

// Pathline worklet
DEFINE_DOADVECT(TemporalGridEvalType,
                vtkm::Particle,
                vtkm::worklet::StreamlineResult,
                vtkm::worklet::Streamline)
DEFINE_DOADVECT(TemporalGridEvalType,
                vtkm::ChargedParticle,
                vtkm::worklet::StreamlineResult,
                vtkm::worklet::Streamline)

}
}
} // namespace vtkm::filter::particleadvection

#endif
