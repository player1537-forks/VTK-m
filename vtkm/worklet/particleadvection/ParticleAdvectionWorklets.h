//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h
#define vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h

#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCast.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/ExecutionObjectBase.h>

#include <vtkm/Particle.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/worklet/particleadvection/Stepper.h>

#ifdef VTKM_CUDA
#include <vtkm/cont/cuda/internal/ScopedCudaStackSize.h>
#endif

namespace vtkm
{
namespace worklet
{
namespace particleadvection
{

class ParticleAdvectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  ParticleAdvectWorklet()
    : ZeroVelocityEps(std::numeric_limits<vtkm::FloatDefault>::epsilon() *
                      static_cast<vtkm::FloatDefault>(10))
  {
  }

  using ControlSignature = void(FieldIn idx,
                                ExecObject integrator,
                                ExecObject integralCurve,
                                FieldIn maxSteps);
  using ExecutionSignature = void(_1 idx, _2 integrator, _3 integralCurve, _4 maxSteps);
  using InputDomain = _1;

  template <typename IntegratorType, typename IntegralCurveType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const IntegratorType& integrator,
                            IntegralCurveType& integralCurve,
                            const vtkm::Id& maxSteps) const
  {
    auto particle = integralCurve.GetParticle(idx);
    vtkm::FloatDefault time = particle.Time;
    bool tookAnySteps = false;

    //the integrator status needs to be more robust:
    // 1. you could have success AND at temporal boundary.
    // 2. could you have success AND at spatial?
    // 3. all three?
    integralCurve.PreStepUpdate(idx);
    do
    {
      particle = integralCurve.GetParticle(idx);
      vtkm::Vec3f outpos;
      auto status = integrator.Step(particle, time, outpos);
      if (status.CheckOk())
      {
        integralCurve.StepUpdate(idx, particle, time, outpos);
        tookAnySteps = true;
      }

      //We can't take a step inside spatial boundary.
      //Try and take a step just past the boundary.
      else if (status.CheckSpatialBounds())
      {
        status = integrator.SmallStep(particle, time, outpos);
        if (status.CheckOk())
        {
          integralCurve.StepUpdate(idx, particle, time, outpos);
          tookAnySteps = true;
        }
      }

      //Check to see if we are in a zero velocity.
      vtkm::FloatDefault dist = vtkm::Magnitude(particle.Pos - outpos);
      if (dist < this->ZeroVelocityEps)
        status.SetZeroVelocity();

      integralCurve.StatusUpdate(idx, status, maxSteps);
    } while (integralCurve.CanContinue(idx));

    //Mark if any steps taken
    integralCurve.UpdateTookSteps(idx, tookAnySteps);
  }
private:
  vtkm::FloatDefault ZeroVelocityEps = std::numeric_limits<vtkm::FloatDefault>::epsilon() * static_cast<vtkm::FloatDefault>(10.0f);
};


template <typename IntegratorType, typename ParticleType>
class ParticleAdvectionWorklet
{
public:
  VTKM_EXEC_CONT ParticleAdvectionWorklet() {}

  ~ParticleAdvectionWorklet() {}

  void Run(const IntegratorType& integrator,
           vtkm::cont::ArrayHandle<ParticleType>& particles,
           vtkm::Id& maxSteps)
  {
    using ParticleArrayType = vtkm::worklet::particleadvection::Particles<ParticleType>;

    vtkm::Id numSeeds = static_cast<vtkm::Id>(particles.GetNumberOfValues());
    //Create and invoke the particle advection.
    vtkm::cont::ArrayHandleConstant<vtkm::Id> maxStepArr(maxSteps, numSeeds);
    vtkm::cont::ArrayHandleIndex idxArray(numSeeds);

    // TODO: The particle advection sometimes behaves incorrectly on CUDA if the stack size
    // is not changed thusly. This is concerning as the compiler should be able to determine
    // staticly the required stack depth. What is even more concerning is that the runtime
    // does not report a stack overflow. Rather, the worklet just silently reports the wrong
    // value. Until we determine the root cause, other problems may pop up.
#ifdef VTKM_CUDA
    // This worklet needs some extra space on CUDA.
    vtkm::cont::cuda::internal::ScopedCudaStackSize stack(16 * 1024);
    (void)stack;
#endif // VTKM_CUDA

    //TODO: max steps in invoker AND particlesObj
    ParticleArrayType particlesObj(particles, maxSteps);

    //Invoke particle advection worklet
    vtkm::cont::Invoker invoker;
    invoker(ParticleAdvectWorklet{}, idxArray, integrator, particlesObj, maxStepArr);
  }
};

namespace detail
{
class GetSteps : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  GetSteps() {}
  using ControlSignature = void(FieldIn, FieldOut);
  using ExecutionSignature = void(_1, _2);

  template <typename ParticleType>
  VTKM_EXEC void operator()(const ParticleType& p, vtkm::Id& numSteps) const
  {
    numSteps = p.NumSteps;
  }
};

class ComputeNumPoints : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  ComputeNumPoints() {}
  using ControlSignature = void(FieldIn, FieldIn, FieldOut);
  using ExecutionSignature = void(_1, _2, _3);

  // Offset is number of points in streamline.
  // 1 (inital point) + number of steps taken (p.NumSteps - initalNumSteps)
  template <typename ParticleType>
  VTKM_EXEC void operator()(const ParticleType& p,
                            const vtkm::Id& initialNumSteps,
                            vtkm::Id& diff) const
  {
    diff = 1 + p.NumSteps - initialNumSteps;
  }
};
} // namespace detail


template <typename IntegratorType, typename ParticleType>
class StreamlineWorklet
{
public:
  template <typename PointStorage, typename PointStorage2>
  void Run(const IntegratorType& it,
           vtkm::cont::ArrayHandle<ParticleType, PointStorage>& particles,
           vtkm::Id& maxSteps,
           vtkm::cont::ArrayHandle<vtkm::Vec3f, PointStorage2>& positions,
           vtkm::cont::CellSetExplicit<>& polyLines)
  {
//    using ParticleWorkletDispatchType = typename vtkm::worklet::DispatcherMapField<
//      vtkm::worklet::particleadvection::ParticleAdvectWorklet>;
    using StreamlineArrayType =
      vtkm::worklet::particleadvection::StateRecordingParticles<ParticleType>;

    vtkm::cont::Invoker invoker;

    vtkm::Id numSeeds = static_cast<vtkm::Id>(particles.GetNumberOfValues());
    vtkm::cont::ArrayHandleIndex idxArray(numSeeds);

    vtkm::cont::ArrayHandle<vtkm::Id> initialStepsTaken;
    invoker(detail::GetSteps{}, particles, initialStepsTaken);


    // This method uses the same workklet as ParticleAdvectionWorklet::Run (and more). Yet for
    // some reason ParticleAdvectionWorklet::Run needs this adjustment while this method does
    // not.
#ifdef VTKM_CUDA
    // // This worklet needs some extra space on CUDA.
    // vtkm::cont::cuda::internal::ScopedCudaStackSize stack(16 * 1024);
    // (void)stack;
#endif // VTKM_CUDA

    //Run streamline worklet
    StreamlineArrayType streamlines(particles, maxSteps);
    //ParticleWorkletDispatchType particleWorkletDispatch;
    vtkm::cont::ArrayHandleConstant<vtkm::Id> maxStepsArr(maxSteps, numSeeds);

    invoker(ParticleAdvectWorklet{}, idxArray, it, streamlines, maxStepsArr);

    //Get the positions
    streamlines.GetCompactedHistory(positions);

    //Create the cells
    vtkm::cont::ArrayHandle<vtkm::Id> numPoints;
    invoker(detail::ComputeNumPoints{}, particles, initialStepsTaken, numPoints);

    vtkm::cont::ArrayHandle<vtkm::Id> cellIndex;
    vtkm::Id connectivityLen = vtkm::cont::Algorithm::ScanExclusive(numPoints, cellIndex);
    vtkm::cont::ArrayHandleCounting<vtkm::Id> connCount(0, 1, connectivityLen);
    vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
    vtkm::cont::ArrayCopy(connCount, connectivity);

    vtkm::cont::ArrayHandle<vtkm::UInt8> cellTypes;
    auto polyLineShape =
      vtkm::cont::make_ArrayHandleConstant<vtkm::UInt8>(vtkm::CELL_SHAPE_POLY_LINE, numSeeds);
    vtkm::cont::ArrayCopy(polyLineShape, cellTypes);

    auto offsets = vtkm::cont::ConvertNumComponentsToOffsets(numPoints);
    polyLines.Fill(positions.GetNumberOfValues(), cellTypes, connectivity, offsets);
  }
};
}
}
} // namespace vtkm::worklet::particleadvection

#endif // vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h
