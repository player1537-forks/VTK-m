//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 Sandia Corporation.
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#ifndef vtk_m_worklet_Gradient_h
#define vtk_m_worklet_Gradient_h

#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

#include <vtkm/worklet/gradient/CellGradient.h>
#include <vtkm/worklet/gradient/Divergence.h>
#include <vtkm/worklet/gradient/GradientOutput.h>
#include <vtkm/worklet/gradient/PointGradient.h>
#include <vtkm/worklet/gradient/QCriterion.h>
#include <vtkm/worklet/gradient/Vorticity.h>

namespace vtkm
{
namespace worklet
{

template <typename T>
struct GradientOutputFields;

namespace gradient
{

//-----------------------------------------------------------------------------
template <typename CoordinateSystem, typename T, typename S, typename Device>
struct DeducedPointGrad
{
  DeducedPointGrad(const CoordinateSystem& coords,
                   const vtkm::cont::ArrayHandle<T, S>& field,
                   GradientOutputFields<T>* result)
    : Points(&coords)
    , Field(&field)
    , Result(result)
  {
  }

  template <typename CellSetType>
  void operator()(const CellSetType& cellset) const
  {
    vtkm::worklet::DispatcherMapTopology<PointGradient<T>, Device> dispatcher;
    dispatcher.Invoke(cellset, //topology to iterate on a per point basis
                      cellset, //whole cellset in
                      *this->Points,
                      *this->Field,
                      *this->Result);
  }

  const CoordinateSystem* const Points;
  const vtkm::cont::ArrayHandle<T, S>* const Field;
  GradientOutputFields<T>* Result;

private:
  void operator=(const DeducedPointGrad<CoordinateSystem, T, S, Device>&) = delete;
};

} //namespace gradient

template <typename T>
struct GradientOutputFields : public vtkm::exec::ExecutionObjectBase
{

  using ValueType = T;
  using BaseTType = typename vtkm::BaseComponent<T>::Type;

  template <typename DeviceAdapter>
  struct ExecutionTypes
  {
    using Portal = vtkm::exec::GradientOutput<T, DeviceAdapter>;
  };

  GradientOutputFields()
    : Gradient()
    , Divergence()
    , Vorticity()
    , QCriterion()
    , StoreGradient(true)
    , ComputeDivergence(false)
    , ComputeVorticity(false)
    , ComputeQCriterion(false)
  {
  }

  GradientOutputFields(bool store, bool divergence, bool vorticity, bool qc)
    : Gradient()
    , Divergence()
    , Vorticity()
    , QCriterion()
    , StoreGradient(store)
    , ComputeDivergence(divergence)
    , ComputeVorticity(vorticity)
    , ComputeQCriterion(qc)
  {
  }

  /// Add divergence field to the output data.
  /// The input array must have 3 components in order to compute this.
  /// The default is off.
  void SetComputeDivergence(bool enable) { ComputeDivergence = enable; }
  bool GetComputeDivergence() const { return ComputeDivergence; }

  /// Add voriticity/curl field to the output data.
  /// The input array must have 3 components in order to compute this.
  /// The default is off.
  void SetComputeVorticity(bool enable) { ComputeVorticity = enable; }
  bool GetComputeVorticity() const { return ComputeVorticity; }

  /// Add Q-criterion field to the output data.
  /// The input array must have 3 components in order to compute this.
  /// The default is off.
  void SetComputeQCriterion(bool enable) { ComputeQCriterion = enable; }
  bool GetComputeQCriterion() const { return ComputeQCriterion; }

  /// Add gradient field to the output data.
  /// The input array must have 3 components in order to disable this.
  /// The default is on.
  void SetComputeGradient(bool enable) { StoreGradient = enable; }
  bool GetComputeGradient() const { return StoreGradient; }

  //todo fix this for scalar
  template <typename DeviceAdapter>
  vtkm::exec::GradientOutput<T, DeviceAdapter> PrepareForOutput(vtkm::Id size, DeviceAdapter)
  {
    vtkm::exec::GradientOutput<T, DeviceAdapter> portal(this->StoreGradient,
                                                        this->ComputeDivergence,
                                                        this->ComputeVorticity,
                                                        this->ComputeQCriterion,
                                                        this->Gradient,
                                                        this->Divergence,
                                                        this->Vorticity,
                                                        this->QCriterion,
                                                        size);
    return portal;
  }

  vtkm::cont::ArrayHandle<vtkm::Vec<T, 3>> Gradient;
  vtkm::cont::ArrayHandle<BaseTType> Divergence;
  vtkm::cont::ArrayHandle<vtkm::Vec<BaseTType, 3>> Vorticity;
  vtkm::cont::ArrayHandle<BaseTType> QCriterion;

private:
  bool StoreGradient;
  bool ComputeDivergence;
  bool ComputeVorticity;
  bool ComputeQCriterion;
};
class PointGradient
{
public:
  template <typename CellSetType,
            typename CoordinateSystem,
            typename T,
            typename S,
            typename DeviceAdapter>
  vtkm::cont::ArrayHandle<vtkm::Vec<T, 3>> Run(const CellSetType& cells,
                                               const CoordinateSystem& coords,
                                               const vtkm::cont::ArrayHandle<T, S>& field,
                                               DeviceAdapter device)
  {
    vtkm::worklet::GradientOutputFields<T> extraOutput(true, false, false, false);
    return this->Run(cells, coords, field, extraOutput, device);
  }

  template <typename CellSetType,
            typename CoordinateSystem,
            typename T,
            typename S,
            typename DeviceAdapter>
  vtkm::cont::ArrayHandle<vtkm::Vec<T, 3>> Run(const CellSetType& cells,
                                               const CoordinateSystem& coords,
                                               const vtkm::cont::ArrayHandle<T, S>& field,
                                               GradientOutputFields<T>& extraOutput,
                                               DeviceAdapter)
  {
    //we are using cast and call here as we pass the cells twice to the invoke
    //and want the type resolved once before hand instead of twice
    //by the dispatcher ( that will cost more in time and binary size )
    gradient::DeducedPointGrad<CoordinateSystem, T, S, DeviceAdapter> func(
      coords, field, &extraOutput);
    vtkm::cont::CastAndCall(cells, func);
    return extraOutput.Gradient;
  }
};

class CellGradient
{
public:
  template <typename CellSetType,
            typename CoordinateSystem,
            typename T,
            typename S,
            typename DeviceAdapter>
  vtkm::cont::ArrayHandle<vtkm::Vec<T, 3>> Run(const CellSetType& cells,
                                               const CoordinateSystem& coords,
                                               const vtkm::cont::ArrayHandle<T, S>& field,
                                               DeviceAdapter device)
  {
    vtkm::worklet::GradientOutputFields<T> extra(true, false, false, false);
    return this->Run(cells, coords, field, extra, device);
  }

  template <typename CellSetType,
            typename CoordinateSystem,
            typename T,
            typename S,
            typename DeviceAdapter>
  vtkm::cont::ArrayHandle<vtkm::Vec<T, 3>> Run(const CellSetType& cells,
                                               const CoordinateSystem& coords,
                                               const vtkm::cont::ArrayHandle<T, S>& field,
                                               GradientOutputFields<T>& extraOutput,
                                               DeviceAdapter)
  {
    using DispatcherType =
      vtkm::worklet::DispatcherMapTopology<vtkm::worklet::gradient::CellGradient<T>, DeviceAdapter>;

    vtkm::worklet::gradient::CellGradient<T> worklet;
    DispatcherType dispatcher(worklet);


    dispatcher.Invoke(cells, coords, field, extraOutput);
    return extraOutput.Gradient;
  }
};
}
} // namespace vtkm::worklet
#endif
