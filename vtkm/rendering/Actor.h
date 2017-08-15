//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2015 Sandia Corporation.
//  Copyright 2015 UT-Battelle, LLC.
//  Copyright 2015 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtk_m_rendering_Actor_h
#define vtk_m_rendering_Actor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

#include <memory>

namespace vtkm
{
namespace rendering
{

class VTKM_RENDERING_EXPORT Actor
{
public:
  Actor(const vtkm::cont::DynamicCellSet& cells,
        const vtkm::cont::CoordinateSystem& coordinates,
        const vtkm::cont::Field& scalarField,
        const vtkm::rendering::ColorTable& colorTable = vtkm::rendering::ColorTable("default"));

  Actor(const vtkm::cont::DynamicCellSet& cells,
        const vtkm::cont::CoordinateSystem& coordinates,
        const vtkm::cont::Field& scalarField,
        const vtkm::rendering::Color& color);

  void Render(vtkm::rendering::Mapper& mapper,
              vtkm::rendering::Canvas& canvas,
              const vtkm::rendering::Camera& camera) const;

  const vtkm::cont::DynamicCellSet& GetCells() const;

  const vtkm::cont::CoordinateSystem& GetCoordinates() const;

  const vtkm::cont::Field& GetScalarField() const;

  const vtkm::rendering::ColorTable& GetColorTable() const;

  const vtkm::Range& GetScalarRange() const;

  const vtkm::Bounds& GetSpatialBounds() const;

  void SetScalarRange(const vtkm::Range& scalarRange);

private:
  struct InternalsType;
  std::shared_ptr<InternalsType> Internals;

  struct RangeFunctor;

  void Init(const vtkm::cont::CoordinateSystem& coordinates, const vtkm::cont::Field& scalarField);
};
}
} //namespace vtkm::rendering

#endif //vtk_m_rendering_Actor_h
