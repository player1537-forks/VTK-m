//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "PointRenderer.hpp"

#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperPoint.h>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>
#include <vtkm/rendering/vtkh/filters/ParticleMerging.hpp>
#include <memory>

namespace vtkh {

PointRenderer::PointRenderer()
{
  typedef vtkm::rendering::MapperPoint TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  this->Mapper = mapper;
}

std::shared_ptr<vtkm::rendering::Canvas>
PointRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

void
PointRenderer::SetBaseRadius(vtkm::Float32 radius)
{
  this->BaseRadius = radius;
  this->RadiusSet = true;
}

void
PointRenderer::PreExecute()
{
  Renderer::PreExecute();

  typedef vtkm::rendering::MapperPoint MapperType;
  std::shared_ptr<MapperType> mesh_mapper =
    std::dynamic_pointer_cast<MapperType>(this->Mapper);

  if(this->UseNodes)
  {
    mesh_mapper->UseNodes();
  }
  else
  {
    mesh_mapper->UseCells();
  }

  vtkm::Float32 radius = this->BaseRadius;
  if(this->RadiusSet)
  {
    mesh_mapper->SetRadius(this->BaseRadius);
  }
  else
  {
    vtkm::Bounds coordBounds = this->Actor.GetDataSet().GetGlobalBounds();
    // set a default radius
    vtkm::Float64 lx = coordBounds.X.Length();
    vtkm::Float64 ly = coordBounds.Y.Length();
    vtkm::Float64 lz = coordBounds.Z.Length();
    vtkm::Float64 mag = vtkm::Sqrt(lx * lx + ly * ly + lz * lz);
    // same as used in vtk ospray
    constexpr vtkm::Float64 heuristic = 1000.;
    radius = static_cast<vtkm::Float32>(mag / heuristic);
    // we likely have a data set with no cells so just set some radius
    if(radius == 0.f)
    {
      radius = 0.00001f;
    }
    mesh_mapper->SetRadius(radius);
  }

  if(!this->UseNodes && vtkh::IsPointMesh(this->Actor.GetDataSet()) && this->UsePointMerging)
  {
    throw vtkm::cont::ErrorBadValue("Need to implement this.");
    /*
    vtkm::Float32 max_radius = radius;
    if(this->UseVariableRadius)
    {
      max_radius = radius + radius * this->DeltaRadius;
    }

    ParticleMerging  merger;
    merger.SetInput(this->Input);
    merger.SetField(this->FieldName);
    merger.SetRadius(max_radius * this->RadiusMult);
    merger.Update();
    this->Input = merger.GetOutput();
    this->DeleteInput = true;
    */
  }

  mesh_mapper->UseVariableRadius(this->UseVariableRadius);
  mesh_mapper->SetRadiusDelta(this->DeltaRadius);

}

void PointRenderer::PostExecute()
{
  throw vtkm::cont::ErrorBadValue("Need to implement this.");
  Renderer::PostExecute();
  if(this->DeleteInput)
  {
//    delete this->Input;
//    this->Input = nullptr;
  }
}

} // namespace vtkh
