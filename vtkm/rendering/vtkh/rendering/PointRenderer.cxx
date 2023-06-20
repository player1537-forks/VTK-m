//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <memory>
#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperPoint.h>
#include <vtkm/rendering/vtkh/rendering/PointRenderer.h>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.h>

namespace vtkh
{

PointRenderer::PointRenderer()
{
  typedef vtkm::rendering::MapperPoint TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  this->Mapper = mapper;
}

std::shared_ptr<vtkm::rendering::Canvas> PointRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

void PointRenderer::SetBaseRadius(vtkm::Float32 radius)
{
  this->BaseRadius = radius;
  this->RadiusSet = true;
}

void PointRenderer::PreExecute(vtkh::Plot& plot)
{
  Renderer::PreExecute(plot);

  typedef vtkm::rendering::MapperPoint MapperType;
  std::shared_ptr<MapperType> meshMapper = std::dynamic_pointer_cast<MapperType>(this->Mapper);

  if (this->UseNodes)
  {
    meshMapper->UseNodes();
  }
  else
  {
    meshMapper->UseCells();
  }

  vtkm::Float32 radius = this->BaseRadius;
  if (this->RadiusSet)
  {
    meshMapper->SetRadius(this->BaseRadius);
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
    if (radius == 0.f)
    {
      radius = 0.00001f;
    }
    meshMapper->SetRadius(radius);
  }

  if (!this->UseNodes && vtkh::IsPointMesh(this->Actor.GetDataSet()) && this->UsePointMerging)
  {
    vtkm::filter::clean_grid::CleanGrid gridCleaner;
    gridCleaner.SetMergePoints(true);
    gridCleaner.SetFastMerge(true);


    vtkm::Float32 maxRadius = radius;
    if (this->UseVariableRadius)
    {
      maxRadius = radius + radius * this->DeltaRadius;
    }
    gridCleaner.SetTolerance(maxRadius * this->RadiusMult);
    gridCleaner.SetActiveField(this->Actor.GetScalarFieldName());
    auto output = gridCleaner.Execute(this->Actor.GetDataSet());
    this->Actor.SetDataSet(output);
    /*

    ParticleMerging  merger;
    merger.SetInput(this->Input);
    merger.SetField(this->FieldName);
    merger.SetRadius(
    merger.Update();
    this->Input = merger.GetOutput();
    this->DeleteInput = true;
*/
  }

  meshMapper->UseVariableRadius(this->UseVariableRadius);
  meshMapper->SetRadiusDelta(this->DeltaRadius);
}

void PointRenderer::PostExecute(vtkh::Plot& plot)
{
  Renderer::PostExecute(plot);
  /*
  if(this->DeleteInput)
  {
//    delete this->Input;
//    this->Input = nullptr;
  }
  */
}

} // namespace vtkh
