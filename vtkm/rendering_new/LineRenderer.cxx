//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering_new/LineRenderer.h>

#include <memory>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperCylinder.h>

namespace vtkm
{
namespace rendering_new
{

LineRenderer::LineRenderer()
  : Radius(0.5f)
  , RadiusSet(false)

{
  typedef vtkm::rendering::MapperCylinder TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  Mapper = mapper;
}

LineRenderer::~LineRenderer() {}

std::shared_ptr<vtkm::rendering::Canvas> LineRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

void LineRenderer::SetRadius(vtkm::Float32 radius)
{
  this->Radius = radius;
  this->RadiusSet = true;
}

void LineRenderer::PreExecute(vtkm::rendering_new::Plot& plot)
{
  Renderer::PreExecute(plot);

  typedef vtkm::rendering::MapperCylinder MapperType;
  std::shared_ptr<MapperType> mapper = std::dynamic_pointer_cast<MapperType>(this->Mapper);

  // allow for the default mapper radius
  if (this->RadiusSet)
  {
    mapper->SetRadius(this->Radius);
  }
}

}
} // namespace vtkm::rendering_new
