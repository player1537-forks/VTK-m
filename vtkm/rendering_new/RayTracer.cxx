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
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering_new/RayTracer.h>

namespace vtkm
{
namespace rendering_new
{

RayTracer::RayTracer()
{
  typedef vtkm::rendering::MapperRayTracer TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  this->Mapper = mapper;
}

RayTracer::~RayTracer() {}

std::shared_ptr<vtkm::rendering::Canvas> RayTracer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

std::string RayTracer::GetName() const
{
  return "vtkm::rendering_new::RayTracer";
}

void RayTracer::SetShadingOn(bool on)
{
  // do nothing by default;
  typedef vtkm::rendering::MapperRayTracer TracerType;
  std::static_pointer_cast<TracerType>(this->Mapper)->SetShadingOn(on);
}


}
} // namespace vtkm::rendering_new
