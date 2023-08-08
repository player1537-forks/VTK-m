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
#include <vtkm/rendering/MapperWireframer.h>
#include <vtkm/rendering_new/MeshRenderer.h>
#include <vtkm/rendering_new/Plot.h>

namespace vtkm
{
namespace rendering_new
{

MeshRenderer::MeshRenderer()
{
  typedef vtkm::rendering::MapperWireframer MapperType;
  auto mapper = std::make_shared<MapperType>();
  this->Mapper = mapper;
}

void MeshRenderer::PreExecute(vtkm::rendering_new::Plot& plot)
{
  Renderer::PreExecute(plot);

  typedef vtkm::rendering::MapperWireframer MapperType;
  std::shared_ptr<MapperType> mesh_mapper = std::dynamic_pointer_cast<MapperType>(this->Mapper);

  mesh_mapper->SetShowInternalZones(this->ShowInternal);
  mesh_mapper->SetIsOverlay(this->IsOverlay);

  if (this->UseForegroundColor)
  {
    vtkm::rendering::Color fg = plot.GetCanvas().GetForegroundColor();
    vtkm::cont::ColorTable single_color;
    vtkm::Vec<vtkm::Float32, 3> fg_vec3_not_4;
    fg_vec3_not_4[0] = fg.Components[0];
    fg_vec3_not_4[1] = fg.Components[1];
    fg_vec3_not_4[2] = fg.Components[2];

    single_color.AddPoint(0.f, fg_vec3_not_4);
    single_color.AddPoint(1.f, fg_vec3_not_4);
    this->SetColorTable(single_color);
    this->HasColorTable = false;
  }
}

std::shared_ptr<vtkm::rendering::Canvas> MeshRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

}
} // namespace vtkm::rendering_new
