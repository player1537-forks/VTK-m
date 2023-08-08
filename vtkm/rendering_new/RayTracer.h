//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_RayTracer_h
#define vtkm_rendering_rendering_RayTracer_h

#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering_new/Renderer.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_EXPORT RayTracer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  RayTracer();
  virtual ~RayTracer();
  std::string GetName() const override;
  void SetShadingOn(bool on) override;
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_rendering_RayTracer_h
