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

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT RayTracer : public Renderer
{
public:
  RayTracer();
  virtual ~RayTracer();
  std::string GetName() const override;
  void SetShadingOn(bool on) override;
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_RayTracer_h
