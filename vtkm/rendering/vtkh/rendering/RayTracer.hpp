//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_RENDERER_RAY_TRACER_HPP
#define VTK_H_RENDERER_RAY_TRACER_HPP

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
#endif
