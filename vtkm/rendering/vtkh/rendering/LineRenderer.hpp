//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_LineRenderer_h
#define vtkm_rendering_rendering_LineRenderer_h

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkh/rendering/Plot.h>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT LineRenderer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  LineRenderer();
  virtual ~LineRenderer();
  std::string GetName() const override {return "vtkh::LineRenderer";}
  void PreExecute(std::vector<vtkh::Plot>& plots) override;
  void SetRadius(vtkm::Float32 radius);
private:
  vtkm::Float32 Radius;
  bool RadiusSet;
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_LineRenderer_h
