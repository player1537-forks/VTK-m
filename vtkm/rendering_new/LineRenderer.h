//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_LineRenderer_h
#define vtkm_rendering_new_LineRenderer_h

#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT LineRenderer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  LineRenderer();
  virtual ~LineRenderer();
  std::string GetName() const override { return "vtkm::rendering_new::LineRenderer"; }
  void PreExecute(vtkm::rendering_new::Plot& plot) override;
  void SetRadius(vtkm::Float32 radius);

private:
  vtkm::Float32 Radius;
  bool RadiusSet;
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_new_LineRenderer_h
