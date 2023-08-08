//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_newAnnotator_h
#define vtkm_rendering_newAnnotator_h

#include <vtkm/rendering/AxisAnnotation3D.h>
#include <vtkm/rendering/BoundingBoxAnnotation.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/ColorBarAnnotation.h>
#include <vtkm/rendering/WorldAnnotator.h>

#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT Annotator
{
public:
  Annotator(vtkm::rendering::Canvas& canvas, vtkm::rendering::Camera& camera, vtkm::Bounds bounds);
  ~Annotator();

  void RenderWorldAnnotations(vtkm::Vec<float, 3> axis_scale);
  void RenderScreenAnnotations(const std::vector<std::string>& field_names,
                               const std::vector<vtkm::Range>& ranges,
                               const std::vector<vtkm::cont::ColorTable>& color_tables);

protected:
  Annotator();

  bool Is3d;
  vtkm::rendering::Canvas& Canvas;
  vtkm::rendering::Camera& Camera;
  vtkm::Bounds Bounds;
  vtkm::rendering::BoundingBoxAnnotation BoxAnnotation;
  vtkm::rendering::AxisAnnotation3D XAxisAnnotation;
  vtkm::rendering::AxisAnnotation3D YAxisAnnotation;
  vtkm::rendering::AxisAnnotation3D ZAxisAnnotation;
  vtkm::rendering::ColorBarAnnotation ColorBarAnnotation;
  vtkm::rendering::WorldAnnotator* WorldAnnotator;
  std::vector<vtkm::Bounds> ColorBarPos;
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_newAnnotator_h
