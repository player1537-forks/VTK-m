//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_RENDERER_POINTS_HPP
#define VTK_H_RENDERER_POINTS_HPP

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>

namespace vtkh {

#include <vtkm/rendering/vtkm_rendering_export.h>

class VTKM_RENDERING_EXPORT PointRenderer : public Renderer
{
public:
  PointRenderer();
  virtual ~PointRenderer();
  std::string GetName() const override;
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);
  void PreExecute() override;
  void PostExecute() override;
  void UseCells();
  void UseNodes();
  void UseVariableRadius(bool useVariableRadius);
  void SetBaseRadius(vtkm::Float32 radius);
  void SetRadiusDelta(vtkm::Float32 delta);
  void UsePointMerging(bool merge);
  // sets the number or radii to merge points
  // defualts to 2 * radius
  void PointMergeRadiusMultiplyer(vtkm::Float32 radius_mult);
private:
  bool m_use_nodes;
  bool m_radius_set;
  bool m_use_variable_radius;
  vtkm::Float32 m_base_radius;
  vtkm::Float32 m_delta_radius;
  bool m_use_point_merging;
  vtkm::Float32 m_radius_mult;
  bool m_delete_input;

};

} // namespace vtkh
#endif
