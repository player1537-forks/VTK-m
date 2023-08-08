//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_PointRenderer_h
#define vtkm_rendering_new_PointRenderer_h

#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT PointRenderer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  PointRenderer();
  virtual ~PointRenderer(){};
  std::string GetName() const override { return "vtkm::rendering_new::PointRenderer"; }
  void PreExecute(vtkm::rendering_new::Plot& plot) override;
  void PostExecute(vtkm::rendering_new::Plot& plot) override;
  void SetUseCells() { this->UseNodes = false; }
  void SetUseNodes() { this->UseNodes = true; }
  void SetUseVariableRadius(bool useVariableRadius) { this->UseVariableRadius = useVariableRadius; }
  void SetBaseRadius(vtkm::Float32 radius);
  void SetRadiusDelta(vtkm::Float32 delta) { this->DeltaRadius = delta; }
  void SetUsePointMerging(bool merge) { this->UsePointMerging = merge; }
  // sets the number or radii to merge points
  // defualts to 2 * radius
  void SetPointMergeRadiusMultiplyer(vtkm::Float32 radius_mult) { this->RadiusMult = radius_mult; }

private:
  bool IsPointMesh() const;

  bool UseNodes = true;
  bool RadiusSet = false;
  bool UseVariableRadius = false;
  vtkm::Float32 BaseRadius = 0.5f;
  vtkm::Float32 DeltaRadius = 0.5f;
  bool UsePointMerging = false;
  vtkm::Float32 RadiusMult = 2.f;
};


}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_new_PointRenderer_h
