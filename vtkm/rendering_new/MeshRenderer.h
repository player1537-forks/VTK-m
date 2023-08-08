//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_newMeshRenderer_h
#define vtkm_rendering_newMeshRenderer_h

#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT MeshRenderer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  MeshRenderer();
  virtual ~MeshRenderer() {}
  std::string GetName() const override { return "vtkm::rendering_new::MeshRenderer"; }


  void SetIsOverlay(bool on) { this->IsOverlay = on; }
  void SetShowInternal(bool on) { this->ShowInternal = on; }
  void SetUseForegroundColor(bool on) { this->UseForegroundColor = on; }
  bool GetIsOverlay() const { return this->IsOverlay; }
  bool GetShowInternal() const { return this->ShowInternal; }

protected:
  void PreExecute(vtkm::rendering_new::Plot& plot) override;

  bool IsOverlay = false;
  bool ShowInternal = false;
  bool UseForegroundColor = false;
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_newMeshRenderer_h
