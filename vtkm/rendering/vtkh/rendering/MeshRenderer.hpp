//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_MeshRenderer_h
#define vtkm_rendering_rendering_MeshRenderer_h

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT MeshRenderer : public Renderer
{
public:
  MeshRenderer();
  virtual ~MeshRenderer() {}
  std::string GetName() const override { return "vtkh::MeshRenderer";}
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);

  void SetIsOverlay(bool on) { this->IsOverlay = on; }
  void SetShowInternal(bool on) { this->ShowInternal = on; }
  void SetUseForegroundColor(bool on) { this->UseForegroundColor = on; }
  bool GetIsOverlay() const { return this->IsOverlay; }
  bool GetShowInternal() const { return this->ShowInternal; }

protected:
  void PreExecute() override;

  bool IsOverlay = false;
  bool ShowInternal = false;
  bool UseForegroundColor = false;
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_MeshRenderer_h
