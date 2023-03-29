//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_RENDERER_MESH_HPP
#define VTK_H_RENDERER_MESH_HPP

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT MeshRenderer : public Renderer
{
public:
  MeshRenderer();
  virtual ~MeshRenderer();
  std::string GetName() const override;
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);

  void SetIsOverlay(bool on);
  void SetShowInternal(bool on);
  void SetUseForegroundColor(bool on);
  bool GetIsOverlay() const;
  bool GetShowInternal() const;
protected:
  void PreExecute() override;
  bool m_use_foreground_color;
  bool m_is_overlay;
  bool m_show_internal;
};

} // namespace vtkh
#endif
