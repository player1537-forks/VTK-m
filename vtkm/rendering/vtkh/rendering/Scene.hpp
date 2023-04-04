//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_Scene_h
#define vtkm_rendering_rendering_Scene_h

#include <vector>
#include <list>
//#include <vtkm/rendering/<vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/rendering/Render.hpp>
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>

namespace vtkh
{

class VTKM_RENDERING_EXPORT Scene
{
private:
  std::list<vtkh::Renderer*>   m_renderers;
  std::vector<vtkh::Render>    m_renders;
  bool                         m_has_volume;
  int                          m_batch_size;
public:
 Scene();
 ~Scene();

  void AddRender(vtkh::Render &render);
  void SetRenders(const std::vector<vtkh::Render> &renders);
  void AddRenderer(vtkh::Renderer *render);
  void Render();
  void Save();
  void SetRenderBatchSize(int batch_size);
  int  GetRenderBatchSize() const;
protected:
  bool IsMesh(vtkh::Renderer *renderer);
  bool IsVolume(vtkh::Renderer *renderer);
  void SynchDepths(std::vector<vtkh::Render> &renders);
}; // class scene

} //namespace  vtkh

#endif //vtkm_rendering_rendering_Scene_h
