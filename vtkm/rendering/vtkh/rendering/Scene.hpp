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
#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/rendering/Plot.h>
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>

namespace vtkh
{

class VTKM_RENDERING_EXPORT Scene
{
public:
  Scene();

  void AddPlot(vtkh::Plot &plot) {  this->Plots.push_back(plot); }
  void SetPlots(const std::vector<vtkh::Plot> &plots) {this->Plots = plots;}
  void AddRenderer(vtkh::Renderer *render);
  void RenderORIG();
  void Render();
  void SetRenderBatchSize(int batch_size);
  int  GetRenderBatchSize() const { return this->BatchSize;}
protected:
  bool IsMesh(vtkh::Renderer *renderer);
  bool IsVolume(vtkh::Renderer *renderer);
  void SynchDepths(std::vector<vtkh::Plot> &plots);

private:
  int                          BatchSize;
  bool                         HasVolume;
  std::list<vtkh::Renderer*>   Renderers;
  std::vector<vtkh::Plot>    Plots;

}; // class scene

} //namespace  vtkh

#endif //vtkm_rendering_rendering_Scene_h
