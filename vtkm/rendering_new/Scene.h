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

#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering_new/Plot.h>

namespace vtkh
{

class VTKM_RENDERING_EXPORT Scene
{
public:
  Scene() {}

  void AddPlot(vtkh::Plot& plot) { this->Plots.push_back(plot); }
  void SetPlots(const std::vector<vtkh::Plot>& plots) { this->Plots = plots; }
  void SetRenderBatchSize(int batch_size);

  void Render();

private:
  int BatchSize = 10;
  bool RenderInBatches = false;
  std::vector<vtkh::Plot> Plots;

}; // class scene

} //namespace  vtkh

#endif //vtkm_rendering_rendering_Scene_h
