//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_newScene_h
#define vtkm_rendering_newScene_h

#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT Scene
{
public:
  Scene() {}

  void AddPlot(vtkm::rendering_new::Plot& plot) { this->Plots.push_back(plot); }
  void SetPlots(const std::vector<vtkm::rendering_new::Plot>& plots) { this->Plots = plots; }
  void SetRenderBatchSize(int batch_size);

  void Render();

private:
  int BatchSize = 10;
  bool RenderInBatches = false;
  std::vector<vtkm::rendering_new::Plot> Plots;

}; // class scene

}
} //namespace  vtkm::rendering_new

#endif //vtkm_rendering_newScene_h
