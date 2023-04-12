//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_VolumeRenderer_h
#define vtkm_rendering_rendering_VolumeRenderer_h

#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/MapperVolume.h>

namespace vtkh {

namespace detail
{
  class VolumeWrapper;
}

class VTKM_RENDERING_EXPORT VolumeRenderer : public Renderer
{
public:
  VolumeRenderer();
  virtual ~VolumeRenderer();
  std::string GetName() const override;
  void SetNumberOfSamples(const int num_samples);
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);

  void Update() override;
  virtual void SetInput(const vtkm::rendering::Actor& actor) override;

protected:
  virtual void Composite() override;
  virtual void PreExecute() override;
  virtual void DoExecute() override;
  virtual void PostExecute() override;

  void RenderOneDomainPerRank();
  void RenderMultipleDomainsPerRank();

  void CorrectOpacity();
  void FindVisibilityOrdering();
  void DepthSort(int num_domains,
                 std::vector<float> &min_depths,
                 std::vector<int> &local_vis_order);
  float FindMinDepth(const vtkm::rendering::Camera &camera,
                     const vtkm::Bounds &bounds) const;

  int NumSamples;
  float SampleDist;
  bool HasUnstructured;
  std::shared_ptr<vtkm::rendering::MapperVolume> Tracer;
  vtkm::cont::ColorTable CorrectedColorTable;
  std::vector<std::vector<int>> VisibilityOrders;

  void ClearWrappers();
  std::vector<detail::VolumeWrapper*> Wrappers;
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_VolumeRenderer_h
