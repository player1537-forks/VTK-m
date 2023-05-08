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

#include <vtkm/rendering/MapperVolume.h>
#include <vtkm/rendering/vtkh/rendering/Renderer.h>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh
{

namespace detail
{
class VolumeWrapper;
}

class VTKM_RENDERING_EXPORT VolumeRenderer : public Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  VolumeRenderer();
  virtual ~VolumeRenderer();
  std::string GetName() const override;
  void SetNumberOfSamples(const int num_samples);


  void Update(vtkh::Plot& plot) override;
  virtual void SetInput(const vtkm::rendering::Actor& actor) override;

protected:
  virtual void Composite(vtkh::Plot& plot) override;
  virtual void PreExecute(vtkh::Plot& plot) override;
  virtual void DoExecute(vtkh::Plot& plot) override;
  virtual void PostExecute(vtkh::Plot& plot) override;

  void RenderOneDomainPerRank(vtkh::Plot& plot);
  void RenderMultipleDomainsPerRank(vtkh::Plot& plot);

  void CorrectOpacity();
  void FindVisibilityOrdering(vtkh::Plot& plot);
  void DepthSort(int num_domains,
                 std::vector<float>& min_depths,
                 std::vector<int>& local_vis_order);
  float FindMinDepth(const vtkm::rendering::Camera& camera, const vtkm::Bounds& bounds) const;

  int NumSamples;
  float SampleDist;
  bool HasUnstructured;
  std::shared_ptr<vtkm::rendering::MapperVolume> Tracer;
  vtkm::cont::ColorTable CorrectedColorTable;
  std::vector<int> VisibilityOrders;

  void ClearWrappers();
  std::vector<detail::VolumeWrapper*> Wrappers;
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_VolumeRenderer_h
