//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_newVolumeRenderer_h
#define vtkm_rendering_newVolumeRenderer_h

#include <vtkm/rendering/MapperVolume.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

namespace detail
{
class VolumeWrapper;
}

class VTKM_RENDERING_NEW_EXPORT VolumeRenderer : public vtkm::rendering_new::Renderer
{
public:
  static std::shared_ptr<vtkm::rendering::Canvas> GetNewCanvas(int width = 1024, int height = 1024);

  VolumeRenderer();
  virtual ~VolumeRenderer();
  std::string GetName() const override;
  void SetNumberOfSamples(const int num_samples);


  void Update(vtkm::rendering_new::Plot& plot) override;
  virtual void SetInput(const vtkm::rendering::Actor& actor) override;

protected:
  virtual void Composite(vtkm::rendering_new::Plot& plot) override;
  virtual void PreExecute(vtkm::rendering_new::Plot& plot) override;
  virtual void DoExecute(vtkm::rendering_new::Plot& plot) override;
  virtual void PostExecute(vtkm::rendering_new::Plot& plot) override;

  void RenderOneDomainPerRank(vtkm::rendering_new::Plot& plot);
  void RenderMultipleDomainsPerRank(vtkm::rendering_new::Plot& plot);

  void CorrectOpacity();
  void FindVisibilityOrdering(vtkm::rendering_new::Plot& plot);
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

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_newVolumeRenderer_h
