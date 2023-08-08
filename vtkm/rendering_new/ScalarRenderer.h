//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_newScalarRenderer_h
#define vtkm_rendering_newScalarRenderer_h

#include <vector>
#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/compositing/PayloadImage.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/ScalarRenderer.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT ScalarRenderer : public vtkm::rendering_new::Renderer
{
public:
  typedef vtkm::rendering::Camera vtkmCamera;
  using Result = vtkm::rendering::ScalarRenderer::Result;

  ScalarRenderer();
  virtual ~ScalarRenderer();
  virtual void Update(vtkm::rendering_new::Plot& plot) override;
  virtual std::string GetName() const override;

  void SetCamera(vtkmCamera& camera);

  int GetNumberOfCameras() const;
  void SetHeight(const int height) { this->Height = height; }
  void SetWidth(const int width) { this->Width = width; }

  //moved from base class.
  vtkm::cont::PartitionedDataSet* GetOutput() { return this->Output; }

protected:
  int Width = 1024;
  int Height = 1024;
  // image related data with cinema support
  vtkmCamera Camera;
  // methods
  virtual void PreExecute(vtkm::rendering_new::Plot& plot) override;
  virtual void PostExecute(vtkm::rendering_new::Plot& plot) override;
  virtual void DoExecute(vtkm::rendering_new::Plot& plot) override;

  PayloadImage* Convert(Result& result);
  ScalarRenderer::Result Convert(PayloadImage& image, std::vector<std::string>& names);
  //void ImageToDataSet(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);

  //moved from base class.
  vtkm::cont::PartitionedDataSet* Output = nullptr;
};


}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_newScalarRenderer_h
