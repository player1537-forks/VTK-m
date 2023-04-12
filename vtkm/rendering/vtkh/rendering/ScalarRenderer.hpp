//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_ScalarRenderer_h
#define vtkm_rendering_rendering_ScalarRenderer_h

#include <vector>
#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkh/rendering/Render.hpp>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/ScalarRenderer.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT ScalarRenderer : public vtkh::Renderer
{
public:
  typedef vtkm::rendering::Camera vtkmCamera;
  using Result = vtkm::rendering::ScalarRenderer::Result;

  ScalarRenderer();
  virtual ~ScalarRenderer();
  virtual void Update() override;
  virtual std::string GetName() const override;

  void SetCamera(vtkmCamera &camera);

  int GetNumberOfCameras() const;
  void SetHeight(const int height) { this->Height = height; }
  void SetWidth(const int width) {this->Width = width; }

  //moved from base class.
  vtkm::cont::PartitionedDataSet*  GetOutput() { return this->Output; }
protected:

  int Width = 1024;
  int Height = 1024;
  // image related data with cinema support
  vtkmCamera  Camera;
  // methods
  virtual void PreExecute() override;
  virtual void PostExecute() override;
  virtual void DoExecute() override;

  PayloadImage * Convert(Result &result);
  ScalarRenderer::Result Convert(PayloadImage &image, std::vector<std::string> &names);
  //void ImageToDataSet(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);

  //moved from base class.
  vtkm::cont::PartitionedDataSet *Output = nullptr;
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_ScalarRenderer_h
