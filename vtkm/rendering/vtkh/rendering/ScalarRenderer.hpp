//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_SCALAR_RENDERER_HPP
#define VTK_H_SCALAR_RENDERER_HPP

#include <vector>
#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/filters/Filter.hpp>
#include <vtkm/rendering/vtkh/rendering/Render.hpp>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/ScalarRenderer.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT ScalarRenderer : public Filter
{
public:
  typedef vtkm::rendering::Camera vtkmCamera;
  using Result = vtkm::rendering::ScalarRenderer::Result;

  ScalarRenderer();
  virtual ~ScalarRenderer();
  virtual void Update();
  virtual std::string GetName() const override;

  void SetCamera(vtkmCamera &camera);

  int GetNumberOfCameras() const;
  vtkm::cont::PartitionedDataSet *GetInput();
  void SetHeight(const int height);
  void SetWidth(const int width);
protected:

  int m_width;
  int m_height;
  // image related data with cinema support
  vtkmCamera  m_camera;
  // methods
  virtual void PreExecute() override;
  virtual void PostExecute() override;
  virtual void DoExecute() override;

  PayloadImage * Convert(Result &result);
  ScalarRenderer::Result Convert(PayloadImage &image, std::vector<std::string> &names);
  //void ImageToDataSet(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
};

} // namespace vtkh
#endif
