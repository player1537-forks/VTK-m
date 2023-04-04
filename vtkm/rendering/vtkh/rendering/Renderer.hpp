//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_Renderer_h
#define vtkm_rendering_rendering_Renderer_h

#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/rendering/Render.hpp>
#include <vtkm/rendering/vtkh/compositing/Image.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

//#include <vector>

namespace vtkh
{

class Compositor;

class VTKM_RENDERING_EXPORT Renderer
{
public:
  typedef std::shared_ptr<vtkm::rendering::Canvas> vtkmCanvasPtr;
  typedef std::shared_ptr<vtkm::rendering::Mapper> vtkmMapperPtr;
  typedef vtkm::rendering::Camera vtkmCamera;

  Renderer();
  virtual ~Renderer();
  virtual std::string GetName() const = 0;
  virtual void SetInput(vtkm::cont::PartitionedDataSet *input) {this->Input = input;}
  virtual void SetShadingOn(bool on);
  virtual void Update();

  void AddRender(vtkh::Render &render);
  void ClearRenders();

  void SetField(const std::string field_name);
  virtual void SetColorTable(const vtkm::cont::ColorTable &color_table);
  void SetDoComposite(bool do_composite);
  void SetRenders(const std::vector<Render> &renders);
  void SetRange(const vtkm::Range &range);
  void DisableColorBar();

  vtkm::cont::ColorTable      GetColorTable() const;
  std::string                 GetFieldName() const;
  int                         GetNumberOfRenders() const;
  std::vector<Render>         GetRenders() const;
  vtkm::cont::PartitionedDataSet*  GetInput() { return this->Input; }
  vtkm::cont::PartitionedDataSet*  GetOutput() { return this->Output; }
  vtkm::Range                 GetRange() const;
  bool                        GetHasColorTable() const;
protected:

  // image related data with cinema support
  vtkm::Bounds                             Bounds;
  vtkm::cont::ColorTable                   ColorTable;
  Compositor                              *Compositor;
  bool                                     DoComposite;
  int                                      FieldIndex;
  std::string                              FieldName;
  bool                                     HasColorTable;
  vtkm::cont::PartitionedDataSet *Input = nullptr;
  vtkm::cont::PartitionedDataSet *Output = nullptr;
  vtkmMapperPtr                            Mapper;
  std::vector<vtkh::Render>                Renders;
  vtkm::Range                              Range;

  // methods
virtual void PreExecute(); // override;
virtual void PostExecute(); // override;
virtual void DoExecute(); // override;

virtual void CheckForRequiredField(const std::string &field_name);

  virtual void Composite(const int &num_images);
  void ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_Renderer_h
