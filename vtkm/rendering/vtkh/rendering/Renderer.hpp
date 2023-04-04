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
  virtual void SetInput(vtkm::cont::PartitionedDataSet *input) {this->m_input = input;}
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
  vtkm::cont::PartitionedDataSet*  GetInput() { return this->m_input; }
  vtkm::cont::PartitionedDataSet*  GetOutput() { return this->m_output; }
  vtkm::Range                 GetRange() const;
  bool                        GetHasColorTable() const;
protected:

  // image related data with cinema support
  vtkm::Bounds                             m_bounds;
  vtkm::cont::ColorTable                   m_color_table;
  Compositor                              *m_compositor;
  bool                                     m_do_composite;
  int                                      m_field_index;
  std::string                              m_field_name;
  bool                                     m_has_color_table;
  vtkm::cont::PartitionedDataSet *m_input = nullptr;
  vtkm::cont::PartitionedDataSet *m_output = nullptr;
  vtkmMapperPtr                            m_mapper;
  std::vector<vtkh::Render>                m_renders;
  vtkm::Range                              m_range;

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
