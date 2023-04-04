//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_Render_h
#define vtkm_rendering_rendering_Render_h

#include <vtkm/cont/PartitionedDataSet.h>
#include <vector>
#include <vtkm/rendering/vtkm_rendering_export.h>
//#include <vtkm/rendering/vtkh/DataSet.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/Mapper.h>

namespace vtkh {
//
// A Render contains the information needed to create a single image.
// There are 'n' canvases that matches the number of domains in the
// data set. It is possible to chain multiple plots together that
// are rendering separate data, i.e. the result of different data
// transformations, to handle this we keep track of the domain ids
// that each canvas is associated with.
//

class VTKM_RENDERING_EXPORT Render
{
public:
  typedef vtkm::rendering::CanvasRayTracer vtkmCanvas;

  Render();
  ~Render();
  Render                          Copy() const;
  vtkmCanvas&                     GetCanvas();
  const vtkm::rendering::Camera&  GetCamera() const;
  std::string                     GetImageName() const;
  std::vector<std::string>        GetComments() const;
  vtkm::Bounds                    GetSceneBounds() const;
  vtkm::Int32                     GetHeight() const;
  vtkm::Int32                     GetWidth() const;
  vtkm::rendering::Color          GetBackgroundColor() const;
bool                            GetShadingOn() const {   return this->DoShading; }
  void                            Print() const;

void                            SetRenderAnnotations(bool on) {   this->DoRenderAnnotations = on; }
void                            SetRenderWorldAnnotations(bool on) {  this->DoRenderWorldAnnotations = on;}
void                            SetRenderScreenAnnotations(bool on){   this->DoRenderScreenAnnotations = on; }
void                            SetRenderBackground(bool on) {  this->DoRenderBackground = on; }
  void                            ScaleWorldAnnotations(float x, float y, float z);
  void                            SetWidth(const vtkm::Int32 width);
  void                            SetHeight(const vtkm::Int32 height);
  void                            SetSceneBounds(const vtkm::Bounds &bounds);
  void                            SetCamera(const vtkm::rendering::Camera &camera);
  void                            SetImageName(const std::string &name);
  void                            SetComments(const std::vector<std::string> &comments);
  void                            SetBackgroundColor(float bg_color[4]);
  void                            SetForegroundColor(float fg_color[4]);
void                            SetShadingOn(bool on) {   this->DoShading = on; }
  void                            RenderWorldAnnotations();
  void                            RenderBackground();
  void                            RenderScreenAnnotations(const std::vector<std::string> &field_names,
                                                          const std::vector<vtkm::Range> &ranges,
                                                          const std::vector<vtkm::cont::ColorTable> &colors);
  void                            Save();
protected:
  vtkm::rendering::Camera      Camera;
  std::string                  ImageName;
  std::vector<std::string>     Comments;
  vtkm::Bounds                 SceneBounds;
  vtkm::Int32                  Width;
  vtkm::Int32                  Height;
  vtkm::rendering::Color       BgColor;
  vtkm::rendering::Color       FgColor;
  vtkmCanvas                   CreateCanvas() const;
  bool                         DoRenderAnnotations;
  bool                         DoRenderWorldAnnotations;
  bool                         DoRenderScreenAnnotations;
  bool                         DoRenderBackground;
  bool                         DoShading;
  vtkmCanvas                   Canvas;
  vtkm::Vec<float,3>           WorldAnnotationScale;
};

static float vtkh_default_bg_color[4] = {0.f, 0.f, 0.f, 1.f};
static float vtkh_default_fg_color[4] = {1.f, 1.f, 1.f, 1.f};

VTKM_RENDERING_EXPORT
vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           const std::string &image_name,
           float bg_color[4] = vtkh_default_bg_color,
           float fg_color[4] = vtkh_default_fg_color);

VTKM_RENDERING_EXPORT
vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           vtkm::rendering::Camera camera,
           const std::string &image_name,
           float bg_color[4] = vtkh_default_bg_color,
           float fg_color[4] = vtkh_default_fg_color);

VTKM_RENDERING_EXPORT
vtkh::Render
MakeRender(int width,
           int height,
           vtkm::rendering::Camera camera,
           vtkm::cont::PartitionedDataSet& data_set,
           const std::string &image_name,
           float bg_color[4] = vtkh_default_bg_color,
           float fg_color[4] = vtkh_default_fg_color);

} // namespace vtkh

#endif //vtkm_rendering_rendering_Render_h
