//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_vtkh_rendering_Plot_h
#define vtkm_rendering_vtkh_rendering_Plot_h

#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/CanvasRayTracer.h>

#include <vector>

namespace vtkm
{
namespace rendering_new
{

//
// A Render contains the information needed to create a single image.
// There are 'n' canvases that matches the number of domains in the
// data set. It is possible to chain multiple plots together that
// are rendering separate data, i.e. the result of different data
// transformations, to handle this we keep track of the domain ids
// that each canvas is associated with.
//

class Renderer;

class VTKM_RENDERING_NEW_EXPORT Plot
{
public:
  typedef vtkm::rendering::CanvasRayTracer vtkmCanvas;

  Plot();
  ~Plot();
  Plot Copy() const;
  vtkmCanvas& GetCanvas() { return this->Canvas; }
  const vtkm::rendering::Camera& GetCamera() const { return this->Camera; }
  std::string GetImageName() const { return this->ImageName; }
  std::vector<std::string> GetComments() const { return this->Comments; }
  vtkm::Bounds GetSceneBounds() const { return this->SceneBounds; }
  vtkm::Int32 GetHeight() const { return this->Height; }
  vtkm::Int32 GetWidth() const { return this->Width; }
  vtkm::rendering::Color GetBackgroundColor() const { return this->BgColor; }
  bool GetShadingOn() const { return this->DoShading; }
  void Print() const;

  void AddRenderer(vtkm::rendering_new::Renderer* renderer);

  void SetPlotAnnotations(bool on) { this->DoPlotAnnotations = on; }
  void SetPlotWorldAnnotations(bool on) { this->DoPlotWorldAnnotations = on; }
  void SetPlotScreenAnnotations(bool on) { this->DoPlotScreenAnnotations = on; }
  void SetPlotBackground(bool on) { this->DoPlotBackground = on; }
  void ScaleWorldAnnotations(float x, float y, float z);
  void SetWidth(const vtkm::Int32 width);
  void SetHeight(const vtkm::Int32 height);
  void SetSceneBounds(const vtkm::Bounds& bounds) { this->SceneBounds = bounds; }
  void SetCamera(const vtkm::rendering::Camera& camera) { this->Camera = camera; }
  void SetImageName(const std::string& name) { this->ImageName = name; }
  void SetComments(const std::vector<std::string>& comments) { this->Comments = comments; }
  void SetBackgroundColor(float bg_color[4]);
  void SetForegroundColor(float fg_color[4]);
  void SetShadingOn(bool on) { this->DoShading = on; }
  void RenderWorldAnnotations();
  void RenderBackground();
  void RenderScreenAnnotations(const std::vector<std::string>& field_names,
                               const std::vector<vtkm::Range>& ranges,
                               const std::vector<vtkm::cont::ColorTable>& colors);
  void SyncDepth();
  void Save();
  void Render();

  std::vector<vtkm::rendering_new::Renderer*> Renderers;

protected:
  vtkm::rendering::Camera Camera;
  std::string ImageName;
  std::vector<std::string> Comments;
  vtkm::Bounds SceneBounds;
  vtkm::Int32 Width;
  vtkm::Int32 Height;
  vtkm::rendering::Color BgColor;
  vtkm::rendering::Color FgColor;
  vtkmCanvas CreateCanvas() const;
  bool DoPlotAnnotations;
  bool DoPlotWorldAnnotations;
  bool DoPlotScreenAnnotations;
  bool DoPlotBackground;
  bool DoShading;
  vtkmCanvas Canvas;
  vtkm::Vec<float, 3> WorldAnnotationScale;
  bool HasVolume = false;
  bool HasMesh = false;
};

static float vtkh_default_bg_color[4] = { 0.f, 0.f, 0.f, 1.f };
static float vtkh_default_fg_color[4] = { 1.f, 1.f, 1.f, 1.f };

VTKM_RENDERING_NEW_EXPORT
vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::Bounds scene_bounds,
                                   const std::string& image_name,
                                   float bg_color[4] = vtkh_default_bg_color,
                                   float fg_color[4] = vtkh_default_fg_color);

VTKM_RENDERING_NEW_EXPORT
vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::Bounds scene_bounds,
                                   vtkm::rendering::Camera camera,
                                   const std::string& image_name,
                                   float bg_color[4] = vtkh_default_bg_color,
                                   float fg_color[4] = vtkh_default_fg_color);

VTKM_RENDERING_NEW_EXPORT
vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::rendering::Camera camera,
                                   vtkm::cont::PartitionedDataSet& data_set,
                                   const std::string& image_name,
                                   float bg_color[4] = vtkh_default_bg_color,
                                   float fg_color[4] = vtkh_default_fg_color);

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_vtkh_rendering_Plot_h
