//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/EnvironmentTracker.h>

#include "Render.hpp"
#include <vtkm/rendering/vtkh/rendering/Annotator.hpp>
#include <vtkm/rendering/vtkh/compositing/PNGEncoder.h>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/View3D.h>

namespace vtkh
{

Render::Render()
  : Width(1024),
    Height(1024),
    DoRenderAnnotations(true),
    DoRenderWorldAnnotations(true),
    DoRenderScreenAnnotations(true),
    DoRenderBackground(true),
    DoShading(true),
    Canvas(this->Width, this->Height)
{
  this->WorldAnnotationScale[0] = 1.f;
  this->WorldAnnotationScale[1] = 1.f;
  this->WorldAnnotationScale[2] = 1.f;
}

Render::~Render()
{
}

Render::vtkmCanvas&
Render::GetCanvas()
{
  return this->Canvas;
}

vtkm::Bounds
Render::GetSceneBounds() const
{
  return this->SceneBounds;
}

void
Render::ScaleWorldAnnotations(float x, float y, float z)
{
  this->WorldAnnotationScale[0] = x;
  this->WorldAnnotationScale[1] = y;
  this->WorldAnnotationScale[2] = z;
}

vtkm::Int32
Render::GetWidth() const
{
  return this->Width;
}

vtkm::Int32
Render::GetHeight() const
{
  return this->Height;
}

void
Render::SetWidth(const vtkm::Int32 width)
{
  if(width == this->Width) return;
  this->Width = width;
  this->Canvas.ResizeBuffers(this->Width, this->Height);
}

void
Render::SetHeight(const vtkm::Int32 height)
{
  if(height == this->Height) return;
  this->Height = height;
  this->Canvas.ResizeBuffers(this->Width, this->Height);
}

void
Render::SetSceneBounds(const vtkm::Bounds &bounds)
{
  this->SceneBounds = bounds;
}

const vtkm::rendering::Camera&
Render::GetCamera() const
{
  return this->Camera;
}

void
Render::SetCamera(const vtkm::rendering::Camera &camera)
{
   this->Camera = camera;
}

void
Render::SetImageName(const std::string &name)
{
  this->ImageName = name;
}

void
Render::SetComments(const std::vector<std::string> &comments)
{
  this->Comments = comments;
}

void
Render::SetBackgroundColor(float bg_color[4])
{
  this->BgColor.Components[0] = bg_color[0];
  this->BgColor.Components[1] = bg_color[1];
  this->BgColor.Components[2] = bg_color[2];
  this->BgColor.Components[3] = bg_color[3];
}

void
Render::SetForegroundColor(float fg_color[4])
{
  this->FgColor.Components[0] = fg_color[0];
  this->FgColor.Components[1] = fg_color[1];
  this->FgColor.Components[2] = fg_color[2];
  this->FgColor.Components[3] = fg_color[3];
}

std::string
Render::GetImageName() const
{
  return this->ImageName;
}

std::vector<std::string>
Render::GetComments() const
{
  return this->Comments;
}

vtkm::rendering::Color
Render::GetBackgroundColor() const
{
  return this->BgColor;
}

void
Render::RenderWorldAnnotations()
{
  if(!this->DoRenderAnnotations) return;
  if(!this->DoRenderWorldAnnotations) return;
#ifdef VTKM_ENABLE_MPI
  //if(vtkh::GetMPIRank() != 0) return;
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() != 0)
    return;
#endif
  this->Canvas.SetBackgroundColor(this->BgColor);
  this->Canvas.SetForegroundColor(this->FgColor);

  Annotator annotator(this->Canvas, this->Camera, this->SceneBounds);
  annotator.RenderWorldAnnotations(this->WorldAnnotationScale);

}

void
Render::RenderScreenAnnotations(const std::vector<std::string> &field_names,
                                const std::vector<vtkm::Range> &ranges,
                                const std::vector<vtkm::cont::ColorTable> &colors)
{
  if(!this->DoRenderAnnotations) return;
  if(!this->DoRenderScreenAnnotations) return;
#ifdef VTKM_ENABLE_MPI
  //if(vtkh::GetMPIRank() != 0) return;
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() != 0)
    return;
#endif
  this->Canvas.SetBackgroundColor(this->BgColor);
  this->Canvas.SetForegroundColor(this->FgColor);
  if(this->DoRenderBackground) this->Canvas.BlendBackground();

  if(!this->DoRenderAnnotations) return;
  Annotator annotator(this->Canvas, this->Camera, this->SceneBounds);
  annotator.RenderScreenAnnotations(field_names, ranges, colors);
}

Render
Render::Copy() const
{
  Render copy;
  copy.Camera = this->Camera;
  copy.ImageName = this->ImageName;
  copy.SceneBounds = this->SceneBounds;
  copy.Width = this->Width;
  copy.Height = this->Height;
  copy.BgColor = this->BgColor;
  copy.FgColor = this->FgColor;
  copy.DoRenderAnnotations = this->DoRenderAnnotations;
  copy.DoRenderBackground = this->DoRenderBackground;
  copy.DoShading = this->DoShading;
  copy.Canvas = this->CreateCanvas();
  copy.WorldAnnotationScale = this->WorldAnnotationScale;
  return copy;
}

void
Render::Print() const
{
  std::cout<<"=== image name  : "<<this->ImageName<<"\n";;
  std::cout<<"=== bounds .... : "<<this->SceneBounds<<"\n";
  std::cout<<"=== width ..... : "<<this->Width<<"\n";
  std::cout<<"=== height .... : "<<this->Height<<"\n";
  std::cout<<"=== bg_color .. : "
            <<this->BgColor.Components[0]<<" "
            <<this->BgColor.Components[1]<<" "
            <<this->BgColor.Components[2]<<" "
            <<this->BgColor.Components[3]<<"\n";
  std::cout<<"=== fg_color .. : "
            <<this->FgColor.Components[0]<<" "
            <<this->FgColor.Components[1]<<" "
            <<this->FgColor.Components[2]<<" "
            <<this->FgColor.Components[3]<<"\n";
  std::cout<<"=== annotations : "
           <<(this->DoRenderAnnotations ? "On" : "Off")
           <<"\n";
  std::cout<<"=== background  : "
           <<(this->DoRenderBackground ? "On" : "Off")
           <<"\n";
  std::cout<<"=== shading ... : "
           <<(this->DoShading ? "On" : "Off")
           <<"\n";
}

void
Render::RenderBackground()
{
  if(this->DoRenderBackground)
    this->Canvas.BlendBackground();
}

Render::vtkmCanvas
Render::CreateCanvas() const
{
  Render::vtkmCanvas canvas(this->Width, this->Height);
  canvas.SetBackgroundColor(this->BgColor);
  canvas.SetForegroundColor(this->FgColor);
  canvas.Clear();
  return canvas;
}

void
Render::Save()
{
  // After rendering and compositing
  // Rank 0 contains the complete image.
#ifdef VTKM_ENABLE_MPI
  //if(vtkh::GetMPIRank() != 0) return;
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() != 0)
    return;
#endif
  float* color_buffer = &GetVTKMPointer(this->Canvas.GetColorBuffer())[0][0];
  int height = this->Canvas.GetHeight();
  int width = this->Canvas.GetWidth();
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(color_buffer, width, height, this->Comments);
  encoder.Save(this->ImageName + ".png");
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(scene_bounds);
  render.SetSceneBounds(scene_bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  camera.SetModeTo3D();

  bool is_2d = scene_bounds.Z.Min == 0. && scene_bounds.Z.Max == 0.;

  if(is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetCamera(camera);
  render.SetImageName(image_name);
  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  return render;
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           vtkm::rendering::Camera camera,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
  render.SetSceneBounds(scene_bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  camera.SetModeTo3D();

  bool is_2d = scene_bounds.Z.Min == 0. && scene_bounds.Z.Max == 0.;

  if(is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetCamera(camera);
  render.SetImageName(image_name);
  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  return render;
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::rendering::Camera camera,
           vtkm::cont::PartitionedDataSet& data_set,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
  render.SetCamera(camera);
  render.SetImageName(image_name);
  vtkm::Bounds bounds = data_set.GetGlobalBounds();
  render.SetSceneBounds(bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  bool is_2d = bounds.Z.Min == 0. && bounds.Z.Max == 0.;

  if(is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  return render;
}
} // namespace vtkh
