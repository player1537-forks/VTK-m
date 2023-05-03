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

#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/vtkh/compositing/PNGEncoder.h>
#include <vtkm/rendering/vtkh/rendering/Annotator.hpp>
#include <vtkm/rendering/vtkh/rendering/MeshRenderer.hpp>
#include <vtkm/rendering/vtkh/rendering/Plot.h>
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkh/rendering/VolumeRenderer.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkh
{

Plot::Plot()
  : Width(1024)
  , Height(1024)
  , DoPlotAnnotations(true)
  , DoPlotWorldAnnotations(true)
  , DoPlotScreenAnnotations(true)
  , DoPlotBackground(true)
  , DoShading(true)
  , Canvas(this->Width, this->Height)
{
  this->WorldAnnotationScale[0] = 1.f;
  this->WorldAnnotationScale[1] = 1.f;
  this->WorldAnnotationScale[2] = 1.f;
}

Plot::~Plot() {}

void Plot::AddRenderer(vtkh::Renderer* renderer)
{
  //Redo this... Implicit ordering is asking for trouble....

  bool isVolume = (dynamic_cast<vtkh::VolumeRenderer*>(renderer) != nullptr);
  bool isMesh = (dynamic_cast<vtkh::MeshRenderer*>(renderer) != nullptr);

  if (isVolume)
  {
    if (this->HasVolume)
    {
      throw vtkm::cont::ErrorBadValue("Scenes only support a single volume plot");
    }
    this->HasVolume = true;

    // make sure that the volume render is last
    this->Renderers.push_back(renderer);
  }
  else if (isMesh)
  {
    // make sure that the mesh plot is last
    // and before the volume plot
    if (this->HasVolume)
    {
      if (this->Renderers.size() == 1)
        this->Renderers.insert(this->Renderers.begin(), renderer);
      else
      {
        auto it = this->Renderers.end();
        it--;
        it--;
        this->Renderers.insert(it, renderer);
      }
    }
    else
    {
      this->Renderers.push_back(renderer);
    }
  }
  else
  {
    this->Renderers.insert(this->Renderers.begin(), renderer);
  }
}

Plot::vtkmCanvas& Plot::GetCanvas()
{
  return this->Canvas;
}

vtkm::Bounds Plot::GetSceneBounds() const
{
  return this->SceneBounds;
}

void Plot::ScaleWorldAnnotations(float x, float y, float z)
{
  this->WorldAnnotationScale[0] = x;
  this->WorldAnnotationScale[1] = y;
  this->WorldAnnotationScale[2] = z;
}

vtkm::Int32 Plot::GetWidth() const
{
  return this->Width;
}

vtkm::Int32 Plot::GetHeight() const
{
  return this->Height;
}

void Plot::SetWidth(const vtkm::Int32 width)
{
  if (width == this->Width)
    return;
  this->Width = width;
  this->Canvas.ResizeBuffers(this->Width, this->Height);
}

void Plot::SetHeight(const vtkm::Int32 height)
{
  if (height == this->Height)
    return;
  this->Height = height;
  this->Canvas.ResizeBuffers(this->Width, this->Height);
}

void Plot::SetSceneBounds(const vtkm::Bounds& bounds)
{
  this->SceneBounds = bounds;
}

const vtkm::rendering::Camera& Plot::GetCamera() const
{
  return this->Camera;
}

void Plot::SetCamera(const vtkm::rendering::Camera& camera)
{
  this->Camera = camera;
}

void Plot::SetImageName(const std::string& name)
{
  this->ImageName = name;
}

void Plot::SetComments(const std::vector<std::string>& comments)
{
  this->Comments = comments;
}

void Plot::SetBackgroundColor(float bg_color[4])
{
  this->BgColor.Components[0] = bg_color[0];
  this->BgColor.Components[1] = bg_color[1];
  this->BgColor.Components[2] = bg_color[2];
  this->BgColor.Components[3] = bg_color[3];
}

void Plot::SetForegroundColor(float fg_color[4])
{
  this->FgColor.Components[0] = fg_color[0];
  this->FgColor.Components[1] = fg_color[1];
  this->FgColor.Components[2] = fg_color[2];
  this->FgColor.Components[3] = fg_color[3];
}

std::string Plot::GetImageName() const
{
  return this->ImageName;
}

std::vector<std::string> Plot::GetComments() const
{
  return this->Comments;
}

vtkm::rendering::Color Plot::GetBackgroundColor() const
{
  return this->BgColor;
}

void Plot::RenderWorldAnnotations()
{
  if (!this->DoPlotAnnotations)
    return;
  if (!this->DoPlotWorldAnnotations)
    return;
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

void Plot::RenderScreenAnnotations(const std::vector<std::string>& field_names,
                                   const std::vector<vtkm::Range>& ranges,
                                   const std::vector<vtkm::cont::ColorTable>& colors)
{
  if (!this->DoPlotAnnotations)
    return;
  if (!this->DoPlotScreenAnnotations)
    return;
#ifdef VTKM_ENABLE_MPI
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() != 0)
    return;
#endif
  this->Canvas.SetBackgroundColor(this->BgColor);
  this->Canvas.SetForegroundColor(this->FgColor);
  if (this->DoPlotBackground)
    this->Canvas.BlendBackground();

  if (!this->DoPlotAnnotations)
    return;
  Annotator annotator(this->Canvas, this->Camera, this->SceneBounds);
  annotator.RenderScreenAnnotations(field_names, ranges, colors);
}

Plot Plot::Copy() const
{
  Plot copy;
  copy.Camera = this->Camera;
  copy.ImageName = this->ImageName;
  copy.SceneBounds = this->SceneBounds;
  copy.Width = this->Width;
  copy.Height = this->Height;
  copy.BgColor = this->BgColor;
  copy.FgColor = this->FgColor;
  copy.DoPlotAnnotations = this->DoPlotAnnotations;
  copy.DoPlotBackground = this->DoPlotBackground;
  copy.DoShading = this->DoShading;
  copy.Canvas = this->CreateCanvas();
  copy.WorldAnnotationScale = this->WorldAnnotationScale;
  return copy;
}

void Plot::Print() const
{
  std::cout << "=== image name  : " << this->ImageName << "\n";
  ;
  std::cout << "=== bounds .... : " << this->SceneBounds << "\n";
  std::cout << "=== width ..... : " << this->Width << "\n";
  std::cout << "=== height .... : " << this->Height << "\n";
  std::cout << "=== bg_color .. : " << this->BgColor.Components[0] << " "
            << this->BgColor.Components[1] << " " << this->BgColor.Components[2] << " "
            << this->BgColor.Components[3] << "\n";
  std::cout << "=== fg_color .. : " << this->FgColor.Components[0] << " "
            << this->FgColor.Components[1] << " " << this->FgColor.Components[2] << " "
            << this->FgColor.Components[3] << "\n";
  std::cout << "=== annotations : " << (this->DoPlotAnnotations ? "On" : "Off") << "\n";
  std::cout << "=== background  : " << (this->DoPlotBackground ? "On" : "Off") << "\n";
  std::cout << "=== shading ... : " << (this->DoShading ? "On" : "Off") << "\n";
}

void Plot::RenderBackground()
{
  if (this->DoPlotBackground)
    this->Canvas.BlendBackground();
}

Plot::vtkmCanvas Plot::CreateCanvas() const
{
  Plot::vtkmCanvas canvas(this->Width, this->Height);
  canvas.SetBackgroundColor(this->BgColor);
  canvas.SetForegroundColor(this->FgColor);
  canvas.Clear();
  return canvas;
}

void Plot::SyncDepth()
{
#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MPI_Comm comm = vtkmdiy::mpi::mpi_cast(diy_comm.handle());

  vtkm::rendering::Canvas& canvas = this->GetCanvas();
  const int image_size = canvas.GetWidth() * canvas.GetHeight();
  float* depth_ptr = GetVTKMPointer(canvas.GetDepthBuffer());
  MPI_Bcast(depth_ptr, image_size, MPI_FLOAT, 0, comm);
#endif
}

void Plot::Save()
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

vtkh::Plot MakePlot(int width,
                    int height,
                    vtkm::Bounds scene_bounds,
                    const std::string& image_name,
                    float bg_color[4],
                    float fg_color[4])
{
  vtkh::Plot render;
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

  if (is_2d)
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

vtkh::Plot MakePlot(int width,
                    int height,
                    vtkm::Bounds scene_bounds,
                    vtkm::rendering::Camera camera,
                    const std::string& image_name,
                    float bg_color[4],
                    float fg_color[4])
{
  vtkh::Plot render;
  render.SetSceneBounds(scene_bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  camera.SetModeTo3D();

  bool is_2d = scene_bounds.Z.Min == 0. && scene_bounds.Z.Max == 0.;

  if (is_2d)
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

vtkh::Plot MakePlot(int width,
                    int height,
                    vtkm::rendering::Camera camera,
                    vtkm::cont::PartitionedDataSet& data_set,
                    const std::string& image_name,
                    float bg_color[4],
                    float fg_color[4])
{
  vtkh::Plot render;
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

  if (is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  return render;
}
} // namespace vtkh
