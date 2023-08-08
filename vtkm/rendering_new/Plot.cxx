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
#include <vtkm/rendering_new/Annotator.h>
#include <vtkm/rendering_new/MeshRenderer.h>
#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/VolumeRenderer.h>
#include <vtkm/rendering_new/compositing/PNGEncoder.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace rendering_new
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

void Plot::AddRenderer(vtkm::rendering_new::Renderer* renderer)
{
  //Redo this... Implicit ordering is asking for trouble....

  bool isVolume = (dynamic_cast<vtkm::rendering_new::VolumeRenderer*>(renderer) != nullptr);
  bool isMesh = (dynamic_cast<vtkm::rendering_new::MeshRenderer*>(renderer) != nullptr);

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

void Plot::ScaleWorldAnnotations(float x, float y, float z)
{
  this->WorldAnnotationScale[0] = x;
  this->WorldAnnotationScale[1] = y;
  this->WorldAnnotationScale[2] = z;
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

void Plot::RenderWorldAnnotations()
{
  if (!this->DoPlotAnnotations)
    return;
  if (!this->DoPlotWorldAnnotations)
    return;
#ifdef VTKM_ENABLE_MPI
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
  copy.Renderers = this->Renderers;
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
  float* depth_ptr = canvas.GetDepthBuffer().WritePortal().GetArray();
  MPI_Bcast(depth_ptr, image_size, MPI_FLOAT, 0, comm);
#endif
}

void Plot::Render()
{
  this->GetCanvas().Clear();
  if (this->Renderers.empty())
    return;

  std::vector<vtkm::Range> ranges;
  std::vector<std::string> fieldNames;
  std::vector<vtkm::cont::ColorTable> colorTables;
  for (const auto renderer : this->Renderers)
  {
    if (renderer->GetHasColorTable())
    {
      ranges.push_back(renderer->GetScalarRange());
      fieldNames.push_back(renderer->GetFieldName());
      colorTables.push_back(renderer->GetColorTable());
    }
  }

  std::size_t totalNum = this->Renderers.size();
  std::size_t numOpaque = totalNum;
  if (this->HasVolume)
    numOpaque--;

  bool syncDepth = (numOpaque > 0);

  //Pass 1: opaque renderers.
  for (std::size_t i = 0; i < numOpaque; i++)
  {
    auto renderer = this->Renderers[i];
    //composite if last.
    renderer->SetDoComposite(i == numOpaque - 1);

    renderer->Update(*this);
  }

  // Pass 2: volume render
  if (this->HasVolume)
  {
    if (syncDepth)
      this->SyncDepth();
    auto renderer = this->Renderers[totalNum - 1];
    renderer->SetDoComposite(true);
    renderer->Update(*this);
  }

  this->RenderWorldAnnotations();
  this->RenderScreenAnnotations(fieldNames, ranges, colorTables);
  this->RenderBackground();
  this->Save();
}

void Plot::Save()
{
  // After rendering and compositing
  // Rank 0 contains the complete image.
#ifdef VTKM_ENABLE_MPI
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() != 0)
    return;
#endif

  float* color_buffer = &(this->Canvas.GetColorBuffer().WritePortal().GetArray()[0][0]);
  int height = this->Canvas.GetHeight();
  int width = this->Canvas.GetWidth();
  vtkm::rendering_new::PNGEncoder encoder;
  encoder.Encode(color_buffer, width, height, this->Comments);
  encoder.Save(this->ImageName + ".png");
}

vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::Bounds scene_bounds,
                                   const std::string& image_name,
                                   float bg_color[4],
                                   float fg_color[4])
{
  vtkm::rendering_new::Plot render;
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

vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::Bounds scene_bounds,
                                   vtkm::rendering::Camera camera,
                                   const std::string& image_name,
                                   float bg_color[4],
                                   float fg_color[4])
{
  vtkm::rendering_new::Plot render;
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

vtkm::rendering_new::Plot MakePlot(int width,
                                   int height,
                                   vtkm::rendering::Camera camera,
                                   vtkm::cont::PartitionedDataSet& data_set,
                                   const std::string& image_name,
                                   float bg_color[4],
                                   float fg_color[4])
{
  vtkm::rendering_new::Plot render;
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

}
} // namespace vtkm::rendering_new
