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
#include <vtkm/cont/FieldRangeGlobalCompute.h>

#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/rendering_new/Renderer.h>
#include <vtkm/rendering_new/compositing/Compositor.h>

namespace vtkm
{
namespace rendering_new
{

Renderer::Renderer()
{
  this->Compositor = new vtkm::rendering_new::Compositor();
}

Renderer::~Renderer()
{
  delete this->Compositor;
}

void Renderer::Composite(vtkm::rendering_new::Plot& plot)
{
  this->Compositor->SetCompositeMode(Compositor::Z_BUFFER_SURFACE);

  float* color_buffer = &(plot.GetCanvas().GetColorBuffer().WritePortal().GetArray()[0][0]);
  float* depth_buffer = plot.GetCanvas().GetDepthBuffer().WritePortal().GetArray();

  int height = plot.GetCanvas().GetHeight();
  int width = plot.GetCanvas().GetWidth();

  this->Compositor->AddImage(color_buffer, depth_buffer, width, height);

  Image result = this->Compositor->Composite();

#ifdef VTKM_ENABLE_MPI
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() == 0)
  {
    ImageToCanvas(result, plot.GetCanvas(), true);
  }
#else
  ImageToCanvas(result, plot.GetCanvas(), true);
#endif
  this->Compositor->ClearImages();
}

void Renderer::PreExecute(vtkm::rendering_new::Plot& vtkmNotUsed(plot))
{
  bool range_set = this->GetScalarRange().IsNonEmpty();
  CheckForRequiredField(this->GetFieldName());

  if (!range_set)
  {
    // we have not been given a range, so ask the data set
    vtkm::cont::ArrayHandle<vtkm::Range> ranges =
      vtkm::cont::FieldRangeGlobalCompute(this->Actor.GetDataSet(), this->GetFieldName());
    int num_components = ranges.GetNumberOfValues();
    //
    // current vtkm renderers only supports single component scalar fields
    //
    if (num_components != 1)
    {
      std::stringstream msg;
      msg << "Renderer '" << this->GetName() << "' cannot render a field with ";
      msg << "'" << num_components << "' components. Field must be a scalar field.";
      throw vtkm::cont::ErrorBadValue(msg.str());
    }

    vtkm::Range global_range = ranges.ReadPortal().Get(0);
    // a min or max may be been set by the user, check to see
    auto currRange = this->Actor.GetScalarRange();
    bool isSet = false;
    if (currRange.Min == vtkm::Infinity64())
    {
      currRange.Min = global_range.Min;
      isSet = true;
    }

    if (currRange.Max == vtkm::NegativeInfinity64())
    {
      currRange.Max = global_range.Max;
      isSet = true;
    }
    if (isSet)
      this->Actor.SetScalarRange(currRange);
  }

  //  this->Bounds = this->Actor.GetDataSet().GetGlobalBounds();
}

void Renderer::Update(vtkm::rendering_new::Plot& plot)
{
  this->PreExecute(plot);
  this->DoExecute(plot);
  this->PostExecute(plot);
}

void Renderer::PostExecute(vtkm::rendering_new::Plot& plot)
{
  if (this->DoComposite)
    this->Composite(plot);
}

void Renderer::DoExecute(vtkm::rendering_new::Plot& plot)
{
  if (this->Mapper.get() == nullptr)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class";
    throw vtkm::cont::ErrorBadValue(msg);
  }

  //  int num_domains = static_cast<int>(this->Input->GetGlobalNumberOfPartitions());
  //  for(int dom = 0; dom < num_domains; ++dom)
  for (const auto& data_set : this->Actor.GetDataSet())
  {
    //    vtkm::cont::DataSet data_set;
    //    vtkm::Id domain_id;
    //    this->Input->GetDomain(dom, data_set, domain_id);
    if (!data_set.HasField(this->GetFieldName()))
    {
      continue;
    }

    const vtkm::cont::UnknownCellSet& cellset = data_set.GetCellSet();
    const vtkm::cont::Field& field = data_set.GetField(this->GetFieldName());
    const vtkm::cont::CoordinateSystem& coords = data_set.GetCoordinateSystem();

    if (cellset.GetNumberOfCells() == 0)
    {
      continue;
    }

    this->SetShadingOn(plot.GetShadingOn());

    this->Mapper->SetActiveColorTable(this->GetColorTable());

    Plot::vtkmCanvas& canvas = plot.GetCanvas();
    const auto& camera = plot.GetCamera();
    this->Mapper->SetCanvas(&canvas);
    this->Mapper->RenderCells(
      cellset, coords, field, this->GetColorTable(), camera, this->GetScalarRange());
  }
}

void Renderer::ImageToCanvas(Image& image, vtkm::rendering::Canvas& canvas, bool get_depth)
{
  const int width = canvas.GetWidth();
  const int height = canvas.GetHeight();
  const int size = width * height;
  const int color_size = size * 4;
  float* color_buffer = &(canvas.GetColorBuffer().WritePortal().GetArray()[0][0]);
  float one_over_255 = 1.f / 255.f;
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int i = 0; i < color_size; ++i)
  {
    color_buffer[i] = static_cast<float>(image.Pixels[i]) * one_over_255;
  }

  float* depth_buffer = canvas.GetDepthBuffer().WritePortal().GetArray();
  if (get_depth)
    memcpy(depth_buffer, &image.Depths[0], sizeof(float) * size);
}

void Renderer::CheckForRequiredField(const std::string& field_name)
{
  if (this->Actor.GetDataSet().GetNumberOfPartitions() == 0)
  {
    std::stringstream msg;
    msg << "Cannot verify required field '" << field_name;
    msg << "' for vkth filter '" << this->GetName() << "' because input is null.";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }

  bool hasFieldLocal = true, hasFieldGlobal = false;
  for (const auto& ds : this->Actor.GetDataSet())
    if (!ds.HasField(field_name))
    {
      hasFieldLocal = false;
      break;
    }
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  vtkmdiy::mpi::all_reduce(comm, hasFieldLocal, hasFieldGlobal, std::logical_and<bool>());
#else
  hasFieldGlobal = hasFieldLocal;
#endif

  if (!hasFieldGlobal)
  {
    std::stringstream msg;
    msg << "Required field '" << field_name;
    msg << "' for vkth filter '" << this->GetName() << "' does not exist";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }
}

}
} // namespace vtkm::rendering_new
