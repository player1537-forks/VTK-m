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
#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkh/compositing/Compositor.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>

//#include <vtkm/rendering/vtkh/Logger.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkm/rendering/raytracing/Logger.h>


namespace vtkh {

Renderer::Renderer()
{
  this->Compositor  = new vtkh::Compositor();
}

Renderer::~Renderer()
{
  delete this->Compositor;
}

void
Renderer::Composite()
{
  //DRP: Logger
  //VTKH_DATA_OPEN("Composite");

  this->Compositor->SetCompositeMode(Compositor::Z_BUFFER_SURFACE);
  std::size_t numImages = this->Plots.size();

  for(std::size_t i = 0; i < numImages; ++i)
  {
    float* color_buffer = &GetVTKMPointer(this->Plots[i].GetCanvas().GetColorBuffer())[0][0];
    float* depth_buffer = GetVTKMPointer(this->Plots[i].GetCanvas().GetDepthBuffer());

    int height = this->Plots[i].GetCanvas().GetHeight();
    int width = this->Plots[i].GetCanvas().GetWidth();

    this->Compositor->AddImage(color_buffer,
                           depth_buffer,
                           width,
                           height);

    Image result = this->Compositor->Composite();

#ifdef VTKM_ENABLE_MPI
    //if(vtkh::GetMPIRank() == 0)
    if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() == 0)
    {
      ImageToCanvas(result, this->Plots[i].GetCanvas(), true);
    }
#else
    ImageToCanvas(result, this->Plots[i].GetCanvas(), true);
#endif
    this->Compositor->ClearImages();
  } // for image
  //DRP: Logger
  //VTKH_DATA_CLOSE();
}

void
Renderer::PreExecute()
{
  bool range_set = this->GetScalarRange().IsNonEmpty();
  CheckForRequiredField(this->GetFieldName());

  if(!range_set)
  {
    // we have not been given a range, so ask the data set
    vtkm::cont::ArrayHandle<vtkm::Range> ranges = vtkm::cont::FieldRangeGlobalCompute(this->Actor.GetDataSet(), this->GetFieldName());
    int num_components = ranges.GetNumberOfValues();
    //
    // current vtkm renderers only supports single component scalar fields
    //
    if(num_components != 1)
    {
      std::stringstream msg;
      msg<<"Renderer '"<<this->GetName()<<"' cannot render a field with ";
      msg<<"'"<<num_components<<"' components. Field must be a scalar field.";
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

    if(currRange.Max == vtkm::NegativeInfinity64())
    {
      currRange.Max = global_range.Max;
      isSet = true;
    }
    if (isSet)
      this->Actor.SetScalarRange(currRange);
  }

//  this->Bounds = this->Actor.GetDataSet().GetGlobalBounds();
}

void
Renderer::Update()
{
  this->PreExecute();
  this->DoExecute();
  this->PostExecute();
}

void
Renderer::PostExecute()
{
  if(this->DoComposite)
    this->Composite();
}

void
Renderer::DoExecute()
{
  if(this->Mapper.get() == nullptr)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class";
    throw vtkm::cont::ErrorBadValue(msg);
  }

  int total_renders = static_cast<int>(this->Plots.size());


//  int num_domains = static_cast<int>(this->Input->GetGlobalNumberOfPartitions());
//  for(int dom = 0; dom < num_domains; ++dom)
  for (const auto& data_set : this->Actor.GetDataSet())
  {
//    vtkm::cont::DataSet data_set;
//    vtkm::Id domain_id;
//    this->Input->GetDomain(dom, data_set, domain_id);
    if(!data_set.HasField(this->GetFieldName()))
    {
      continue;
    }

    const vtkm::cont::UnknownCellSet &cellset = data_set.GetCellSet();
    const vtkm::cont::Field &field = data_set.GetField(this->GetFieldName());
    const vtkm::cont::CoordinateSystem &coords = data_set.GetCoordinateSystem();

    if(cellset.GetNumberOfCells() == 0)
    {
      continue;
    }

    for(int i = 0; i < total_renders; ++i)
    {
      if(this->Plots[i].GetShadingOn())
      {
        this->SetShadingOn(true);
      }
      else
      {
        this->SetShadingOn(false);
      }

      this->Mapper->SetActiveColorTable(this->GetColorTable());

      Plot::vtkmCanvas &canvas = this->Plots[i].GetCanvas();
      const auto& camera = this->Plots[i].GetCamera();
      this->Mapper->SetCanvas(&canvas);
      this->Mapper->RenderCells(cellset,
                                coords,
                                field,
                                this->GetColorTable(),
                                camera,
                                this->GetScalarRange());
    }
  }


}

void
Renderer::ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth)
{
  const int width = canvas.GetWidth();
  const int height = canvas.GetHeight();
  const int size = width * height;
  const int color_size = size * 4;
  float* color_buffer = &GetVTKMPointer(canvas.GetColorBuffer())[0][0];
  float one_over_255 = 1.f / 255.f;
#ifdef VTKH_OPENMP_ENABLED
  #pragma omp parallel for
#endif
  for(int i = 0; i < color_size; ++i)
  {
    color_buffer[i] = static_cast<float>(image.Pixels[i]) * one_over_255;
  }

  float* depth_buffer = GetVTKMPointer(canvas.GetDepthBuffer());
  if(get_depth) memcpy(depth_buffer, &image.Depths[0], sizeof(float) * size);
}

void
Renderer::CheckForRequiredField(const std::string &field_name)
{
  if (this->Actor.GetDataSet().GetNumberOfPartitions() == 0)
  {
    std::stringstream msg;
    msg<<"Cannot verify required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' because input is null.";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }

  if (!vtkh::GlobalHasField(this->Actor.GetDataSet(), field_name))
  {
    std::stringstream msg;
    msg<<"Required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' does not exist";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }
}

} // namespace vtkh
