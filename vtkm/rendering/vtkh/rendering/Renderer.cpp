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
#include "Renderer.hpp"
#include <vtkm/rendering/vtkh/compositing/Compositor.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>

//#include <vtkm/rendering/vtkh/Logger.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkm/rendering/raytracing/Logger.h>


namespace vtkh {

Renderer::Renderer()
  : ColorTable("Cool to Warm"),
   DoComposite(true),
   FieldIndex(0),
    HasColorTable(true)
{
  this->Compositor  = new vtkh::Compositor();
}

Renderer::~Renderer()
{
  delete this->Compositor;
}

void
Renderer::SetShadingOn(bool)

{
  // do nothing by default;
}

void Renderer::DisableColorBar()
{
  // not all plots have color bars, so
  // we only give the option to turn it off
  this->HasColorTable = false;
}

void
Renderer::SetField(const std::string field_name)
{
  this->FieldName = field_name;
}

std::string
Renderer::GetFieldName() const
{
  return this->FieldName;
}

bool
Renderer::GetHasColorTable() const
{
  return this->HasColorTable;
}

void
Renderer::SetDoComposite(bool do_composite)
{
  this->DoComposite = do_composite;
}

void
Renderer::AddRender(vtkh::Render &render)
{
  this->Renders.push_back(render);
}

void
Renderer::SetRenders(const std::vector<vtkh::Render> &renders)
{
  this->Renders = renders;
}

int
Renderer::GetNumberOfRenders() const
{
  return static_cast<int>(this->Renders.size());
}

void
Renderer::ClearRenders()
{
  this->Renders.clear();
}

void Renderer::SetColorTable(const vtkm::cont::ColorTable &color_table)
{
  this->ColorTable = color_table;
}

vtkm::cont::ColorTable Renderer::GetColorTable() const
{
  return this->ColorTable;
}

void
Renderer::Composite(const int &num_images)
{
  //DRP: Logger
  //VTKH_DATA_OPEN("Composite");
  this->Compositor->SetCompositeMode(Compositor::Z_BUFFER_SURFACE);
  for(int i = 0; i < num_images; ++i)
  {
    float* color_buffer = &GetVTKMPointer(this->Renders[i].GetCanvas().GetColorBuffer())[0][0];
    float* depth_buffer = GetVTKMPointer(this->Renders[i].GetCanvas().GetDepthBuffer());

    int height = this->Renders[i].GetCanvas().GetHeight();
    int width = this->Renders[i].GetCanvas().GetWidth();

    this->Compositor->AddImage(color_buffer,
                           depth_buffer,
                           width,
                           height);

    Image result = this->Compositor->Composite();

#ifdef VTKM_ENABLE_MPI
    //if(vtkh::GetMPIRank() == 0)
    if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() == 0)
    {
      ImageToCanvas(result, this->Renders[i].GetCanvas(), true);
    }
#else
    ImageToCanvas(result, this->Renders[i].GetCanvas(), true);
#endif
    this->Compositor->ClearImages();
  } // for image
  //DRP: Logger
  //VTKH_DATA_CLOSE();
}

void
Renderer::PreExecute()
{
  bool range_set = this->Range.IsNonEmpty();
  CheckForRequiredField(this->FieldName);

  if(!range_set)
  {
    // we have not been given a range, so ask the data set
    vtkm::cont::ArrayHandle<vtkm::Range> ranges = vtkm::cont::FieldRangeGlobalCompute(*this->Input, this->FieldName);
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
    if(this->Range.Min == vtkm::Infinity64())
    {
      this->Range.Min = global_range.Min;
    }

    if(this->Range.Max == vtkm::NegativeInfinity64())
    {
      this->Range.Max = global_range.Max;
    }
  }

  this->Bounds = this->Input->GetGlobalBounds();
}

void
Renderer::Update()
{
  //DRP: Logger
  //VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  long long int in_cells = this->Input->GetNumberOfCells();
  //DRP: Logger
  //VTKH_DATA_ADD("input_cells", in_cells);
#endif
  PreExecute();
  DoExecute();
  PostExecute();
  //DRP: Logger
  //VTKH_DATA_CLOSE();
}

void
Renderer::PostExecute()
{
  int total_renders = static_cast<int>(this->Renders.size());
  if(this->DoComposite)
  {
    this->Composite(total_renders);
  }
}

void
Renderer::DoExecute()
{
  if(this->Mapper.get() == 0)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class";
    throw vtkm::cont::ErrorBadValue(msg);
  }

  int total_renders = static_cast<int>(this->Renders.size());


//  int num_domains = static_cast<int>(this->Input->GetGlobalNumberOfPartitions());
//  for(int dom = 0; dom < num_domains; ++dom)
  for (const auto& data_set : this->Input->GetPartitions())
  {
//    vtkm::cont::DataSet data_set;
//    vtkm::Id domain_id;
//    this->Input->GetDomain(dom, data_set, domain_id);
    if(!data_set.HasField(this->FieldName))
    {
      continue;
    }

    const vtkm::cont::UnknownCellSet &cellset = data_set.GetCellSet();
    const vtkm::cont::Field &field = data_set.GetField(this->FieldName);
    const vtkm::cont::CoordinateSystem &coords = data_set.GetCoordinateSystem();

    if(cellset.GetNumberOfCells() == 0)
    {
      continue;
    }

    for(int i = 0; i < total_renders; ++i)
    {
      if(this->Renders[i].GetShadingOn())
      {
        this->SetShadingOn(true);
      }
      else
      {
        this->SetShadingOn(false);
      }

      this->Mapper->SetActiveColorTable(this->ColorTable);

      Render::vtkmCanvas &canvas = this->Renders[i].GetCanvas();
      const vtkmCamera &camera = this->Renders[i].GetCamera();
      this->Mapper->SetCanvas(&canvas);
      this->Mapper->RenderCells(cellset,
                            coords,
                            field,
                            this->ColorTable,
                            camera,
                            this->Range);
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

std::vector<Render>
Renderer::GetRenders() const
{
  return this->Renders;
}

vtkm::Range
Renderer::GetRange() const
{
  return this->Range;
}

void
Renderer::SetRange(const vtkm::Range &range)
{
  this->Range = range;
}

void
Renderer::CheckForRequiredField(const std::string &field_name)
{
  if(this->Input == nullptr)
  {
    std::stringstream msg;
    msg<<"Cannot verify required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' because input is null.";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }

  if (!vtkh::GlobalHasField(*this->Input, field_name))
  {
    std::stringstream msg;
    msg<<"Required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' does not exist";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }
}


} // namespace vtkh
