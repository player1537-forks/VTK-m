//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <memory>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperPoint.h>
#include <vtkm/rendering_new/PointRenderer.h>

namespace vtkm
{
namespace rendering_new
{

//this needs to be moved someplace more general...
namespace detail
{
bool IsStructured(const vtkm::cont::UnknownCellSet& cell_set, int& topo_dims)
{
  bool is_structured = false;
  topo_dims = -1;

  if (cell_set.IsType<vtkm::cont::CellSetStructured<1>>())
  {
    is_structured = true;
    topo_dims = 1;
  }
  else if (cell_set.IsType<vtkm::cont::CellSetStructured<2>>())
  {
    is_structured = true;
    topo_dims = 2;
  }
  else if (cell_set.IsType<vtkm::cont::CellSetStructured<3>>())
  {
    is_structured = true;
    topo_dims = 3;
  }

  return is_structured;
}

bool IsSingleCellShape(const vtkm::cont::UnknownCellSet& cell_set, vtkm::UInt8& shape_id)
{
  int dims;
  shape_id = 0;
  bool is_single_shape = false;
  if (IsStructured(cell_set, dims))
  {
    is_single_shape = true;
    shape_id = 12;
  }
  else
  {
    // we have an explicit cell set so we have to look deeper
    if (cell_set.IsType<vtkm::cont::CellSetSingleType<>>())
    {
      vtkm::cont::CellSetSingleType<> single =
        cell_set.AsCellSet<vtkm::cont::CellSetSingleType<>>();
      is_single_shape = true;
      shape_id = single.GetCellShape(0);
    }
    else if (cell_set.IsType<vtkm::cont::CellSetExplicit<>>())
    {
      vtkm::cont::CellSetExplicit<> exp = cell_set.AsCellSet<vtkm::cont::CellSetExplicit<>>();
      const vtkm::cont::ArrayHandle<vtkm::UInt8> shapes =
        exp.GetShapesArray(vtkm::TopologyElementTagCell(), vtkm::TopologyElementTagPoint());

      vtkm::UInt8 init_min = 255;
      vtkm::UInt8 min = vtkm::cont::Algorithm::Reduce(shapes, init_min, vtkm::Minimum());

      vtkm::UInt8 init_max = 0;
      vtkm::UInt8 max = vtkm::cont::Algorithm::Reduce(shapes, init_max, vtkm::Maximum());
      if (min == max)
      {
        is_single_shape = true;
        shape_id = max;
      }
    }
  }
  return is_single_shape;
}
}

PointRenderer::PointRenderer()
{
  typedef vtkm::rendering::MapperPoint TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  this->Mapper = mapper;
}

bool PointRenderer::IsPointMesh() const
{
  const auto& pds = this->Actor.GetDataSet();
  vtkm::Id totNumCells = pds.GetGlobalNumberOfCells();
  if (totNumCells == 0)
    return false;

  bool isPoints = true;
  for (const auto& ds : pds.GetPartitions())
  {
    vtkm::UInt8 shape_type;
    bool single_type = vtkm::rendering_new::detail::IsSingleCellShape(ds.GetCellSet(), shape_type);

    if (ds.GetCellSet().GetNumberOfCells() > 0)
    {
      isPoints = (single_type && (shape_type == 1)) && isPoints;
    }
  }
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  bool globalIsPoints;
  vtkmdiy::mpi::all_reduce(comm, isPoints, globalIsPoints, std::logical_and<bool>());
  return globalIsPoints;
}

std::shared_ptr<vtkm::rendering::Canvas> PointRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

void PointRenderer::SetBaseRadius(vtkm::Float32 radius)
{
  this->BaseRadius = radius;
  this->RadiusSet = true;
}

void PointRenderer::PreExecute(vtkm::rendering_new::Plot& plot)
{
  Renderer::PreExecute(plot);

  typedef vtkm::rendering::MapperPoint MapperType;
  std::shared_ptr<MapperType> meshMapper = std::dynamic_pointer_cast<MapperType>(this->Mapper);

  if (this->UseNodes)
  {
    meshMapper->UseNodes();
  }
  else
  {
    meshMapper->UseCells();
  }

  vtkm::Float32 radius = this->BaseRadius;
  if (this->RadiusSet)
  {
    meshMapper->SetRadius(this->BaseRadius);
  }
  else
  {
    vtkm::Bounds coordBounds = this->Actor.GetDataSet().GetGlobalBounds();
    // set a default radius
    vtkm::Float64 lx = coordBounds.X.Length();
    vtkm::Float64 ly = coordBounds.Y.Length();
    vtkm::Float64 lz = coordBounds.Z.Length();
    vtkm::Float64 mag = vtkm::Sqrt(lx * lx + ly * ly + lz * lz);
    // same as used in vtk ospray
    constexpr vtkm::Float64 heuristic = 1000.;
    radius = static_cast<vtkm::Float32>(mag / heuristic);
    // we likely have a data set with no cells so just set some radius
    if (radius == 0.f)
    {
      radius = 0.00001f;
    }
    meshMapper->SetRadius(radius);
  }

  if (!this->UseNodes && this->IsPointMesh() && this->UsePointMerging)
  {
    vtkm::filter::clean_grid::CleanGrid gridCleaner;
    gridCleaner.SetMergePoints(true);
    gridCleaner.SetFastMerge(true);


    vtkm::Float32 maxRadius = radius;
    if (this->UseVariableRadius)
    {
      maxRadius = radius + radius * this->DeltaRadius;
    }
    gridCleaner.SetTolerance(maxRadius * this->RadiusMult);
    gridCleaner.SetActiveField(this->Actor.GetScalarFieldName());
    auto output = gridCleaner.Execute(this->Actor.GetDataSet());
    this->Actor.SetDataSet(output);
    /*
    ParticleMerging  merger;
    merger.SetInput(this->Input);
    merger.SetField(this->FieldName);
    merger.SetRadius(
    merger.Update();
    this->Input = merger.GetOutput();
    this->DeleteInput = true;
*/
  }

  meshMapper->UseVariableRadius(this->UseVariableRadius);
  meshMapper->SetRadiusDelta(this->DeltaRadius);
}

void PointRenderer::PostExecute(vtkm::rendering_new::Plot& plot)
{
  Renderer::PostExecute(plot);
  /*
  if(this->DeleteInput)
  {
//    delete this->Input;
//    this->Input = nullptr;
  }
  */
}

}
} // namespace vtkm::rendering_new
