//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/Actor.h>

#include <vtkm/Assert.h>
#include <vtkm/cont/FieldRangeGlobalCompute.h>
#include <vtkm/cont/TryExecute.h>
#include <vtkm/cont/UnknownCellSet.h>

namespace vtkm
{
namespace rendering
{

struct Actor::InternalsType
{
  vtkm::cont::ColorTable ColorTable;
  vtkm::cont::PartitionedDataSet DataSet;
  vtkm::Bounds SpatialBounds;
  std::string ScalarFieldName;
  vtkm::Range ScalarRange;

  VTKM_CONT
  InternalsType(const vtkm::cont::PartitionedDataSet& dataSet,
                const std::string& scalarFieldName,
                const vtkm::rendering::Color& color)
    : ColorTable(vtkm::Range{ 0, 1 }, color.Components, color.Components)
    , DataSet(dataSet)
    , ScalarFieldName(scalarFieldName)
  {
    this->Init();
  }

  VTKM_CONT
  InternalsType(const vtkm::cont::PartitionedDataSet& dataSet,
                const std::string& scalarFieldName,
                const vtkm::cont::ColorTable& colorTable = vtkm::cont::ColorTable::Preset::Default)
    : ColorTable(colorTable)
    , DataSet(dataSet)
    , ScalarFieldName(scalarFieldName)
  {
    this->Init();
  }

  void Init()
  {
    std::cout << __FILE__ << " " << __LINE__ << " Should this be global bounds/range??"
              << std::endl;
    this->SpatialBounds = this->DataSet.GetGlobalBounds();

    for (const auto& ds : this->DataSet.GetPartitions())
    {
      if (!ds.HasField(this->ScalarFieldName))
        throw vtkm::cont::ErrorBadValue("Field name not on data set.");
    }

    auto ranges = vtkm::cont::FieldRangeGlobalCompute(this->DataSet, this->ScalarFieldName);
    if (ranges.GetNumberOfValues() != 1)
      throw vtkm::cont::ErrorBadValue("Can only render scalar fields");

    this->ScalarRange = ranges.ReadPortal().Get(0);
  }
};

Actor::Actor(const vtkm::cont::PartitionedDataSet& dataSet,
             const std::string& scalarFieldName,
             const vtkm::rendering::Color& color)
  : Internals(new InternalsType(dataSet, scalarFieldName, color))
{
}

Actor::Actor(const vtkm::cont::DataSet& dataSet,
             const std::string& scalarFieldName,
             const vtkm::rendering::Color& color)
  : Internals(new InternalsType({ dataSet }, scalarFieldName, color))
{
}

Actor::Actor(const vtkm::cont::PartitionedDataSet& dataSet,
             const std::string& scalarFieldName,
             const vtkm::cont::ColorTable& colorTable)
  : Internals(new InternalsType(dataSet, scalarFieldName, colorTable))
{
}

Actor::Actor(const vtkm::cont::DataSet& dataSet,
             const std::string& scalarFieldName,
             const vtkm::cont::ColorTable& colorTable)
  : Internals(new InternalsType({ dataSet }, scalarFieldName, colorTable))
{
}

void Actor::Render(vtkm::rendering::Mapper& mapper,
                   vtkm::rendering::Canvas& canvas,
                   const vtkm::rendering::Camera& camera) const
{
  mapper.SetCanvas(&canvas);
  mapper.SetActiveColorTable(this->Internals->ColorTable);
  for (const auto& ds : this->Internals->DataSet)
  {
    const auto& cells = ds.GetCellSet();
    const auto& coords = ds.GetCoordinateSystem();
    const auto& field = ds.GetField(this->Internals->ScalarFieldName);
    mapper.RenderCells(
      cells, coords, field, this->Internals->ColorTable, camera, this->Internals->ScalarRange);
  }
}

const vtkm::cont::PartitionedDataSet& Actor::GetDataSet() const
{
  return this->Internals->DataSet;
}

void Actor::SetDataSet(const vtkm::cont::PartitionedDataSet& pds)
{
  this->Internals->DataSet = pds;
}

std::string Actor::GetScalarFieldName() const
{
  return this->Internals->ScalarFieldName;
}

void Actor::SetScalarFieldName(const std::string& fieldName)
{
  this->Internals->ScalarFieldName = fieldName;
}

const vtkm::cont::ColorTable& Actor::GetColorTable() const
{
  return this->Internals->ColorTable;
}

void Actor::SetColorTable(const vtkm::cont::ColorTable& colorTable)
{
  this->Internals->ColorTable = colorTable;
}

const vtkm::Range& Actor::GetScalarRange() const
{
  return this->Internals->ScalarRange;
}

const vtkm::Bounds& Actor::GetSpatialBounds() const
{
  return this->Internals->SpatialBounds;
}

void Actor::SetScalarRange(const vtkm::Range& scalarRange)
{
  this->Internals->ScalarRange = scalarRange;
}

}
} // namespace vtkm::rendering
