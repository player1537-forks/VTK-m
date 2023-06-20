//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_rendering_Actor_h
#define vtk_m_rendering_Actor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

#include <memory>

namespace vtkm
{
namespace rendering
{

class VTKM_RENDERING_EXPORT Actor
{
public:
  Actor() {} //REMOVE ME.
  Actor(const vtkm::cont::PartitionedDataSet& dataSet,
        const std::string& scalarFieldName,
        const vtkm::cont::ColorTable& colorTable);

  Actor(const vtkm::cont::DataSet& dataSet,
        const std::string& scalarFieldName,
        const vtkm::cont::ColorTable& colorTable);

  Actor(const vtkm::cont::PartitionedDataSet& dataSet,
        const std::string& scalarFieldName,
        const vtkm::rendering::Color& color);

  Actor(const vtkm::cont::DataSet& dataSet,
        const std::string& scalarFieldName,
        const vtkm::rendering::Color& color);

  void Render(vtkm::rendering::Mapper& mapper,
              vtkm::rendering::Canvas& canvas,
              const vtkm::rendering::Camera& camera) const;

  const vtkm::cont::PartitionedDataSet& GetDataSet() const;
  void SetDataSet(const vtkm::cont::PartitionedDataSet& pds);

  void SetScalarFieldName(const std::string& fieldName);
  std::string GetScalarFieldName() const;

  const vtkm::cont::ColorTable& GetColorTable() const;
  void SetColorTable(const vtkm::cont::ColorTable& colorTable);

  const vtkm::Range& GetScalarRange() const;

  const vtkm::Bounds& GetSpatialBounds() const;

  void SetScalarRange(const vtkm::Range& scalarRange);

private:
  struct InternalsType;
  std::shared_ptr<InternalsType> Internals;
};

}
} //namespace vtkm::rendering

#endif //vtk_m_rendering_Actor_h
