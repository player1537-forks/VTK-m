//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_cont_CellLocatorUniformBins_h
#define vtk_m_cont_CellLocatorUniformBins_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellSetList.h>

#include <vtkm/cont/internal/CellLocatorBase.h>

#include <vtkm/exec/CellLocatorMultiplexer.h>
#include <vtkm/exec/CellLocatorUniformBins.h>


namespace vtkm
{
namespace cont
{

class VTKM_CONT_EXPORT CellLocatorUniformBins
  : public vtkm::cont::internal::CellLocatorBase<CellLocatorUniformBins>
{
  using Superclass = vtkm::cont::internal::CellLocatorBase<CellLocatorUniformBins>;

  template <typename CellSetCont>
  using CellSetContToExec =
    typename CellSetCont::template ExecConnectivityType<vtkm::TopologyElementTagCell,
                                                        vtkm::TopologyElementTagPoint>;

public:
  using SupportedCellSets = VTKM_DEFAULT_CELL_SET_LIST;

  using CellExecObjectList = vtkm::ListTransform<SupportedCellSets, CellSetContToExec>;
  using CellLocatorExecList =
    vtkm::ListTransform<CellExecObjectList, vtkm::exec::CellLocatorUniformBins>;

  using ExecObjType = vtkm::ListApply<CellLocatorExecList, vtkm::exec::CellLocatorMultiplexer>;
  using LastCell = typename ExecObjType::LastCell;

  CellLocatorUniformBins() {}
  void SetDims(const vtkm::Id3& dims) { this->UniformDims = dims; }
  vtkm::Id3 GetDims() const { return this->UniformDims; }

  void PrintSummary(std::ostream& out) const;

public:
  ExecObjType PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                  vtkm::cont::Token& token) const;

private:
  friend Superclass;
  VTKM_CONT void Build();

  vtkm::Vec3f InvSpacing;
  vtkm::Vec3f MaxPoint;
  vtkm::Vec3f Origin;
  vtkm::Id3 UniformDims;
  vtkm::Id3 MaxCellIds;

  vtkm::cont::ArrayHandle<vtkm::Id> CellCount;
  vtkm::cont::ArrayHandle<vtkm::Id> CellIds;
  vtkm::cont::ArrayHandle<vtkm::Id> CellStartIdx;

  struct MakeExecObject;
};

}
} // vtkm::cont

#endif // vtk_m_cont_CellLocatorUniformBins_h
