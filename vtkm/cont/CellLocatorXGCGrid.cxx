//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleStride.h>
#include <vtkm/cont/ArrayHandleXGCCoordinates.h>
#include <vtkm/cont/CellLocatorTwoLevel.h>
#include <vtkm/cont/CellLocatorXGCGrid.h>
#include <vtkm/cont/CellSetExtrude.h>
#include <vtkm/cont/CellSetSingleType.h>
#include <vtkm/exec/ConnectivityExtrude.h>

namespace vtkm
{
namespace cont
{

using XGCType = vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault>;
using ExtrudedCell = vtkm::cont::CellSetExtrude;

void CellLocatorXGCGrid::Build()
{
  vtkm::cont::CoordinateSystem coords = this->GetCoordinates();
  vtkm::cont::DynamicCellSet cellSet = this->GetCellSet();


  if (!coords.GetData().IsType<XGCType>())
    throw vtkm::cont::ErrorBadType("Coordinates are not XGC type.");
  if (!cellSet.IsSameType(ExtrudedCell()))
    throw vtkm::cont::ErrorBadType("Cells are not Extruded type.");

  //Create the coords/cellset for the plane.

  XGCType xgcCoordSys;
  coords.GetData().AsArrayHandle(xgcCoordSys);
  this->IsCylindrical = xgcCoordSys.GetUseCylindrical();

  auto xgcCellSet = cellSet.Cast<vtkm::cont::CellSetExtrude>();
  this->NumPlanes = static_cast<vtkm::Id>(xgcCellSet.GetNumberOfPlanes());
  this->CellsPerPlane = static_cast<vtkm::Id>(xgcCellSet.GetNumberOfCellsPerPlane());

  //If cylindrical, create a TwoLevelLocator on the 3D points.
  if (this->IsCylindrical)
  {
    this->CellLocator.SetCoordinates(coords);

    //The cell locator needs to work on a non-periodic connectivity.
    //Create a non-periodic cellset and set the periodic flag.
    this->IsPeriodic = xgcCellSet.GetIsPeriodic();
    if (this->IsPeriodic)
    {
      vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault> xgcCoords;
      this->GetCoordinates().GetData().AsArrayHandle(xgcCoords);
      auto nonPeriodicCellSet = vtkm::cont::make_CellSetExtrude(
        xgcCellSet.GetConnectivityArray(), xgcCoords, xgcCellSet.GetNextNodeArray(), false);
      this->CellLocator.SetCellSet(nonPeriodicCellSet);
    }
    else
      this->CellLocator.SetCellSet(cellSet);
  }
  else
  {
    //For cartesian, create a locator on the RZ plane points.
    auto rzPts = xgcCoordSys.GetArray();

    //Create the RZ 2D plane.
    vtkm::Id num = rzPts.GetNumberOfValues();
    ArrayHandleStride<vtkm::FloatDefault> RCoords(rzPts, num / 2, 2, 0);
    ArrayHandleStride<vtkm::FloatDefault> ZCoords(rzPts, num / 2, 2, 1);
    auto ZeroArray = vtkm::cont::make_ArrayHandleConstant(vtkm::FloatDefault(0), num / 2);

    vtkm::cont::ArrayHandle<vtkm::Vec3f> planePts;
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleCompositeVector(RCoords, ZCoords, ZeroArray),
                          planePts);

    vtkm::cont::ArrayHandle<vtkm::Id> conn;
    vtkm::cont::ArrayCopy(
      vtkm::cont::make_ArrayHandleCast<vtkm::Id>(xgcCellSet.GetConnectivityArray()), conn);
    vtkm::cont::CellSetSingleType<> planeCells;
    vtkm::Id ptsPerPlane = static_cast<vtkm::Id>(xgcCellSet.GetNumberOfPointsPerPlane());
    planeCells.Fill(ptsPerPlane, vtkm::CELL_SHAPE_TRIANGLE, 3, conn);
    this->CellsPerPlane = planeCells.GetNumberOfCells();

    //Use RZ plane for the locator
    this->CellLocator.SetCellSet(planeCells);
    this->CellLocator.SetCoordinates(vtkm::cont::CoordinateSystem("coords", planePts));
  }
}

vtkm::exec::CellLocatorXGCGrid CellLocatorXGCGrid::PrepareForExecution(
  vtkm::cont::DeviceAdapterId device,
  vtkm::cont::Token& token) const
{
  this->CellLocator.Update();
  this->Update();

  auto locMux = this->CellLocator.PrepareForExecution(device, token);

  vtkm::cont::DynamicCellSet cellSet = this->GetCellSet();
  auto xgcCellSet = cellSet.Cast<vtkm::cont::CellSetExtrude>();
  vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault> xgcCoords;
  this->GetCoordinates().GetData().AsArrayHandle(xgcCoords);

  auto cellSetExec = xgcCellSet.PrepareForInput(
    device, vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, token);

  auto coordsExec = xgcCoords.PrepareForInput(device, token);

  if (this->IsCylindrical)
  {
    using ExtrudeConnectivityType = vtkm::exec::ConnectivityExtrude;

    auto locator = locMux.Locators.Get<vtkm::exec::CellLocatorTwoLevel<ExtrudeConnectivityType>>();
    return vtkm::exec::CellLocatorXGCGrid(
      cellSetExec, coordsExec, locator, this->NumPlanes, this->CellsPerPlane, this->IsPeriodic);
  }
  else
  {
    using SingleConnectivityType =
      vtkm::cont::CellSetSingleType<>::ExecConnectivityType<vtkm::TopologyElementTagCell,
                                                            vtkm::TopologyElementTagPoint>;
    auto locator = locMux.Locators.Get<vtkm::exec::CellLocatorTwoLevel<SingleConnectivityType>>();
    return vtkm::exec::CellLocatorXGCGrid(
      cellSetExec, coordsExec, locator, this->NumPlanes, this->CellsPerPlane);
  }
}

} //namespace cont
} //namespace vtkm
