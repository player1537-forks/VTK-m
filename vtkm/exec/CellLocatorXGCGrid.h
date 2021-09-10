//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtkm_exec_celllocatorxgcgrid_h
#define vtkm_exec_celllocatorxgcgrid_h

#include <vtkm/Bounds.h>
#include <vtkm/Math.h>
#include <vtkm/TopologyElementTag.h>
#include <vtkm/Types.h>
#include <vtkm/VecFromPortalPermute.h>

#include <vtkm/cont/ArrayHandleXGCCoordinates.h>
#include <vtkm/cont/CellSetSingleType.h>
#include <vtkm/cont/CellSetExtrude.h>

#include <vtkm/exec/CellInside.h>
#include <vtkm/exec/ConnectivityExtrude.h>
#include <vtkm/exec/ParametricCoordinates.h>
#include <vtkm/exec/CellLocatorTwoLevel.h>

namespace vtkm
{

namespace exec
{

using CellLocatorType = vtkm::cont::CellSetSingleType<>::ExecConnectivityType<vtkm::TopologyElementTagCell, vtkm::TopologyElementTagPoint>;

/*
using CoordsPortalType = vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault>::ReadPortalType;
//using ExtrudeCellSetPortalType = vtkm::cont::CellSetExtrude::ReadPortalType;

template <typename T>
using ReadPortal = typename vtkm::cont::ArrayHandle<T>::ReadPortalType;
*/

using CoordsPortalType = vtkm::internal::ArrayPortalXGCCoordinates<vtkm::internal::ArrayPortalBasicRead<vtkm::FloatDefault>>;

class VTKM_ALWAYS_EXPORT CellLocatorXGCGrid
{
public:
  VTKM_CONT
CellLocatorXGCGrid(const vtkm::exec::ConnectivityExtrude& conn,
                   const CoordsPortalType& coords,
                   const vtkm::exec::CellLocatorTwoLevel<CellLocatorType>& planeLocator,
                   const vtkm::Id& numPlanes,
                   const vtkm::Id& cellsPerPlane)
: CellsPerPlane(cellsPerPlane)
, Coords(coords)
, Connectivity(conn)
, NumPlanes(numPlanes)
, PlaneLocator(planeLocator)
, ThetaSpacing(vtkm::TwoPi() / numPlanes)
{
}

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                           vtkm::Id& cellId,
                           vtkm::Vec3f& parametric) const
{
    vtkm::FloatDefault x = point[0];
    vtkm::FloatDefault y = point[1];
    vtkm::FloatDefault z = point[2];

    vtkm::FloatDefault r = vtkm::Sqrt(x*x + y*y);
    vtkm::FloatDefault theta = vtkm::ATan2(y, x);
    if (theta < 0)
      theta += vtkm::TwoPi();

    vtkm::Vec3f cylPt(r, z, theta);

    std::cout<<"FindCell: "<<point<<" --> "<<cylPt<<std::endl;
    cylPt[2] = 0;
    auto res = this->PlaneLocator.FindCell(cylPt, cellId, parametric);
    if (res == vtkm::ErrorCode::Success)
    {
      vtkm::Id planeIdx = vtkm::Floor(theta / this->ThetaSpacing);

      if (planeIdx > 0)
        cellId += (planeIdx*this->CellsPerPlane);

      auto indices = this->Connectivity.GetIndices(cellId);
      auto pts = vtkm::make_VecFromPortalPermute(&indices, this->Coords);
      for (int i = 0; i < 6; i++)
        std::cout<<"Pt_"<<i<<" "<<indices[i]<<std::endl;

      /*
      vtkm::Id cid = this->CellIds.Get(i);
      auto indices = this->CellSet.GetIndices(cid);
      */

      std::cout<<" *** cellId= "<<cellId<<" p: "<<parametric<<std::endl;
      std::cout<<"   planeIdx= "<<planeIdx<<" cellId= "<<cellId<<std::endl;
      std::cout<<"   parametric= "<<parametric<<std::endl;
    }

    return res;
}

private:
  vtkm::exec::ConnectivityExtrude Connectivity;
  CoordsPortalType Coords;
  vtkm::exec::CellLocatorTwoLevel<CellLocatorType> PlaneLocator;
  bool IsCylindrical;
  vtkm::Id NumPlanes;
  vtkm::Id CellsPerPlane;
  vtkm::FloatDefault ThetaSpacing;

#if 0
    CellLocatorXGCGrid(const vtkm::cont::DynamicCellSet& cellSet,
                       const vtkm::exec::CellLocatorTwoLevel<CellLocatorType>& planeLocator,
                       bool isCyl,
                       vtkm::Id cellsPerPlane,
                       vtkm::Id numPlanes)
    : PlaneLocator(planeLocator)
    , IsCylindrical(isCyl)
    , NumPlanes(numPlanes)
    , CellsPerPlane(cellsPerPlane)
    , ThetaSpacing(vtkm::TwoPi() / numPlanes)
    //, CellSet(cellSet)
    //.PrepareForInput(device,
    //vtkm::TopologyElementTagCell{},
    //vtkm::TopologyElementTagPoint{},
    //token))
  {
    std::cout<<"THETA SPACING: "<<this->ThetaSpacing<<std::endl;
    //auto cellSetExec = cellSet.PrepareForInput(device, vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, token);
  }
/*
  VTKM_CONT
    CellLocatorXGCGrid(const vtkm::cont::CellSetExtrude& cellSet,
                       const vtkm::exec::CellLocatorTwoLevel<CellLocatorType>& planeLocator,
                       bool isCyl,
                       vtkm::Id cellsPerPlane,
                       vtkm::Id numPlanes,
                       vtkm::cont::DeviceAdapterId device,
                       vtkm::cont::Token& token)
//                       const vtkm::exec::ConnectivityExtrude& bumBum)
    : PlaneLocator(planeLocator)
    , IsCylindrical(isCyl)
    , NumPlanes(numPlanes)
    , CellsPerPlane(cellsPerPlane)
    , ThetaSpacing(vtkm::TwoPi() / numPlanes)
    , CellSet(cellSet)
    //.PrepareForInput(device,
    //vtkm::TopologyElementTagCell{},
    //vtkm::TopologyElementTagPoint{},
    //token))
  {
    std::cout<<"THETA SPACING: "<<this->ThetaSpacing<<std::endl;
    auto cellSetExec = cellSet.PrepareForInput(device, vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, token);
  }
*/

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                           vtkm::Id& cellId,
                           vtkm::Vec3f& parametric) const
  {
    vtkm::FloatDefault x = point[0];
    vtkm::FloatDefault y = point[1];
    vtkm::FloatDefault z = point[2];

    vtkm::FloatDefault r = vtkm::Sqrt(x*x + y*y);
    vtkm::FloatDefault theta = vtkm::ATan2(y, x);
    if (theta < 0)
      theta += vtkm::TwoPi();

    vtkm::Vec3f cylPt(r, z, theta);

    std::cout<<"FindCell: "<<point<<" --> "<<cylPt<<std::endl;
    cylPt[2] = 0;
    auto res = this->PlaneLocator.FindCell(cylPt, cellId, parametric);
    if (res == vtkm::ErrorCode::Success)
    {
      vtkm::Id planeIdx = vtkm::Floor(theta / this->ThetaSpacing);

      if (planeIdx > 0)
        cellId += (planeIdx*this->CellsPerPlane);

      /*
      vtkm::Id cid = this->CellIds.Get(i);
      auto indices = this->CellSet.GetIndices(cid);
      */

      std::cout<<" *** cellId= "<<cellId<<" p: "<<parametric<<std::endl;
      std::cout<<"   planeIdx= "<<planeIdx<<" cellId= "<<cellId<<std::endl;
      std::cout<<"   parametric= "<<parametric<<std::endl;
    }

    return res;
  }

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorXGCGrid* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorXGCGrid* operator->() const { return this; }

private:
  vtkm::cont::CellSetExtrude CellSet;
  vtkm::exec::CellLocatorTwoLevel<CellLocatorType> PlaneLocator;
  bool IsCylindrical;
  vtkm::Id NumPlanes;
  vtkm::Id CellsPerPlane;
  vtkm::FloatDefault ThetaSpacing;

  //CoordsPortalType Coords;
  //ExtrudeCellSetPortalType CellSet2;

  //ReadPortal<vtkm::Id> Connectivity;
#endif
};
}
}

#endif //vtkm_exec_celllocatorxgcgrid_h
