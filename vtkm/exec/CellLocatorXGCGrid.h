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
#include <vtkm/cont/CellSetExtrude.h>
#include <vtkm/cont/CellSetSingleType.h>

#include <vtkm/exec/CellInside.h>
#include <vtkm/exec/CellLocatorMultiplexer.h>
#include <vtkm/exec/CellLocatorTwoLevel.h>
#include <vtkm/exec/ConnectivityExtrude.h>
#include <vtkm/exec/ParametricCoordinates.h>

namespace vtkm
{

namespace exec
{

using FloatVec3 = vtkm::Vec3f;
using ConnSingleType =
  vtkm::cont::CellSetSingleType<>::ExecConnectivityType<vtkm::TopologyElementTagCell,
                                                        vtkm::TopologyElementTagPoint>;

using CoordsPortalType = vtkm::internal::ArrayPortalXGCCoordinates<
  vtkm::internal::ArrayPortalBasicRead<vtkm::FloatDefault>>;

using ConnExtrudeType = vtkm::exec::ConnectivityExtrude;

class VTKM_ALWAYS_EXPORT CellLocatorXGCGrid
{
public:
  VTKM_CONT
  CellLocatorXGCGrid(const vtkm::exec::ConnectivityExtrude& conn,
                     const CoordsPortalType& coords,
                     const vtkm::exec::CellLocatorTwoLevel<ConnSingleType>& locator,
                     const vtkm::Id& numPlanes,
                     const vtkm::Id& cellsPerPlane)
    : CellsPerPlane(cellsPerPlane)
    , Connectivity(conn)
    , Coords(coords)
    , LocatorMux(locator)
    , Periodic(true)
    , ThetaSpacing(static_cast<vtkm::FloatDefault>(vtkm::TwoPi() / numPlanes))
    , UseCylindrical(false)
  {
    this->LastPlaneTheta = this->ThetaSpacing * static_cast<vtkm::FloatDefault>(numPlanes - 1);
  }

  VTKM_CONT
  CellLocatorXGCGrid(const vtkm::exec::ConnectivityExtrude& conn,
                     const CoordsPortalType& coords,
                     const vtkm::exec::CellLocatorTwoLevel<ConnExtrudeType>& locator,
                     const vtkm::Id& numPlanes,
                     const vtkm::Id& cellsPerPlane,
                     bool periodic)
    : CellsPerPlane(cellsPerPlane)
    , Connectivity(conn)
    , Coords(coords)
    , LocatorMux(locator)
    , Periodic(periodic)
    , ThetaSpacing(static_cast<vtkm::FloatDefault>(vtkm::TwoPi() / numPlanes))
    , UseCylindrical(true)
  {
    this->LastPlaneTheta = this->ThetaSpacing * static_cast<vtkm::FloatDefault>(numPlanes - 1);

    //Non-periodic case not supported.
    //VTKM_ASSERT(this->Connectivity.GetPeriodic());
  }

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                           vtkm::Id& cellId,
                           vtkm::Vec3f& parametric) const
  {
    if (this->UseCylindrical)
      return this->FindCellCylindrical(point, cellId, parametric);
    else
      return this->FindCellCartesian(point, cellId, parametric);
  }

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorXGCGrid* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorXGCGrid* operator->() const { return this; }

private:
  VTKM_EXEC
  vtkm::ErrorCode FindCellCylindrical(const vtkm::Vec3f& point,
                                      vtkm::Id& cellId,
                                      vtkm::Vec3f& parametric) const
  {
    using LocType = vtkm::exec::CellLocatorTwoLevel<ConnExtrudeType>;
    auto loc = this->LocatorMux.Locators.Get<LocType>();

    //Case where point is *INSIDE* the wrap around cell.
    if (this->Periodic && point[1] > this->LastPlaneTheta)
    {
      //Shift the point BACK to the cell that isn't a wrap around.
      vtkm::Vec3f shiftPt(point[0], point[1] - this->ThetaSpacing, point[2]);
      //std::cout<<"****** SHIFT: "<<point<<" ---> "<<shiftPt<<std::endl;
      auto status = loc.FindCell(shiftPt, cellId, parametric);

      //Found the cell.  Now, the *RIGHT* cell will be offset by numCellsPerPlane
      if (status == vtkm::ErrorCode::Success)
        cellId += this->CellsPerPlane;
      return status;
    }
    else
      return loc.FindCell(point, cellId, parametric);
  }

  VTKM_EXEC
  vtkm::ErrorCode FindCellCartesian(const vtkm::Vec3f& point,
                                    vtkm::Id& cellId,
                                    vtkm::Vec3f& parametric) const
  {
    using LocType = vtkm::exec::CellLocatorTwoLevel<ConnSingleType>;
    LocType locator = this->LocatorMux.Locators.Get<LocType>();

    //Convert to cylindrical
    vtkm::FloatDefault x = point[0];
    vtkm::FloatDefault y = point[1];
    vtkm::FloatDefault z = point[2];

    vtkm::FloatDefault r = vtkm::Sqrt(x * x + y * y);
    vtkm::FloatDefault theta = vtkm::ATan2(y, x);
    //std::cout<<"  FindCellCart: "<<point<<" R= "<<r<<" Z= "<<z<<" theta= "<<theta<<" "<<(theta*57.2958)<<std::endl;

    if (theta < 0)
      theta += static_cast<vtkm::FloatDefault>(vtkm::TwoPi());

    //Point in cylindrical R,Z space for plane locator.
    vtkm::Vec3f cylPtRZ(r, z, 0);

    vtkm::Id cid = -1;
    auto res = locator.FindCell(cylPtRZ, cid, parametric);
    if (res != vtkm::ErrorCode::Success)
    {
      /*
      std::cout << "*******************************************************" << std::endl;
      std::cout << "  Pt NOT in 2D plane. " << std::endl;
      std::cout << "   cylPtRZ= " << cylPtRZ << std::endl;
      std::cout << "*******************************************************" << std::endl;
      */
      return res;
    }

    vtkm::Id planeIdx = static_cast<vtkm::Id>(vtkm::Floor(theta / this->ThetaSpacing));
    if (planeIdx > 0)
      cid += (planeIdx * this->CellsPerPlane);
    auto indices = this->Connectivity.GetIndices(cid);
    auto pts = vtkm::make_VecFromPortalPermute(&indices, this->Coords);

    vtkm::Vec3f pc;
    bool inside;
    VTKM_RETURN_ON_ERROR(this->PointInsideCell(point, pts, pc, inside));
    if (inside)
    {
      VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
        pts, point, vtkm::CellShapeTagWedge{}, parametric));
      inside = vtkm::exec::CellInside(parametric, vtkm::CellShapeTagWedge{});

      cellId = cid;
      return vtkm::ErrorCode::Success;
    }
    return vtkm::ErrorCode::CellNotFound;
  }

  template <typename PointsVecType>
  VTKM_EXEC vtkm::Bounds ComputeCellBounds(const PointsVecType& points) const
  {
    using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
    CoordsType minp = points[0], maxp = points[0];

    for (vtkm::IdComponent i = 1; i < 6; i++)
    {
      minp = vtkm::Min(minp, points[i]);
      maxp = vtkm::Max(maxp, points[i]);
    }

    return { FloatVec3(minp), FloatVec3(maxp) };
  }

  template <typename CoordsType>
  VTKM_EXEC vtkm::ErrorCode PointInsideCell(FloatVec3 point,
                                            CoordsType cellPoints,
                                            FloatVec3& parametricCoordinates,
                                            bool& inside) const
  {
    vtkm::Bounds bounds = this->ComputeCellBounds(cellPoints);

    inside = false;
    if (bounds.Contains(point))
    {
      VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
        cellPoints, point, vtkm::CellShapeTagWedge{}, parametricCoordinates));
      //      std::cout << "   PARAMETRIC: " << parametricCoordinates << std::endl;
      inside = vtkm::exec::CellInside(parametricCoordinates, vtkm::CellShapeTagWedge{});
    }

    // Return success error code even point is not inside this cell
    return vtkm::ErrorCode::Success;
  }

  vtkm::Id CellsPerPlane;
  vtkm::exec::ConnectivityExtrude Connectivity;
  CoordsPortalType Coords;
  vtkm::FloatDefault LastPlaneTheta;
  vtkm::exec::CellLocatorMultiplexer<vtkm::exec::CellLocatorTwoLevel<ConnSingleType>,
                                     vtkm::exec::CellLocatorTwoLevel<ConnExtrudeType>>
    LocatorMux;
  bool Periodic;
  vtkm::FloatDefault ThetaSpacing;
  bool UseCylindrical;
};
}
}

#endif //vtkm_exec_celllocatorxgcgrid_h
