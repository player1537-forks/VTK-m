//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_exec_CellLocatorTwoLevel_h
#define vtk_m_exec_CellLocatorTwoLevel_h

#include <vtkm/exec/CellInside.h>
#include <vtkm/exec/ParametricCoordinates.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>

#include <vtkm/Math.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/TopologyElementTag.h>
#include <vtkm/Types.h>
#include <vtkm/VecFromPortalPermute.h>
#include <vtkm/VecTraits.h>

#include <vtkm/exec/LocatorTimer.h>

#include <iomanip>

extern vtkm::exec::LocatorTimer locTimer;

namespace vtkm
{
namespace internal
{
namespace cl_uniform_bins
{

using DimensionType = vtkm::Int16;
using DimVec3 = vtkm::Vec<DimensionType, 3>;
using FloatVec3 = vtkm::Vec3f;

struct Grid
{
  DimVec3 Dimensions;
  // Bug in CUDA 9.2 where having this gap for alignment was for some reason setting garbage
  // in a union with other cell locators (or perhaps not properly copying data). This appears
  // to be fixed by CUDA 10.2.
  DimensionType Padding;
  FloatVec3 Origin;
  FloatVec3 BinSize;
};

struct Bounds
{
  FloatVec3 Min;
  FloatVec3 Max;
};

VTKM_EXEC inline vtkm::Id ComputeFlatIndex(const DimVec3& idx, const DimVec3 dim)
{
  return idx[0] + (dim[0] * (idx[1] + (dim[1] * idx[2])));
}

VTKM_EXEC inline Grid ComputeLeafGrid(const DimVec3& idx, const DimVec3& dim, const Grid& l1Grid)
{
  return { dim,
           0,
           l1Grid.Origin + (static_cast<FloatVec3>(idx) * l1Grid.BinSize),
           l1Grid.BinSize / static_cast<FloatVec3>(dim) };
}

template <typename PointsVecType>
VTKM_EXEC inline Bounds ComputeCellBounds(const PointsVecType& points)
{
  using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
  auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

  CoordsType minp = points[0], maxp = points[0];
  for (vtkm::IdComponent i = 1; i < numPoints; ++i)
  {
    minp = vtkm::Min(minp, points[i]);
    maxp = vtkm::Max(maxp, points[i]);
  }

  return { FloatVec3(minp), FloatVec3(maxp) };
}
}
}
} // vtkm::internal::cl_uniform_bins

namespace vtkm
{
namespace exec
{

//--------------------------------------------------------------------
template <typename CellStructureType>
class VTKM_ALWAYS_EXPORT CellLocatorTwoLevel
{
private:
  using DimVec3 = vtkm::internal::cl_uniform_bins::DimVec3;
  using FloatVec3 = vtkm::internal::cl_uniform_bins::FloatVec3;

  template <typename T>
  using ReadPortal = typename vtkm::cont::ArrayHandle<T>::ReadPortalType;

  using CoordsPortalType =
    typename vtkm::cont::CoordinateSystem::MultiplexerArrayType::ReadPortalType;

  VTKM_EXEC
  static vtkm::FloatDefault Area(const vtkm::FloatDefault& x1, const vtkm::FloatDefault& y1,
                                 const vtkm::FloatDefault& x2, const vtkm::FloatDefault& y2,
                                 const vtkm::FloatDefault& x3, const vtkm::FloatDefault& y3)
  {
    return vtkm::Abs((x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2.0);
    //return vtkm::Abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
  }

  template <typename CoordsType>
  VTKM_EXEC
  static vtkm::Vec3f PCoords(const vtkm::Vec3f& pt,
                             const CoordsType& cellPoints)
  {
    vtkm::Vec3f Norm(0,0,1), pc(0,0,0);

    for (int i = 0; i < 2; i++)
    {
      const auto& p0 = cellPoints[0];
      const auto& p1 = cellPoints[i+1];
      const auto& p2 = cellPoints[2-i];
      const auto planeNorm = vtkm::Cross(Norm, (p2-p0));
      pc[i] = vtkm::Dot((pt-p0), planeNorm) / vtkm::Dot((p1-p0), planeNorm);
    }

    return pc;
  }

  template <typename CoordsType>
  VTKM_EXEC
  static bool PtInCell(const vtkm::Vec3f& point,
                       const CoordsType& cellPoints,
                       vtkm::Vec3f& pcoords)
  {
    auto bounds = vtkm::internal::cl_uniform_bins::ComputeCellBounds(cellPoints);

    auto xMin = vtkm::Min(cellPoints[0][0], vtkm::Min(cellPoints[1][0], cellPoints[2][0]));
    auto xMax = vtkm::Max(cellPoints[0][0], vtkm::Max(cellPoints[1][0], cellPoints[2][0]));
    auto yMin = vtkm::Min(cellPoints[0][1], vtkm::Min(cellPoints[1][1], cellPoints[2][1]));
    auto yMax = vtkm::Max(cellPoints[0][1], vtkm::Max(cellPoints[1][1], cellPoints[2][1]));
    if (point[0] < xMin || point[0] > xMax || point[1] < yMin || point[1] > yMax)
      return false;

    pcoords = PCoords(point, cellPoints);
    return (pcoords[0] >=0 && pcoords[0] <= 1 && pcoords[1] >= 0 && pcoords[1] <= 1 && pcoords[0]+pcoords[1] <= 1);
  }


  // TODO: This function may return false positives for non 3D cells as the
  // tests are done on the projection of the point on the cell. Extra checks
  // should be added to test if the point actually falls on the cell.
  template <typename CellShapeTag, typename CoordsType>
  VTKM_EXEC static vtkm::ErrorCode PointInsideCell(FloatVec3 point,
                                                   CellShapeTag cellShape,
                                                   CoordsType cellPoints,
                                                   FloatVec3& pcoords,
                                                   bool& inside)
  {
#if 1
    inside = PtInCell(point, cellPoints, pcoords);

#if 0
    auto xMin = vtkm::Min(cellPoints[0][0], vtkm::Min(cellPoints[1][0], cellPoints[2][0]));
    auto xMax = vtkm::Max(cellPoints[0][0], vtkm::Max(cellPoints[1][0], cellPoints[2][0]));
    auto yMin = vtkm::Min(cellPoints[0][1], vtkm::Min(cellPoints[1][1], cellPoints[2][1]));
    auto yMax = vtkm::Max(cellPoints[0][1], vtkm::Max(cellPoints[1][1], cellPoints[2][1]));
    if (point[0] < xMin || point[0] > xMax || point[1] < yMin || point[1] > yMax)
      inside = false;
    else
    {
      pcoords = PCoords(point, cellPoints);
      inside = (pcoords[0] >=0 && pcoords[0] <= 1 && pcoords[1] >= 0 && pcoords[1] <= 1 && pcoords[0]+pcoords[1] <= 1);
    }

    /*
    const auto& p1 = cellPoints[0];
    const auto& p2 = cellPoints[1];
    const auto& p3 = cellPoints[2];

    auto A = Area(    p1[0],    p1[1],    p2[0],    p2[1],    p3[0],    p3[1]);
    auto A1 = Area(point[0], point[1],    p2[0],    p2[1],    p3[0],    p3[1]);
    auto A2 = Area(   p1[0],    p1[1], point[0], point[1],    p3[0],    p3[1]);
    auto A3 = Area(   p1[0],    p1[1],    p2[0],    p2[1], point[0], point[1]);
    auto dA = A-(A1+A2+A3));
    inside = vtkm::Abs(dA) <= 1e-6;
    */
#endif
#else

    auto bounds = vtkm::internal::cl_uniform_bins::ComputeCellBounds(cellPoints);
    if (point[0] >= bounds.Min[0] && point[0] <= bounds.Max[0] && point[1] >= bounds.Min[1] &&
        point[1] <= bounds.Max[1] && point[2] >= bounds.Min[2] && point[2] <= bounds.Max[2])
    {
      VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
        cellPoints, point, cellShape, pcoords));

      inside = vtkm::exec::CellInside(pcoords, cellShape);
    }
    else
    {
      inside = false;
    }

    vtkm::Vec3f pcoords2;
    bool inside2 = PtInCell(point, cellPoints, pcoords2);
    if (inside != inside2)
    {
      std::cout<<"Wrong answer."<<std::endl;
    }

    if (inside && vtkm::Magnitude(pcoords-pcoords2) > 1e-14)
    {
      std::cout<<"***** "<<std::endl;
      std::cout<<std::setprecision(20);
      std::cout<<pcoords<<" "<<pcoords2<<" :: "<<vtkm::Magnitude(pcoords-pcoords2)<<std::endl;
      inside2 = PtInCell(point, cellPoints, pcoords2);
      inside2 = PtInCell(point, cellPoints, pcoords2);
      inside2 = PtInCell(point, cellPoints, pcoords2);
    }
#endif

    /*
   if (!inside && omp_get_thread_num() == 0)
   {
     auto p1 = cellPoints[0];
     auto p2 = cellPoints[1];
     auto p3 = cellPoints[2];

//     p1[0] = -2; p1[1] =  3;
//     p2[0] = -3; p2[1] = -1;
//     p3[0] =  3; p3[1] = -2;

     auto A = Area(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
     auto A1 = Area(point[0], point[1], p2[0], p2[1], p3[0], p3[1]);
     auto A2 = Area(p1[0], p1[1], point[0], point[1], p3[0], p3[1]);
     auto A3 = Area(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
     auto dA = vtkm::Abs(A-(A1+A2+A3));
     if (dA > 0)
     {
       std::cout<<std::setprecision(20)<<" ********ERROR:  "<<dA<<" "<<A<<" "<<A1<<" "<<A2<<" "<<A3<<std::endl;
       std::cout<<p1<<" "<<p2<<" "<<p3<<std::endl;
     }

     auto pc = PCoords(point, cellPoints);

     auto b1 = A1/A;
     auto b2 = A2/A;
     auto b3 = A3/A;
     std::cout<<" Param: "<<parametricCoordinates<<" me: "<<pc<<std::endl;
   }
    */

    return vtkm::ErrorCode::Success;
  }

public:
  template <typename CellSetType>
  VTKM_CONT CellLocatorTwoLevel(const vtkm::internal::cl_uniform_bins::Grid& topLevelGrid,
                                const vtkm::cont::ArrayHandle<DimVec3>& leafDimensions,
                                const vtkm::cont::ArrayHandle<vtkm::Id>& leafStartIndex,
                                const vtkm::cont::ArrayHandle<vtkm::Id>& cellStartIndex,
                                const vtkm::cont::ArrayHandle<vtkm::Id>& cellCount,
                                const vtkm::cont::ArrayHandle<vtkm::Id>& cellIds,
                                const CellSetType& cellSet,
                                const vtkm::cont::CoordinateSystem& coords,
                                vtkm::cont::DeviceAdapterId device,
                                vtkm::cont::Token& token)
    : TopLevel(topLevelGrid)
    , LeafDimensions(leafDimensions.PrepareForInput(device, token))
    , LeafStartIndex(leafStartIndex.PrepareForInput(device, token))
    , CellStartIndex(cellStartIndex.PrepareForInput(device, token))
    , CellCount(cellCount.PrepareForInput(device, token))
    , CellIds(cellIds.PrepareForInput(device, token))
    , CellSet(cellSet.PrepareForInput(device,
                                      vtkm::TopologyElementTagCell{},
                                      vtkm::TopologyElementTagPoint{},
                                      token))
    , Coords(coords.GetDataAsMultiplexer().PrepareForInput(device, token))
  {
  }

  struct LastCell
  {
    vtkm::Id CellId = -1;
    vtkm::Id LeafIdx = -1;
  };

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const FloatVec3& point, vtkm::Id& cellId, FloatVec3& parametric) const
  {
    LastCell lastCell;
    return this->FindCellImpl(point, cellId, parametric, lastCell);
  }

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const FloatVec3& point,
                           vtkm::Id& cellId,
                           FloatVec3& parametric,
                           LastCell& lastCell) const
  {
    locTimer.Start();
    vtkm::Vec3f pc;
    //See if point is inside the last cell.
    if ((lastCell.CellId >= 0) && (lastCell.CellId < this->CellSet.GetNumberOfElements()) &&
        this->PointInCell(point, lastCell.CellId, pc) == vtkm::ErrorCode::Success)
    {
      parametric = pc;
      cellId = lastCell.CellId;
      locTimer.Update(LocatorTimer::FindCell0_IDX);

      return vtkm::ErrorCode::Success;
    }

    //See if it's in the last leaf.
    if ((lastCell.LeafIdx >= 0) && (lastCell.LeafIdx < this->CellCount.GetNumberOfValues()) &&
        this->PointInLeaf(point, lastCell.LeafIdx, cellId, pc) == vtkm::ErrorCode::Success)
    {
      parametric = pc;
      lastCell.CellId = cellId;
      locTimer.Update(LocatorTimer::FindCell1_IDX);
      return vtkm::ErrorCode::Success;
    }
    locTimer.Update(LocatorTimer::FindCell2_IDX);


    //Call the full point search.
    return this->FindCellImpl(point, cellId, parametric, lastCell);
  }

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorTwoLevel* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorTwoLevel* operator->() const { return this; }

private:
  VTKM_EXEC
  vtkm::ErrorCode PointInCell(const vtkm::Vec3f& point,
                              const vtkm::Id& cid,
                              vtkm::Vec3f& parametric) const
  {
    locTimer.Start();
    auto indices = this->CellSet.GetIndices(cid);
    auto pts = vtkm::make_VecFromPortalPermute(&indices, this->Coords);
    vtkm::Vec3f pc;
    bool inside;
    auto status = PointInsideCell(point, this->CellSet.GetCellShape(cid), pts, pc, inside);
    if (status == vtkm::ErrorCode::Success && inside)
    {
      parametric = pc;
      locTimer.Update(LocatorTimer::PointInCell_IDX);

      return vtkm::ErrorCode::Success;
    }
    locTimer.Update(LocatorTimer::PointInCell_IDX);
    return vtkm::ErrorCode::CellNotFound;
  }

  VTKM_EXEC
  vtkm::ErrorCode PointInLeaf(const FloatVec3& point,
                              const vtkm::Id& leafIdx,
                              vtkm::Id& cellId,
                              FloatVec3& parametric) const
  {
    vtkm::Id start = this->CellStartIndex.Get(leafIdx);
    vtkm::Id end = start + this->CellCount.Get(leafIdx);

    for (vtkm::Id i = start; i < end; ++i)
    {
      vtkm::Vec3f pc;

      vtkm::Id cid = this->CellIds.Get(i);
      if (this->PointInCell(point, cid, pc) == vtkm::ErrorCode::Success)
      {
        cellId = cid;
        parametric = pc;
        return vtkm::ErrorCode::Success;
      }
    }

    return vtkm::ErrorCode::CellNotFound;
  }


  VTKM_EXEC
  vtkm::ErrorCode FindCellImpl(const FloatVec3& point,
                               vtkm::Id& cellId,
                               FloatVec3& parametric,
                               LastCell& lastCell) const
  {
    using namespace vtkm::internal::cl_uniform_bins;
    locTimer.Start();

    cellId = -1;
    lastCell.CellId = -1;
    lastCell.LeafIdx = -1;

    DimVec3 binId3 = static_cast<DimVec3>((point - this->TopLevel.Origin) / this->TopLevel.BinSize);
    if (binId3[0] >= 0 && binId3[0] < this->TopLevel.Dimensions[0] && binId3[1] >= 0 &&
        binId3[1] < this->TopLevel.Dimensions[1] && binId3[2] >= 0 &&
        binId3[2] < this->TopLevel.Dimensions[2])
    {
      vtkm::Id binId = ComputeFlatIndex(binId3, this->TopLevel.Dimensions);

      auto ldim = this->LeafDimensions.Get(binId);
      if (!ldim[0] || !ldim[1] || !ldim[2])
      {
        locTimer.Update(LocatorTimer::FindCellImpl0_IDX);
        return vtkm::ErrorCode::CellNotFound;
      }

      auto leafGrid = ComputeLeafGrid(binId3, ldim, this->TopLevel);

      DimVec3 leafId3 = static_cast<DimVec3>((point - leafGrid.Origin) / leafGrid.BinSize);
      // precision issues may cause leafId3 to be out of range so clamp it
      leafId3 = vtkm::Max(DimVec3(0), vtkm::Min(ldim - DimVec3(1), leafId3));

      vtkm::Id leafStart = this->LeafStartIndex.Get(binId);
      vtkm::Id leafIdx = leafStart + ComputeFlatIndex(leafId3, leafGrid.Dimensions);

      if (this->PointInLeaf(point, leafIdx, cellId, parametric) == vtkm::ErrorCode::Success)
      {
        lastCell.CellId = cellId;
        lastCell.LeafIdx = leafIdx;
        locTimer.Update(LocatorTimer::FindCellImpl1_IDX);
        return vtkm::ErrorCode::Success;
      }
    }
    locTimer.Update(LocatorTimer::FindCellImpl2_IDX);
    return vtkm::ErrorCode::CellNotFound;
  }

  vtkm::internal::cl_uniform_bins::Grid TopLevel;

  ReadPortal<DimVec3> LeafDimensions;
  ReadPortal<vtkm::Id> LeafStartIndex;

  ReadPortal<vtkm::Id> CellStartIndex;
  ReadPortal<vtkm::Id> CellCount;
  ReadPortal<vtkm::Id> CellIds;

  CellStructureType CellSet;
  CoordsPortalType Coords;

//  vtkm::FloatDefault FindCellTime = 0.0;
//  vtkm::FloatDefault TriIntersectType = 0.0;
};
}
} // vtkm::exec

#endif //vtk_m_exec_CellLocatorTwoLevel_h
