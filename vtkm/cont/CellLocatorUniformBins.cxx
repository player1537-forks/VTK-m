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
<<<<<<< HEAD
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/CellLocatorUniformBins.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
=======
#include <vtkm/cont/CellLocatorUniformBins.h>
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
#include <vtkm/cont/Invoker.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

namespace vtkm
{
  namespace cont
  {

<<<<<<< HEAD
  namespace
  {

  VTKM_EXEC
  inline vtkm::Id ComputeFlatIndex(const vtkm::Id3& idx, const vtkm::Id3& dims)
=======
  namespace detail
  {

  VTKM_EXEC
  static inline vtkm::Id ComputeFlatIndex(const vtkm::Id3& idx, const vtkm::Id3& dims)
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  {
    return idx[0] + (dims[0] * (idx[1] + (dims[1] * idx[2])));
  }

  VTKM_EXEC
<<<<<<< HEAD
  inline vtkm::Id3 GetFlatIndex(const vtkm::Vec3f& pt,
                                const vtkm::Vec3f& origin,
                                const vtkm::Vec3f& invSpacing,
                                const vtkm::Id3& maxCellIds)
=======
  static inline vtkm::Id3 GetFlatIndex(const vtkm::Vec3f& pt,
                                       const vtkm::Vec3f& origin,
                                       const vtkm::Vec3f& invSpacing,
                                       const vtkm::Id3& maxCellIds)
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  {
    auto temp = pt - origin;
    temp = temp * invSpacing;

    auto logicalIdx = vtkm::Min(vtkm::Id3(temp), maxCellIds);
    return logicalIdx;
  }

<<<<<<< HEAD
  template <typename PointsVecType>
  VTKM_EXEC inline void MinMaxIndicesForCellPoints(const PointsVecType& points,
                                                   const vtkm::Vec3f& origin,
                                                   const vtkm::Vec3f& invSpacing,
                                                   const vtkm::Id3& maxCellIds,
                                                   vtkm::Id3& minIdx,
                                                   vtkm::Id3& maxIdx)
  {
    auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

    vtkm::Bounds bounds;
    for (vtkm::IdComponent i = 0; i < numPoints; ++i)
      bounds.Include(points[i]);

    //Get 8 corners of bbox.
    minIdx = GetFlatIndex(bounds.MinCorner(), origin, invSpacing, maxCellIds);
    maxIdx = GetFlatIndex(bounds.MaxCorner(), origin, invSpacing, maxCellIds);
  }

=======
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  class CountCellBins : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cellset, FieldInPoint coords, FieldOutCell bincount);
    using ExecutionSignature = void(_2, _3);
    using InputDomain = _1;

    CountCellBins(const vtkm::Vec3f& origin,
                  const vtkm::Vec3f& invSpacing,
                  const vtkm::Id3& maxCellIds)
      : InvSpacing(invSpacing)
      , MaxCellIds(maxCellIds)
      , Origin(origin)
    {
    }

    template <typename PointsVecType>
    VTKM_EXEC void operator()(const PointsVecType& points, vtkm::Id& numBins) const
    {
<<<<<<< HEAD
      vtkm::Id3 idx000, idx111;
      MinMaxIndicesForCellPoints(
        points, this->Origin, this->InvSpacing, this->MaxCellIds, idx000, idx111);

      //Count the number of bins for this cell
      numBins =
        (idx111[0] - idx000[0] + 1) * (idx111[1] - idx000[1] + 1) * (idx111[2] - idx000[2] + 1);
    }

  private:
=======
      auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

      vtkm::Vec3f bbox[2] = { points[0], points[0] };
      for (vtkm::IdComponent i = 1; i < numPoints; ++i)
      {
        for (vtkm::IdComponent j = 0; j < 3; j++)
        {
          bbox[0][j] = vtkm::Min(bbox[0][j], static_cast<vtkm::FloatDefault>(points[i][j]));
          bbox[1][j] = vtkm::Max(bbox[1][j], static_cast<vtkm::FloatDefault>(points[i][j]));
        }
      }

      //Get 8 corners of bbox.
      auto idx000 = GetFlatIndex(bbox[0], this->Origin, this->InvSpacing, this->MaxCellIds);
      auto idx111 = GetFlatIndex(bbox[1], this->Origin, this->InvSpacing, this->MaxCellIds);

      //Count the number of bins
      numBins = 0;
      for (vtkm::Id i = idx000[0]; i <= idx111[0]; i++)
        for (vtkm::Id j = idx000[1]; j <= idx111[1]; j++)
          for (vtkm::Id k = idx000[2]; k <= idx111[2]; k++)
          {
            numBins++;
          }
    }

  private:
    template <typename PointsVecType>
    VTKM_EXEC vtkm::Bounds ComputeCellBounds(const PointsVecType& points)
    {
      using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
      auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

      CoordsType minp = points[0], maxp = points[0];
      for (vtkm::IdComponent i = 1; i < numPoints; ++i)
      {
        minp = vtkm::Min(minp, points[i]);
        maxp = vtkm::Max(maxp, points[i]);
      }

      return { vtkm::Vec3f(minp), vtkm::Vec3f(maxp) };
    }

>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    vtkm::Vec3f InvSpacing;
    vtkm::Id3 MaxCellIds;
    vtkm::Vec3f Origin;
  };

  class RecordBinsPerCell : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cellset,
                                  FieldInPoint coords,
                                  FieldInCell start,
                                  WholeArrayInOut binsPerCell,
<<<<<<< HEAD
                                  WholeArrayInOut cellIds,
                                  AtomicArrayInOut cellCounts);
    using ExecutionSignature = void(InputIndex, _2, _3, _4, _5, _6);
=======
                                  WholeArrayInOut cellIds);
    using ExecutionSignature = void(InputIndex, _2, _3, _4, _5);
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    using InputDomain = _1;

    RecordBinsPerCell(const vtkm::Vec3f& origin,
                      const vtkm::Vec3f& invSpacing,
                      const vtkm::Id3& dims,
                      const vtkm::Id3& maxCellIds)
      : Dims(dims)
      , InvSpacing(invSpacing)
      , MaxCellIds(maxCellIds)
      , Origin(origin)
    {
    }


<<<<<<< HEAD
    template <typename PointsVecType, typename ResultArrayType, typename CellCountType>
=======
    template <typename PointsVecType, typename ResultArrayType>
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    VTKM_EXEC void operator()(const vtkm::Id& cellIdx,
                              const PointsVecType& points,
                              const vtkm::Id& start,
                              ResultArrayType& binsPerCell,
<<<<<<< HEAD
                              ResultArrayType& cellIds,
                              CellCountType cellCounts) const
    {
      vtkm::Id3 idx000, idx111;
      MinMaxIndicesForCellPoints(
        points, this->Origin, this->InvSpacing, this->MaxCellIds, idx000, idx111);

      vtkm::Id cnt = 0;
      vtkm::Id sliceStart = ComputeFlatIndex(idx000, this->Dims);
      for (vtkm::Id k = idx000[2]; k <= idx111[2]; k++)
      {
        vtkm::Id shaftStart = sliceStart;
        for (vtkm::Id j = idx000[1]; j <= idx111[1]; j++)
        {
          vtkm::Id flatIdx = shaftStart;
          for (vtkm::Id i = idx000[0]; i <= idx111[0]; i++)
          {
            // set portals and increment cnt...
            binsPerCell.Set(start + cnt, flatIdx);
            cellIds.Set(start + cnt, cellIdx);
            cellCounts.Add(flatIdx, 1);
            ++flatIdx;
            ++cnt;
          }
          shaftStart += this->Dims[0];
        }
        sliceStart += this->Dims[0] * this->Dims[1];
      }
=======
                              ResultArrayType& cellIds) const
    {
      auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

      vtkm::Vec3f bbox[2] = { points[0], points[0] };
      for (vtkm::IdComponent i = 1; i < numPoints; ++i)
      {
        for (vtkm::IdComponent j = 0; j < 3; j++)
        {
          bbox[0][j] = vtkm::Min(bbox[0][j], static_cast<vtkm::FloatDefault>(points[i][j]));
          bbox[1][j] = vtkm::Max(bbox[1][j], static_cast<vtkm::FloatDefault>(points[i][j]));
        }
      }

      //Get 8 corners of bbox.
      auto idx000 = GetFlatIndex(bbox[0], this->Origin, this->InvSpacing, this->MaxCellIds);
      auto idx111 = GetFlatIndex(bbox[1], this->Origin, this->InvSpacing, this->MaxCellIds);

      //Set the indices and counts for each bin
      vtkm::Id cnt = 0;
      for (vtkm::Id i = idx000[0]; i <= idx111[0]; i++)
        for (vtkm::Id j = idx000[1]; j <= idx111[1]; j++)
          for (vtkm::Id k = idx000[2]; k <= idx111[2]; k++)
          {
            vtkm::Id flatIdx = ComputeFlatIndex(vtkm::Id3(i, j, k), this->Dims);
            binsPerCell.Set(start + cnt, flatIdx);
            cellIds.Set(start + cnt, cellIdx);
            cnt++;
          }
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    }

  private:
    vtkm::Id3 Dims;
    vtkm::Vec3f InvSpacing;
    vtkm::Id3 MaxCellIds;
    vtkm::Vec3f Origin;
  };

<<<<<<< HEAD
=======
  class CountBins : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn input, AtomicArrayInOut out);
    using ExecutionSignature = void(_1, _2);
    using InputDomain = _1;

    template <typename AtomicArrayType>
    VTKM_EXEC void operator()(vtkm::Id value, const AtomicArrayType& bins) const
    {
      bins.Add(value, 1);
    }
  };

>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  } //namespace detail


  //----------------------------------------------------------------------------
  /// Builds the cell locator lookup structure
  ///
  VTKM_CONT void CellLocatorUniformBins::Build()
  {
    if (this->UniformDims[0] <= 0 || this->UniformDims[1] <= 0 || this->UniformDims[2] <= 0)
      throw vtkm::cont::ErrorBadValue("Grid dimensions of CellLocatorUniformBins must be > 0");

    VTKM_LOG_SCOPE(vtkm::cont::LogLevel::Perf, "CellLocatorUniformBins::Build");

    this->MaxCellIds = (vtkm::Max(this->UniformDims - vtkm::Id3(1), vtkm::Id3(0)));
    vtkm::Id totalNumBins = this->UniformDims[0] * this->UniformDims[1] * this->UniformDims[2];

    vtkm::cont::Invoker invoker;

    auto cellset = this->GetCellSet();
    const auto& coords = this->GetCoordinates();

    auto bounds = coords.GetBounds();
<<<<<<< HEAD
    this->Origin = vtkm::Vec3f(bounds.MinCorner());
    this->MaxPoint = vtkm::Vec3f(bounds.MaxCorner());

=======
    this->Origin = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Min),
                               static_cast<vtkm::FloatDefault>(bounds.Y.Min),
                               static_cast<vtkm::FloatDefault>(bounds.Z.Min));
    this->MaxPoint = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Max),
                                 static_cast<vtkm::FloatDefault>(bounds.Y.Max),
                                 static_cast<vtkm::FloatDefault>(bounds.Z.Max));
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    auto size = this->MaxPoint - this->Origin;
    vtkm::Vec3f spacing(size[0] / this->UniformDims[0],
                        size[1] / this->UniformDims[1],
                        size[2] / this->UniformDims[2]);

    for (vtkm::IdComponent i = 0; i < 3; i++)
    {
      if (vtkm::Abs(spacing[i]) > 0)
        this->InvSpacing[i] = 1.0f / spacing[i];
      else
        this->InvSpacing[i] = 0;
    }

<<<<<<< HEAD
    // The following example will be used in the explanation below
    // Dataset with 3 cells: c0, c1, c2
    // 2x2 uniform grid: b0, b1, b2, b3
    // Assume that the bounding box for each cell overlaps as follows:
    // c0: b0, b1, b2
    // c1: b1
    // c2: b2
    //
    // The acceleration structure is an array of cell ids that are grouped
    // by the overlapping bin. This information can be represented using
    // vtkm::cont::ArrayHandleGroupVecVariable
    // In the example above:
    // CellIds = {c0,  c0,c1,  c0,c2,  }
    //            b0    b1       b2   b3
    //
    // The algorithm runs as follows:
    //  Given a point p, find the bin (b) that contains p.
    //  Do a point-in-cell test for each cell in bin b.
    //
    // Example:  point p is in b=1
    // vtkm::cont::ArrayHandleGroupVecVariable provides the offset and number
    // cells in bin b=1. The offset is 1 and the count is 2.
    //  Do a point-in-cell test on the 2 cells that start at offset 1
    //    CellIds[ 1 + 0], which is c0
    //    CellIds[ 1 + 1], which is c1


    //Computing this involves several steps which are described below.
    //Step 1:
    // For each cell in cellSet, count the number of bins that overlap with the cell bbox.
    // For the example above
    // binCountsPerCell = {3, 1, 1}
    // cell0 overlaps with 3 bins
    // cell1 overlaps with 1 bin
    // cell2 overlaps with 1 bin
    vtkm::cont::ArrayHandle<vtkm::Id> binCountsPerCell;
    CountCellBins countCellBins(this->Origin, this->InvSpacing, this->MaxCellIds);
    invoker(countCellBins, cellset, coords, binCountsPerCell);

    //2: Compute number of unique cell/bin pairs and start indices.

    //Step 2:
    // Given the number of bins for each cell, we can compute the offset for each cell.
    // For the example above, binOffset is:
    // {0, 3, 4}
    // and the total number, num is 5
    // c0 begins at idx=0
    // c1 begins at idx=3
    // c2 begins at idx=4
    vtkm::cont::ArrayHandle<vtkm::Id> binOffset;
    auto num = vtkm::cont::Algorithm::ScanExclusive(binCountsPerCell, binOffset);

    //Step 3:
    // Now that we know the start index and numbers, we can fill an array of bin ids.
    // binsPerCell is the list of binIds for each cell. In the example above,
    // binsPerCell = {b0,b1,b2,   b2,       b2}
    //               \ cell0 /    cell1    cell2
    // We can also compute the cellIds and number of cells per bin, cellCount
    // cids      = {c0,c0,c0, c1, c2}
    // cellCount = {3, 1, 1, 0}
    // These are set using RecordBinsPerCell worklet, which does the following
    // for each cell
    //   compute cell bbox and list of overlaping bins
    //   for each overlaping bin
    //     add the bin id to binsPerCell starting at binOffset
    //     add the cell id to the CellIds starting at binOffset
    //     increment CellCount for the bin (uses an atomic for thread safety).

    vtkm::cont::ArrayHandle<vtkm::Id> binsPerCell, cids, cellCount;
    binsPerCell.AllocateAndFill(num, 0);
    cids.Allocate(num);
    cellCount.AllocateAndFill(totalNumBins, 0);
    RecordBinsPerCell recordBinsPerCell(
      this->Origin, this->InvSpacing, this->UniformDims, this->MaxCellIds);
    invoker(recordBinsPerCell, cellset, coords, binOffset, binsPerCell, cids, cellCount);

    //Step 4:
    // binsPerCell is the overlapping bins for each cell.
    // We want to sort CellIds by the bin ID.  SortByKey does this.
    vtkm::cont::Algorithm::SortByKey(binsPerCell, cids);

    // Convert the cell counts to offsets using the helper function
    // vtkm::cont::ConvertNumComponentsToOffsets, and create the
    // CellIds that are used as the acceleration structure.
    this->CellIds = vtkm::cont::make_ArrayHandleGroupVecVariable(
      cids, vtkm::cont::ConvertNumComponentsToOffsets(cellCount));
=======
    //1: Count number of (cell,bin) pairs.
    vtkm::cont::ArrayHandle<vtkm::Id> binCounts;
    detail::CountCellBins countCellBins(this->Origin, this->InvSpacing, this->MaxCellIds);
    invoker(countCellBins, cellset, coords, binCounts);

    //2: Compute number of unique cell/bin pairs and start indices.
    vtkm::cont::ArrayHandle<vtkm::Id> binStartIdx;
    auto num = vtkm::cont::Algorithm::ScanExclusive(binCounts, binStartIdx);

    //3: Allocate and set the cellids for each bin.
    //   uses an atomic to ensure each bin is updated correctly.
    vtkm::cont::ArrayHandle<vtkm::Id> binsPerCell;
    binsPerCell.AllocateAndFill(num, 0);
    this->CellIds.Allocate(num);
    detail::RecordBinsPerCell recordBinsPerCell(
      this->Origin, this->InvSpacing, this->UniformDims, this->MaxCellIds);
    invoker(recordBinsPerCell, cellset, coords, binStartIdx, binsPerCell, this->CellIds);

    //4: Sort CellIds by the bins in each cell.
    vtkm::cont::Algorithm::SortByKey(binsPerCell, this->CellIds);


    //5: Calculate counts for each bin and the start index.
    this->CellCount.AllocateAndFill(totalNumBins, 0);
    invoker(detail::CountBins{}, binsPerCell, this->CellCount);
    vtkm::cont::Algorithm::ScanExclusive(this->CellCount, this->CellStartIdx);
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  }

  //----------------------------------------------------------------------------
  struct CellLocatorUniformBins::MakeExecObject
  {
    template <typename CellSetType>
    VTKM_CONT void operator()(const CellSetType& cellSet,
                              vtkm::cont::DeviceAdapterId device,
                              vtkm::cont::Token& token,
                              const CellLocatorUniformBins& self,
                              ExecObjType& execObject) const
    {
<<<<<<< HEAD
      using CellStructureType = CellSetContToExec<CellSetType>;

      execObject = vtkm::exec::CellLocatorUniformBins<CellStructureType>(self.UniformDims,
                                                                         self.Origin,
                                                                         self.MaxPoint,
                                                                         self.InvSpacing,
                                                                         self.MaxCellIds,
                                                                         self.CellIds,
                                                                         cellSet,
                                                                         self.GetCoordinates(),
                                                                         device,
                                                                         token);
=======
      using CellStructuredType = CellSetContToExec<CellSetType>;

      execObject = vtkm::exec::CellLocatorUniformBins<CellStructuredType>(self.UniformDims,
                                                                          self.Origin,
                                                                          self.MaxPoint,
                                                                          self.InvSpacing,
                                                                          self.MaxCellIds,
                                                                          self.CellCount,
                                                                          self.CellStartIdx,
                                                                          self.CellIds,
                                                                          cellSet,
                                                                          self.GetCoordinates(),
                                                                          device,
                                                                          token);
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
    }
  };

  CellLocatorUniformBins::ExecObjType CellLocatorUniformBins::PrepareForExecution(
    vtkm::cont::DeviceAdapterId device,
    vtkm::cont::Token& token) const
  {
    this->Update();
    ExecObjType execObject;
    vtkm::cont::CastAndCall(this->GetCellSet(), MakeExecObject{}, device, token, *this, execObject);
    return execObject;
  }

  //----------------------------------------------------------------------------
  void CellLocatorUniformBins::PrintSummary(std::ostream& out) const
  {
    out << std::endl;
    out << "CellLocatorUniformBins" << std::endl;
    out << " UniformDims: " << this->UniformDims << std::endl;
    out << " Origin: " << this->Origin << std::endl;
    out << " MaxPoint: " << this->MaxPoint << std::endl;
    out << " InvSpacing: " << this->InvSpacing << std::endl;
    out << " MaxCellIds: " << this->MaxCellIds << std::endl;

    out << "Input CellSet: \n";
    this->GetCellSet().PrintSummary(out);
    out << "Input Coordinates: \n";
    this->GetCoordinates().PrintSummary(out);

<<<<<<< HEAD
    /*
=======
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  out << "  CellStartIdx:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellStartIdx, out);
  out << "  CellCount:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellCount, out);
  out << "  CellIds:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellIds, out);
<<<<<<< HEAD
  */
=======
>>>>>>> df93dd221... Uniform bins cell locator for unstructured grids.
  }
  }
} // vtkm::cont
