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
#include <vtkm/cont/CellLocatorUniformBins.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

namespace vtkm
{
namespace cont
{

namespace detail
{

VTKM_EXEC
static inline vtkm::Id ComputeFlatIndex(const vtkm::Id3& idx, const vtkm::Id3& dims)
{
  return idx[0] + (dims[0] * (idx[1] + (dims[1] * idx[2])));
}

VTKM_EXEC
static inline vtkm::Id3 GetFlatIndex(const vtkm::Vec3f& pt,
                                     const vtkm::Vec3f& origin,
                                     const vtkm::Vec3f& invSpacing,
                                     const vtkm::Id3& maxCellIds)
{
  auto temp = pt - origin;
  temp = temp * invSpacing;

  auto logicalIdx = vtkm::Min(vtkm::Id3(temp), maxCellIds);
  return logicalIdx;
}

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
    auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

    vtkm::Bounds bounds;
    for (vtkm::IdComponent i = 0; i < numPoints; ++i)
    {
      bounds.Include(points[i]);
    }

    //Get 8 corners of bbox.
    auto idx000 = GetFlatIndex({ bounds.Min.X, bounds.Min.Y, bounds.Min.Z },
                               this->Origin,
                               this->InvSpacing,
                               this->MaxCellIds);
    auto idx111 = GetFlatIndex({ bounds.Max.X, bounds.Max.Y, bounds.Max.Z },
                               this->Origin,
                               this->InvSpacing,
                               this->MaxCellIds);

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
                                WholeArrayInOut cellIds);
  using ExecutionSignature = void(InputIndex, _2, _3, _4, _5);
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


  template <typename PointsVecType, typename ResultArrayType>
  VTKM_EXEC void operator()(const vtkm::Id& cellIdx,
                            const PointsVecType& points,
                            const vtkm::Id& start,
                            ResultArrayType& binsPerCell,
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
  }

private:
  vtkm::Id3 Dims;
  vtkm::Vec3f InvSpacing;
  vtkm::Id3 MaxCellIds;
  vtkm::Vec3f Origin;
};

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
  this->Origin = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Min),
                             static_cast<vtkm::FloatDefault>(bounds.Y.Min),
                             static_cast<vtkm::FloatDefault>(bounds.Z.Min));
  this->MaxPoint = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Max),
                               static_cast<vtkm::FloatDefault>(bounds.Y.Max),
                               static_cast<vtkm::FloatDefault>(bounds.Z.Max));
  auto size = this->MaxPoint - this->Origin;
  vtkm::Vec3f spacing(
    size[0] / this->UniformDims[0], size[1] / this->UniformDims[1], size[2] / this->UniformDims[2]);

  for (vtkm::IdComponent i = 0; i < 3; i++)
  {
    if (vtkm::Abs(spacing[i]) > 0)
      this->InvSpacing[i] = 1.0f / spacing[i];
    else
      this->InvSpacing[i] = 0;
  }

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

  out << "  CellStartIdx:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellStartIdx, out);
  out << "  CellCount:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellCount, out);
  out << "  CellIds:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellIds, out);
}
}
} // vtkm::cont
