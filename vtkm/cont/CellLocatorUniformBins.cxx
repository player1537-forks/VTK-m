//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#include <vtkm/cont/CellLocatorUniformBins.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

namespace vtkm
{
namespace cont
{

class CountCellBins : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn cellset, FieldInPoint coords, FieldOutCell bincount);
  using ExecutionSignature = void(InputIndex, _2, _3);
  using InputDomain = _1;

  CountCellBins(const vtkm::Vec3f& origin,
                const vtkm::Vec3f& invSpacing,
                const vtkm::Id3& dims,
                const vtkm::Id3& maxCellIds)
    : Dims(dims)
    , InvSpacing(invSpacing)
    , MaxCellIds(maxCellIds)
    , Origin(origin)
  {
  }

  template <typename PointsVecType>
  VTKM_EXEC void operator()(const vtkm::Id& /*cellIdx*/, const PointsVecType& points, vtkm::Id& numBins) const
  {
    //using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
    auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

    vtkm::Vec3f bbox[2] = {points[0], points[0]};
    for (vtkm::IdComponent i = 1; i < numPoints; ++i)
    {
      for (vtkm::IdComponent j = 0; j < 3; j++)
      {
        bbox[0][j] = vtkm::Min(bbox[0][j], static_cast<vtkm::FloatDefault>(points[i][j]));
        bbox[1][j] = vtkm::Max(bbox[1][j], static_cast<vtkm::FloatDefault>(points[i][j]));
      }
    }

    //std::cout<<"CellId: "<<cellIdx<<" bbox: "<<bbox[0]<<" "<<bbox[1]<<std::endl;

    //Get 8 corners of bbox.
    auto idx000 = this->GetFlatIndex(bbox[0]);
    auto idx111 = this->GetFlatIndex(bbox[1]);

    //avoid the loops and do some math..
    numBins = 0;
    for (vtkm::Id i = idx000[0]; i <= idx111[0]; i++)
      for (vtkm::Id j = idx000[1]; j <= idx111[1]; j++)
        for (vtkm::Id k = idx000[2]; k <= idx111[2]; k++)
        {
          numBins++;
          //vtkm::Id flatIdx = this->ComputeFlatIndex(vtkm::Id3(i,j,k));
          //std::cout<<"  flatIdx= "<<flatIdx<<std::endl;
        }
    //std::cout<<"Cell: "<<cellIdx<<" nb= "<<numBins<<std::endl;
  }

  VTKM_EXEC vtkm::Id ComputeFlatIndex(const vtkm::Id3& idx) const
  {
    return idx[0] + (this->Dims[0] * (idx[1] + (this->Dims[1] * idx[2])));
  }

private:

  vtkm::Id3 GetFlatIndex(const vtkm::Vec3f& pt) const
  {
    vtkm::Id3 logicalIdx(0, 0, 0);
    auto temp = pt - this->Origin;
    temp = temp * this->InvSpacing;

    logicalIdx = vtkm::Min(vtkm::Id3(temp), this->MaxCellIds);
    //std::cout<<"GetFlatIndex: "<<pt<<" --> "<<logicalIdx<<std::endl;
    return logicalIdx;
  }

  struct BinsBBox
  {
    vtkm::Id3 Min;
    vtkm::Id3 Max;

    VTKM_EXEC
    bool Empty() const
    {
      return (this->Max[0] < this->Min[0]) || (this->Max[1] < this->Min[1]) ||
        (this->Max[2] < this->Min[2]);
    }
  };

  template <typename PointsVecType>
  VTKM_EXEC
  vtkm::Bounds ComputeCellBounds(const PointsVecType& points)
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

  VTKM_EXEC BinsBBox ComputeIntersectingBins(const Bounds cellBounds)
  {
    vtkm::Vec3f cbMin(cellBounds.X.Min, cellBounds.Y.Min, cellBounds.Z.Min);
    vtkm::Vec3f cbMax(cellBounds.X.Max, cellBounds.Y.Max, cellBounds.Z.Max);
    auto minb = static_cast<vtkm::Id3>((cbMin - this->Origin) * this->InvSpacing);
    auto maxb = static_cast<vtkm::Id3>((cbMax - this->Origin) * this->InvSpacing);

    return { vtkm::Max(vtkm::Id3(0), minb), vtkm::Min(this->Dims - vtkm::Id3(1), maxb) };
  }

  vtkm::Id3 Dims;
  vtkm::Vec3f InvSpacing;
  vtkm::Id3 MaxCellIds;
  vtkm::Vec3f Origin;
};

class RecordBinsPerCell : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn cellset, FieldInPoint coords, FieldInCell count, FieldInCell start, WholeArrayInOut result, WholeArrayInOut result2);
  using ExecutionSignature = void(InputIndex, _2, _3, _4, _5, _6);
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
                            const vtkm::Id& count,
                            const vtkm::Id& start,
                            ResultArrayType& result,
                            ResultArrayType& result2) const
  {
    //using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
    auto numPoints = vtkm::VecTraits<PointsVecType>::GetNumberOfComponents(points);

    vtkm::Vec3f bbox[2] = {points[0], points[0]};
    for (vtkm::IdComponent i = 1; i < numPoints; ++i)
    {
      for (vtkm::IdComponent j = 0; j < 3; j++)
      {
        bbox[0][j] = vtkm::Min(bbox[0][j], static_cast<vtkm::FloatDefault>(points[i][j]));
        bbox[1][j] = vtkm::Max(bbox[1][j], static_cast<vtkm::FloatDefault>(points[i][j]));
      }
    }

    //std::cout<<"CellId: "<<cellIdx<<" bbox: "<<bbox[0]<<" "<<bbox[1]<<std::endl;

    //Get 8 corners of bbox.
    auto idx000 = this->GetFlatIndex(bbox[0]);
    auto idx111 = this->GetFlatIndex(bbox[1]);

    //avoid the loops and do some math..
    vtkm::Id cnt = 0;
    for (vtkm::Id i = idx000[0]; i <= idx111[0]; i++)
      for (vtkm::Id j = idx000[1]; j <= idx111[1]; j++)
        for (vtkm::Id k = idx000[2]; k <= idx111[2]; k++)
        {
          vtkm::Id flatIdx = this->ComputeFlatIndex(vtkm::Id3(i,j,k));
          result.Set(start + cnt, flatIdx);
          result2.Set(start + cnt, cellIdx);
          cnt++;
          if (cnt > count) std::cout<<"Overflow!!!!!!"<<std::endl;
        }
  }

  VTKM_EXEC vtkm::Id3 GetFlatIndex(const vtkm::Vec3f& pt) const
  {
    vtkm::Id3 logicalIdx(0, 0, 0);
    auto temp = pt - this->Origin;
    temp = temp * this->InvSpacing;

    logicalIdx = vtkm::Min(vtkm::Id3(temp), this->MaxCellIds);
    //std::cout<<"GetFlatIndex: "<<pt<<" --> "<<logicalIdx<<std::endl;
    return logicalIdx;
  }

  VTKM_EXEC vtkm::Id ComputeFlatIndex(const vtkm::Id3& idx) const
  {
    return idx[0] + (this->Dims[0] * (idx[1] + (this->Dims[1] * idx[2])));
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
  VTKM_EXEC void operator() (vtkm::Id value, const AtomicArrayType& bins) const
  {
    bins.Add(value, 1);
  }
};


//----------------------------------------------------------------------------
/// Builds the cell locator lookup structure
///
VTKM_CONT void CellLocatorUniformBins::Build()
{
  std::cout<<"CellLocatorUniformBins::Build()"<<std::endl;

  VTKM_LOG_SCOPE(vtkm::cont::LogLevel::Perf, "CellLocatorUniformBins::Build");

  this->MaxCellIds = (vtkm::Max(this->UniformCellDims - vtkm::Id3(1), vtkm::Id3(0)));
  vtkm::Id totalNumBins = this->UniformCellDims[0] * this->UniformCellDims[1] * this->UniformCellDims[2];

  vtkm::cont::Invoker invoker;

  auto cellset = this->GetCellSet();
  const auto& coords = this->GetCoordinates();

  // 1: Compute the top level grid
  auto bounds = coords.GetBounds();
  this->Origin = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Min),
                             static_cast<vtkm::FloatDefault>(bounds.Y.Min),
                             static_cast<vtkm::FloatDefault>(bounds.Z.Min));
  this->MaxPoint = vtkm::Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Max),
                               static_cast<vtkm::FloatDefault>(bounds.Y.Max),
                               static_cast<vtkm::FloatDefault>(bounds.Z.Max));
  auto size = this->MaxPoint - this->Origin;
  vtkm::Vec3f spacing(size[0] / this->UniformCellDims[0],
                      size[1] / this->UniformCellDims[1],
                      size[2] / this->UniformCellDims[2]);

  this->InvSpacing = 1.0f / spacing;
  /*
  this->InvSpacing[2] = 0.0; std::cout<<"FIX ME ****"<<std::endl;
  */

  //1: Count number of (cell,bin) pairs.
  std::cout<<"countBins: "<<this->Origin<<" "<<this->MaxPoint<<" "<<this->InvSpacing<<" "<<this->UniformCellDims<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> binCounts;
  CountCellBins countCellBins(this->Origin, this->InvSpacing, this->UniformCellDims, this->MaxCellIds);
  invoker(countCellBins, cellset, coords, binCounts);

  /*
  std::cout<<"*************************************"<<std::endl;
  std::cout<<"Cell_i is in this many bins"<<std::endl;
  std::cout<<"binCounts: ";
  vtkm::cont::printSummary_ArrayHandle(binCounts, std::cout);
  */

  //2: number of unique cell/bin pairs.
  vtkm::cont::ArrayHandle<vtkm::Id> binStartIdx;
  auto num = vtkm::cont::Algorithm::ScanExclusive(binCounts, binStartIdx);
  /*
  std::cout<<"start index for each cell"<<std::endl;
  std::cout<<"binStartIdx: ";
  vtkm::cont::printSummary_ArrayHandle(binStartIdx, std::cout);
  std::cout<<"Num= "<<num<<std::endl;
  */

  //3: Allocate and set the cellids for each bin.
  vtkm::cont::ArrayHandle<vtkm::Id> binsPerCell;
  binsPerCell.AllocateAndFill(num, 0);
  this->CellIds.Allocate(num);
  RecordBinsPerCell recordBinsPerCell(this->Origin, this->InvSpacing, this->UniformCellDims, this->MaxCellIds);
  invoker(recordBinsPerCell, cellset, coords, binCounts, binStartIdx, binsPerCell, this->CellIds);
  /*
  std::cout<<"List of bins for each cell"<<std::endl;
  std::cout<<"binsPerCell"<<std::endl;
  vtkm::cont::printSummary_ArrayHandle(binsPerCell, std::cout);
  std::cout<<"List of cellIds for each entry above."<<std::endl;
  std::cout<<"this->CellIds"<<std::endl;
  vtkm::cont::printSummary_ArrayHandle(this->CellIds, std::cout);
  */

  //4: Sort CellIds by the bins in each cell.
  vtkm::cont::Algorithm::SortByKey(binsPerCell, this->CellIds);
  /*
  std::cout<<"***** The elusive array of cells!"<<std::endl;
  vtkm::cont::printSummary_ArrayHandle(this->CellIds, std::cout);
  */


  //5: Calculate counts for each bin and the start index.
  //vtkm::cont::ArrayHandle<vtkm::Id> binCellCount;
  //binCellCount.AllocateAndFill(totalNumBins, 0);
  this->CellCount.AllocateAndFill(totalNumBins, 0);
  invoker(CountBins{}, binsPerCell, this->CellCount);
  /*
  std::cout<<"this->CellCount"<<std::endl;
  std::cout<<" number of cells in each bin"<<std::endl;
  vtkm::cont::printSummary_ArrayHandle(this->CellCount, std::cout);
  */

  vtkm::cont::Algorithm::ScanExclusive(this->CellCount, this->CellStartIdx);
  /*
  std::cout<<"this->CellStartIdx"<<std::endl;
  vtkm::cont::printSummary_ArrayHandle(this->CellStartIdx, std::cout);
  */

  auto minCount = vtkm::cont::Algorithm::Reduce(this->CellCount, vtkm::Id(this->CellCount.ReadPortal().Get(0)), vtkm::Minimum());
  auto maxCount = vtkm::cont::Algorithm::Reduce(this->CellCount, vtkm::Id(0), vtkm::Maximum());
  std::cout<<"CellLocatorUniformBins:: Cell counts: "<<minCount<<" "<<maxCount<<std::endl;
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

    execObject = vtkm::exec::CellLocatorUniformBins<CellStructuredType>(self.UniformCellDims,
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
  out<<std::endl;

  /*
  out << "DensityL1: " << this->DensityL1 << "\n";
  out << "DensityL2: " << this->DensityL2 << "\n";
  out << "Input CellSet: \n";
  this->GetCellSet().PrintSummary(out);
  out << "Input Coordinates: \n";
  this->GetCoordinates().PrintSummary(out);
  out << "LookupStructure:\n";
  out << "  TopLevelGrid\n";
  out << "    Dimensions: " << this->TopLevel.Dimensions << "\n";
  out << "    Origin: " << this->TopLevel.Origin << "\n";
  out << "    BinSize: " << this->TopLevel.BinSize << "\n";
  out << "  LeafDimensions:\n";
  vtkm::cont::printSummary_ArrayHandle(this->LeafDimensions, out);
  out << "  LeafStartIndex:\n";
  vtkm::cont::printSummary_ArrayHandle(this->LeafStartIndex, out);
  out << "  CellStartIndex:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellStartIndex, out);
  out << "  CellCount:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellCount, out);
  out << "  CellIds:\n";
  vtkm::cont::printSummary_ArrayHandle(this->CellIds, out);
  */
}
}
} // vtkm::cont
