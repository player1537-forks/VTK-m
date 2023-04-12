//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/StaticAssert.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/BoundsCompute.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/internal/Configure.h>

// clang-format off
VTKM_THIRDPARTY_PRE_INCLUDE
#include <vtkm/thirdparty/diy/diy.h>
VTKM_THIRDPARTY_POST_INCLUDE
// clang-format on

namespace vtkm
{
namespace cont
{

VTKM_CONT
PartitionedDataSet::PartitionedDataSet(const vtkm::cont::DataSet& ds)
{
  this->Partitions.insert(this->Partitions.end(), ds);
}

VTKM_CONT
PartitionedDataSet::PartitionedDataSet(const std::vector<vtkm::cont::DataSet>& partitions)
{
  this->Partitions = partitions;
}

VTKM_CONT
PartitionedDataSet::PartitionedDataSet(vtkm::Id size)
{
  this->Partitions.reserve(static_cast<std::size_t>(size));
}

VTKM_CONT
vtkm::cont::Field PartitionedDataSet::GetFieldFromPartition(const std::string& field_name,
                                                            int partition_index) const
{
  assert(partition_index >= 0);
  assert(static_cast<std::size_t>(partition_index) < this->Partitions.size());
  return this->Partitions[static_cast<std::size_t>(partition_index)].GetField(field_name);
}

VTKM_CONT
vtkm::Id PartitionedDataSet::GetNumberOfPartitions() const
{
  return static_cast<vtkm::Id>(this->Partitions.size());
}

VTKM_CONT
vtkm::Id PartitionedDataSet::GetGlobalNumberOfPartitions() const
{
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  vtkm::Id globalSize = 0;
  vtkmdiy::mpi::all_reduce(comm, this->GetNumberOfPartitions(), globalSize, std::plus<vtkm::Id>{});
  return globalSize;
#else
  return this->GetNumberOfPartitions();
#endif
}

VTKM_CONT
const vtkm::cont::DataSet& PartitionedDataSet::GetPartition(vtkm::Id blockId) const
{
  return this->Partitions[static_cast<std::size_t>(blockId)];
}

VTKM_CONT
const std::vector<vtkm::cont::DataSet>& PartitionedDataSet::GetPartitions() const
{
  return this->Partitions;
}

VTKM_CONT
void PartitionedDataSet::AppendPartition(const vtkm::cont::DataSet& ds)
{
  this->Partitions.insert(this->Partitions.end(), ds);
}

VTKM_CONT
void PartitionedDataSet::AppendPartitions(const std::vector<vtkm::cont::DataSet>& partitions)
{
  this->Partitions.insert(this->Partitions.end(), partitions.begin(), partitions.end());
}

VTKM_CONT
void PartitionedDataSet::InsertPartition(vtkm::Id index, const vtkm::cont::DataSet& ds)
{
  if (index <= static_cast<vtkm::Id>(this->Partitions.size()))
  {
    this->Partitions.insert(this->Partitions.begin() + index, ds);
  }
  else
  {
    std::string msg = "invalid insert position\n ";
    throw ErrorBadValue(msg);
  }
}

VTKM_CONT
void PartitionedDataSet::ReplacePartition(vtkm::Id index, const vtkm::cont::DataSet& ds)
{
  if (index < static_cast<vtkm::Id>(this->Partitions.size()))
    this->Partitions.at(static_cast<std::size_t>(index)) = ds;
  else
  {
    std::string msg = "invalid replace position\n ";
    throw ErrorBadValue(msg);
  }
}

VTKM_CONT
void PartitionedDataSet::CopyPartitions(const vtkm::cont::PartitionedDataSet& source)
{
  this->Partitions = source.Partitions;
}

VTKM_CONT
void PartitionedDataSet::PrintSummary(std::ostream& stream) const
{
  stream << "PartitionedDataSet [" << this->Partitions.size() << " partitions]:\n";

  for (size_t part = 0; part < this->Partitions.size(); ++part)
  {
    stream << "Partition " << part << ":\n";
    this->Partitions[part].PrintSummary(stream);
  }

  stream << "  Fields[" << this->GetNumberOfFields() << "]\n";
  for (vtkm::Id index = 0; index < this->GetNumberOfFields(); index++)
  {
    this->GetField(index).PrintSummary(stream);
  }
}

VTKM_CONT
vtkm::Id PartitionedDataSet::GetNumberOfCells() const
{
  vtkm::Id numCells = 0;
  for (const auto& ds : this->Partitions)
    numCells += ds.GetNumberOfCells();

  return numCells;
}

VTKM_CONT
vtkm::Id PartitionedDataSet::GetGlobalNumberOfCells() const
{
  vtkm::Id numCells = this->GetNumberOfCells();
#ifdef VTKM_ENABLE_MPI
  vtkm::Id globalNumCells = 0;
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  vtkmdiy::mpi::all_reduce(comm, numCells, globalNumCells, std::plus<vtkm::Id>{});
  numCells = globalNumCells;
#endif

  return numCells;
}

VTKM_CONT
vtkm::Bounds PartitionedDataSet::GetBounds(vtkm::Id coordinateIndex) const
{
  return vtkm::cont::BoundsCompute(*this, coordinateIndex);
}

VTKM_CONT
vtkm::Bounds PartitionedDataSet::GetGlobalBounds(vtkm::Id coordinateIndex) const
{
  vtkm::Bounds localBounds = this->GetBounds(coordinateIndex);

#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.size() == 1)
    return localBounds;

  std::vector<vtkm::Float64> vals(3), globalMinVals(3, 0), globalMaxVals(3, 0);
  vals[0] = localBounds.X.Min;
  vals[1] = localBounds.Y.Min;
  vals[2] = localBounds.Z.Min;
  vtkmdiy::mpi::all_reduce(comm, vals, globalMinVals, vtkmdiy::mpi::minimum<vtkm::Float64>{});

  vals[0] = localBounds.X.Max;
  vals[1] = localBounds.Y.Max;
  vals[2] = localBounds.Z.Max;
  vtkmdiy::mpi::all_reduce(comm, vals, globalMaxVals, vtkmdiy::mpi::maximum<vtkm::Float64>{});

  return vtkm::Bounds(globalMinVals[0],
                      globalMaxVals[0],
                      globalMinVals[1],
                      globalMaxVals[1],
                      globalMinVals[2],
                      globalMaxVals[2]);
#else

  return localBounds;
#endif
}

VTKM_CONT
vtkm::Range PartitionedDataSet::GetScalarFieldRange(const std::string& fieldName) const
{
  vtkm::Range range;
  for (const auto& ds : this->Partitions)
  {
    if (ds.HasField(fieldName))
    {
      auto rangeArray = ds.GetField(fieldName).GetRange();
      auto portal = rangeArray.ReadPortal();
      if (portal.GetNumberOfValues() != 1)
        throw vtkm::cont::ErrorBadValue("Error in range computation: Field is not a scalar");

      range.Include(portal.Get(0));
    }
  }

  return range;
}

VTKM_CONT
vtkm::Range PartitionedDataSet::GetGlobalScalarFieldRange(const std::string& fieldName) const
{
  vtkm::Range localRange = this->GetScalarFieldRange(fieldName);
#ifdef VTKM_ENABLE_MPI
  vtkm::Float64 globalMin, globalMax;

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.size() == 0)
    return localRange;

  vtkmdiy::mpi::all_reduce(comm, localRange.Min, globalMin, vtkmdiy::mpi::minimum<vtkm::Float64>{});
  vtkmdiy::mpi::all_reduce(comm, localRange.Max, globalMax, vtkmdiy::mpi::maximum<vtkm::Float64>{});

  return vtkm::Range(globalMin, globalMax);
#else
  return localRange;
#endif
}


}
} // namespace vtkm::cont
