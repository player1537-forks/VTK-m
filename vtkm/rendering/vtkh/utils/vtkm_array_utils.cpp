//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>

#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/rendering/vtkh/utils/vtkm_dataset_info.hpp>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
#include <vtkm/thirdparty/diy/mpi-cast.h>
#include <mpi.h>
#endif

namespace vtkh
{

namespace detail_utils
{
bool GlobalAgreement(bool local)
{
  bool agreement = local;
#ifdef VTKM_ENABLE_MPI
  int local_boolean = local ? 1 : 0;
  int global_boolean;
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MPI_Comm mpi_comm = vtkmdiy::mpi::mpi_cast(diy_comm.handle());
  MPI_Allreduce((void *)(&local_boolean),
                (void *)(&global_boolean),
                1,
                MPI_INT,
                MPI_SUM,
                mpi_comm);

  if(global_boolean != diy_comm.size())
  {
    agreement = false;
  }
#endif
  return agreement;
}
} //detail


bool
IsEmpty(const vtkm::cont::PartitionedDataSet& pds)
{
  for (const auto& ds : pds.GetPartitions())
  {
    auto cellset = ds.GetCellSet();
    if (cellset.GetNumberOfCells() > 0)
      return false;
  }

  return true;
}


bool
GlobalHasOnePartition(const vtkm::cont::PartitionedDataSet& pds)
{
  bool one = pds.GetNumberOfPartitions() == 1;
  return detail_utils::GlobalAgreement(one);
}

bool
GlobalHasField(const vtkm::cont::PartitionedDataSet& pds, const std::string& fieldName)
{
  bool haveLocal = true;
  for (const auto& ds : pds.GetPartitions())
    if (!ds.HasField(fieldName))
    {
      haveLocal = false;
      break;
    }

  return detail_utils::GlobalAgreement(haveLocal);
}

bool
GlobalIsEmpty(const vtkm::cont::PartitionedDataSet& pds)
{
  return detail_utils::GlobalAgreement(vtkh::IsEmpty(pds));
}

bool
IsPointMesh(const vtkm::cont::PartitionedDataSet& pds)
{
  const bool is_empty = GlobalIsEmpty(pds);
  if (is_empty)
    return false;

  // since we are not empty, start with the affirmative is_points.
  // if someone is not points, the we will figure it out here
  bool is_points = true;

  for (const auto& ds : pds.GetPartitions())
  {
    vtkm::UInt8 shape_type;
    bool single_type = vtkh::VTKMDataSetInfo::IsSingleCellShape(ds.GetCellSet(), shape_type);

    if (ds.GetCellSet().GetNumberOfCells() > 0)
    {
      is_points = (single_type && (shape_type == 1)) && is_points;
    }
  }

  is_points = detail_utils::GlobalAgreement(is_points);
  return is_points;
}

} //namespace vtkh
