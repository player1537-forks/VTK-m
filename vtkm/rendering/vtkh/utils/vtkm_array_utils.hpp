//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTKH_VTKM_ARRAY_UTILS_HPP
#define VTKH_VTKM_ARRAY_UTILS_HPP

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/PartitionedDataSet.h>

namespace vtkh
{

template<typename T>
T *
GetVTKMPointer(vtkm::cont::ArrayHandle<T> &handle)
{
  return handle.WritePortal().GetArray();
}

bool
GlobalHasOnePartition(const vtkm::cont::PartitionedDataSet& pds);

bool
IsEmpty(const vtkm::cont::PartitionedDataSet& pds);

bool
GlobalIsEmpty(const vtkm::cont::PartitionedDataSet& pds);


bool
GlobalHasField(const vtkm::cont::PartitionedDataSet& pds, const std::string& fieldName);

bool
IsPointMesh(const vtkm::cont::PartitionedDataSet& pds);

}//namespace vtkh
#endif
