//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_PartialCompoistor_h
#define vtkm_rendering_new_PartialCompoistor_h

#include <iostream>
#include <vector>
#include <vtkm/Types.h>
#include <vtkm/rendering_new/compositing/AbsorptionPartial.h>
#include <vtkm/rendering_new/compositing/EmissionPartial.h>
#include <vtkm/rendering_new/compositing/VolumePartial.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#endif


namespace vtkm
{
namespace rendering_new
{

template <typename PartialType>
class VTKM_RENDERING_NEW_EXPORT PartialCompositor
{
public:
  PartialCompositor();
  ~PartialCompositor();
  void composite(std::vector<std::vector<PartialType>>& partial_images,
                 std::vector<PartialType>& output_partials);
  void set_background(std::vector<vtkm::Float32>& background_values);
  void set_background(std::vector<vtkm::Float64>& background_values);
#ifdef VTKM_ENABLE_MPI
  //  void set_comm_handle(int mpi_comm_id);
  void set_comm_handle(MPI_Comm mpi_comm) { this->mpiComm = mpi_comm; }
#endif
protected:
  void merge(const std::vector<std::vector<PartialType>>& in_partials,
             std::vector<PartialType>& partials,
             int& global_min_pixel,
             int& global_max_pixel);

  void composite_partials(std::vector<PartialType>& partials,
                          std::vector<PartialType>& output_partials);

  std::vector<typename PartialType::ValueType> BackgroundValues;
#ifdef VTKM_ENABLE_MPI
  MPI_Comm mpiComm;
#endif
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_new_PartialCompoistor_h
