//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_partial_compositor_h
#define vtkm_rendering_compositing_partial_compositor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <iostream>
#include <vector>
#include <vtkm/Types.h>
#include <vtkm/rendering/compositing/AbsorptionPartial.h>
#include <vtkm/rendering/compositing/EmissionPartial.h>
#include <vtkm/rendering/compositing/VolumePartial.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

template <typename PartialType>
class VTKM_RENDERING_EXPORT PartialCompositor
{
public:
  PartialCompositor();
  ~PartialCompositor();
  void composite(std::vector<std::vector<PartialType>>& partial_images,
                 std::vector<PartialType>& output_partials);
  void set_background(std::vector<vtkm::Float32>& background_values);
  void set_background(std::vector<vtkm::Float64>& background_values);
  void set_comm_handle(int mpi_comm_id);

protected:
  void merge(const std::vector<std::vector<PartialType>>& in_partials,
             std::vector<PartialType>& partials,
             int& global_min_pixel,
             int& global_max_pixel);

  void composite_partials(std::vector<PartialType>& partials,
                          std::vector<PartialType>& output_partials);

  std::vector<typename PartialType::ValueType> m_background_values;
  int m_mpi_comm_id;
};


}
}
} //vtkm::rendering::compositing

#endif //vtkm_rendering_compositing_partial_compositor_h
