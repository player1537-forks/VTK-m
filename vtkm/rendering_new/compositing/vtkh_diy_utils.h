//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_vtkh_diy_utils_h
#define vtkm_rendering_new_vtkh_diy_utils_h

//#include <vtkmdiy/decomposition.hpp>
#include <vtkm/Bounds.h>

namespace vtkm
{
namespace rendering_new
{

static vtkm::Bounds DIYBoundsToVTKM(const vtkmdiy::DiscreteBounds& bounds)
{
  vtkm::Bounds vtkm_bounds;

  vtkm_bounds.X.Min = bounds.min[0];
  vtkm_bounds.Y.Min = bounds.min[1];
  vtkm_bounds.Z.Min = bounds.min[2];

  vtkm_bounds.X.Max = bounds.max[0];
  vtkm_bounds.Y.Max = bounds.max[1];
  vtkm_bounds.Z.Max = bounds.max[2];
  return vtkm_bounds;
}

static vtkmdiy::DiscreteBounds VTKMBoundsToDIY(const vtkm::Bounds& bounds)
{
  vtkmdiy::DiscreteBounds diy_bounds(3);

  diy_bounds.min[0] = bounds.X.Min;
  diy_bounds.min[1] = bounds.Y.Min;

  diy_bounds.max[0] = bounds.X.Max;
  diy_bounds.max[1] = bounds.Y.Max;

  if (bounds.Z.IsNonEmpty())
  {
    diy_bounds.min[2] = bounds.Z.Min;
    diy_bounds.max[2] = bounds.Z.Max;
  }
  else
  {
    diy_bounds.min[2] = 0;
    diy_bounds.max[2] = 0;
  }
  return diy_bounds;
}

}
} //namespace vtkm::rendering_new

#endif //vtkm_rendering_new_vtkh_diy_utils_h
