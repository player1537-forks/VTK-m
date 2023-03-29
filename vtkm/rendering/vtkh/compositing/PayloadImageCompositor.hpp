//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTKH_DIY_PAYLOAD_IMAGE_COMPOSITOR_HPP
#define VTKH_DIY_PAYLOAD_IMAGE_COMPOSITOR_HPP

#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>
#include <algorithm>
#include <cmath>

#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh
{

class VTKM_RENDERING_EXPORT PayloadImageCompositor
{
public:

void ZBufferComposite(vtkh::PayloadImage &front, const vtkh::PayloadImage &image)
{
  if(front.m_payload_bytes != image.m_payload_bytes)
  {
    std::cout<<"very bad\n";
  }
  assert(front.m_depths.size() == front.m_payloads.size() / front.m_payload_bytes);
  assert(front.m_bounds.X.Min == image.m_bounds.X.Min);
  assert(front.m_bounds.Y.Min == image.m_bounds.Y.Min);
  assert(front.m_bounds.X.Max == image.m_bounds.X.Max);
  assert(front.m_bounds.Y.Max == image.m_bounds.Y.Max);

  const int size = static_cast<int>(front.m_depths.size());
  const bool nan_check = image.m_default_value != image.m_default_value;
#ifdef VTKH_OPENMP_ENABLED
  #pragma omp parallel for
#endif
  for(int i = 0; i < size; ++i)
  {
    const float depth = image.m_depths[i];
    const float fdepth = front.m_depths[i];
    // this should handle NaNs correctly
    const bool take_back = fmin(depth, fdepth) == depth;

    if(take_back)
    {
      const int offset = i * 4;
      front.m_depths[i] = depth;
      const size_t p_offset = i * front.m_payload_bytes;
      std::copy(&image.m_payloads[p_offset],
                &image.m_payloads[p_offset] + front.m_payload_bytes,
                &front.m_payloads[p_offset]);
    }
  }
}


};

} // namespace vtkh
#endif
