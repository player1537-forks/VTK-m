//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_PayloadImageCompositor_h
#define vtkm_rendering_new_PayloadImageCompositor_h

#include <algorithm>
#include <cmath>
#include <vtkm/rendering_new/compositing/PayloadImage.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT PayloadImageCompositor
{
public:
  void ZBufferComposite(vtkm::rendering_new::PayloadImage& front,
                        const vtkm::rendering_new::PayloadImage& image)
  {
    if (front.PayloadBytes != image.PayloadBytes)
    {
      std::cout << "very bad\n";
    }
    assert(front.Depths.size() == front.Payloads.size() / front.PayloadBytes);
    assert(front.Bounds.X.Min == image.Bounds.X.Min);
    assert(front.Bounds.Y.Min == image.Bounds.Y.Min);
    assert(front.Bounds.X.Max == image.Bounds.X.Max);
    assert(front.Bounds.Y.Max == image.Bounds.Y.Max);

    const int size = static_cast<int>(front.Depths.size());
    //const bool nan_check = image.m_default_value != image.m_default_value;
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      const float depth = image.Depths[i];
      const float fdepth = front.Depths[i];
      // this should handle NaNs correctly
      const bool take_back = fmin(depth, fdepth) == depth;

      if (take_back)
      {
        //const int offset = i * 4;
        front.Depths[i] = depth;
        const size_t p_offset = i * front.PayloadBytes;
        std::copy(&image.Payloads[p_offset],
                  &image.Payloads[p_offset] + front.PayloadBytes,
                  &front.Payloads[p_offset]);
      }
    }
  }
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_new_PayloadImageCompositor_h
