//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_VolumePartial_h
#define vtkm_rendering_compositing_VolumePartial_h

#include <limits>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

template <typename FloatType>
struct VolumePartial
{
  typedef FloatType ValueType;
  int m_pixel_id;
  float m_depth;
  float m_pixel[3];
  float m_alpha;

  VolumePartial()
    : m_pixel_id(0)
    , m_depth(0.f)
    , m_alpha(0.f)
  {
    m_pixel[0] = 0;
    m_pixel[1] = 0;
    m_pixel[2] = 0;
  }

  void print() const
  {
    std::cout << "[id : " << m_pixel_id << ", red : " << m_pixel[0] << ","
              << " green : " << m_pixel[1] << ", blue : " << m_pixel[2] << ", alpha " << m_alpha
              << ", depth : " << m_depth << "]\n";
  }

  bool operator<(const VolumePartial& other) const
  {
    if (m_pixel_id != other.m_pixel_id)
    {
      return m_pixel_id < other.m_pixel_id;
    }
    else
    {
      return m_depth < other.m_depth;
    }
  }

  inline void blend(const VolumePartial& other)
  {
    if (m_alpha >= 1.f || other.m_alpha == 0.f)
      return;
    const float opacity = (1.f - m_alpha);
    m_pixel[0] += opacity * other.m_pixel[0];
    m_pixel[1] += opacity * other.m_pixel[1];
    m_pixel[2] += opacity * other.m_pixel[2];
    m_alpha += opacity * other.m_alpha;
    m_alpha = m_alpha > 1.f ? 1.f : m_alpha;
  }

  static void composite_background(std::vector<VolumePartial>& partials,
                                   const std::vector<FloatType>& background)
  {
    VolumePartial bg_color;
    bg_color.m_pixel[0] = static_cast<float>(background[0]);
    bg_color.m_pixel[1] = static_cast<float>(background[1]);
    bg_color.m_pixel[2] = static_cast<float>(background[2]);
    bg_color.m_alpha = static_cast<float>(background[3]);
    //
    // Gather the unique pixels into the output
    //
    const int total_pixels = static_cast<int>(partials.size());
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < total_pixels; ++i)
    {
      partials[i].blend(bg_color);
    }
  }
};

}
}
} //vtkm::rendering::compositing

#endif //vtkm_rendering_compositing_VolumePartial_h
