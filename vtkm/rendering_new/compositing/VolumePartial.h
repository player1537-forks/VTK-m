//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_VolumePartial_h
#define vtkm_rendering_new_VolumePartial_h

#include <iostream>
#include <limits>
#include <vector>

namespace vtkm
{
namespace rendering_new
{

template <typename FloatType>
struct VolumePartial
{
  typedef FloatType ValueType;
  int PixelId;
  float Depth;
  float Pixel[3];
  float Alpha;

  VolumePartial()
    : PixelId(0)
    , Depth(0.f)
    , Alpha(0.f)
  {
    Pixel[0] = 0;
    Pixel[1] = 0;
    Pixel[2] = 0;
  }

  void print() const
  {
    std::cout << "[id : " << this->PixelId << ", red : " << this->Pixel[0] << ","
              << " green : " << this->Pixel[1] << ", blue : " << this->Pixel[2] << ", alpha "
              << this->Alpha << ", depth : " << this->Depth << "]\n";
  }

  bool operator<(const VolumePartial& other) const
  {
    if (this->PixelId != other.PixelId)
    {
      return this->PixelId < other.PixelId;
    }
    else
    {
      return this->Depth < other.Depth;
    }
  }

  inline void blend(const VolumePartial& other)
  {
    if (this->Alpha >= 1.f || other.Alpha == 0.f)
      return;
    const float opacity = (1.f - this->Alpha);
    this->Pixel[0] += opacity * other.Pixel[0];
    this->Pixel[1] += opacity * other.Pixel[1];
    this->Pixel[2] += opacity * other.Pixel[2];
    this->Alpha += opacity * other.Alpha;
    this->Alpha = this->Alpha > 1.f ? 1.f : this->Alpha;
  }

  static void composite_background(std::vector<VolumePartial>& partials,
                                   const std::vector<FloatType>& background)
  {
    VolumePartial bg_color;
    bg_color.Pixel[0] = static_cast<float>(background[0]);
    bg_color.Pixel[1] = static_cast<float>(background[1]);
    bg_color.Pixel[2] = static_cast<float>(background[2]);
    bg_color.Alpha = static_cast<float>(background[3]);
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
} //vtkm::rendering_new

#endif //vtkm_rendering_new_VolumePartial_h
