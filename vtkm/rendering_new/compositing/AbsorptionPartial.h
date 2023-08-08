//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_AbsorptionPartial_h
#define vtkm_rendering_new_AbsorptionPartial_h

#include <assert.h>

namespace vtkm
{
namespace rendering_new
{

template <typename FloatType>
struct AbsorptionPartial
{
  typedef FloatType ValueType;
  int PixelId;
  double Depth;
  std::vector<FloatType> Bins;

  AbsorptionPartial()
    : PixelId(0)
    , Depth(0.f)
  {
  }

  void print() {}

  bool operator<(const AbsorptionPartial<FloatType>& other) const
  {
    //
    // In absorption only we can blend the same
    // pixel ids in any order
    //
    return this->PixelId < other.PixelId;
  }

  inline void blend(const AbsorptionPartial<FloatType>& other)
  {
    const int num_bins = static_cast<int>(this->Bins.size());
    assert(num_bins == (int)other.Bins.size());
    for (int i = 0; i < num_bins; ++i)
    {
      this->Bins[i] *= other.Bins[i];
    }
  }

  static void composite_background(std::vector<AbsorptionPartial>& partials,
                                   const std::vector<FloatType>& background)
  {
    const int size = static_cast<int>(partials.size());
    AbsorptionPartial<FloatType> bg;
    bg.Bins = background;
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      partials[i].blend(bg);
    }
  }
};

}
} // vtkm::rendering_new


#endif //vtkm_rendering_new_AbsorptionPartial_h
