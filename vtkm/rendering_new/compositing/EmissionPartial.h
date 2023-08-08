//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_EmissionPartial_h
#define vtkm_rendering_new_EmissionPartial_h

#include <assert.h>

namespace vtkm
{
namespace rendering_new
{

template <typename FloatType>
struct EmissionPartial
{
  typedef FloatType ValueType;

  int PixelId;
  double Depth;
  std::vector<FloatType> Bins;
  std::vector<FloatType> EmissionBins;

  EmissionPartial()
    : PixelId(0)
    , Depth(0.f)
  {
  }

  void alter_bin(int bin, FloatType value)
  {
    this->Bins[bin] = value;
    this->EmissionBins[bin] = value;
  }

  void print()
  {
    std::cout << "Partial id " << this->PixelId << "\n";
    std::cout << "Absorption : ";
    for (int i = 0; i < this->Bins.size(); ++i)
    {
      std::cout << this->Bins[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Emission: ";
    for (int i = 0; i < this->Bins.size(); ++i)
    {
      std::cout << this->EmissionBins[i] << " ";
    }
    std::cout << "\n";
  }

  bool operator<(const EmissionPartial<FloatType>& other) const
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

  inline void blend_absorption(const EmissionPartial<FloatType>& other)
  {
    const int num_bins = static_cast<int>(this->Bins.size());
    assert(num_bins == (int)other.Bins.size());
    for (int i = 0; i < num_bins; ++i)
    {
      this->Bins[i] *= other.Bins[i];
    }
  }

  inline void blend_emission(EmissionPartial<FloatType>& other)
  {
    const int num_bins = static_cast<int>(this->Bins.size());
    assert(num_bins == (int)other.Bins.size());
    for (int i = 0; i < num_bins; ++i)
    {
      this->EmissionBins[i] *= other.Bins[i];
    }
  }

  inline void add_emission(EmissionPartial<FloatType>& other)
  {
    const int num_bins = static_cast<int>(this->Bins.size());
    assert(num_bins == (int)other.Bins.size());
    for (int i = 0; i < num_bins; ++i)
    {
      this->EmissionBins[i] += other.EmissionBins[i];
    }
  }

  static void composite_background(std::vector<EmissionPartial>& /*partials*/,
                                   const std::vector<FloatType>& /*background*/)
  {
    //for(
  }
};

}
} // vtkm::rendering_new


#endif //vtkm_rendering_new_EmissionPartial_h
