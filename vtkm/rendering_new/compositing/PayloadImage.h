//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_PayloadImage_h
#define vtkm_rendering_new_PayloadImage_h

#include <sstream>
#include <vector>
#include <vtkm/Bounds.h>

#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

struct VTKM_RENDERING_NEW_EXPORT PayloadImage
{
  // The image bounds are indicated by a grid starting at
  // 1-width and 1-height. Actual width would be calculated
  // this->Bounds.X.Max - this->Bounds.X.Min + 1
  // 1024 - 1 + 1 = 1024
  vtkm::Bounds OrigBounds;
  vtkm::Bounds Bounds;
  std::vector<unsigned char> Payloads;
  std::vector<float> Depths;
  int OrigRank;
  int PayloadBytes; // Size of the payload in bytes
  float DefaultValue;

  PayloadImage() {}

  PayloadImage(const vtkm::Bounds& bounds, const int payload_bytes)
    : OrigBounds(bounds)
    , Bounds(bounds)
    , OrigRank(-1)
    , PayloadBytes(payload_bytes)
  {
    this->DefaultValue = vtkm::Nan32();
    const int dx = bounds.X.Max - bounds.X.Min + 1;
    const int dy = bounds.Y.Max - bounds.Y.Min + 1;
    this->Payloads.resize(dx * dy * this->PayloadBytes);
    this->Depths.resize(dx * dy);
  }

  void InitOriginal(const PayloadImage& other)
  {
    this->OrigBounds = other.OrigBounds;
    this->Bounds = other.OrigBounds;
    this->PayloadBytes = other.PayloadBytes;
    this->DefaultValue = other.DefaultValue;

    const int dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;
    this->Payloads.resize(dx * dy * this->PayloadBytes);
    this->Depths.resize(dx * dy);

    this->OrigRank = -1;
  }

  int GetNumberOfPixels() const { return static_cast<int>(this->Depths.size()); }

  void Init(const unsigned char* payload_buffer, const float* depth_buffer, int width, int height)
  {
    this->Bounds.X.Min = 1;
    this->Bounds.Y.Min = 1;
    this->Bounds.X.Max = width;
    this->Bounds.Y.Max = height;
    this->OrigBounds = this->Bounds;
    const int size = width * height;
    this->Payloads.resize(size * this->PayloadBytes);
    this->Depths.resize(size);

    std::copy(payload_buffer, payload_buffer + size * this->PayloadBytes, &this->Payloads[0]);

    std::copy(depth_buffer, depth_buffer + size, &this->Depths[0]);
  }

  //
  // Fill this image with a sub-region of another image
  //
  void SubsetFrom(const PayloadImage& image, const vtkm::Bounds& sub_region)
  {
    this->OrigBounds = image.OrigBounds;
    this->Bounds = sub_region;
    this->OrigRank = image.OrigRank;
    this->PayloadBytes = image.PayloadBytes;

    assert(sub_region.X.Min >= image.Bounds.X.Min);
    assert(sub_region.Y.Min >= image.Bounds.Y.Min);
    assert(sub_region.X.Max <= image.Bounds.X.Max);
    assert(sub_region.Y.Max <= image.Bounds.Y.Max);

    const int s_dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int s_dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;

    const int dx = image.Bounds.X.Max - image.Bounds.X.Min + 1;
    //const int dy  = image.Bounds.Y.Max - image.Bounds.Y.Min + 1;

    const int start_x = this->Bounds.X.Min - image.Bounds.X.Min;
    const int start_y = this->Bounds.Y.Min - image.Bounds.Y.Min;
    const int end_y = start_y + s_dy;

    size_t buffer_size = s_dx * s_dy * this->PayloadBytes;

    this->Payloads.resize(buffer_size);
    this->Depths.resize(s_dx * s_dy);


#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = start_y; y < end_y; ++y)
    {
      const int copy_to = (y - start_y) * s_dx;
      const int copy_from = y * dx + start_x;

      std::copy(&image.Payloads[copy_from * this->PayloadBytes],
                &image.Payloads[copy_from * this->PayloadBytes] + s_dx * this->PayloadBytes,
                &this->Payloads[copy_to * this->PayloadBytes]);
      std::copy(&image.Depths[copy_from], &image.Depths[copy_from] + s_dx, &this->Depths[copy_to]);
    }
  }

  //
  // Fills the passed in image with the contents of this image
  //
  void SubsetTo(PayloadImage& image) const
  {
    assert(this->Bounds.X.Min >= image.Bounds.X.Min);
    assert(this->Bounds.Y.Min >= image.Bounds.Y.Min);
    assert(this->Bounds.X.Max <= image.Bounds.X.Max);
    assert(this->Bounds.Y.Max <= image.Bounds.Y.Max);

    const int s_dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int s_dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;

    const int dx = image.Bounds.X.Max - image.Bounds.X.Min + 1;
    //const int dy  = image.Bounds.Y.Max - image.Bounds.Y.Min + 1;

    const int start_x = this->Bounds.X.Min - image.Bounds.X.Min;
    const int start_y = this->Bounds.Y.Min - image.Bounds.Y.Min;

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = 0; y < s_dy; ++y)
    {
      const int copy_to = (y + start_y) * dx + start_x;
      const int copy_from = y * s_dx;

      std::copy(&this->Payloads[copy_from * this->PayloadBytes],
                &this->Payloads[copy_from * this->PayloadBytes] + s_dx * this->PayloadBytes,
                &image.Payloads[copy_to * this->PayloadBytes]);

      std::copy(&this->Depths[copy_from], &this->Depths[copy_from] + s_dx, &image.Depths[copy_to]);
    }
  }

  void Swap(PayloadImage& other)
  {
    vtkm::Bounds orig = this->OrigBounds;
    vtkm::Bounds bounds = this->Bounds;

    this->OrigBounds = other.OrigBounds;
    this->Bounds = other.Bounds;

    other.OrigBounds = orig;
    other.Bounds = bounds;

    this->Payloads.swap(other.Payloads);
    this->Depths.swap(other.Depths);
  }

  void Clear()
  {
    vtkm::Bounds empty;
    this->OrigBounds = empty;
    this->Bounds = empty;
    this->Payloads.clear();
    this->Depths.clear();
  }

  std::string ToString() const
  {
    std::stringstream ss;
    ss << "Total size pixels " << (int)this->Depths.size();
    ss << " tile dims: {" << this->Bounds.X.Min << "," << this->Bounds.Y.Min << "} - ";
    ss << "{" << this->Bounds.X.Max << "," << this->Bounds.Y.Max << "}\n";
    ;
    return ss.str();
  }

  void Save(const std::string& name, const std::vector<std::string>& comments);
};

}
} //namespace  vtkm::rendering_new

#endif //vtkm_rendering_new_PayloadImage_h
