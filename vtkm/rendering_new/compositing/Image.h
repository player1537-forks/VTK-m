//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_Image_h
#define vtkm_rendering_new_Image_h

#include <vtkm/Bounds.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

#include <sstream>
#include <vector>

namespace vtkm
{
namespace rendering_new
{

struct VTKM_RENDERING_NEW_EXPORT Image
{
  // The image bounds are indicated by a grid starting at
  // 1-width and 1-height. Actual width would be calculated
  // Bounds.X.Max - Bounds.X.Min + 1
  // 1024 - 1 + 1 = 1024
  vtkm::Bounds OrigBounds;
  vtkm::Bounds Bounds;
  std::vector<unsigned char> Pixels;
  std::vector<float> Depths;
  int OrigRank;
  bool HasTransparency;
  int CompositeOrder;

  Image()
    : OrigRank(-1)
    , HasTransparency(false)
    , CompositeOrder(-1)
  {
  }


  Image(const vtkm::Bounds& bounds)
    : OrigBounds(bounds)
    , Bounds(bounds)
    , OrigRank(-1)
    , HasTransparency(false)
    , CompositeOrder(-1)

  {
    const int dx = bounds.X.Max - bounds.X.Min + 1;
    const int dy = bounds.Y.Max - bounds.Y.Min + 1;
    this->Pixels.resize(dx * dy * 4);
    this->Depths.resize(dx * dy);
  }

  // init this image based on the original bounds
  // of the other image
  void InitOriginal(const Image& other)
  {
    this->OrigBounds = other.OrigBounds;
    this->Bounds = other.OrigBounds;

    const int dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;
    this->Pixels.resize(dx * dy * 4);
    this->Depths.resize(dx * dy);

    this->OrigRank = -1;
    this->HasTransparency = false;
    this->CompositeOrder = -1;
  }

  int GetNumberOfPixels() const { return static_cast<int>(this->Pixels.size() / 4); }

  void SetHasTransparency(bool has_transparency) { this->HasTransparency = has_transparency; }

  bool GetHasTransparency() { return this->HasTransparency; }

  void Init(const float* color_buffer,
            const float* depth_buffer,
            int width,
            int height,
            int composite_order = -1)
  {
    this->CompositeOrder = composite_order;
    this->Bounds.X.Min = 1;
    this->Bounds.Y.Min = 1;
    this->Bounds.X.Max = width;
    this->Bounds.Y.Max = height;
    this->OrigBounds = this->Bounds;
    const int size = width * height;
    this->Pixels.resize(size * 4);
    this->Depths.resize(size);

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      const int offset = i * 4;
      this->Pixels[offset + 0] = static_cast<unsigned char>(color_buffer[offset + 0] * 255.f);
      this->Pixels[offset + 1] = static_cast<unsigned char>(color_buffer[offset + 1] * 255.f);
      this->Pixels[offset + 2] = static_cast<unsigned char>(color_buffer[offset + 2] * 255.f);
      this->Pixels[offset + 3] = static_cast<unsigned char>(color_buffer[offset + 3] * 255.f);
      float depth = depth_buffer[i];
      //make sure we can do a single comparison on depth
      //deal with negative depth values
      //TODO: This may not be the best way
      depth = depth < 0 ? abs(depth) : depth;
      this->Depths[i] = depth;
    }
  }

  void Init(const unsigned char* color_buffer,
            const float* depth_buffer,
            int width,
            int height,
            int composite_order = -1)
  {
    this->CompositeOrder = composite_order;
    this->Bounds.X.Min = 1;
    this->Bounds.Y.Min = 1;
    this->Bounds.X.Max = width;
    this->Bounds.Y.Max = height;
    this->OrigBounds = this->Bounds;

    const int size = width * height;
    this->Pixels.resize(size * 4);
    this->Depths.resize(size);

    std::copy(color_buffer, color_buffer + size * 4, &this->Pixels[0]);

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      float depth = depth_buffer[i];
      //make sure we can do a single comparison on depth
      depth = depth < 0 ? 2.f : depth;
      this->Depths[i] = depth;
    } // for
  }


  void CompositeBackground(const float* color)
  {

    const int size = static_cast<int>(this->Pixels.size() / 4);
    unsigned char bg_color[4];
    for (int i = 0; i < 4; ++i)
    {
      bg_color[i] = static_cast<unsigned char>(color[i] * 255.f);
    }

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      const int offset = i * 4;
      unsigned int alpha = static_cast<unsigned int>(this->Pixels[offset + 3]);
      const float opacity = (255 - alpha);
      this->Pixels[offset + 0] += static_cast<unsigned char>(opacity * bg_color[0] / 255);
      this->Pixels[offset + 1] += static_cast<unsigned char>(opacity * bg_color[1] / 255);
      this->Pixels[offset + 2] += static_cast<unsigned char>(opacity * bg_color[2] / 255);
      this->Pixels[offset + 3] += static_cast<unsigned char>(opacity * bg_color[3] / 255);
    }
  }
  //
  // Fill this image with a sub-region of another image
  //
  void SubsetFrom(const Image& image, const vtkm::Bounds& sub_region)
  {
    this->OrigBounds = image.OrigBounds;
    this->Bounds = sub_region;
    this->OrigRank = image.OrigRank;
    this->CompositeOrder = image.CompositeOrder;

    assert(sub_region.X.Min >= image.Bounds.X.Min);
    assert(sub_region.Y.Min >= image.Bounds.Y.Min);
    assert(sub_region.X.Max <= image.Bounds.X.Max);
    assert(sub_region.Y.Max <= image.Bounds.Y.Max);

    const int s_dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int s_dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;

    const int dx = image.Bounds.X.Max - image.Bounds.X.Min + 1;
    //const int dy  = image.this->Bounds.Y.Max - image.this->Bounds.Y.Min + 1;

    const int start_x = this->Bounds.X.Min - image.Bounds.X.Min;
    const int start_y = this->Bounds.Y.Min - image.Bounds.Y.Min;
    const int end_y = start_y + s_dy;

    this->Pixels.resize(s_dx * s_dy * 4);
    this->Depths.resize(s_dx * s_dy);



#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = start_y; y < end_y; ++y)
    {
      const int copy_to = (y - start_y) * s_dx;
      const int copy_from = y * dx + start_x;

      std::copy(&image.Pixels[copy_from * 4],
                &image.Pixels[copy_from * 4] + s_dx * 4,
                &this->Pixels[copy_to * 4]);
      std::copy(&image.Depths[copy_from], &image.Depths[copy_from] + s_dx, &this->Depths[copy_to]);
    }
  }

  void Color(int color)
  {
    unsigned char c[4];
    c[3] = 255;

    c[0] = 0;
    c[1] = 0;
    c[2] = 0;
    int index = color % 3;
    c[index] = 255 - color * 11;
    ;
    const int size = static_cast<int>(this->Pixels.size());
    for (int i = 0; i < size; ++i)
    {
      float d = this->Depths[i / 4];
      if (d > 0 && d < 1)
      {
        this->Pixels[i] = c[i % 4];
      }
      else
      {
        this->Pixels[i] = 155;
      }
    }
  }
  //
  // Fills the passed in image with the contents of this image
  //
  void SubsetTo(Image& image) const
  {
    image.CompositeOrder = this->CompositeOrder;
    assert(this->Bounds.X.Min >= image.Bounds.X.Min);
    assert(this->Bounds.Y.Min >= image.Bounds.Y.Min);
    assert(this->Bounds.X.Max <= image.Bounds.X.Max);
    assert(this->Bounds.Y.Max <= image.Bounds.Y.Max);

    const int s_dx = this->Bounds.X.Max - this->Bounds.X.Min + 1;
    const int s_dy = this->Bounds.Y.Max - this->Bounds.Y.Min + 1;

    const int dx = image.Bounds.X.Max - image.Bounds.X.Min + 1;
    //const int dy  = image.this->Bounds.Y.Max - image.this->Bounds.Y.Min + 1;

    const int start_x = this->Bounds.X.Min - image.Bounds.X.Min;
    const int start_y = this->Bounds.Y.Min - image.Bounds.Y.Min;

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = 0; y < s_dy; ++y)
    {
      const int copy_to = (y + start_y) * dx + start_x;
      const int copy_from = y * s_dx;

      std::copy(&this->Pixels[copy_from * 4],
                &this->Pixels[copy_from * 4] + s_dx * 4,
                &image.Pixels[copy_to * 4]);

      std::copy(&this->Depths[copy_from], &this->Depths[copy_from] + s_dx, &image.Depths[copy_to]);
    }
  }

  void Swap(Image& other)
  {
    vtkm::Bounds orig = this->OrigBounds;
    vtkm::Bounds bounds = this->Bounds;

    this->OrigBounds = other.OrigBounds;
    this->Bounds = other.Bounds;

    other.OrigBounds = orig;
    other.Bounds = bounds;

    this->Pixels.swap(other.Pixels);
    this->Depths.swap(other.Depths);
  }

  void Clear()
  {
    vtkm::Bounds empty;
    this->OrigBounds = empty;
    this->Bounds = empty;
    this->Pixels.clear();
    this->Depths.clear();
  }

  std::string ToString() const
  {
    std::stringstream ss;
    ss << "Total size pixels " << (int)this->Pixels.size() / 4;
    ss << " tile dims: {" << this->Bounds.X.Min << "," << this->Bounds.Y.Min << "} - ";
    ss << "{" << this->Bounds.X.Max << "," << this->Bounds.Y.Max << "}\n";
    ;
    return ss.str();
  }

  void Save(const std::string& name, const std::vector<std::string>& comments) const;
};

struct CompositeOrderSort
{
  inline bool operator()(const Image& lhs, const Image& rhs) const
  {
    return lhs.CompositeOrder < rhs.CompositeOrder;
  }
};

}
} //namespace  vtkm::rendering_new


#endif //vtkm_rendering_new_Image_h
