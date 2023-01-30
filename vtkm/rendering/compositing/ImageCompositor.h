//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_ImageCompositor_h
#define vtk_m_rendering_compositing_ImageCompositor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <algorithm>
#include <vtkm/rendering/compositing/Image.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

class VTKM_RENDERING_EXPORT ImageCompositor
{
public:
  void Blend(vtkm::rendering::compositing::Image& front, vtkm::rendering::compositing::Image& back)
  {
    assert(front.Bounds.X.Min == back.Bounds.X.Min);
    assert(front.Bounds.Y.Min == back.Bounds.Y.Min);
    assert(front.Bounds.X.Max == back.Bounds.X.Max);
    assert(front.Bounds.Y.Max == back.Bounds.Y.Max);
    const int size = static_cast<int>(front.Pixels.size() / 4);

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      const int offset = i * 4;
      unsigned int alpha = front.Pixels[offset + 3];
      const unsigned int opacity = 255 - alpha;

      front.Pixels[offset + 0] +=
        static_cast<unsigned char>(opacity * back.Pixels[offset + 0] / 255);
      front.Pixels[offset + 1] +=
        static_cast<unsigned char>(opacity * back.Pixels[offset + 1] / 255);
      front.Pixels[offset + 2] +=
        static_cast<unsigned char>(opacity * back.Pixels[offset + 2] / 255);
      front.Pixels[offset + 3] +=
        static_cast<unsigned char>(opacity * back.Pixels[offset + 3] / 255);

      float d1 = std::min(front.Depths[i], 1.001f);
      float d2 = std::min(back.Depths[i], 1.001f);
      float depth = std::min(d1, d2);
      front.Depths[i] = depth;
    }
  }

  void ZBufferComposite(vtkm::rendering::compositing::Image& front,
                        const vtkm::rendering::compositing::Image& image)
  {
    assert(front.Depths.size() == front.Pixels.size() / 4);
    assert(front.Bounds.X.Min == image.Bounds.X.Min);
    assert(front.Bounds.Y.Min == image.Bounds.Y.Min);
    assert(front.Bounds.X.Max == image.Bounds.X.Max);
    assert(front.Bounds.Y.Max == image.Bounds.Y.Max);

    const int size = static_cast<int>(front.Depths.size());

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      const float depth = image.Depths[i];
      if (depth > 1.f || front.Depths[i] < depth)
      {
        continue;
      }
      const int offset = i * 4;
      front.Depths[i] = abs(depth);
      front.Pixels[offset + 0] = image.Pixels[offset + 0];
      front.Pixels[offset + 1] = image.Pixels[offset + 1];
      front.Pixels[offset + 2] = image.Pixels[offset + 2];
      front.Pixels[offset + 3] = image.Pixels[offset + 3];
    }
  }

  void OrderedComposite(std::vector<vtkm::rendering::compositing::Image>& images)
  {
    const int total_images = images.size();
    std::sort(images.begin(), images.end(), CompositeOrderSort());
    for (int i = 1; i < total_images; ++i)
    {
      Blend(images[0], images[i]);
    }
  }

  void ZBufferComposite(std::vector<vtkm::rendering::compositing::Image>& images)
  {
    const int total_images = images.size();
    for (int i = 1; i < total_images; ++i)
    {
      ZBufferComposite(images[0], images[i]);
    }
  }

  struct Pixel
  {
    unsigned char Color[4];
    float Depth;
    int PixelId; // local (sub-image) pixels id

    bool operator<(const Pixel& other) const
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
  };

  void CombineImages(const std::vector<vtkm::rendering::compositing::Image>& images,
                     std::vector<Pixel>& pixels)
  {

    const int num_images = static_cast<int>(images.size());
    for (int i = 0; i < num_images; ++i)
    {
      //
      //  Extract the partial composites into a contiguous array
      //

      const int image_size = images[i].GetNumberOfPixels();
      const int offset = i * image_size;
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
      for (int j = 0; j < image_size; ++j)
      {
        const int image_offset = j * 4;
        pixels[offset + j].Color[0] = images[i].Pixels[image_offset + 0];
        pixels[offset + j].Color[1] = images[i].Pixels[image_offset + 1];
        pixels[offset + j].Color[2] = images[i].Pixels[image_offset + 2];
        pixels[offset + j].Color[3] = images[i].Pixels[image_offset + 3];
        pixels[offset + j].Depth = images[i].Depths[j];
        pixels[offset + j].PixelId = j;
      } // for pixels
    }   // for images
  }

  void ZBufferBlend(std::vector<vtkm::rendering::compositing::Image>& images)
  {
    const int image_pixels = images[0].GetNumberOfPixels();
    const int num_images = static_cast<int>(images.size());
    std::vector<Pixel> pixels;
    CombineImages(images, pixels);
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < image_pixels; ++i)
    {
      const int begin = image_pixels * i;
      const int end = image_pixels * i - 1;
      std::sort(pixels.begin() + begin, pixels.begin() + end);
    }

    // check to see if that worked
    int pixel_id_0 = pixels[0].PixelId;
    for (int i = 1; i < num_images; ++i)
    {
      assert(pixel_id_0 == pixels[i].PixelId);
    }


#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < image_pixels; ++i)
    {
      const int index = i * num_images;
      Pixel pixel = pixels[index];
      for (int j = 1; j < num_images; ++j)
      {
        if (pixel.Color[3] == 255 || pixel.Depth > 1.f)
        {
          break;
        }
        unsigned int alpha = pixel.Color[3];
        const unsigned int opacity = 255 - alpha;
        pixel.Color[0] += static_cast<unsigned char>(opacity * pixels[index + j].Color[0] / 255);
        pixel.Color[1] += static_cast<unsigned char>(opacity * pixels[index + j].Color[1] / 255);
        pixel.Color[2] += static_cast<unsigned char>(opacity * pixels[index + j].Color[2] / 255);
        pixel.Color[3] += static_cast<unsigned char>(opacity * pixels[index + j].Color[3] / 255);
        pixel.Depth = pixels[index + j].Depth;
      } // for each image
      images[0].Pixels[i * 4 + 0] = pixel.Color[0];
      images[0].Pixels[i * 4 + 1] = pixel.Color[1];
      images[0].Pixels[i * 4 + 2] = pixel.Color[2];
      images[0].Pixels[i * 4 + 3] = pixel.Color[3];
      images[0].Depths[i] = pixel.Depth;
    } // for each pixel
  }
};

}
}
} //namespace vtkm::rendering::compositing

#endif //vtk_m_rendering_compositing_ImageComposititing_h
