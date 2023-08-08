//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_Compositor_h
#define vtkm_rendering_new_Compositor_h

#include <sstream>
#include <vtkm/rendering_new/compositing/Image.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT Compositor
{
public:
  enum CompositeMode
  {
    Z_BUFFER_SURFACE, // zbuffer composite no transparency
    Z_BUFFER_BLEND,   // zbuffer composite with transparency
    VIS_ORDER_BLEND   // blend images in a specific order
  };
  Compositor();

  virtual ~Compositor();

  void SetCompositeMode(enum CompositeMode composite_mode);

  void ClearImages();

  void AddImage(const unsigned char* color_buffer,
                const float* depth_buffer,
                const int width,
                const int height);

  void AddImage(const float* color_buffer,
                const float* depth_buffer,
                const int width,
                const int height);

  void AddImage(const unsigned char* color_buffer,
                const float* depth_buffer,
                const int width,
                const int height,
                const int vis_order);

  void AddImage(const float* color_buffer,
                const float* depth_buffer,
                const int width,
                const int height,
                const int vis_order);

  Image Composite();

  virtual void Cleanup();

  std::string GetLogString();

  unsigned char* ConvertBuffer(const float* buffer, const int size)
  {
    unsigned char* ubytes = new unsigned char[size];

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i)
    {
      ubytes[i] = static_cast<unsigned char>(buffer[i] * 255.f);
    }

    return ubytes;
  }

protected:
  virtual void CompositeZBufferSurface();
  virtual void CompositeZBufferBlend();
  virtual void CompositeVisOrder();

  std::stringstream m_log_stream;
  CompositeMode CompositeMode;
  std::vector<Image> Images;
};

}
} //vtkm::rendering_new

#endif //vtkm_rendering_new_Compositor_h
