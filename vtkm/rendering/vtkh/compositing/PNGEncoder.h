//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_PNGEncoder_h
#define vtk_m_rendering_compositing_PNGEncoder_h

#include <string>
#include <vector>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

class PNGEncoder
{
public:
  PNGEncoder();
  ~PNGEncoder();

  void Encode(const unsigned char* rgba_in, const int width, const int height);
  void Encode(const float* rgba_in, const int width, const int height);
  void Encode(const unsigned char* rgba_in,
              const int width,
              const int height,
              const std::vector<std::string>& comments);
  void Encode(const float* rgba_in,
              const int width,
              const int height,
              const std::vector<std::string>& comments);
  void Save(const std::string& filename);

  void* PngBuffer();
  size_t PngBufferSize();

  void Cleanup();

private:
  unsigned char* m_buffer;
  size_t m_buffer_size;
};


}
}
} //namespace vtkm::rendering::compositing

#endif //vtk_m_rendering_compositing_PNGEncoder_h
