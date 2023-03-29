//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/vtkh/compositing/Image.hpp>
//#include <png_utils/ascent_png_encoder.hpp>
#include <vtkm/rendering/vtkh/compositing/PNGEncoder.h>

namespace vtkh
{

void Image::Save(const std::string &name,
                 const std::vector<std::string> &comments)
{
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&m_pixels[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

void Image::Save(const std::string &name,
                 const std::vector<std::string> &comments) const
{
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&m_pixels[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}
} // namespace vtkh
