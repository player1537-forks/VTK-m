//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/compositing/PNGEncoder.h>
#include <vtkm/rendering/compositing/PayloadImage.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

void PayloadImage::Save(const std::string& name, const std::vector<std::string>& comments)
{
  PNGEncoder encoder;
  encoder.Encode(&m_payloads[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

}
}
} //namespace vtkm::rendering::compositing
