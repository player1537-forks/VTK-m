//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering_new/compositing/Image.h>
#include <vtkm/rendering_new/compositing/PNGEncoder.h>

namespace vtkm
{
namespace rendering_new
{

void Image::Save(const std::string& name, const std::vector<std::string>& comments) const
{

  vtkm::rendering_new::PNGEncoder encoder;
  encoder.Encode(&this->Pixels[0],
                 this->Bounds.X.Max - this->Bounds.X.Min + 1,
                 this->Bounds.Y.Max - this->Bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

}
} // vtkm::rendering_new
