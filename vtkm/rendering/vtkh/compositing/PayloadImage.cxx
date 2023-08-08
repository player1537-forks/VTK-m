//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/vtkh/compositing/PNGEncoder.h>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.h>

namespace vtkh
{

void PayloadImage::Save(const std::string& name, const std::vector<std::string>& comments)
{
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&this->Payloads[0],
                 this->Bounds.X.Max - this->Bounds.X.Min + 1,
                 this->Bounds.Y.Max - this->Bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

} // namespace vtkh