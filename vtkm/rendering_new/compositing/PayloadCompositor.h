//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_PayloadCompositor_h
#define vtkm_rendering_new_PayloadCompositor_h

#include <sstream>
#include <vtkm/rendering_new/compositing/PayloadImage.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

namespace vtkm
{
namespace rendering_new
{

class VTKM_RENDERING_NEW_EXPORT PayloadCompositor
{
public:
  PayloadCompositor();

  void ClearImages();

  void AddImage(PayloadImage& image);

  PayloadImage Composite();

protected:
  std::vector<PayloadImage> m_images;
};

}
} //vtkm::rendering_new

#endif //vtkm_rendering_new_PayloadCompoistor_h
