//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_PayloadCompositor_h
#define vtk_m_rendering_compositing_PayloadCompositor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/compositing/PayloadImage.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

class VTKM_RENDERING_EXPORT PayloadCompositor
{
public:
  PayloadCompositor();

  void ClearImages();

  void AddImage(vtkm::rendering::compositing::PayloadImage& image);

  vtkm::rendering::compositing::PayloadImage Composite();

protected:
  std::vector<vtkm::rendering::compositing::PayloadImage> m_images;
};

}
}
} // namespace vtkm:rendering::compositing


#endif //vtk_m_rendering_compositing_PayloadCompositor_h
