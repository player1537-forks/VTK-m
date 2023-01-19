#ifndef VTKH_PAYLOAD_COMPOSITOR_HPP
#define VTKH_PAYLOAD_COMPOSITOR_HPP

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/compositing/PayloadImage.h>

namespace vtkh
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

};

#endif
