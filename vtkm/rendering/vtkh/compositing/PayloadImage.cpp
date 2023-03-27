// See License.txt

#include "PayloadImage.hpp"
#include <vtkm/rendering/vtkh/compositing/PNGEncoder.h>

namespace vtkh
{

void PayloadImage::Save(const std::string &name,
                 const std::vector<std::string> &comments)
{
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&m_payloads[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

} // namespace vtkh
