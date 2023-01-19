// See License.txt

#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/compositing/PNGEncoder.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

void Image::Save(const std::string& name, const std::vector<std::string>& comments)
{
  PNGEncoder encoder;
  encoder.Encode(&m_pixels[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

void Image::Save(const std::string& name, const std::vector<std::string>& comments) const
{
  PNGEncoder encoder;
  encoder.Encode(&m_pixels[0],
                 m_bounds.X.Max - m_bounds.X.Min + 1,
                 m_bounds.Y.Max - m_bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

}
}
} //namespace vtkm::rendering::compositing
