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
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&this->Pixels[0],
                 this->Bounds.X.Max - this->Bounds.X.Min + 1,
                 this->Bounds.Y.Max - this->Bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

void Image::Save(const std::string& name, const std::vector<std::string>& comments) const
{
  vtkm::rendering::compositing::PNGEncoder encoder;
  encoder.Encode(&this->Pixels[0],
                 this->Bounds.X.Max - this->Bounds.X.Min + 1,
                 this->Bounds.Y.Max - this->Bounds.Y.Min + 1,
                 comments);
  encoder.Save(name);
}

}
}
} //namespace vtkm::rendering::compositing
