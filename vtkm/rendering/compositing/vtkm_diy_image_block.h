#ifndef VTKH_DIY_IMAGE_BLOCK_HPP
#define VTKH_DIY_IMAGE_BLOCK_HPP

#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/compositing/PayloadImage.h>
#include <vtkmdiy/master.hpp>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

template <typename ImageType>
struct ImageBlock
{
  ImageType& m_image;
  ImageBlock(ImageType& image)
    : m_image(image)
  {
  }
};

struct MultiImageBlock
{
  std::vector<vtkm::rendering::compositing::Image>& m_images;
  vtkm::rendering::compositing::Image& m_output;
  MultiImageBlock(std::vector<vtkm::rendering::compositing::Image>& images,
                  vtkm::rendering::compositing::Image& output)
    : m_images(images)
    , m_output(output)
  {
  }
};

template <typename ImageType>
struct AddImageBlock
{
  ImageType& m_image;
  const vtkmdiy::Master& m_master;

  AddImageBlock(vtkmdiy::Master& master, ImageType& image)
    : m_image(image)
    , m_master(master)
  {
  }
  template <typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType&, // local_bounds
                  const BoundsType&, // local_with_ghost_bounds
                  const BoundsType&, // domain_bounds
                  const LinkType& link) const
  {
    ImageBlock<ImageType>* block = new ImageBlock<ImageType>(m_image);
    LinkType* linked = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(m_master);
    master.add(gid, block, linked);
  }
};

struct AddMultiImageBlock
{
  const vtkmdiy::Master& m_master;
  std::vector<vtkm::rendering::compositing::Image>& m_images;
  vtkm::rendering::compositing::Image& m_output;

  AddMultiImageBlock(vtkmdiy::Master& master,
                     std::vector<vtkm::rendering::compositing::Image>& images,
                     vtkm::rendering::compositing::Image& output)
    : m_master(master)
    , m_images(images)
    , m_output(output)
  {
  }
  template <typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType&, // local_bounds
                  const BoundsType&, // local_with_ghost_bounds
                  const BoundsType&, // domain_bounds
                  const LinkType& link) const
  {
    MultiImageBlock* block = new MultiImageBlock(m_images, m_output);
    LinkType* linked = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(m_master);
    int lid = master.add(gid, block, linked);
  }
};

}
}
} //namespace vtkm::rendering::compositing

namespace vtkmdiy
{

template <>
struct Serialization<vtkm::rendering::compositing::PayloadImage>
{
  static void save(BinaryBuffer& bb, const vtkm::rendering::compositing::PayloadImage& image)
  {
    vtkmdiy::save(bb, image.OrigBounds.X.Min);
    vtkmdiy::save(bb, image.OrigBounds.Y.Min);
    vtkmdiy::save(bb, image.OrigBounds.Z.Min);
    vtkmdiy::save(bb, image.OrigBounds.X.Max);
    vtkmdiy::save(bb, image.OrigBounds.Y.Max);
    vtkmdiy::save(bb, image.OrigBounds.Z.Max);

    vtkmdiy::save(bb, image.Bounds.X.Min);
    vtkmdiy::save(bb, image.Bounds.Y.Min);
    vtkmdiy::save(bb, image.Bounds.Z.Min);
    vtkmdiy::save(bb, image.Bounds.X.Max);
    vtkmdiy::save(bb, image.Bounds.Y.Max);
    vtkmdiy::save(bb, image.Bounds.Z.Max);

    vtkmdiy::save(bb, image.Payloads);
    vtkmdiy::save(bb, image.PayloadBytes);
    vtkmdiy::save(bb, image.Depths);
    vtkmdiy::save(bb, image.OrigRank);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering::compositing::PayloadImage& image)
  {
    vtkmdiy::load(bb, image.OrigBounds.X.Min);
    vtkmdiy::load(bb, image.OrigBounds.Y.Min);
    vtkmdiy::load(bb, image.OrigBounds.Z.Min);
    vtkmdiy::load(bb, image.OrigBounds.X.Max);
    vtkmdiy::load(bb, image.OrigBounds.Y.Max);
    vtkmdiy::load(bb, image.OrigBounds.Z.Max);

    vtkmdiy::load(bb, image.Bounds.X.Min);
    vtkmdiy::load(bb, image.Bounds.Y.Min);
    vtkmdiy::load(bb, image.Bounds.Z.Min);
    vtkmdiy::load(bb, image.Bounds.X.Max);
    vtkmdiy::load(bb, image.Bounds.Y.Max);
    vtkmdiy::load(bb, image.Bounds.Z.Max);

    vtkmdiy::load(bb, image.Payloads);
    vtkmdiy::load(bb, image.PayloadBytes);
    vtkmdiy::load(bb, image.Depths);
    vtkmdiy::load(bb, image.OrigRank);
  }
};

template <>
struct Serialization<vtkm::rendering::compositing::Image>
{
  static void save(BinaryBuffer& bb, const vtkm::rendering::compositing::Image& image)
  {
    vtkmdiy::save(bb, image.OrigBounds.X.Min);
    vtkmdiy::save(bb, image.OrigBounds.Y.Min);
    vtkmdiy::save(bb, image.OrigBounds.Z.Min);
    vtkmdiy::save(bb, image.OrigBounds.X.Max);
    vtkmdiy::save(bb, image.OrigBounds.Y.Max);
    vtkmdiy::save(bb, image.OrigBounds.Z.Max);

    vtkmdiy::save(bb, image.Bounds.X.Min);
    vtkmdiy::save(bb, image.Bounds.Y.Min);
    vtkmdiy::save(bb, image.Bounds.Z.Min);
    vtkmdiy::save(bb, image.Bounds.X.Max);
    vtkmdiy::save(bb, image.Bounds.Y.Max);
    vtkmdiy::save(bb, image.Bounds.Z.Max);

    vtkmdiy::save(bb, image.Pixels);
    vtkmdiy::save(bb, image.Depths);
    vtkmdiy::save(bb, image.OrigRank);
    vtkmdiy::save(bb, image.CompositeOrder);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering::compositing::Image& image)
  {
    vtkmdiy::load(bb, image.OrigBounds.X.Min);
    vtkmdiy::load(bb, image.OrigBounds.Y.Min);
    vtkmdiy::load(bb, image.OrigBounds.Z.Min);
    vtkmdiy::load(bb, image.OrigBounds.X.Max);
    vtkmdiy::load(bb, image.OrigBounds.Y.Max);
    vtkmdiy::load(bb, image.OrigBounds.Z.Max);

    vtkmdiy::load(bb, image.Bounds.X.Min);
    vtkmdiy::load(bb, image.Bounds.Y.Min);
    vtkmdiy::load(bb, image.Bounds.Z.Min);
    vtkmdiy::load(bb, image.Bounds.X.Max);
    vtkmdiy::load(bb, image.Bounds.Y.Max);
    vtkmdiy::load(bb, image.Bounds.Z.Max);

    vtkmdiy::load(bb, image.Pixels);
    vtkmdiy::load(bb, image.Depths);
    vtkmdiy::load(bb, image.OrigRank);
    vtkmdiy::load(bb, image.CompositeOrder);
  }
};
} //namespace vtkmdiy

#endif
