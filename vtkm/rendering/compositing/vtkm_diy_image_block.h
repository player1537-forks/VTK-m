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
  std::vector<vtkm::rendering::compositing::Image>& m_images;
  vtkm::rendering::compositing::Image& m_output;
  const vtkmdiy::Master& m_master;

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
    vtkmdiy::save(bb, image.m_orig_bounds.X.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.Y.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.Z.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.X.Max);
    vtkmdiy::save(bb, image.m_orig_bounds.Y.Max);
    vtkmdiy::save(bb, image.m_orig_bounds.Z.Max);

    vtkmdiy::save(bb, image.m_bounds.X.Min);
    vtkmdiy::save(bb, image.m_bounds.Y.Min);
    vtkmdiy::save(bb, image.m_bounds.Z.Min);
    vtkmdiy::save(bb, image.m_bounds.X.Max);
    vtkmdiy::save(bb, image.m_bounds.Y.Max);
    vtkmdiy::save(bb, image.m_bounds.Z.Max);

    vtkmdiy::save(bb, image.m_payloads);
    vtkmdiy::save(bb, image.m_payload_bytes);
    vtkmdiy::save(bb, image.m_depths);
    vtkmdiy::save(bb, image.m_orig_rank);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering::compositing::PayloadImage& image)
  {
    vtkmdiy::load(bb, image.m_orig_bounds.X.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.Y.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.Z.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.X.Max);
    vtkmdiy::load(bb, image.m_orig_bounds.Y.Max);
    vtkmdiy::load(bb, image.m_orig_bounds.Z.Max);

    vtkmdiy::load(bb, image.m_bounds.X.Min);
    vtkmdiy::load(bb, image.m_bounds.Y.Min);
    vtkmdiy::load(bb, image.m_bounds.Z.Min);
    vtkmdiy::load(bb, image.m_bounds.X.Max);
    vtkmdiy::load(bb, image.m_bounds.Y.Max);
    vtkmdiy::load(bb, image.m_bounds.Z.Max);

    vtkmdiy::load(bb, image.m_payloads);
    vtkmdiy::load(bb, image.m_payload_bytes);
    vtkmdiy::load(bb, image.m_depths);
    vtkmdiy::load(bb, image.m_orig_rank);
  }
};

template <>
struct Serialization<vtkm::rendering::compositing::Image>
{
  static void save(BinaryBuffer& bb, const vtkm::rendering::compositing::Image& image)
  {
    vtkmdiy::save(bb, image.m_orig_bounds.X.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.Y.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.Z.Min);
    vtkmdiy::save(bb, image.m_orig_bounds.X.Max);
    vtkmdiy::save(bb, image.m_orig_bounds.Y.Max);
    vtkmdiy::save(bb, image.m_orig_bounds.Z.Max);

    vtkmdiy::save(bb, image.m_bounds.X.Min);
    vtkmdiy::save(bb, image.m_bounds.Y.Min);
    vtkmdiy::save(bb, image.m_bounds.Z.Min);
    vtkmdiy::save(bb, image.m_bounds.X.Max);
    vtkmdiy::save(bb, image.m_bounds.Y.Max);
    vtkmdiy::save(bb, image.m_bounds.Z.Max);

    vtkmdiy::save(bb, image.m_pixels);
    vtkmdiy::save(bb, image.m_depths);
    vtkmdiy::save(bb, image.m_orig_rank);
    vtkmdiy::save(bb, image.m_composite_order);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering::compositing::Image& image)
  {
    vtkmdiy::load(bb, image.m_orig_bounds.X.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.Y.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.Z.Min);
    vtkmdiy::load(bb, image.m_orig_bounds.X.Max);
    vtkmdiy::load(bb, image.m_orig_bounds.Y.Max);
    vtkmdiy::load(bb, image.m_orig_bounds.Z.Max);

    vtkmdiy::load(bb, image.m_bounds.X.Min);
    vtkmdiy::load(bb, image.m_bounds.Y.Min);
    vtkmdiy::load(bb, image.m_bounds.Z.Min);
    vtkmdiy::load(bb, image.m_bounds.X.Max);
    vtkmdiy::load(bb, image.m_bounds.Y.Max);
    vtkmdiy::load(bb, image.m_bounds.Z.Max);

    vtkmdiy::load(bb, image.m_pixels);
    vtkmdiy::load(bb, image.m_depths);
    vtkmdiy::load(bb, image.m_orig_rank);
    vtkmdiy::load(bb, image.m_composite_order);
  }
};
} //namespace vtkmdiy

#endif
