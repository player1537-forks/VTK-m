//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_vtkh_diy_image_block_h
#define vtkm_rendering_compositing_vtkh_diy_image_block_h

#include <vtkm/rendering/vtkh/compositing/Image.h>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.h>

namespace vtkh
{

template <typename ImageType>
struct ImageBlock
{
  ImageType& Image;
  ImageBlock(ImageType& image)
    : Image(image)
  {
  }
};

struct MultiImageBlock
{
  std::vector<Image>& Images;
  Image& Output;
  MultiImageBlock(std::vector<Image>& images, Image& output)
    : Images(images)
    , Output(output)
  {
  }
};

template <typename ImageType>
struct AddImageBlock
{
  ImageType& Image;
  const vtkmdiy::Master& Master;

  AddImageBlock(vtkmdiy::Master& master, ImageType& image)
    : Image(image)
    , Master(master)
  {
  }
  template <typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType&, // local_bounds
                  const BoundsType&, // local_with_ghost_bounds
                  const BoundsType&, // domain_bounds
                  const LinkType& link) const
  {
    ImageBlock<ImageType>* block = new ImageBlock<ImageType>(this->Image);
    LinkType* linked = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(this->Master);
    master.add(gid, block, linked);
  }
};

struct AddMultiImageBlock
{
  std::vector<Image>& Images;
  const vtkmdiy::Master& Master;
  Image& Output;

  AddMultiImageBlock(vtkmdiy::Master& master, std::vector<Image>& images, Image& output)
    : Images(images)
    , Master(master)
    , Output(output)
  {
  }
  template <typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType&, // local_bounds
                  const BoundsType&, // local_with_ghost_bounds
                  const BoundsType&, // domain_bounds
                  const LinkType& link) const
  {
    MultiImageBlock* block = new MultiImageBlock(Images, Output);
    LinkType* linked = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(this->Master);
    master.add(gid, block, linked);
  }
};

} //namespace  vtkh

namespace vtkmdiy
{

template <>
struct Serialization<vtkh::PayloadImage>
{
  static void save(BinaryBuffer& bb, const vtkh::PayloadImage& image)
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

  static void load(BinaryBuffer& bb, vtkh::PayloadImage& image)
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
struct Serialization<vtkh::Image>
{
  static void save(BinaryBuffer& bb, const vtkh::Image& image)
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

  static void load(BinaryBuffer& bb, vtkh::Image& image)
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

} // namespace diy

#endif //vtkm_rendering_compositing_vtkh_diy_image_block_h
