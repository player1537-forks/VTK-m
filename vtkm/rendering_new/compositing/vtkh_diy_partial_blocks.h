//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_vtkh_diy_partial_blocks_h
#define vtkm_rendering_new_vtkh_diy_partial_blocks_h

#include <vtkm/thirdparty/diy/master.h>

#include <vtkm/rendering_new/compositing/AbsorptionPartial.h>
#include <vtkm/rendering_new/compositing/EmissionPartial.h>
#include <vtkm/rendering_new/compositing/VolumePartial.h>

namespace vtkm
{
namespace rendering_new
{

//--------------------------------------Volume Block Structure-----------------------------------
template <typename FloatType>
struct VolumeBlock
{
  typedef vtkmdiy::DiscreteBounds Bounds;
  typedef VolumePartial<FloatType> PartialType;
  std::vector<VolumePartial<FloatType>>& Partials;
  VolumeBlock(std::vector<VolumePartial<FloatType>>& partials)
    : Partials(partials)
  {
  }
};


//--------------------------------------Absorption Block Structure------------------------------
template <typename FloatType>
struct AbsorptionBlock
{
  typedef vtkmdiy::DiscreteBounds Bounds;
  typedef AbsorptionPartial<FloatType> PartialType;
  std::vector<AbsorptionPartial<FloatType>>& Partials;

  AbsorptionBlock(std::vector<AbsorptionPartial<FloatType>>& partials)
    : Partials(partials)
  {
  }
};

//--------------------------------------Emission Block Structure------------------------------
template <typename FloatType>
struct EmissionBlock
{
  typedef vtkmdiy::DiscreteBounds Bounds;
  typedef EmissionPartial<FloatType> PartialType;
  std::vector<EmissionPartial<FloatType>>& Partials;

  EmissionBlock(std::vector<EmissionPartial<FloatType>>& partials)
    : Partials(partials)
  {
  }
};

//--------------------------------------Add Block Template-----------------------------------
template <typename BlockType>
struct AddBlock
{
  typedef typename BlockType::PartialType PartialType;
  typedef BlockType Block;
  const vtkmdiy::Master& Master;
  std::vector<PartialType>& Partials;

  AddBlock(vtkmdiy::Master& master, std::vector<PartialType>& partials)
    : Master(master)
    , Partials(partials)
  {
  }
  template <typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType& local_bounds,
                  const BoundsType& local_with_ghost_bounds,
                  const BoundsType& domain_bounds,
                  const LinkType& link) const
  {
    (void)local_bounds;
    (void)domain_bounds;
    (void)local_with_ghost_bounds;
    Block* block = new Block(this->Partials);
    LinkType* rg_link = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(this->Master);
    int lid = master.add(gid, block, rg_link);
    (void)lid;
  }
};

}
} //namespace vtkm::rendering_new

//-------------------------------Serialization Specializations--------------------------------
namespace vtkmdiy
{

template <>
struct Serialization<vtkm::rendering_new::AbsorptionPartial<double>>
{

  static void save(BinaryBuffer& bb, const vtkm::rendering_new::AbsorptionPartial<double>& partial)
  {
    vtkmdiy::save(bb, partial.Bins);
    vtkmdiy::save(bb, partial.PixelId);
    vtkmdiy::save(bb, partial.Depth);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering_new::AbsorptionPartial<double>& partial)
  {
    vtkmdiy::load(bb, partial.Bins);
    vtkmdiy::load(bb, partial.PixelId);
    vtkmdiy::load(bb, partial.Depth);
  }
};

template <>
struct Serialization<vtkm::rendering_new::AbsorptionPartial<float>>
{

  static void save(BinaryBuffer& bb, const vtkm::rendering_new::AbsorptionPartial<float>& partial)
  {
    vtkmdiy::save(bb, partial.Bins);
    vtkmdiy::save(bb, partial.PixelId);
    vtkmdiy::save(bb, partial.Depth);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering_new::AbsorptionPartial<float>& partial)
  {
    vtkmdiy::load(bb, partial.Bins);
    vtkmdiy::load(bb, partial.PixelId);
    vtkmdiy::load(bb, partial.Depth);
  }
};

template <>
struct Serialization<vtkm::rendering_new::EmissionPartial<double>>
{

  static void save(BinaryBuffer& bb, const vtkm::rendering_new::EmissionPartial<double>& partial)
  {
    vtkmdiy::save(bb, partial.Bins);
    vtkmdiy::save(bb, partial.EmissionBins);
    vtkmdiy::save(bb, partial.PixelId);
    vtkmdiy::save(bb, partial.Depth);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering_new::EmissionPartial<double>& partial)
  {
    vtkmdiy::load(bb, partial.Bins);
    vtkmdiy::load(bb, partial.EmissionBins);
    vtkmdiy::load(bb, partial.PixelId);
    vtkmdiy::load(bb, partial.Depth);
  }
};

template <>
struct Serialization<vtkm::rendering_new::EmissionPartial<float>>
{

  static void save(BinaryBuffer& bb, const vtkm::rendering_new::EmissionPartial<float>& partial)
  {
    vtkmdiy::save(bb, partial.Bins);
    vtkmdiy::save(bb, partial.EmissionBins);
    vtkmdiy::save(bb, partial.PixelId);
    vtkmdiy::save(bb, partial.Depth);
  }

  static void load(BinaryBuffer& bb, vtkm::rendering_new::EmissionPartial<float>& partial)
  {
    vtkmdiy::load(bb, partial.Bins);
    vtkmdiy::load(bb, partial.EmissionBins);
    vtkmdiy::load(bb, partial.PixelId);
    vtkmdiy::load(bb, partial.Depth);
  }
};

} //namespace vtkmdiy

#endif //vtkm_rendering_new_vtkh_diy_partial_blocks_h
