//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef rover_blocks_h
#define rover_blocks_h

#include <vtkm/thirdparty/diy/master.h>

#include "AbsorptionPartial.hpp"
#include "EmissionPartial.hpp"
#include "VolumePartial.hpp"

namespace vtkh {

//--------------------------------------Volume Block Structure-----------------------------------
template<typename FloatType>
struct VolumeBlock
{
  typedef vtkmdiy::DiscreteBounds            Bounds;
  typedef VolumePartial<FloatType>       PartialType;
  std::vector<VolumePartial<FloatType>> &m_partials;
  VolumeBlock(std::vector<VolumePartial<FloatType>> &partials)
    : m_partials(partials)
  {}
};


//--------------------------------------Absorption Block Structure------------------------------
template<typename FloatType>
struct AbsorptionBlock
{
  typedef vtkmdiy::DiscreteBounds Bounds;
  typedef AbsorptionPartial<FloatType> PartialType;
  std::vector<AbsorptionPartial<FloatType>>   &m_partials;

  AbsorptionBlock(std::vector<AbsorptionPartial<FloatType>> &partials)
    : m_partials(partials)
  {
  }
};

//--------------------------------------Emission Block Structure------------------------------
template<typename FloatType>
struct EmissionBlock
{
  typedef vtkmdiy::DiscreteBounds Bounds;
  typedef EmissionPartial<FloatType> PartialType;
  std::vector<EmissionPartial<FloatType>>   &m_partials;

  EmissionBlock(std::vector<EmissionPartial<FloatType>> &partials)
    : m_partials(partials)
  {}
};

//--------------------------------------Add Block Template-----------------------------------
template<typename BlockType>
struct AddBlock
{
  typedef typename BlockType::PartialType PartialType;
  typedef BlockType                       Block;
  const vtkmdiy::Master &m_master;
  std::vector<PartialType> &m_partials;

  AddBlock(vtkmdiy::Master &master,std::vector<PartialType> &partials)
    : m_master(master), m_partials(partials)
  {
  }
  template<typename BoundsType, typename LinkType>
  void operator()(int gid,
                  const BoundsType &local_bounds,
                  const BoundsType &local_with_ghost_bounds,
                  const BoundsType &domain_bounds,
                  const LinkType &link) const
  {
    (void) local_bounds;
    (void) domain_bounds;
    (void) local_with_ghost_bounds;
    Block *block = new Block(m_partials);
    LinkType *rg_link = new LinkType(link);
    vtkmdiy::Master& master = const_cast<vtkmdiy::Master&>(m_master);
    int lid = master.add(gid, block, rg_link);
    (void) lid;
  }
};

} //namespace vtkh

//-------------------------------Serialization Specializations--------------------------------
namespace vtkmdiy {

template<>
struct Serialization<vtkh::AbsorptionPartial<double>>
{

  static void save(BinaryBuffer& bb, const vtkh::AbsorptionPartial<double> &partial)
  {
    vtkmdiy::save(bb, partial.m_bins);
    vtkmdiy::save(bb, partial.m_pixel_id);
    vtkmdiy::save(bb, partial.m_depth);
  }

  static void load(BinaryBuffer& bb, vtkh::AbsorptionPartial<double> &partial)
  {
    vtkmdiy::load(bb, partial.m_bins);
    vtkmdiy::load(bb, partial.m_pixel_id);
    vtkmdiy::load(bb, partial.m_depth);
  }
};

template<>
struct Serialization<vtkh::AbsorptionPartial<float>>
{

  static void save(BinaryBuffer& bb, const vtkh::AbsorptionPartial<float> &partial)
  {
    vtkmdiy::save(bb, partial.m_bins);
    vtkmdiy::save(bb, partial.m_pixel_id);
    vtkmdiy::save(bb, partial.m_depth);
  }

  static void load(BinaryBuffer& bb, vtkh::AbsorptionPartial<float> &partial)
  {
    vtkmdiy::load(bb, partial.m_bins);
    vtkmdiy::load(bb, partial.m_pixel_id);
    vtkmdiy::load(bb, partial.m_depth);
  }
};

template<>
struct Serialization<vtkh::EmissionPartial<double>>
{

  static void save(BinaryBuffer& bb, const vtkh::EmissionPartial<double> &partial)
  {
    vtkmdiy::save(bb, partial.m_bins);
    vtkmdiy::save(bb, partial.m_emission_bins);
    vtkmdiy::save(bb, partial.m_pixel_id);
    vtkmdiy::save(bb, partial.m_depth);
  }

  static void load(BinaryBuffer& bb, vtkh::EmissionPartial<double> &partial)
  {
    vtkmdiy::load(bb, partial.m_bins);
    vtkmdiy::load(bb, partial.m_emission_bins);
    vtkmdiy::load(bb, partial.m_pixel_id);
    vtkmdiy::load(bb, partial.m_depth);
  }
};

template<>
struct Serialization<vtkh::EmissionPartial<float>>
{

  static void save(BinaryBuffer& bb, const vtkh::EmissionPartial<float> &partial)
  {
    vtkmdiy::save(bb, partial.m_bins);
    vtkmdiy::save(bb, partial.m_emission_bins);
    vtkmdiy::save(bb, partial.m_pixel_id);
    vtkmdiy::save(bb, partial.m_depth);
  }

  static void load(BinaryBuffer& bb, vtkh::EmissionPartial<float> &partial)
  {
    vtkmdiy::load(bb, partial.m_bins);
    vtkmdiy::load(bb, partial.m_emission_bins);
    vtkmdiy::load(bb, partial.m_pixel_id);
    vtkmdiy::load(bb, partial.m_depth);
  }
};

} // namespace diy

#endif
