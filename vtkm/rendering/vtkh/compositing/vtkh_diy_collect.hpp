//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_vtkh_diy_collect_h
#define vtkm_rendering_compositing_vtkh_diy_collect_h

#include <vtkm/rendering/vtkh/compositing/Image.hpp>
#include <vtkm/rendering/vtkh/compositing/vtkh_diy_image_block.hpp>

namespace vtkh
{

template<typename ImageType>
struct CollectImages
{
  const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> &m_decomposer;

  CollectImages(const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> &decomposer)
    : m_decomposer(decomposer)
  {}

  void operator()(void *b, const vtkmdiy::ReduceProxy &proxy) const
  {
    ImageBlock<ImageType> *block = reinterpret_cast<ImageBlock<ImageType>*>(b);
    //
    // first round we have no incoming. Take the images we have
    // and sent them to to the right rank
    //
    const int collection_rank = 0;
    if(proxy.in_link().size() == 0)
    {

      if(proxy.gid() != collection_rank)
      {
        int dest_gid = collection_rank;
        vtkmdiy::BlockID dest = proxy.out_link().target(dest_gid);

        proxy.enqueue(dest, block->m_image);
        block->m_image.Clear();
      }
    } // if
    else if(proxy.gid() == collection_rank)
    {
      ImageType final_image;
      final_image.InitOriginal(block->m_image);
      block->m_image.SubsetTo(final_image);

      for(int i = 0; i < proxy.in_link().size(); ++i)
      {
        int gid = proxy.in_link().target(i).gid;

        if(gid == collection_rank)
        {
          continue;
        }
        ImageType incoming;
        proxy.dequeue(gid, incoming);
        incoming.SubsetTo(final_image);
      } // for
      block->m_image.Swap(final_image);
    } // else

  } // operator
};

} // namespace vtkh

#endif //vtkm_rendering_compositing_vtkh_diy_collect_h
