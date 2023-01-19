//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_vtkm_diy_collect_h
#define vtk_m_rendering_compositing_vtkm_diy_collect_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/compositing/vtkm_diy_image_block.h>
#include <vtkmdiy/master.hpp>
#include <vtkmdiy/partners/swap.hpp>
#include <vtkmdiy/reduce-operations.hpp>
#include <vtkmdiy/reduce.hpp>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

template <typename ImageType>
struct CollectImages
{
  const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>& m_decomposer;

  CollectImages(const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>& decomposer)
    : m_decomposer(decomposer)
  {
  }

  void operator()(void* b, const vtkmdiy::ReduceProxy& proxy) const
  {
    ImageBlock<ImageType>* block = reinterpret_cast<ImageBlock<ImageType>*>(b);
    //
    // first round we have no incoming. Take the images we have
    // and sent them to to the right rank
    //
    const int collection_rank = 0;
    if (proxy.in_link().size() == 0)
    {

      if (proxy.gid() != collection_rank)
      {
        int dest_gid = collection_rank;
        vtkmdiy::BlockID dest = proxy.out_link().target(dest_gid);

        proxy.enqueue(dest, block->m_image);
        block->m_image.Clear();
      }
    } // if
    else if (proxy.gid() == collection_rank)
    {
      ImageType final_image;
      final_image.InitOriginal(block->m_image);
      block->m_image.SubsetTo(final_image);

      for (int i = 0; i < proxy.in_link().size(); ++i)
      {
        int gid = proxy.in_link().target(i).gid;

        if (gid == collection_rank)
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

}
}
} //namespace vtkm::rendering::compositing

#endif //vtk_m_rendering_compositing_vtkm_diy_collect_h
