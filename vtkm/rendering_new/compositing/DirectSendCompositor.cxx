//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering_new/compositing/DirectSendCompositor.h>
#include <vtkm/rendering_new/compositing/ImageCompositor.h>
#include <vtkm/rendering_new/compositing/vtkh_diy_collect.h>
#include <vtkm/rendering_new/compositing/vtkh_diy_utils.h>

namespace vtkm
{
namespace rendering_new
{

struct Redistribute
{
  typedef vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> Decomposer;
  const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>& m_decomposer;
  Redistribute(const Decomposer& decomposer)
    : m_decomposer(decomposer)
  {
  }

  void operator()(void* v_block, const vtkmdiy::ReduceProxy& proxy) const
  {
    MultiImageBlock* block = static_cast<MultiImageBlock*>(v_block);
    //
    // first round we have no incoming. Take the image we have,
    // chop it up into pieces, and send it to the domain resposible
    // for that portion
    //
    const int world_size = m_decomposer.nblocks;
    const int local_images = block->Images.size();
    if (proxy.in_link().size() == 0)
    {
      std::map<vtkmdiy::BlockID, std::vector<Image>> outgoing;

      for (int i = 0; i < world_size; ++i)
      {
        vtkmdiy::DiscreteBounds sub_image_bounds(3);
        m_decomposer.fill_bounds(sub_image_bounds, i);
        vtkm::Bounds vtkm_sub_bounds = DIYBoundsToVTKM(sub_image_bounds);

        vtkmdiy::BlockID dest = proxy.out_link().target(i);
        outgoing[dest].resize(local_images);

        for (int img = 0; img < local_images; ++img)
        {
          outgoing[dest][img].SubsetFrom(block->Images[img], vtkm_sub_bounds);
        }
      } //for

      typename std::map<vtkmdiy::BlockID, std::vector<Image>>::iterator it;
      for (it = outgoing.begin(); it != outgoing.end(); ++it)
      {
        proxy.enqueue(it->first, it->second);
      }
    } // if
    else if (block->Images.at(0).CompositeOrder != -1)
    {
      // blend images according to vis order
      std::vector<Image> images;
      for (int i = 0; i < proxy.in_link().size(); ++i)
      {

        std::vector<Image> incoming;
        int gid = proxy.in_link().target(i).gid;
        proxy.dequeue(gid, incoming);
        const int in_size = incoming.size();
        for (int img = 0; img < in_size; ++img)
        {
          images.emplace_back(incoming[img]);
          //std::cout<<"rank "<<rank<<" rec "<<incoming[img].ToString()<<"\n";
        }
      } // for

      ImageCompositor compositor;
      compositor.OrderedComposite(images);

      block->Output.Swap(images[0]);
    } // else if
    else if (block->Images.at(0).CompositeOrder == -1 && block->Images.at(0).GetHasTransparency())
    {
      std::vector<Image> images;
      for (int i = 0; i < proxy.in_link().size(); ++i)
      {

        std::vector<Image> incoming;
        int gid = proxy.in_link().target(i).gid;
        proxy.dequeue(gid, incoming);
        const int in_size = incoming.size();
        for (int img = 0; img < in_size; ++img)
        {
          images.emplace_back(incoming[img]);
          //std::cout<<"rank "<<rank<<" rec "<<incoming[img].ToString()<<"\n";
        }
      } // for

      //
      // we have images with a depth buffer and transparency
      //
      ImageCompositor compositor;
      compositor.ZBufferBlend(images);
    }

  } // operator
};

DirectSendCompositor::DirectSendCompositor() {}

DirectSendCompositor::~DirectSendCompositor() {}

void DirectSendCompositor::CompositeVolume(vtkmdiy::mpi::communicator& diy_comm,
                                           std::vector<Image>& images)
{
  vtkmdiy::DiscreteBounds global_bounds = VTKMBoundsToDIY(images.at(0).OrigBounds);

  const int num_threads = 1;
  const int num_blocks = diy_comm.size();
  const int magic_k = 8;
  Image sub_image;
  //
  // DIY does not seem to like being called with different block types
  // so we isolate them within separate blocks
  //
  {
    vtkmdiy::Master master(diy_comm, num_threads, -1, 0, [](void* b) {
      ImageBlock<Image>* block = reinterpret_cast<ImageBlock<Image>*>(b);
      delete block;
    });

    // create an assigner with one block per rank
    vtkmdiy::ContiguousAssigner assigner(num_blocks, num_blocks);

    AddMultiImageBlock create(master, images, sub_image);

    const int dims = 2;
    vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> decomposer(dims, global_bounds, num_blocks);
    decomposer.decompose(diy_comm.rank(), assigner, create);

    vtkmdiy::all_to_all(master, assigner, Redistribute(decomposer), magic_k);
  }

  {
    vtkmdiy::Master master(diy_comm, num_threads, -1, 0, [](void* b) {
      ImageBlock<Image>* block = reinterpret_cast<ImageBlock<Image>*>(b);
      delete block;
    });
    vtkmdiy::ContiguousAssigner assigner(num_blocks, num_blocks);

    const int dims = 2;
    vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> decomposer(dims, global_bounds, num_blocks);
    AddImageBlock<Image> all_create(master, sub_image);
    decomposer.decompose(diy_comm.rank(), assigner, all_create);
    diy_comm.barrier();
    //MPI_Barrier(diy_comm);

    //MPICollect(sub_image,diy_comm);
    vtkmdiy::all_to_all(master, assigner, CollectImages<Image>(decomposer), magic_k);
  }

  images.at(0).Swap(sub_image);
}

std::string DirectSendCompositor::GetTimingString()
{
  std::string res(m_timing_log.str());
  m_timing_log.str("");
  return res;
}

}
} //vtkm::rendering_new
