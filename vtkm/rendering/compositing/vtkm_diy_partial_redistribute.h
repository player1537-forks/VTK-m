//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_vtkm_diy_partial_redistribute_h
#define vtkm_rendering_compositing_vtkm_diy_partial_redistribute_h

#include <map>
#include <vtkm/rendering/compositing/vtkm_diy_partial_blocks.h>
#include <vtkm/thirdparty/diy/assigner.h>
#include <vtkm/thirdparty/diy/decomposition.h>
#include <vtkm/thirdparty/diy/master.h>
#include <vtkm/thirdparty/diy/point.h>
#include <vtkm/thirdparty/diy/reduce-operations.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

//
// Redistributes partial composites to the ranks that owns
// that sectoon of the image. Currently, the domain is decomposed
// in 1-D from min_pixel to max_pixel.
//
template <typename BlockType>
struct Redistribute
{
  const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>& m_decomposer;

  Redistribute(const vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds>& decomposer)
    : m_decomposer(decomposer)
  {
  }

  void operator()(void* v_block, const vtkmdiy::ReduceProxy& proxy) const
  {
    BlockType* block = static_cast<BlockType*>(v_block);
    //
    // first round we have no incoming. Take the partials we have
    // and sent them to to the right rank
    //
    if (proxy.in_link().size() == 0)
    {
      const int size = block->m_partials.size();
      std::map<vtkmdiy::BlockID, std::vector<typename BlockType::PartialType>> outgoing;

      for (int i = 0; i < size; ++i)
      {
        vtkmdiy::Point<int, VTKMDIY_MAX_DIM> point;
        point[0] = block->m_partials[i].m_pixel_id;
        int dest_gid = m_decomposer.point_to_gid(point);
        vtkmdiy::BlockID dest = proxy.out_link().target(dest_gid);
        outgoing[dest].push_back(block->m_partials[i]);
      } //for

      block->m_partials.clear();


      for (int i = 0; i < proxy.out_link().size(); ++i)
      {
        int dest_gid = proxy.out_link().target(i).gid;
        vtkmdiy::BlockID dest = proxy.out_link().target(dest_gid);
        proxy.enqueue(dest, outgoing[dest]);
        //outgoing[dest].clear();
      }

    } // if
    else
    {
      for (int i = 0; i < proxy.in_link().size(); ++i)
      {
        int gid = proxy.in_link().target(i).gid;
        std::vector<typename BlockType::PartialType> incoming_partials;
        proxy.dequeue(gid, incoming_partials);
        const int incoming_size = incoming_partials.size();
        // TODO: make this a std::copy
        for (int j = 0; j < incoming_size; ++j)
        {
          block->m_partials.push_back(incoming_partials[j]);
        }
      } // for

    }                            // else
    MPI_Barrier(MPI_COMM_WORLD); //HACK
  }                              // operator
};


template <typename AddBlockType>
void redistribute_detail(std::vector<typename AddBlockType::PartialType>& partials,
                         MPI_Comm comm,
                         const int& domain_min_pixel,
                         const int& domain_max_pixel)
{
  typedef typename AddBlockType::Block Block;

  vtkmdiy::mpi::communicator world(vtkmdiy::mpi::make_DIY_MPI_Comm(comm));
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "             DRP: Is this the right dimension???" << std::endl;
  vtkmdiy::DiscreteBounds global_bounds(1); //DRP???
  global_bounds.min[0] = domain_min_pixel;
  global_bounds.max[0] = domain_max_pixel;

  // tells diy to use all availible threads
  const int num_threads = 1;
  const int num_blocks = world.size();
  const int magic_k = 2;

  vtkmdiy::Master master(world, num_threads);

  // create an assigner with one block per rank
  vtkmdiy::ContiguousAssigner assigner(num_blocks, num_blocks);
  AddBlockType create(master, partials);

  const int dims = 1;
  vtkmdiy::RegularDecomposer<vtkmdiy::DiscreteBounds> decomposer(dims, global_bounds, num_blocks);
  decomposer.decompose(world.rank(), assigner, create);
  vtkmdiy::all_to_all(master, assigner, Redistribute<Block>(decomposer), magic_k);
}

//
// Define a default template that cannot be instantiated
//
template <typename T>
void redistribute(std::vector<T>& partials,
                  MPI_Comm comm,
                  const int& domain_min_pixel,
                  const int& domain_max_pixel);
// ----------------------------- VolumePartial Specialization------------------------------------------
template <>
void redistribute<VolumePartial<float>>(std::vector<VolumePartial<float>>& partials,
                                        MPI_Comm comm,
                                        const int& domain_min_pixel,
                                        const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<VolumeBlock<float>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

template <>
void redistribute<VolumePartial<double>>(std::vector<VolumePartial<double>>& partials,
                                         MPI_Comm comm,
                                         const int& domain_min_pixel,
                                         const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<VolumeBlock<double>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

// ----------------------------- AbsorpPartial Specialization------------------------------------------
template <>
void redistribute<AbsorptionPartial<double>>(std::vector<AbsorptionPartial<double>>& partials,
                                             MPI_Comm comm,
                                             const int& domain_min_pixel,
                                             const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<AbsorptionBlock<double>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

template <>
void redistribute<AbsorptionPartial<float>>(std::vector<AbsorptionPartial<float>>& partials,
                                            MPI_Comm comm,
                                            const int& domain_min_pixel,
                                            const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<AbsorptionBlock<float>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

// ----------------------------- EmissPartial Specialization------------------------------------------
template <>
void redistribute<EmissionPartial<double>>(std::vector<EmissionPartial<double>>& partials,
                                           MPI_Comm comm,
                                           const int& domain_min_pixel,
                                           const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<EmissionBlock<double>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

template <>
void redistribute<EmissionPartial<float>>(std::vector<EmissionPartial<float>>& partials,
                                          MPI_Comm comm,
                                          const int& domain_min_pixel,
                                          const int& domain_max_pixel)
{
  redistribute_detail<AddBlock<EmissionBlock<float>>>(
    partials, comm, domain_min_pixel, domain_max_pixel);
}

}
}
} //vtkm::rendering::compositing

#endif //vtkm_rendering_compositing_vtkm_diy_partial_redistribute_h
