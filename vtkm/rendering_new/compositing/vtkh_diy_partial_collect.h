//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_vtkh_diy_partial_collect_h
#define vtkm_rendering_new_vtkh_diy_partial_collect_h

#include <vtkm/rendering_new/compositing/AbsorptionPartial.h>
#include <vtkm/rendering_new/compositing/EmissionPartial.h>
#include <vtkm/rendering_new/compositing/VolumePartial.h>
#include <vtkm/rendering_new/compositing/vtkh_diy_partial_blocks.h>

#include <vtkm/thirdparty/diy/assigner.h>
#include <vtkm/thirdparty/diy/decomposition.h>
#include <vtkm/thirdparty/diy/master.h>
#include <vtkm/thirdparty/diy/reduce-operations.h>

namespace vtkm
{
namespace rendering_new
{

//
// Collect struct sends all data to a single node.
//
template <typename BlockType>
struct Collect
{
  const vtkmdiy::RegularDecomposer<vtkmdiy::ContinuousBounds>& m_decomposer;

  Collect(const vtkmdiy::RegularDecomposer<vtkmdiy::ContinuousBounds>& decomposer)
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
    const int collection_rank = 0;
    if (proxy.in_link().size() == 0 && proxy.gid() != collection_rank)
    {
      int dest_gid = collection_rank;
      vtkmdiy::BlockID dest = proxy.out_link().target(dest_gid);
      proxy.enqueue(dest, block->Partials);

      block->Partials.clear();

    } // if
    else if (proxy.gid() == collection_rank)
    {

      for (int i = 0; i < proxy.in_link().size(); ++i)
      {
        int gid = proxy.in_link().target(i).gid;
        if (gid == collection_rank)
        {
          continue;
        }
        //TODO: leave the paritals that start here, here
        std::vector<typename BlockType::PartialType> incoming_partials;
        proxy.dequeue(gid, incoming_partials);
        const int incoming_size = incoming_partials.size();
        // TODO: make this a std::copy
        for (int j = 0; j < incoming_size; ++j)
        {
          block->Partials.push_back(incoming_partials[j]);
        }
      } // for
    }   // else

  } // operator
};

//
// collect uses the all-to-all construct to perform a gather to
// the root rank. All other ranks will have no data
//
template <typename AddBlockType>
void collect_detail(std::vector<typename AddBlockType::PartialType>& partials, MPI_Comm comm)
{
  typedef typename AddBlockType::Block Block;

  vtkmdiy::mpi::communicator world(vtkmdiy::mpi::make_DIY_MPI_Comm(comm));
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "             DRP: Is this the right dimension???" << std::endl;
  vtkmdiy::ContinuousBounds global_bounds(1);
  global_bounds.min[0] = 0;
  global_bounds.max[0] = 1;

  // tells diy to use all availible threads
  const int num_threads = -1;
  const int num_blocks = world.size();
  const int magic_k = 2;

  vtkmdiy::Master master(world, num_threads);

  // create an assigner with one block per rank
  vtkmdiy::ContiguousAssigner assigner(num_blocks, num_blocks);
  AddBlockType create(master, partials);

  const int dims = 1;
  vtkmdiy::RegularDecomposer<vtkmdiy::ContinuousBounds> decomposer(dims, global_bounds, num_blocks);
  decomposer.decompose(world.rank(), assigner, create);

  vtkmdiy::all_to_all(master, assigner, Collect<Block>(decomposer), magic_k);
}

template <typename T>
void collect(std::vector<T>& partials, MPI_Comm comm);

template <>
void collect<VolumePartial<float>>(std::vector<VolumePartial<float>>& partials, MPI_Comm comm)
{
  collect_detail<AddBlock<VolumeBlock<float>>>(partials, comm);
}

template <>
void collect<VolumePartial<double>>(std::vector<VolumePartial<double>>& partials, MPI_Comm comm)
{
  collect_detail<AddBlock<VolumeBlock<double>>>(partials, comm);
}

template <>
void collect<AbsorptionPartial<double>>(std::vector<AbsorptionPartial<double>>& partials,
                                        MPI_Comm comm)
{
  collect_detail<AddBlock<AbsorptionBlock<double>>>(partials, comm);
}

template <>
void collect<AbsorptionPartial<float>>(std::vector<AbsorptionPartial<float>>& partials,
                                       MPI_Comm comm)
{
  collect_detail<AddBlock<AbsorptionBlock<float>>>(partials, comm);
}

template <>
void collect<EmissionPartial<double>>(std::vector<EmissionPartial<double>>& partials, MPI_Comm comm)
{
  collect_detail<AddBlock<EmissionBlock<double>>>(partials, comm);
}

template <>
void collect<EmissionPartial<float>>(std::vector<EmissionPartial<float>>& partials, MPI_Comm comm)
{
  collect_detail<AddBlock<EmissionBlock<float>>>(partials, comm);
}

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_new_vtkh_diy_partial_collect_h
