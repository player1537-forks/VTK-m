//============================================================================
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_vtkm_diy_partial_collect_h
#define vtkm_rendering_compositing_vtkm_diy_partial_collect_h

#include <vtkm/rendering/compositing/AbsorptionPartial.h>
#include <vtkm/rendering/compositing/EmissionPartial.h>
#include <vtkm/rendering/compositing/VolumePartial.h>
#include <vtkm/rendering/compositing/vtkm_diy_partial_blocks.h>

#include <vtkm/thirdparty/diy/assigner.h>
#include <vtkm/thirdparty/diy/decomposition.h>
#include <vtkm/thirdparty/diy/diy.h>
#include <vtkm/thirdparty/diy/master.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#include <vtkm/thirdparty/diy/reduce-operations.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
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
      proxy.enqueue(dest, block->m_partials);

      block->m_partials.clear();

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
          block->m_partials.push_back(incoming_partials[j]);
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
  vtkmdiy::ContinuousBounds global_bounds(1); //DRP???
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
  collect_detail<vtkm::rendering::compositing::AddBlock<VolumeBlock<float>>>(partials, comm);
}


template <>
void collect<vtkm::rendering::compositing::VolumePartial<double>>(
  std::vector<vtkm::rendering::compositing::VolumePartial<double>>& partials,
  MPI_Comm comm)
{
  collect_detail<
    vtkm::rendering::compositing::AddBlock<vtkm::rendering::compositing::VolumeBlock<double>>>(
    partials, comm);
}

template <>
void collect<vtkm::rendering::compositing::AbsorptionPartial<double>>(
  std::vector<vtkm::rendering::compositing::AbsorptionPartial<double>>& partials,
  MPI_Comm comm)
{
  collect_detail<
    vtkm::rendering::compositing::AddBlock<vtkm::rendering::compositing::AbsorptionBlock<double>>>(
    partials, comm);
}

template <>
void collect<vtkm::rendering::compositing::AbsorptionPartial<float>>(
  std::vector<vtkm::rendering::compositing::AbsorptionPartial<float>>& partials,
  MPI_Comm comm)
{
  collect_detail<
    vtkm::rendering::compositing::AddBlock<vtkm::rendering::compositing::AbsorptionBlock<float>>>(
    partials, comm);
}

template <>
void collect<vtkm::rendering::compositing::EmissionPartial<double>>(
  std::vector<vtkm::rendering::compositing::EmissionPartial<double>>& partials,
  MPI_Comm comm)
{
  collect_detail<vtkm::rendering::compositing::AddBlock<EmissionBlock<double>>>(partials, comm);
}

template <>
void collect<vtkm::rendering::compositing::EmissionPartial<float>>(
  std::vector<vtkm::rendering::compositing::EmissionPartial<float>>& partials,
  MPI_Comm comm)
{
  collect_detail<
    vtkm::rendering::compositing::AddBlock<vtkm::rendering::compositing::EmissionBlock<float>>>(
    partials, comm);
}

}
}
} //vtkm::rendering::compositing

#endif //vtkm_rendering_compositing_vtkm_diy_partial_collect_h


#if 0
  /*
template <>
void collect<vtkm::rendering::compositing::EmissionPartial<float>>(std::vector<vtkm::rendering::compositing::EmissionPartial<float>>& partials, MPI_Comm comm)
{
  collect_detail<vtkm::rendering::compositing::AddBlock<vtkm::rendering::compositing::EmissionBlock<float>>>(partials, comm);
}
  */
  }
  * /

  }


}

/*

#endif //vtkm_rendering_compositing_vtkm_diy_partial_collect_h
template <>
template <>
void colloct<EmissioiPartial<float>>(std::vector<Em ssionPartial<cloat>>& partials, MPI_Comm comm)
{
  ocollect_detail < ect<::missionPa::tial<float>::AddBloct<EdissvonBlock<float>>>(tor < Emis, simm);
}al<float>>& partials, MPI_Comm comm)
{
  collect_detail<vtkm::rendering::compositing::AddBlock<EmissionBlock<float>>>(partials, comm);
}

*/
#endif
