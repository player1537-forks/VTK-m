//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_flow_internal_AdvectAlgorithm_h
#define vtk_m_filter_flow_internal_AdvectAlgorithm_h

#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/filter/flow/internal/BoundsMap.h>
#include <vtkm/filter/flow/internal/DataSetIntegrator.h>
#include <vtkm/filter/flow/internal/ParticleMessenger.h>

#ifdef VTKM_ENABLE_MPI
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace filter
{
namespace flow
{
namespace internal
{

template <typename DSIType, template <typename> class ResultType, typename ParticleType>
class AdvectAlgorithm
{
public:
  AdvectAlgorithm(const vtkm::filter::flow::internal::BoundsMap& bm, std::vector<DSIType>& blocks)
    : Blocks(blocks)
    , BoundsMap(bm)
    , NumRanks(this->Comm.size())
    , Rank(this->Comm.rank())
  {
  }

  void Execute(vtkm::Id numSteps,
               vtkm::FloatDefault stepSize,
               const vtkm::cont::ArrayHandle<ParticleType>& seeds)
  {
    this->SetNumberOfSteps(numSteps);
    this->SetStepSize(stepSize);
    this->SetSeeds(seeds);
    this->Go();
  }

  vtkm::cont::PartitionedDataSet GetOutput() const
  {
    vtkm::cont::PartitionedDataSet output;

    for (const auto& b : this->Blocks)
    {
      vtkm::cont::DataSet ds;
      if (b.template GetOutput<ParticleType>(ds))
        output.AppendPartition(ds);
    }

    return output;
  }

  void SetStepSize(vtkm::FloatDefault stepSize) { this->StepSize = stepSize; }
  void SetNumberOfSteps(vtkm::Id numSteps) { this->NumberOfSteps = numSteps; }
  void SetSeeds(const vtkm::cont::ArrayHandle<ParticleType>& seeds)
  {
    this->ClearParticles();

    vtkm::Id n = seeds.GetNumberOfValues();
    auto portal = seeds.ReadPortal();

    std::vector<std::vector<vtkm::Id>> blockIDs;
    std::vector<ParticleType> particles;
    for (vtkm::Id i = 0; i < n; i++)
    {
      const ParticleType p = portal.Get(i);
      std::vector<vtkm::Id> ids = this->BoundsMap.FindBlocks(p.Pos);

      if (!ids.empty() && this->BoundsMap.FindRank(ids[0]) == this->Rank)
      {
        particles.emplace_back(p);
        blockIDs.emplace_back(ids);
      }
    }

    this->ValidateSeeds(particles);

    this->SetSeedArray(particles, blockIDs);
  }

  //Advect all the particles.
  virtual void Go()
  {
    vtkm::filter::flow::internal::ParticleMessenger<ParticleType> messenger(
      this->Comm, this->BoundsMap, 1, 128);

    vtkm::Id nLocal = static_cast<vtkm::Id>(this->Active.size() + this->Inactive.size());
    this->ComputeTotalNumParticles(nLocal);

    while (this->TotalNumTerminatedParticles < this->TotalNumParticles)
    {
      std::vector<ParticleType> v;
      vtkm::Id numTerm = 0, blockId = -1;
      if (this->GetActiveParticles(v, blockId))
      {
        //make this a pointer to avoid the copy?
        auto& block = this->GetDataSet(blockId);
        DSIHelperInfoType bb =
          DSIHelperInfo<ParticleType>(v, this->BoundsMap, this->ParticleBlockIDsMap);
        block.Advect(bb, this->StepSize, this->NumberOfSteps);
        numTerm = this->UpdateResult(bb.Get<DSIHelperInfo<ParticleType>>());
      }

      vtkm::Id numTermMessages = 0;
      this->Communicate(messenger, numTerm, numTermMessages);

      this->TotalNumTerminatedParticles += (numTerm + numTermMessages);
      if (this->TotalNumTerminatedParticles > this->TotalNumParticles)
        throw vtkm::cont::ErrorFilterExecution("Particle count error");
    }
  }

protected:
  virtual void ClearParticles()
  {
    this->Active.clear();
    this->Inactive.clear();
    this->ParticleBlockIDsMap.clear();
  }

  virtual bool GetBlockAndWait(const vtkm::Id& numLocalTerm)
  {
    //There are only two cases where blocking would deadlock.
    //1. There are active particles.
    //2. numLocalTerm + this->TotalNumberOfTerminatedParticles == this->TotalNumberOfParticles
    //So, if neither are true, we can safely block and wait for communication to come in.

    if (this->Active.empty() && this->Inactive.empty() &&
        (numLocalTerm + this->TotalNumTerminatedParticles < this->TotalNumParticles))
    {
      return true;
    }

    return false;
  }

  virtual void SetSeedArray(const std::vector<ParticleType>& particles,
                            const std::vector<std::vector<vtkm::Id>>& blockIds)
  {
    VTKM_ASSERT(particles.size() == blockIds.size());

    auto pit = particles.begin();
    auto bit = blockIds.begin();
    while (pit != particles.end() && bit != blockIds.end())
    {
      this->ParticleBlockIDsMap[pit->ID] = *bit;
      pit++;
      bit++;
    }

    this->Active.insert(this->Active.end(), particles.begin(), particles.end());
  }

  virtual bool GetActiveParticles(std::vector<ParticleType>& particles, vtkm::Id& blockId)
  {
    particles.clear();
    blockId = -1;
    if (this->Active.empty())
      return false;

    blockId = this->ParticleBlockIDsMap[this->Active.front().ID][0];
    auto it = this->Active.begin();
    while (it != this->Active.end())
    {
      auto p = *it;
      if (blockId == this->ParticleBlockIDsMap[p.ID][0])
      {
        particles.emplace_back(p);
        it = this->Active.erase(it);
      }
      else
        it++;
    }

    return !particles.empty();
  }

  virtual void Communicate(vtkm::filter::flow::internal::ParticleMessenger<ParticleType>& messenger,
                           vtkm::Id numLocalTerminations,
                           vtkm::Id& numTermMessages)
  {
    std::vector<ParticleType> incoming;
    std::unordered_map<vtkm::Id, std::vector<vtkm::Id>> incomingIDs;
    numTermMessages = 0;
    messenger.Exchange(this->Inactive,
                       this->ParticleBlockIDsMap,
                       numLocalTerminations,
                       incoming,
                       incomingIDs,
                       numTermMessages,
                       this->GetBlockAndWait(numLocalTerminations));

    this->Inactive.clear();
    this->UpdateActive(incoming, incomingIDs);
  }

  virtual void UpdateActive(const std::vector<ParticleType>& particles,
                            const std::unordered_map<vtkm::Id, std::vector<vtkm::Id>>& idsMap)
  {
    this->Update(this->Active, particles, idsMap);
  }

  virtual void UpdateInactive(const std::vector<ParticleType>& particles,
                              const std::unordered_map<vtkm::Id, std::vector<vtkm::Id>>& idsMap)
  {
    this->Update(this->Inactive, particles, idsMap);
  }

  void ComputeTotalNumParticles(const vtkm::Id& numLocal)
  {
    long long total = static_cast<long long>(numLocal);
#ifdef VTKM_ENABLE_MPI
    MPI_Comm mpiComm = vtkmdiy::mpi::mpi_cast(this->Comm.handle());
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, mpiComm);
#endif
    this->TotalNumParticles = static_cast<vtkm::Id>(total);
  }

  DataSetIntegrator<DSIType>& GetDataSet(vtkm::Id id)
  {
    for (auto& it : this->Blocks)
      if (it.GetID() == id)
        return it;

    throw vtkm::cont::ErrorFilterExecution("Bad block");
  }

  vtkm::Id UpdateResult(const DSIHelperInfo<ParticleType>& stuff)
  {
    this->UpdateActive(stuff.A, stuff.IdMapA);
    this->UpdateInactive(stuff.I, stuff.IdMapI);

    vtkm::Id numTerm = static_cast<vtkm::Id>(stuff.TermID.size());
    //Update terminated particles.
    if (numTerm > 0)
    {
      for (const auto& id : stuff.TermID)
        this->ParticleBlockIDsMap.erase(id);
    }

    return numTerm;
  }

private:
  void ValidateSeeds(const std::vector<ParticleType>& particles)
  {
    long numLocal = static_cast<long>(particles.size());

    std::vector<long> pids;
    pids.reserve(numLocal);
    for (const auto& p : particles)
      pids.emplace_back(p.ID);

    //First, do a local check of my particles.
    std::set<long> pidSet(pids.begin(), pids.end());
    if (pidSet.size() != pids.size())
      throw vtkm::cont::ErrorFilterExecution("Duplicate IDs detected");

#ifdef VTKM_ENABLE_MPI
    //Need to do a global check for parallel case.
    auto MPIComm = vtkmdiy::mpi::mpi_cast(this->Comm.handle());
    //Make sure we do not have duplicate IDs.
    std::vector<long> recvData;
    if (this->Rank == 0)
      recvData.resize(this->NumRanks);

    MPI_Gather(&numLocal, 1, MPI_LONG_INT, recvData.data(), 1, MPI_LONG_INT, 0, MPIComm);
    if (this->Rank == 0)
    {
      long totalNum = std::accumulate(recvData.begin(), recvData.end(), 0);
      recvData.resize(totalNum);
    }
    std::vector<int> counts, offsets;
    if (this->Rank == 0)
    {
      counts.resize(this->NumRanks);
      offsets.resize(this->NumRanks);
    }

    MPI_Gatherv(pids.data(), numLocal, MPI_LONG_INT, recvData.data(), counts.data(), offsets.data(), MPI_LONG_INT, 0, MPIComm);

    if (this->Rank == 0)
    {
      //Putting the ids in a set will detect if we have duplicate IDs.
      std::set<long> pidSet2(pids.begin(), pids.end());
      if (pidSet2.size() != pids.size())
        throw vtkm::cont::ErrorFilterExecution("Duplicate IDs detected");
    }
#endif
  }

  void Update(std::vector<ParticleType>& arr,
              const std::vector<ParticleType>& particles,
              const std::unordered_map<vtkm::Id, std::vector<vtkm::Id>>& idsMap)
  {
    VTKM_ASSERT(particles.size() == idsMap.size());

    arr.insert(arr.end(), particles.begin(), particles.end());
    for (const auto& it : idsMap)
      this->ParticleBlockIDsMap[it.first] = it.second;
  }

protected:
  //Member data
  std::vector<ParticleType> Active;
  std::vector<DSIType> Blocks;
  vtkm::filter::flow::internal::BoundsMap BoundsMap;
  vtkmdiy::mpi::communicator Comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  std::vector<ParticleType> Inactive;
  vtkm::Id NumberOfSteps;
  vtkm::Id NumRanks;
  std::unordered_map<vtkm::Id, std::vector<vtkm::Id>> ParticleBlockIDsMap;
  vtkm::Id Rank;
  vtkm::FloatDefault StepSize;
  vtkm::Id TotalNumParticles = 0;
  vtkm::Id TotalNumTerminatedParticles = 0;
};

}
}
}
} //vtkm::filter::flow::internal

#endif //vtk_m_filter_flow_internal_AdvectAlgorithm_h
