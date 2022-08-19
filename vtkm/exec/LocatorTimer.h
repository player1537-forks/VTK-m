#ifndef vtk_m_exec_LocatorTimer_h
#define vtk_m_exec_LocatorTimer_h

#include <vector>
#include <numeric>
#include <functional>
#include <omp.h>

namespace vtkm
{
namespace exec
{

class LocatorTimer
{
public:
  static const std::size_t FindCell0_IDX = 0;
  static const std::size_t FindCell1_IDX = 1;
  static const std::size_t FindCell2_IDX = 2;
  static const std::size_t FindCellImpl0_IDX = 3;
  static const std::size_t FindCellImpl1_IDX = 4;
  static const std::size_t FindCellImpl2_IDX = 5;
  static const std::size_t PointInCell_IDX = 6;
  static const std::size_t NUM_IDX = 7;

  const std::vector<std::string> Names =
  {
    "FindCell_CId      ",
    "FindCell_Leaf     ",
    "FindCell_None     ",
    "FindCellImpl_OOB  ",
    "FindCellImpl_Good ",
    "FindCellImpl_Fail ",
    "PointInCell       "
  };

  LocatorTimer(int nThreads)
  {
    this->Counters.resize(this->NUM_IDX);
    this->Timers.resize(this->NUM_IDX);

    for (std::size_t i = 0; i < this->NUM_IDX; i++)
    {
      this->Counters[i].resize(nThreads, 0);
      this->Timers[i].resize(nThreads, 0.0);
    }

    this->Starts.resize(nThreads);
  }

  void Reset()
  {
    for (auto& c : this->Counters) for (auto& v : c) v = 0;
    for (auto& t : this->Timers) for (auto& v : t) v = 0.0;
  }

  void Start()
  {
#if !defined(VTKM_CUDA)
    int tid = omp_get_thread_num();
    std::chrono::time_point<std::chrono::steady_clock> s = std::chrono::steady_clock::now();
    //x = this->Starts[tid];
    this->Starts[tid] = s;
#endif
  }

  void Update(std::size_t idx)
  {
#if !defined(VTKM_CUDA)
    int tid = omp_get_thread_num();
    std::chrono::duration<double> dt = std::chrono::steady_clock::now() - this->Starts[tid];
    auto T = dt.count();
    this->Timers[idx][tid] += T;

    this->Counters[idx][tid]++;
#endif
  }

  std::vector<double> GetTime() const
  {
    std::vector<double> vals;

    for (const auto& v : this->Timers)
    {
      auto sum = std::accumulate(v.begin(), v.end(), 0.0);
      vals.push_back(sum);
    }

    return vals;
  }

  std::vector<vtkm::Id> GetCounters() const
  {
    std::vector<vtkm::Id> vals;

    for (const auto& v : this->Counters)
    {
      auto sum = std::accumulate(v.begin(), v.end(), 0);
      vals.push_back(sum);
    }

    return vals;
  }

  std::vector<std::vector<vtkm::Id>> Counters;
  std::vector<std::vector<double>> Timers;
  std::vector<std::chrono::time_point<std::chrono::steady_clock>> Starts;
};

}
}

#endif //vtk_m_exec_LocatorTimer_h
