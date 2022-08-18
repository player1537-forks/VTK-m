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
  static const std::size_t FindCell_IDX = 0;
  static const std::size_t FindCellImpl_IDX = 1;
  static const std::size_t PointInCell_IDX = 2;

  LocatorTimer(int nThreads)
  {
    this->Timers.resize(3);
    this->Counters.resize(3);
    for (int i = 0; i < 3; i++)
    {
      this->Timers[i].resize(nThreads, 0.0);
      this->Counters[i].resize(nThreads, 0);
    }
  }

  void Reset()
  {
    for (auto& t : this->Timers)
      for (auto& v : t)
        v = 0.0;

    for (auto& c : this->Counters)
      for (auto& v : c)
        v = 0;
  }

  void Update(std::size_t idx, const std::chrono::time_point<std::chrono::steady_clock>& start)
  {
    int tid = omp_get_thread_num();
    std::chrono::duration<double> dt = std::chrono::steady_clock::now()-start;
    auto T = dt.count();
    this->Timers[idx][tid] += T;

    this->Counters[idx][tid]++;
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

  std::vector<std::vector<double>> Timers;
  std::vector<std::vector<vtkm::Id>> Counters;
};

}
}

#endif //vtk_m_exec_LocatorTimer_h
