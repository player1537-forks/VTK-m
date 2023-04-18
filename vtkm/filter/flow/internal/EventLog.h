//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_flow_internal_EventLog_h
#define vtk_m_filter_flow_internal_EventLog_h

#include <vtkm/cont/Timer.h>

#include <vtkm/filter/flow/vtkm_filter_flow_export.h>

#include <sstream>
#include <string>
#include <vector>

//#define VTKM_ENABLE_EVENTLOG 1

namespace vtkm
{
namespace filter
{
namespace flow
{
namespace internal
{

class VTKM_FILTER_FLOW_EXPORT EventLog
{
public:
  VTKM_CONT EventLog() {}
  VTKM_CONT void Init();
  VTKM_CONT void Output() const;
  VTKM_CONT void SetTimeBufferOutputFile(const std::string& fname);
  VTKM_CONT void SetEventBufferOutputFile(const std::string& fname);
  VTKM_CONT void SetStep(const vtkm::Id& s);
  VTKM_CONT void StartTimer();

  VTKM_CONT void RecordTime(const std::string& eventName);

  template <typename T>
  inline VTKM_CONT void RecordEvent(const std::string& eventName, const T& val);
  template <typename T>
  inline VTKM_CONT void RecordEvent(const std::string& eventName, const std::vector<T>& vals);

private:

#ifdef VTKM_ENABLE_EVENTLOG
  VTKM_CONT void Clear();
  VTKM_CONT vtkm::Float64 GetElapsedTime() const;

  std::stringstream TimeBuffer;
  std::string TimeBufferOutputName = "timer";
  std::stringstream EventBuffer;
  std::string EventBufferOutputName = "event";
  vtkm::Id Step = 0;
  vtkm::cont::Timer Timer;
#endif
};

template <typename T>
void EventLog::RecordEvent(const std::string& eventName, const T& val)
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->EventBuffer<<eventName<<"_"<<this->Step<<" "<<val<<std::endl;
#else
  (void) eventName;
  (void) val;
#endif
}

template <typename T>
void EventLog::RecordEvent(const std::string& eventName, const std::vector<T>& vals)
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->EventBuffer<<eventName<<"_"<<this->Step<<" ";
  for (auto it = vals.begin(); it != vals.end(); it++)
  {
    if (it != vals.begin())
      this->EventBuffer << " ";
    this->EventBuffer << *it;
  }
  this->EventBuffer<<std::endl;
#else
  (void) eventName;
  (void) vals;
#endif
}

}
}
}
} //namespace vtkm::filter::flow::internal


#endif //vtk_m_filter_flow_internal_EventLog_h
