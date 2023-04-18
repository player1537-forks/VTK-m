//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/ErrorExecution.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/filter/flow/internal/EventLog.h>

//#define VTKM_ENABLE_EVENTLOG 1

namespace vtkm
{
namespace filter
{
namespace flow
{
namespace internal
{

void
EventLog::Init()
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->Clear();
#endif
}

void EventLog::SetTimeBufferOutputFile(const std::string& fname)
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->TimeBufferOutputName = fname;
#else
  (void)fname;
  shit();
#endif
}

void EventLog::SetEventBufferOutputFile(const std::string& fname)
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->EventBufferOutputName = fname;
#else
  (void)fname;
  shit();
#endif
}

void EventLog::SetStep(const vtkm::Id& s)
{
#ifdef VTKM_ENABLE_EVENTLOG
  this->Step = s;
#else
  (void)s;
#endif
}

void EventLog::RecordTime(const std::string& eventName)
{
#ifdef VTKM_ENABLE_EVENTLOG
  if (!this->Timer.Started())
    throw vtkm::cont::ErrorExecution("Error: Timer not started.");
  this->TimeBuffer<<eventName<<"_"<<this->Step<<" "<<this->GetElapsedTime()<<std::endl;
#else
  (void) eventName;
#endif
}

void EventLog::StartTimer()
{
#ifdef VTKM_ENABLE_EVENTLOG
  if (this->Timer.Started())
    throw vtkm::cont::ErrorExecution("Error: Timer already started.");
  this->Timer.Start();
#endif
}

void EventLog::Output() const
{
#ifdef VTKM_ENABLE_EVENTLOG
  std::stringstream timerOutName, eventOutName;

  timerOutName << this->TimeBufferOutputName;
  eventOutName << this->EventBufferOutputName;
#ifdef VTKM_ENABLE_MPI
  int rank = vtkm::cont::EnvironmentTracker::GetCommunicator().rank();
  timerOutName<<"."<<rank;
  eventOutName<<"."<<rank;
#endif
  timerOutName<<".log";
  eventOutName<<".log";

  std::ofstream timerOut(timerOutName.str(), std::ofstream::out);
  std::ofstream eventOut(eventOutName.str(), std::ofstream::out);

  timerOut<<this->TimeBuffer.rdbuf();
  timerOut.close();
  eventOut<<this->EventBuffer.rdbuf();
  eventOut.close();
#endif
}

#ifdef VTKM_ENABLE_EVENTLOG
void
EventLog::Clear()
{
  this->Step = 0;
  this->TimeBuffer.clear();
  this->EventBuffer.clear();
  if (this->Timer.Started())
    this->Timer.Stop();
}

vtkm::Float64 EventLog::GetElapsedTime() const
{
  //convert to miliseconds.
  return (this->Timer.GetElapsedTime() * 1000);
}
#endif

}
}
}
} //namespace vtkm::filter::flow::internal
