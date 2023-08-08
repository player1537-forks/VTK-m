//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_DirectSendCompositor_h
#define vtkm_rendering_new_DirectSendCompositor_h

#include <sstream>
#include <vtkm/rendering_new/compositing/Image.h>
#include <vtkm/thirdparty/diy/diy.h>

namespace vtkm
{
namespace rendering_new
{

class DirectSendCompositor
{
public:
  DirectSendCompositor();
  ~DirectSendCompositor();
  void CompositeVolume(vtkmdiy::mpi::communicator& diy_comm, std::vector<Image>& images);
  std::string GetTimingString();

private:
  std::stringstream m_timing_log;
};

}
} // namespace vtkm::rendering_new


#endif //vtkm_rendering_new_DirectSendCompositor_h
