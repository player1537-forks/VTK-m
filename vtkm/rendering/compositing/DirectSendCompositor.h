//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_DirectSendCompositor_h
#define vtk_m_rendering_compositing_DirectSendCompositor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/compositing/Image.h>

#include <vtkm/thirdparty/diy/diy.h>

namespace vtkm
{
namespace rendering
{
namespace compositing
{

class VTKM_RENDERING_EXPORT DirectSendCompositor
{
public:
  DirectSendCompositor();
  ~DirectSendCompositor();
  void CompositeVolume(vtkmdiy::mpi::communicator& diy_comm,
                       std::vector<vtkm::rendering::compositing::Image>& images);
  std::string GetTimingString();

private:
  std::stringstream m_timing_log;
};

}
}
} //namespace vtkm::rendering::compositing

#endif //vtk_m_rendering_compositing_DirectSendCompositor_h
