//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_rendering_compositing_RadixKCompositor_h
#define vtk_m_rendering_compositing_RadixKCompositor_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/compositing/PayloadImage.h>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
//#include <mpi.h>
//#include <vtkm/thirdparty/diy/mpi.h>
//#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace rendering
{
namespace compositing
{

class VTKM_RENDERING_EXPORT RadixKCompositor
{
public:
  RadixKCompositor();
  ~RadixKCompositor();
  void CompositeSurface(vtkmdiy::mpi::communicator& diy_comm, Image& image);
  void CompositeSurface(vtkmdiy::mpi::communicator& diy_comm, PayloadImage& image);

  template <typename ImageType>
  void CompositeImpl(vtkmdiy::mpi::communicator& diy_comm, ImageType& image);

  std::string GetTimingString();

private:
  std::stringstream m_timing_log;
};

}
}
} //namespace vtkm::rendering::compositing

#endif //vtk_m_rendering_compositing_RadixKCompositor_h
