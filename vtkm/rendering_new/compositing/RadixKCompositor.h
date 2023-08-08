//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_new_RadixKCompositor_h
#define vtkm_rendering_new_RadixKCompositor_h

#include <vtkm/rendering_new/compositing/Image.h>
#include <vtkm/rendering_new/compositing/PayloadImage.h>
#include <vtkm/thirdparty/diy/diy.h>

namespace vtkm
{
namespace rendering_new
{

class RadixKCompositor
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
} // vtkm::rendering_new

#endif //vtkm_rendering_new_RadixKCompositor_h
