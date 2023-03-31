//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTKH_DIY_RADIX_K_HPP
#define VTKH_DIY_RADIX_K_HPP

#include <vtkm/rendering/vtkh/compositing/Image.hpp>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>
#include <vtkm/thirdparty/diy/diy.h>

namespace vtkh
{

class RadixKCompositor
{
public:
  RadixKCompositor();
  ~RadixKCompositor();
  void CompositeSurface(vtkmdiy::mpi::communicator &diy_comm, Image &image);
  void CompositeSurface(vtkmdiy::mpi::communicator &diy_comm, PayloadImage &image);

  template<typename ImageType>
  void CompositeImpl(vtkmdiy::mpi::communicator &diy_comm, ImageType &image);

  std::string GetTimingString();
private:
  std::stringstream m_timing_log;
};

} // namspace vtkh

#endif
