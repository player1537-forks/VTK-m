#ifndef VTKH_DIY_RADIX_K_HPP
#define VTKH_DIY_RADIX_K_HPP

#include <vtkm/rendering/vtkh/compositing/Image.hpp>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>
#include <vtkm/thirdparty/diy/diy.h>

/*
#include <diy/mpi.hpp>
#include <sstream>
*/

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
