#ifndef VTKH_DIY_DIRECT_SEND_HPP
#define VTKH_DIY_DIRECT_SEND_HPP

#include <vtkm/rendering/vtkh/compositing/Image.hpp>
#include <vtkm/thirdparty/diy/diy.h>
#include <sstream>

namespace vtkh
{

class DirectSendCompositor
{
public:
  DirectSendCompositor();
  ~DirectSendCompositor();
  void CompositeVolume(vtkmdiy::mpi::communicator &diy_comm,
                       std::vector<Image>     &images);
  std::string GetTimingString();
private:
  std::stringstream m_timing_log;
};

} // namespace vtkh
#endif
