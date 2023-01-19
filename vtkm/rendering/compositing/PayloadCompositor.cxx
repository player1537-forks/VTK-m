#include <vtkm/rendering/compositing/PayloadCompositor.h>
#include <vtkm/rendering/compositing/PayloadImageCompositor.h>

#include <algorithm>
#include <assert.h>

#ifdef VTKH_PARALLEL
#include <diy/mpi.hpp>
#include <mpi.h>
#include <vtkh/compositing/RadixKCompositor.hpp>
#include <vtkh/vtkh.hpp>
#endif

using namespace vtkm::rendering::compositing;

namespace vtkh
{

PayloadCompositor::PayloadCompositor() {}

void PayloadCompositor::ClearImages()
{
  m_images.clear();
}

void PayloadCompositor::AddImage(PayloadImage& image)
{
  assert(image.GetNumberOfPixels() != 0);

  if (m_images.size() == 0)
  {
    m_images.push_back(image);
  }
  else
  {
    //
    // Do local composite and keep a single image
    //
    PayloadImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
}

PayloadImage PayloadCompositor::Composite()
{
  assert(m_images.size() != 0);
  // nothing to do here in serial. Images were composited as
  // they were added to the compositor
#ifdef VTKH_PARALLEL
  vtkhdiy::mpi::communicator diy_comm;
  diy_comm = vtkhdiy::mpi::communicator(MPI_Comm_f2c(GetMPICommHandle()));

  assert(m_images.size() == 1);
  RadixKCompositor compositor;
  compositor.CompositeSurface(diy_comm, this->m_images[0]);
#endif
  // Make this a param to avoid the copy?
  return m_images[0];
}


} // namespace vtkh
