//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/rendering/compositing/PayloadCompositor.h>
#include <vtkm/rendering/compositing/PayloadImageCompositor.h>
#include <vtkm/rendering/compositing/RadixKCompositor.h>

#include <algorithm>
#include <assert.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
//#include <vtkh/vtkh.hpp>
#include <vtkm/rendering/compositing/RadixKCompositor.h>
//#include <diy/mpi.hpp>
#endif


namespace vtkm
{
namespace rendering
{
namespace compositing
{

PayloadCompositor::PayloadCompositor() {}

void PayloadCompositor::ClearImages()
{
  m_images.clear();
}

void PayloadCompositor::AddImage(vtkm::rendering::compositing::PayloadImage& image)
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
    vtkm::rendering::compositing::PayloadImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
}

vtkm::rendering::compositing::PayloadImage PayloadCompositor::Composite()
{
  assert(m_images.size() != 0);
  // nothing to do here in serial. Images were composited as
  // they were added to the compositor
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  //  vtkmdiy::mpi::communicator diy_comm;
  //  diy_comm = vtkmdiy::mpi::communicator(MPI_Comm_f2c(GetMPICommHandle()));

  assert(m_images.size() == 1);
  vtkm::rendering::compositing::RadixKCompositor compositor;
  compositor.CompositeSurface(comm, this->m_images[0]);
#endif
  // Make this a param to avoid the copy?
  return m_images[0];
}


}
}
} // namespace vtkm:rendering::compositing
