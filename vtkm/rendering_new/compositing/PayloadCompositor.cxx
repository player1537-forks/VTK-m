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
#include <vtkm/rendering_new/compositing/PayloadCompositor.h>
#include <vtkm/rendering_new/compositing/PayloadImageCompositor.h>

#include <algorithm>
#include <assert.h>
#include <vtkm/thirdparty/diy/diy.h>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
#include <vtkm/rendering_new/compositing/RadixKCompositor.h>
#endif

namespace vtkm
{
namespace rendering_new
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
    vtkm::rendering_new::PayloadImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
}

PayloadImage PayloadCompositor::Composite()
{
  assert(m_images.size() != 0);
  // nothing to do here in serial. Images were composited as
  // they were added to the compositor
#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  //vtkmdiy::mpi::communicator diy_comm;
  //diy_comm = vtkmdiy::mpi::communicator(MPI_Comm_f2c(GetMPICommHandle()));

  assert(m_images.size() == 1);
  RadixKCompositor compositor;
  compositor.CompositeSurface(diy_comm, this->m_images[0]);
#endif
  // Make this a param to avoid the copy?
  return m_images[0];
}


}
} //vtkm::rendering_new
