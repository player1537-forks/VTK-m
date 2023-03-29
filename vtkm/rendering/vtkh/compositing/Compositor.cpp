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
#include <vtkm/rendering/vtkh/compositing/Compositor.hpp>
#include <vtkm/rendering/vtkh/compositing/ImageCompositor.hpp>

#include <assert.h>
#include <algorithm>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
#include <vtkm/rendering/vtkh/compositing/DirectSendCompositor.hpp>
#include <vtkm/rendering/vtkh/compositing/RadixKCompositor.hpp>
#endif

namespace vtkh
{

Compositor::Compositor()
  : m_composite_mode(Z_BUFFER_SURFACE)
{

}

Compositor::~Compositor()
{

}

void
Compositor::SetCompositeMode(CompositeMode composite_mode)
{
  // assure we don't have mixed image types
  assert(m_images.size() == 0);
  m_composite_mode = composite_mode;
}

void
Compositor::ClearImages()
{
  m_images.clear();
}

void
Compositor::AddImage(const unsigned char *color_buffer,
                     const float *        depth_buffer,
                     const int            width,
                     const int            height)
{
  assert(m_composite_mode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if(m_images.size() == 0)
  {
    m_images.push_back(image);
    m_images[0].Init(color_buffer,
                     depth_buffer,
                     width,
                     height);
    //m_images[0].Save("first.png");
  }
  else if(m_composite_mode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer,
               depth_buffer,
               width,
               height);
    vtkh::ImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0],image);
  }
  else
  {
    const size_t image_index = m_images.size();
    m_images.push_back(image);
    m_images[image_index].Init(color_buffer,
                               depth_buffer,
                               width,
                               height);
  }

}

void
Compositor::AddImage(const float *color_buffer,
                     const float *depth_buffer,
                     const int    width,
                     const int    height)
{
  assert(m_composite_mode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if(m_images.size() == 0)
  {
    m_images.push_back(image);
    m_images[0].Init(color_buffer,
                     depth_buffer,
                     width,
                     height);
  }
  else if(m_composite_mode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer,
               depth_buffer,
               width,
               height);

    vtkh::ImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0],image);
  }
  else
  {
    const size_t image_index = m_images.size();
    m_images.push_back(image);
    m_images[image_index].Init(color_buffer,
                               depth_buffer,
                               width,
                               height);
  }

}

void
Compositor::AddImage(const unsigned char *color_buffer,
                     const float         *depth_buffer,
                     const int            width,
                     const int            height,
                     const int            vis_order)
{
  assert(m_composite_mode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = m_images.size();
  m_images.push_back(image);
  m_images[image_index].Init(color_buffer,
                             depth_buffer,
                             width,
                             height,
                             vis_order);
}

void
Compositor::AddImage(const float *color_buffer,
                     const float *depth_buffer,
                     const int    width,
                     const int    height,
                     const int    vis_order)
{
  assert(m_composite_mode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = m_images.size();
  m_images.push_back(image);

  m_images[image_index].Init(color_buffer,
                             depth_buffer,
                             width,
                             height,
                             vis_order);
}

Image
Compositor::Composite()
{
  assert(m_images.size() != 0);

  if(m_composite_mode == Z_BUFFER_SURFACE)
  {
    CompositeZBufferSurface();
  }
  else if(m_composite_mode == Z_BUFFER_BLEND)
  {
    CompositeZBufferBlend();
  }
  else if(m_composite_mode == VIS_ORDER_BLEND)
  {
    CompositeVisOrder();
  }
  // Make this a param to avoid the copy?
  return m_images[0];
}

void
Compositor::Cleanup()
{

}

std::string
Compositor::GetLogString()
{
  std::string res = m_log_stream.str();
  m_log_stream.str("");
  return res;
}

void
Compositor::CompositeZBufferSurface()
{
  // nothing to do here in serial. Images were composited as
  // they were added to the compositor
#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(m_images.size() == 1);
  RadixKCompositor compositor;
  compositor.CompositeSurface(diy_comm, this->m_images[0]);
  m_log_stream<<compositor.GetTimingString();
#endif
}

void
Compositor::CompositeZBufferBlend()
{
  assert("this is not implemented yet" == "error");
}

void
Compositor::CompositeVisOrder()
{

#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(m_images.size() != 0);
  DirectSendCompositor compositor;
  compositor.CompositeVolume(diy_comm, this->m_images);
#else
  vtkh::ImageCompositor compositor;
  compositor.OrderedComposite(m_images);
#endif
}

} // namespace vtkh
