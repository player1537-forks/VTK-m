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
#include <vtkm/cont/ErrorBadValue.h>

#include <algorithm>
#include <vtkm/rendering/compositing/Compositor.h>
#include <vtkm/rendering/compositing/ImageCompositor.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/rendering/compositing/DirectSendCompositor.h>
#include <vtkm/rendering/compositing/RadixKCompositor.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace rendering
{
namespace compositing
{

Compositor::Compositor()
  : m_composite_mode(Z_BUFFER_SURFACE)
{
}

Compositor::~Compositor() {}

void Compositor::SetCompositeMode(CompositeMode composite_mode)
{
  // assure we don't have mixed image types
  assert(m_images.size() == 0);
  m_composite_mode = composite_mode;
}

void Compositor::ClearImages()
{
  m_images.clear();
}

void Compositor::AddImage(vtkm::rendering::Canvas& canvas)
{
  auto colors = &(canvas.GetColorBuffer().ReadPortal().GetArray()[0][0]);
  auto depths = canvas.GetDepthBuffer().ReadPortal().GetArray();
  vtkm::Id width = canvas.GetWidth();
  vtkm::Id height = canvas.GetHeight();

  assert(m_composite_mode != VIS_ORDER_BLEND);
  assert(depths != NULL);
  Image image;
  if (m_images.size() == 0)
  {
    m_images.push_back(image);
    m_images[0].Init(colors, depths, width, height);
    //m_images[0].Save("first.png");
  }
  else if (m_composite_mode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(colors, depths, width, height);
    vtkm::rendering::compositing::ImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
  else
  {
    const size_t image_index = m_images.size();
    m_images.push_back(image);
    m_images[image_index].Init(colors, depths, width, height);
  }
}


/*
void Compositor::AddImage(const vtkm::cont::ArrayHandle<vtkm::Vec4<T>>& colors,
                          const vtkm::cont::ArrayHandle<T>& depths,
                          vtkm::Id width,
                          vtkm::Id height)
{
  auto c = colors.WritePortal().GetArray();
  auto d = depths.WritePortal().GetArray();
  this->AddImage(c, d, width, height);
}

void Compositor::AddImage(const unsigned char* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height)
{
  assert(m_composite_mode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (m_images.size() == 0)
  {
    m_images.push_back(image);
    m_images[0].Init(color_buffer, depth_buffer, width, height);
    //m_images[0].Save("first.png");
  }
  else if (m_composite_mode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);
    vtkm::rendering::compositing::ImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
  else
  {
    const size_t image_index = m_images.size();
    m_images.push_back(image);
    m_images[image_index].Init(color_buffer, depth_buffer, width, height);
  }
}

void Compositor::AddImage(const float* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height)
{
  assert(m_composite_mode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (m_images.size() == 0)
  {
    m_images.push_back(image);
    m_images[0].Init(color_buffer, depth_buffer, width, height);
  }
  else if (m_composite_mode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);

    vtkm::rendering::compositing::ImageCompositor compositor;
    compositor.ZBufferComposite(m_images[0], image);
  }
  else
  {
    const size_t image_index = m_images.size();
    m_images.push_back(image);
    m_images[image_index].Init(color_buffer, depth_buffer, width, height);
  }
}

void Compositor::AddImage(const unsigned char* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height,
                          const int vis_order)
{
  assert(m_composite_mode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = m_images.size();
  m_images.push_back(image);
  m_images[image_index].Init(color_buffer, depth_buffer, width, height, vis_order);
}


void Compositor::AddImage(const float* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height,
                          const int vis_order)
{
  assert(m_composite_mode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = m_images.size();
  m_images.push_back(image);

  m_images[image_index].Init(color_buffer, depth_buffer, width, height, vis_order);
}
*/

Image Compositor::Composite()
{
  assert(m_images.size() != 0);

  if (m_composite_mode == Z_BUFFER_SURFACE)
  {
    CompositeZBufferSurface();
  }
  else if (m_composite_mode == Z_BUFFER_BLEND)
  {
    CompositeZBufferBlend();
  }
  else if (m_composite_mode == VIS_ORDER_BLEND)
  {
    CompositeVisOrder();
  }
  // Make this a param to avoid the copy?
  return m_images[0];
}

void Compositor::Cleanup() {}

std::string Compositor::GetLogString()
{
  std::string res = m_log_stream.str();
  m_log_stream.str("");
  return res;
}

void Compositor::CompositeZBufferSurface()
{
  // nothing to do here in serial. Images were composited as
  // they were added to the compositor
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(m_images.size() == 1);
  RadixKCompositor compositor;
  compositor.CompositeSurface(comm, this->m_images[0]);
  m_log_stream << compositor.GetTimingString();
#endif
}

void Compositor::CompositeZBufferBlend()
{
  throw vtkm::cont::ErrorBadValue("Not implemented");
}

void Compositor::CompositeVisOrder()
{

#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  assert(m_images.size() != 0);
  vtkm::rendering::compositing::DirectSendCompositor compositor;
  compositor.CompositeVolume(comm, this->m_images);
#else
  vtkm::rendering::compositing::ImageCompositor compositor;
  compositor.OrderedComposite(m_images);
#endif
}

}
}
} //namespace vtkm::rendering::compositing
