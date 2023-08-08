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
#include <vtkm/rendering_new/compositing/Compositor.h>
#include <vtkm/rendering_new/compositing/ImageCompositor.h>

#include <algorithm>
#include <assert.h>

#include <vtkm/thirdparty/diy/diy.h>
#ifdef VTKM_ENABLE_MPI
#include <vtkm/rendering_new/compositing/DirectSendCompositor.h>
#include <vtkm/rendering_new/compositing/RadixKCompositor.h>
#endif

namespace vtkm
{
namespace rendering_new
{

Compositor::Compositor()
  : CompositeMode(Z_BUFFER_SURFACE)
{
}

Compositor::~Compositor() {}

void Compositor::SetCompositeMode(enum CompositeMode composite_mode)
{
  // assure we don't have mixed image types
  assert(this->Images.size() == 0);
  this->CompositeMode = composite_mode;
}

void Compositor::ClearImages()
{
  this->Images.clear();
}

void Compositor::AddImage(const unsigned char* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height)
{
  assert(this->CompositeMode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (this->Images.size() == 0)
  {
    this->Images.push_back(image);
    this->Images[0].Init(color_buffer, depth_buffer, width, height);
    //m_images[0].Save("first.png");
  }
  else if (this->CompositeMode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);
    vtkm::rendering_new::ImageCompositor compositor;
    compositor.ZBufferComposite(this->Images[0], image);
  }
  else
  {
    const size_t image_index = this->Images.size();
    this->Images.push_back(image);
    this->Images[image_index].Init(color_buffer, depth_buffer, width, height);
  }
}

void Compositor::AddImage(const float* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height)
{
  assert(this->CompositeMode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (this->Images.size() == 0)
  {
    this->Images.push_back(image);
    this->Images[0].Init(color_buffer, depth_buffer, width, height);
  }
  else if (this->CompositeMode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);

    vtkm::rendering_new::ImageCompositor compositor;
    compositor.ZBufferComposite(this->Images[0], image);
  }
  else
  {
    const size_t image_index = this->Images.size();
    this->Images.push_back(image);
    this->Images[image_index].Init(color_buffer, depth_buffer, width, height);
  }
}

void Compositor::AddImage(const unsigned char* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height,
                          const int vis_order)
{
  assert(this->CompositeMode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = this->Images.size();
  this->Images.push_back(image);
  this->Images[image_index].Init(color_buffer, depth_buffer, width, height, vis_order);
}

void Compositor::AddImage(const float* color_buffer,
                          const float* depth_buffer,
                          const int width,
                          const int height,
                          const int vis_order)
{
  assert(this->CompositeMode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = this->Images.size();
  this->Images.push_back(image);

  this->Images[image_index].Init(color_buffer, depth_buffer, width, height, vis_order);
}

Image Compositor::Composite()
{
  assert(this->Images.size() != 0);

  if (this->CompositeMode == Z_BUFFER_SURFACE)
  {
    CompositeZBufferSurface();
  }
  else if (this->CompositeMode == Z_BUFFER_BLEND)
  {
    CompositeZBufferBlend();
  }
  else if (this->CompositeMode == VIS_ORDER_BLEND)
  {
    CompositeVisOrder();
  }
  // Make this a param to avoid the copy?
  return this->Images[0];
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
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(this->Images.size() == 1);
  if (diy_comm.size() > 1)
  {
    RadixKCompositor compositor;
    compositor.CompositeSurface(diy_comm, this->Images[0]);
    m_log_stream << compositor.GetTimingString();
  }
#endif
}

void Compositor::CompositeZBufferBlend()
{
  throw vtkm::cont::ErrorBadValue("Compoistor::CompositeZBufferBlend not implemented");
}

void Compositor::CompositeVisOrder()
{

#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(this->Images.size() != 0);
  DirectSendCompositor compositor;
  compositor.CompositeVolume(diy_comm, this->Images);
#else
  vtkm::rendering_new::ImageCompositor compositor;
  compositor.OrderedComposite(this->Images);
#endif
}

}
} // namespace vtkm::rendering_new
