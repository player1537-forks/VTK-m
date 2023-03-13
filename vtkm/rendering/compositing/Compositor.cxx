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
  : CompositingMode(Z_BUFFER_SURFACE)
{
}

Compositor::~Compositor() {}

void Compositor::SetCompositeMode(CompositeMode composite_mode)
{
  // assure we don't have mixed image types
  assert(this->Images.size() == 0);
  this->CompositingMode = composite_mode;
}

void Compositor::ClearImages()
{
  this->Images.clear();
}

void Compositor::AddImage(vtkm::rendering::Canvas& canvas)
{
  auto colors = &(canvas.GetColorBuffer().ReadPortal().GetArray()[0][0]);
  auto depths = canvas.GetDepthBuffer().ReadPortal().GetArray();
  vtkm::Id width = canvas.GetWidth();
  vtkm::Id height = canvas.GetHeight();

  // assert(this->CompositingMode != VIS_ORDER_BLEND);
  assert(depths != NULL);
  Image image;
  if (this->Images.size() == 0)
  {
    this->Images.push_back(image);
    this->Images[0].Init(colors, depths, width, height);
    //this->Images[0].Save("first.png");
  }
  else if (this->CompositingMode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(colors, depths, width, height);
    vtkm::rendering::compositing::ImageCompositor compositor;
    compositor.ZBufferComposite(this->Images[0], image);
  }
  else
  {
    const size_t image_index = this->Images.size();
    this->Images.push_back(image);
    this->Images[image_index].Init(colors, depths, width, height);
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
  assert(this->CompositingMode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (this->Images.size() == 0)
  {
    this->Images.push_back(image);
    this->Images[0].Init(color_buffer, depth_buffer, width, height);
    //this->Images[0].Save("first.png");
  }
  else if (this->CompositingMode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);
    vtkm::rendering::compositing::ImageCompositor compositor;
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
  assert(this->CompositingMode != VIS_ORDER_BLEND);
  assert(depth_buffer != NULL);
  Image image;
  if (this->Images.size() == 0)
  {
    this->Images.push_back(image);
    this->Images[0].Init(color_buffer, depth_buffer, width, height);
  }
  else if (this->CompositingMode == Z_BUFFER_SURFACE)
  {
    //
    // Do local composite and keep a single image
    //
    image.Init(color_buffer, depth_buffer, width, height);

    vtkm::rendering::compositing::ImageCompositor compositor;
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
  assert(this->CompositingMode == VIS_ORDER_BLEND);
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
  assert(this->CompositingMode == VIS_ORDER_BLEND);
  Image image;
  const size_t image_index = this->Images.size();
  this->Images.push_back(image);

  this->Images[image_index].Init(color_buffer, depth_buffer, width, height, vis_order);
}
*/

Image Compositor::Composite()
{
  assert(this->Images.size() != 0);

  if (this->CompositingMode == Z_BUFFER_SURFACE)
  {
    CompositeZBufferSurface();
  }
  else if (this->CompositingMode == Z_BUFFER_BLEND)
  {
    CompositeZBufferBlend();
  }
  else if (this->CompositingMode == VIS_ORDER_BLEND)
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
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  assert(this->Images.size() == 1);
  RadixKCompositor compositor;
  compositor.CompositeSurface(comm, this->Images[0]);
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
  assert(this->Images.size() != 0);
  vtkm::rendering::compositing::DirectSendCompositor compositor;
  compositor.CompositeVolume(comm, this->Images);
#else
  vtkm::rendering::compositing::ImageCompositor compositor;
  compositor.OrderedComposite(this->Images);
#endif
}

}
}
} //namespace vtkm::rendering::compositing
