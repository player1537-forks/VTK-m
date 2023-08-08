##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

set(headers
  ## Add VTK-h headers
  compositing/AbsorptionPartial.h
  compositing/Compositor.h
  compositing/DirectSendCompositor.h
  compositing/EmissionPartial.h
  compositing/ImageCompositor.h
  compositing/Image.h
  compositing/PartialCompositor.h
  compositing/PayloadCompositor.h
  compositing/PayloadImageCompositor.h
  compositing/PayloadImage.h
  compositing/PNGEncoder.h
  compositing/RadixKCompositor.h
  compositing/VolumePartial.h
  compositing/vtkh_diy_collect.h
  compositing/vtkh_diy_image_block.h
  compositing/vtkh_diy_partial_blocks.h
  compositing/vtkh_diy_partial_collect.h
  compositing/vtkh_diy_partial_redistribute.h
  compositing/vtkh_diy_utils.h

  Annotator.h
  LineRenderer.h
  MeshRenderer.h
  PointRenderer.h
  RayTracer.h
  Renderer.h
  Plot.h
  ScalarRenderer.h
  Scene.h
  VolumeRenderer.h
  )

set(sources
  ##vtk-h sources
  compositing/Compositor.cxx
  compositing/DirectSendCompositor.cxx
  compositing/Image.cxx
  compositing/PartialCompositor.cxx
  compositing/PayloadCompositor.cxx
  compositing/PayloadImage.cxx
  compositing/PNGEncoder.cxx
  compositing/RadixKCompositor.cxx

  Annotator.cxx
  LineRenderer.cxx
  MeshRenderer.cxx
  PointRenderer.cxx
  RayTracer.cxx
  Plot.cxx
  Renderer.cxx
  ScalarRenderer.cxx
  Scene.cxx
  )

# This list of sources has code that uses devices and so might need to be
# compiled with a device-specific compiler (like CUDA).
set(device_sources
  VolumeRenderer.cxx
  )

# the None backend supports not building the opengl version

#-----------------------------------------------------------------------------
vtkm_library(
  NAME vtkm_rendering_new
  SOURCES ${sources}
  HEADERS ${headers}
  DEVICE_SOURCES ${device_sources}
  )

#-----------------------------------------------------------------------------
if(UNIX AND NOT APPLE)
  target_link_libraries(vtkm_rendering_new PRIVATE rt)
endif()

#-----------------------------------------------------------------------------
#add_subdirectory(internal)
#add_subdirectory(raytracing)

#-----------------------------------------------------------------------------
# Pretty much obsolete, but might be expected to be exported.
include(VTKmRenderingContexts)
