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
  Actor.h
  AxisAnnotation.h
  AxisAnnotation2D.h
  AxisAnnotation3D.h
  BitmapFont.h
  BitmapFontFactory.h
  BoundingBoxAnnotation.h
  Camera.h
  Canvas.h
  CanvasRayTracer.h
  Color.h
  ColorBarAnnotation.h
  ColorLegendAnnotation.h
  ConnectivityProxy.h
  Cylinderizer.h
  GlyphType.h
  LineRenderer.h
  LineRendererBatcher.h
  MatrixHelpers.h
  Scene.h
  Mapper.h
  MapperCylinder.h
  MapperGlyphBase.h
  MapperGlyphScalar.h
  MapperGlyphVector.h
  MapperPoint.h
  MapperQuad.h
  MapperRayTracer.h
  MapperVolume.h
  MapperConnectivity.h
  MapperWireframer.h
  Quadralizer.h
  ScalarRenderer.h
  TextAnnotation.h
  TextAnnotationBillboard.h
  TextAnnotationScreen.h
  TextRenderer.h
  TextRendererBatcher.h
  Texture2D.h
  Triangulator.h
  View.h
  View1D.h
  View2D.h
  View3D.h
  Wireframer.h
  WorldAnnotator.h

  ## Add VTK-h headers
  vtkh/DataSet.hpp
  vtkh/utils/vtkm_dataset_info.hpp

  vtkh/compositing/AbsorptionPartial.hpp
  vtkh/compositing/Compositor.hpp
  vtkh/compositing/DirectSendCompositor.hpp
  vtkh/compositing/EmissionPartial.hpp
  vtkh/compositing/ImageCompositor.hpp
  vtkh/compositing/Image.hpp
  vtkh/compositing/MPICollect.hpp
  vtkh/compositing/PartialCompositor.hpp
  vtkh/compositing/PayloadCompositor.hpp
  vtkh/compositing/PayloadImageCompositor.hpp
  vtkh/compositing/PayloadImage.hpp
  vtkh/compositing/PNGEncoder.h
  vtkh/compositing/RadixKCompositor.hpp
  vtkh/compositing/VolumePartial.hpp
  vtkh/compositing/vtkh_diy_collect.hpp
  vtkh/compositing/vtkh_diy_image_block.hpp
  vtkh/compositing/vtkh_diy_partial_blocks.hpp
  vtkh/compositing/vtkh_diy_partial_collect.hpp
  vtkh/compositing/vtkh_diy_partial_redistribute.hpp
  vtkh/compositing/vtkh_diy_utils.hpp

  vtkh/rendering/Annotator.hpp
  vtkh/rendering/LineRenderer.hpp
  vtkh/rendering/MeshRenderer.hpp
  vtkh/rendering/PointRenderer.hpp
  vtkh/rendering/RayTracer.hpp
  vtkh/rendering/Renderer.hpp
  vtkh/rendering/Render.hpp
  vtkh/rendering/ScalarRenderer.hpp
  vtkh/rendering/Scene.hpp
  vtkh/rendering/VolumeRenderer.hpp

  vtkh/filters/Filter.hpp
  vtkh/filters/ParticleMerging.hpp
  )

set(sources
  Actor.cxx
  AxisAnnotation.cxx
  AxisAnnotation2D.cxx
  AxisAnnotation3D.cxx
  BitmapFont.cxx
  BitmapFontFactory.cxx
  BoundingBoxAnnotation.cxx
  Camera.cxx
  Color.cxx
  ColorBarAnnotation.cxx
  ColorLegendAnnotation.cxx
  LineRenderer.cxx
  Mapper.cxx
  MapperConnectivity.cxx
  Scene.cxx
  TextAnnotation.cxx
  TextAnnotationBillboard.cxx
  TextAnnotationScreen.cxx
  TextRenderer.cxx
  View.cxx
  View1D.cxx
  View2D.cxx
  View3D.cxx
  WorldAnnotator.cxx

  raytracing/Logger.cxx
  raytracing/MeshConnectivityContainers.cxx
  raytracing/TriangleExtractor.cxx

  ##vtk-h sources


  vtkh/compositing/Compositor.cpp
  vtkh/compositing/DirectSendCompositor.cpp
  vtkh/compositing/Image.cpp
  vtkh/compositing/PartialCompositor.cpp
  vtkh/compositing/PayloadCompositor.cpp
  vtkh/compositing/PayloadImage.cpp
  vtkh/compositing/PNGEncoder.cxx
  vtkh/compositing/RadixKCompositor.cpp

  vtkh/rendering/Annotator.cpp
  vtkh/rendering/LineRenderer.cpp
  vtkh/rendering/MeshRenderer.cpp
  vtkh/rendering/PointRenderer.cpp
  vtkh/rendering/RayTracer.cpp
  vtkh/rendering/Render.cpp
  vtkh/rendering/Renderer.cpp
  vtkh/rendering/ScalarRenderer.cpp
  vtkh/rendering/Scene.cpp

  vtkh/filters/Filter.cpp
  )

# This list of sources has code that uses devices and so might need to be
# compiled with a device-specific compiler (like CUDA).
set(device_sources
  Canvas.cxx
  CanvasRayTracer.cxx
  ConnectivityProxy.cxx
  LineRendererBatcher.cxx
  MapperCylinder.cxx
  MapperGlyphBase.cxx
  MapperGlyphScalar.cxx
  MapperGlyphVector.cxx
  MapperPoint.cxx
  MapperQuad.cxx
  MapperRayTracer.cxx
  MapperVolume.cxx
  MapperWireframer.cxx
  ScalarRenderer.cxx
  TextRendererBatcher.cxx

  internal/RunTriangulator.cxx

  raytracing/BoundingVolumeHierarchy.cxx
  raytracing/Camera.cxx
  raytracing/ChannelBuffer.cxx
  raytracing/ConnectivityTracer.cxx
  raytracing/CylinderExtractor.cxx
  raytracing/CylinderIntersector.cxx
  raytracing/GlyphExtractor.cxx
  raytracing/GlyphExtractorVector.cxx
  raytracing/GlyphIntersector.cxx
  raytracing/GlyphIntersectorVector.cxx
  raytracing/MeshConnectivityBuilder.cxx
  raytracing/QuadExtractor.cxx
  raytracing/QuadIntersector.cxx
  raytracing/RayOperations.cxx
  raytracing/RayTracer.cxx
  raytracing/ScalarRenderer.cxx
  raytracing/ShapeIntersector.cxx
  raytracing/SphereExtractor.cxx
  raytracing/SphereIntersector.cxx
  raytracing/TriangleIntersector.cxx
  raytracing/VolumeRendererStructured.cxx


  vtkh/DataSet.cpp

  vtkh/utils/vtkm_dataset_info.cpp

  vtkh/rendering/VolumeRenderer.cpp

  vtkh/filters/ParticleMerging.cpp
  )

# the None backend supports not building the opengl version

#-----------------------------------------------------------------------------
vtkm_library(
  NAME vtkm_rendering
  SOURCES ${sources}
  HEADERS ${headers}
  DEVICE_SOURCES ${device_sources}
  )

#-----------------------------------------------------------------------------
if(UNIX AND NOT APPLE)
  target_link_libraries(vtkm_rendering PRIVATE rt)
endif()

#-----------------------------------------------------------------------------
add_subdirectory(internal)
add_subdirectory(raytracing)

#-----------------------------------------------------------------------------
# Pretty much obsolete, but might be expected to be exported.
include(VTKmRenderingContexts)
