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
    Annotator.hpp
    LineRenderer.hpp
    MeshRenderer.hpp
    RayTracer.hpp
    Render.hpp
    Renderer.hpp
    PointRenderer.hpp
    ScalarRenderer.hpp
    Scene.hpp
    VolumeRenderer.hpp
  )

set_source_files_properties(CellTables.h
                            PROPERTIES VTKm_CANT_BE_HEADER_TESTED TRUE)

vtkm_declare_headers(${headers})
