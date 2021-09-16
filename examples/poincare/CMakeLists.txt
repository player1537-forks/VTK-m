##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(Poincare CXX)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
find_package(ADIOS2 REQUIRED)
find_package(Fides REQUIRED)

add_executable(Poincare Poincare.cxx io.cxx io.h)
target_link_libraries(Poincare PRIVATE vtkm_filter vtkm_io fides)
vtkm_add_target_information(Poincare
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES Poincare.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(Poincare PRIVATE BUILDING_TBB_VERSION)
endif()
