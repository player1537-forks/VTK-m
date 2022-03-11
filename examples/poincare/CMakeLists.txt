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

add_executable(Simple2.3 Simple2.3.cxx)
target_link_libraries(Simple2.3 PRIVATE vtkm_filter vtkm_io adios2::adios2)

vtkm_add_target_information(Simple2.3
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES Simple2.3.cxx)

if(TARGET vtkm::tbb)
  target_compile_definitions(Poincare PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple2 PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple3 PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple4 PRIVATE BUILDING_TBB_VERSION)
endif()
