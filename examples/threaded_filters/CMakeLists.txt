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
project(ThreadedFilters CXX)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)

add_executable(ThreadedFilters ThreadedFilters.cxx)
target_link_libraries(ThreadedFilters PRIVATE vtkm_filter vtkm_io)
vtkm_add_target_information(ThreadedFilters
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES ThreadedFilters.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(ThreadedFilters PRIVATE BUILDING_TBB_VERSION)
endif()
