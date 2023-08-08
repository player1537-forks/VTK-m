##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

vtkm_declare_headers(${headers})

#set(unit_tests
#)

#vtkm_unit_tests(SOURCES ${unit_tests})

#add distributed tests i.e.test to run with MPI
#if MPI is enabled.
if (VTKm_ENABLE_MPI)
  set(mpi_unit_tests
    UnitTestRendering_par.cxx
  )
  vtkm_unit_tests(
    MPI
    DEVICE_SOURCES ${mpi_unit_tests}
    LIBRARIES vtkm_source vtkm_filter_contour
    USE_VTKM_JOB_POOL
  )
endif()
