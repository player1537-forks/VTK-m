##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

set(unit_tests
  UnitTestSerialArrayHandle.cxx
  UnitTestSerialArrayHandleFancy.cxx
  UnitTestSerialArrayHandleMultiplexer.cxx
  UnitTestSerialBitField.cxx
  UnitTestSerialCellLocatorRectilinearGrid.cxx
  UnitTestSerialCellLocatorTwoLevel.cxx
  UnitTestSerialCellLocatorUniformGrid.cxx
  UnitTestSerialComputeRange.cxx
  UnitTestSerialColorTable.cxx
  UnitTestSerialDataSetExplicit.cxx
  UnitTestSerialDataSetSingleType.cxx
  UnitTestSerialDeviceAdapter.cxx
  UnitTestSerialImplicitFunction.cxx
  UnitTestSerialPointLocatorSparseGrid.cxx
  )

if (NOT VTKm_NO_DEPRECATED_VIRTUAL)
  set(unit_tests ${unit_tests}
    UnitTestSerialVirtualObjectHandle.cxx
    )
endif()

# these cannot be compiled with HIPCC
if(HIP IN_LIST Kokkos_DEVICES)
  list(REMOVE_ITEM unit_tests
       UnitTestSerialArrayHandleFancy.cxx
       UnitTestSerialCellLocatorRectilinearGrid.cxx
       UnitTestSerialCellLocatorTwoLevel.cxx
       UnitTestSerialCellLocatorUniformGrid.cxx
       UnitTestSerialDataSetExplicit.cxx
       UnitTestSerialImplicitFunction.cxx
       UnitTestSerialPointLocatorSparseGrid.cxx
       )
endif()


vtkm_unit_tests(SOURCES ${unit_tests}
                LABEL "SERIAL"
                DEFINES VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG
                LIBRARIES vtkm_worklet)
