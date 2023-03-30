//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_HPP
#define VTK_H_HPP

#include <vtkm/rendering/vtkm_rendering_export.h>
//#include <vtkh/vtkh_config.h>
#include <string>

namespace vtkh
{
#if 0
  VTKM_RENDERING_EXPORT std::string AboutVTKH();
  VTKM_RENDERING_EXPORT void        Initialize();

  // is backend support compiled in
  VTKM_RENDERING_EXPORT bool        IsSerialAvailable();
  VTKM_RENDERING_EXPORT bool        IsOpenMPAvailable();
  VTKM_RENDERING_EXPORT bool        IsCUDAAvailable();
  VTKM_RENDERING_EXPORT bool        IsKokkosAvailable();

  // is backend enabled (e.g., ForceX)
  VTKM_RENDERING_EXPORT bool        IsSerialEnabled();
  VTKM_RENDERING_EXPORT bool        IsOpenMPEnabled();
  VTKM_RENDERING_EXPORT bool        IsCUDAEnabled();
  VTKM_RENDERING_EXPORT bool        IsKokkosEnabled();


  VTKM_RENDERING_EXPORT bool        IsMPIEnabled();

  VTKM_RENDERING_EXPORT int         CUDADeviceCount();
  VTKM_RENDERING_EXPORT void        SelectCUDADevice(int device_index);
  VTKM_RENDERING_EXPORT void        InitializeKokkos();
  VTKM_RENDERING_EXPORT void        SelectKokkosDevice(int device_index);

  VTKM_RENDERING_EXPORT void        ForceSerial();
  VTKM_RENDERING_EXPORT void        ForceOpenMP();
  VTKM_RENDERING_EXPORT void        ForceCUDA();
  VTKM_RENDERING_EXPORT void        ForceKokkos();
  VTKM_RENDERING_EXPORT void        ResetDevices();
  VTKM_RENDERING_EXPORT std::string GetCurrentDevice();

  VTKM_RENDERING_EXPORT int         GetMPIRank();
  VTKM_RENDERING_EXPORT int         GetMPISize();

  VTKM_RENDERING_EXPORT void        SetMPICommHandle(int mpi_comm_id);
  VTKM_RENDERING_EXPORT int         GetMPICommHandle();
#endif
}
#endif
