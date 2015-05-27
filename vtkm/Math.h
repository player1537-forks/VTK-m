//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2015 Sandia Corporation.
//  Copyright 2015 UT-Battelle, LLC.
//  Copyright 2015 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//
//=============================================================================
// **** DO NOT EDIT THIS FILE!!! ****
// This file is automatically generated by FunctionInterfaceDetailPre.h.in

#ifndef vtk_m_Math_h
#define vtk_m_Math_h

#include <vtkm/Types.h>

#ifndef VTKM_CUDA
#include <math.h>
#endif

#define VTKM_SYS_MATH_FUNCTION_32(func) func ## f
#define VTKM_SYS_MATH_FUNCTION_64(func) func




namespace vtkm {

/// Compute the square root of \p x.
///
VTKM_EXEC_CONT_EXPORT
vtkm::Float32 Sqrt(vtkm::Float32 x) {
  return VTKM_SYS_MATH_FUNCTION_32(sqrt)(x);
}
VTKM_EXEC_CONT_EXPORT
vtkm::Float64 Sqrt(vtkm::Float64 x) {
  return VTKM_SYS_MATH_FUNCTION_64(sqrt)(x);
}
template<typename T, vtkm::IdComponent N>
vtkm::Vec<T,N> Sqrt(const vtkm::Vec<T,N> x) {
  vtkm::Vec<T,N> result;
  for (vtkm::IdComponent index = 0; index < N; index++)
  {
    result[index] = vtkm::Sqrt(x[index]);
  }
  return result;
}
template<typename T>
VTKM_EXEC_CONT_EXPORT
vtkm::Vec<T,4> Sqrt(const vtkm::Vec<T,4> x) {
  return vtkm::Vec<T,4>(vtkm::Sqrt(x[0]),
                        vtkm::Sqrt(x[1]),
                        vtkm::Sqrt(x[2]),
                        vtkm::Sqrt(x[3]));
}
template<typename T>
VTKM_EXEC_CONT_EXPORT
vtkm::Vec<T,3> Sqrt(const vtkm::Vec<T,3> x) {
  return vtkm::Vec<T,3>(vtkm::Sqrt(x[0]),
                        vtkm::Sqrt(x[1]),
                        vtkm::Sqrt(x[2]));
}
template<typename T>
VTKM_EXEC_CONT_EXPORT
vtkm::Vec<T,2> Sqrt(const vtkm::Vec<T,2> x) {
  return vtkm::Vec<T,2>(vtkm::Sqrt(x[0]),
                        vtkm::Sqrt(x[1]));
}


} // namespace vtkm

#endif //vtk_m_Math_h
