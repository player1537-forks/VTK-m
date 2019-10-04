//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.md for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_c_internal_Config_h
#define vtk_c_internal_Config_h

#include <cstdint>
#include <type_traits>

#ifdef __CUDACC__
# define VTKC_EXEC __device__ __host__
#else
# define VTKC_EXEC
#endif

namespace vtkc
{

namespace internal
{
template <typename T>
using ClosestFloatType =
  typename std::enable_if<std::is_arithmetic<T>::value,
                          typename std::conditional<sizeof(T) <= 4, float, double>::type>::type;
}

using IdShape = std::int8_t;
using IdComponent = std::int32_t;

} // namespace vtkc

#endif // vtk_c_internal_Config_h