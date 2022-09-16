//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_flow_Lagrangian_deprecated_h
#define vtk_m_filter_flow_Lagrangian_deprecated_h

#include <vtkm/filter/flow/Lagrangian.h>
#include <vtkm/filter/flow/vtkm_filter_flow_export.h>

namespace vtkm
{
namespace filter
{

class VTKM_FILTER_FLOW_EXPORT Lagrangian_deprecated : public vtkm::filter::flow::Lagrangian
{
private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& inData) override;
};

}
} //vtkm::filter

#endif // #define vtk_m_filter_Lagrangian_deprecated_h
