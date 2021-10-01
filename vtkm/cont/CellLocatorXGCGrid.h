//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtkm_cont_CellLocatorXGCGrid_h
#define vtkm_cont_CellLocatorXGCGrid_h

#include <vtkm/cont/CellLocatorTwoLevel.h>
#include <vtkm/cont/internal/CellLocatorBase.h>

#include <vtkm/exec/CellLocatorXGCGrid.h>

namespace vtkm
{
namespace cont
{

class VTKM_CONT_EXPORT CellLocatorXGCGrid
  : public vtkm::cont::internal::CellLocatorBase<CellLocatorXGCGrid>
{
  using Superclass = vtkm::cont::internal::CellLocatorBase<CellLocatorXGCGrid>;

public:
  VTKM_CONT vtkm::exec::CellLocatorXGCGrid PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                                               vtkm::cont::Token& token) const;

private:
  vtkm::cont::CellLocatorTwoLevel CellLocator;
  vtkm::Id CellsPerPlane = 0;
  bool IsCylindrical = false;
  bool IsPeriodic = false;
  vtkm::Id NumPlanes = 0;

  friend Superclass;
  VTKM_CONT void Build();
};
}
} // vtkm::cont

#endif //vtkm_cont_CellLocatorXGCGrid_h
