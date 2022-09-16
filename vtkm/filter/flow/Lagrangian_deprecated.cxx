//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/filter/flow/Lagrangian_deprecated.h>

namespace vtkm
{
namespace filter
{

static vtkm::Id global_Cycle = 0;
static vtkm::cont::ArrayHandle<vtkm::Particle> global_BasisParticles;
static vtkm::cont::ArrayHandle<vtkm::Particle> global_BasisParticlesOriginal;
static vtkm::cont::ArrayHandle<vtkm::Id> global_BasisParticlesValidity;


VTKM_CONT vtkm::cont::DataSet Lagrangian_deprecated::DoExecute(const vtkm::cont::DataSet& input)
{
  //Initialize the filter with the static variables
  this->SetCycle(global_Cycle);
  this->SetBasisParticles(global_BasisParticles);
  this->SetBasisParticlesOriginal(global_BasisParticlesOriginal);
  this->SetBasisParticleValidity(global_BasisParticlesValidity);

  auto output = vtkm::filter::flow::Lagrangian::DoExecute(input);

  //Set the static variables with the current values.
  global_Cycle = this->GetCycle();
  global_BasisParticles = this->GetBasisParticles();
  global_BasisParticlesOriginal = this->GetBasisParticlesOriginal();
  global_BasisParticlesValidity = this->GetBasisParticleValidity();

  return output;
}

}
} // namespace vtkm::filter
