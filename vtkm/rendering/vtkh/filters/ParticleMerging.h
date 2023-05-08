//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_PARTICLE_MERGING_HPP
#define VTK_H_PARTICLE_MERGING_HPP

#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh
{

class VTKM_RENDERING_EXPORT ParticleMerging
{
public:
  ParticleMerging();
  virtual ~ParticleMerging();
  std::string GetName() const { return "vtkh::ParticleMerging"; }
  void SetField(const std::string& field_name);
  void SetRadius(const vtkm::Float64 radius);

  virtual void SetInput(vtkm::cont::PartitionedDataSet* input) { this->m_input = input; }
  vtkm::cont::PartitionedDataSet* GetOutput() { return this->m_output; }

  vtkm::cont::PartitionedDataSet* Update()
  {
    this->PreExecute();
    this->DoExecute();
    this->PostExecute();
    return this->m_output;
  }

protected:
  void PreExecute();
  void PostExecute();
  void DoExecute();

  std::string m_field_name;
  vtkm::Float64 m_radius;

  vtkm::cont::PartitionedDataSet* m_input = nullptr;
  vtkm::cont::PartitionedDataSet* m_output = nullptr;
};

} //namespace vtkh
#endif
