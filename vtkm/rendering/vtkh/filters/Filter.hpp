//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef VTK_H_FILTER_HPP
#define VTK_H_FILTER_HPP


#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/filter/FieldSelection.h>

namespace vtkh
{

class VTKM_RENDERING_EXPORT Filter
{
public:
virtual ~Filter() {}

  virtual std::string GetName() const = 0;
  virtual void SetInput(vtkm::cont::PartitionedDataSet *input) {this->m_input = input;}
  vtkm::cont::PartitionedDataSet* GetOutput() { return this->m_output; }

  vtkm::cont::PartitionedDataSet* Update();

  void AddMapField(const std::string &field_name) {this->m_map_fields.push_back(field_name);}
  void ClearMapFields() { this->m_map_fields.clear(); }

protected:
  virtual void DoExecute() = 0;
  virtual void PreExecute();
  virtual void PostExecute();

  //@{
  /// These are all temporary methods added to gets things building again
  /// while we totally deprecate vtk-h compnents
  ///
  vtkm::filter::FieldSelection GetFieldSelection() const;
  //@}

  std::vector<std::string> m_map_fields;

  vtkm::cont::PartitionedDataSet *m_input = nullptr;
  vtkm::cont::PartitionedDataSet *m_output = nullptr;

  void MapAllFields();

  void PropagateMetadata();

  void CheckForRequiredField(const std::string &field_name);
};

} //namespace vtkh
#endif
