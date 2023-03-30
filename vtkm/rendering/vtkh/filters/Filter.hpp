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
//#include <vtkh/vtkh.hpp>
#include <vtkm/rendering/vtkh/DataSet.hpp>
#include <vtkm/filter/FieldSelection.h>

namespace vtkh
{

class VTKM_RENDERING_EXPORT Filter
{
public:
  Filter();
  virtual ~Filter();
  virtual void SetInput(DataSet *input);
  virtual std::string GetName() const = 0;

  DataSet* GetOutput();
  DataSet* Update();

  void AddMapField(const std::string &field_name);

  void ClearMapFields();

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

  DataSet *m_input;
  DataSet *m_output;

  void MapAllFields();

  void PropagateMetadata();

  void CheckForRequiredField(const std::string &field_name);
};

} //namespace vtkh
#endif
