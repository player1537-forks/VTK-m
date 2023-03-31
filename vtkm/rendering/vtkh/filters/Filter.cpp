//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/vtkh/filters/Filter.hpp>
#include <vtkm/rendering/vtkh/utils/vtkm_array_utils.hpp>

namespace vtkh
{

vtkm::cont::PartitionedDataSet*
Filter::Update()
{
  //DRP: Logger
  //VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  //DRP: Logger
  //VTKH_DATA_ADD("device", GetCurrentDevice());
  long long int in_cells = this->m_input->GetNumberOfCells();
  //DRP: Logger
  //VTKH_DATA_ADD("input_cells", in_cells);
  //VTKH_DATA_ADD("input_domains", this->m_input->GetNumberOfDomains());
  int in_topo_dims;
  bool in_structured = this->m_input->IsStructured(in_topo_dims);
  if(in_structured)
  {
    //DRP: Logger
    //VTKH_DATA_ADD("in_topology", "structured");
  }
  else
  {
    //DRP: Logger
    //VTKH_DATA_ADD("in_topology", "unstructured");
  }
#endif
  PreExecute();
  DoExecute();
  PostExecute();
#ifdef VTKH_ENABLE_LOGGING
  long long int out_cells = this->m_output->GetNumberOfCells();
  //DRP: Logger
  //VTKH_DATA_ADD("output_cells", out_cells);
  //VTKH_DATA_ADD("output_domains", this->m_output->GetNumberOfDomains());
  int out_topo_dims;
  bool out_structured = this->m_output->IsStructured(out_topo_dims);

  if(out_structured)
  {
    //DRP: Logger
    //VTKH_DATA_ADD("output_topology", "structured");
  }
  else
  {
    //DRP: Logger
    //VTKH_DATA_ADD("output_topology", "unstructured");
  }
#endif
  //DRP: Logger
  //VTKH_DATA_CLOSE();
  return m_output;
}

void
Filter::PreExecute()
{
  if(this->m_input == nullptr)
  {
    std::stringstream msg;
    msg<<"Input for vtkh filter '"<<this->GetName()<<"' is null.";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }

  if(m_map_fields.size() == 0)
  {
    this->MapAllFields();
  }

};

void
Filter::PostExecute()
{
  this->PropagateMetadata();
};

void
Filter::MapAllFields()
{
  if (this->m_input->GetNumberOfPartitions() > 0)
  {
    vtkm::cont::DataSet dom = this->m_input->GetPartition(0);
    vtkm::IdComponent num_fields = dom.GetNumberOfFields();
    for(vtkm::IdComponent i = 0; i < num_fields; ++i)
    {
      std::string field_name = dom.GetField(i).GetName();
      m_map_fields.push_back(field_name);
    }
  }
}

void
Filter::CheckForRequiredField(const std::string &field_name)
{
  if(this->m_input == nullptr)
  {
    std::stringstream msg;
    msg<<"Cannot verify required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' because input is null.";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }

  if (!vtkh::GlobalHasField(*this->m_input, field_name))
  {
    std::stringstream msg;
    msg<<"Required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' does not exist";
    throw vtkm::cont::ErrorBadValue(msg.str());
  }
}

void
Filter::PropagateMetadata()
{
  vtkm::Id numFields = this->m_input->GetNumberOfFields();
  for (vtkm::Id i = 0; i < numFields; i++)
    this->m_output->AddField(this->m_input->GetField(i));
}

vtkm::filter::FieldSelection
Filter::GetFieldSelection() const
{
  vtkm::filter::FieldSelection sel;
  for (const auto& str : this->m_map_fields)
  {
    sel.AddField(str);
  }
  return sel;
}


} //namespace vtkh
