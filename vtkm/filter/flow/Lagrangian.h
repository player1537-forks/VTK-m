//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_flow_Lagrangian_h
#define vtk_m_filter_flow_Lagrangian_h

#include <vtkm/filter/NewFilterField.h>
#include <vtkm/filter/flow/vtkm_filter_flow_export.h>

namespace vtkm
{
namespace filter
{
namespace flow
{

class VTKM_FILTER_FLOW_EXPORT Lagrangian : public vtkm::filter::NewFilterField
{
public:
  VTKM_CONT
  Lagrangian();

  VTKM_CONT
  bool CanThread() const override { return false; }

  VTKM_CONT
  void SetInitFlag(bool val) { this->initFlag = val; }

  VTKM_CONT
  void SetExtractFlows(bool val) { this->extractFlows = val; }

  VTKM_CONT
  void SetResetParticles(bool val) { this->resetParticles = val; }

  VTKM_CONT
  void SetStepSize(vtkm::Float32 val) { this->stepSize = val; }

  VTKM_CONT
  void SetWriteFrequency(vtkm::Id val) { this->writeFrequency = val; }

  VTKM_CONT
  void SetSeedResolutionInX(vtkm::Id val) { this->x_res = val; }

  VTKM_CONT
  void SetSeedResolutionInY(vtkm::Id val) { this->y_res = val; }

  VTKM_CONT
  void SetSeedResolutionInZ(vtkm::Id val) { this->z_res = val; }

  VTKM_CONT
  void SetCustomSeedResolution(vtkm::Id val) { this->cust_res = val; }

  VTKM_CONT
  void SetSeedingResolution(vtkm::Id3 val) { this->SeedRes = val; }

  VTKM_CONT
  void UpdateSeedResolution(vtkm::cont::DataSet input);

  VTKM_CONT
  void InitializeSeedPositions(const vtkm::cont::DataSet& input);

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& inData) override;

  VTKM_CONT
  void InitializeCoordinates(const vtkm::cont::DataSet& input,
                             std::vector<Float64>& xC,
                             std::vector<Float64>& yC,
                             std::vector<Float64>& zC);

  bool initFlag = true;
  bool extractFlows = false;
  bool resetParticles = true;
  vtkm::Float32 stepSize;
  vtkm::Id x_res = 0, y_res = 0, z_res = 0;
  vtkm::Id cust_res = 0;
  vtkm::Id3 SeedRes = { 1, 1, 1 };
  vtkm::Id writeFrequency = 0;
};

}
}
} //vtkm::filter::flow

#endif // #define vtk_m_filter_flow_Lagrangian_h
