//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/Types.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ErrorFilterExecution.h>

#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/filter/flow/worklet/Field.h>
#include <vtkm/filter/flow/worklet/GridEvaluators.h>
#include <vtkm/filter/flow/worklet/ParticleAdvection.h>
#include <vtkm/filter/flow/worklet/RK4Integrator.h>
#include <vtkm/filter/flow/worklet/Stepper.h>

#include <vtkm/filter/flow/Lagrangian_deprecated.h>

#include <cstring>
#include <sstream>
#include <string.h>

static vtkm::Id cycle = 0;
static vtkm::cont::ArrayHandle<vtkm::Particle> BasisParticles;
static vtkm::cont::ArrayHandle<vtkm::Particle> BasisParticlesOriginal;
static vtkm::cont::ArrayHandle<vtkm::Id> BasisParticlesValidity;

namespace
{
class ValidityCheck : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn end_point, FieldInOut output);
  using ExecutionSignature = void(_1, _2);
  using InputDomain = _1;

  ValidityCheck(vtkm::Bounds b)
    : bounds(b)
  {
  }

  template <typename ValidityType>
  VTKM_EXEC void operator()(const vtkm::Particle& end_point, ValidityType& res) const
  {
    vtkm::Id steps = end_point.NumSteps;
    if (steps > 0 && res == 1)
    {
      if (bounds.Contains(end_point.Pos))
      {
        res = 1;
      }
      else
      {
        res = 0;
      }
    }
    else
    {
      res = 0;
    }
  }

private:
  vtkm::Bounds bounds;
};

class DisplacementCalculation : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn end_point, FieldIn start_point, FieldInOut output);
  using ExecutionSignature = void(_1, _2, _3);
  using InputDomain = _1;

  template <typename DisplacementType>
  VTKM_EXEC void operator()(const vtkm::Particle& end_point,
                            const vtkm::Particle& start_point,
                            DisplacementType& res) const
  {
    res[0] = end_point.Pos[0] - start_point.Pos[0];
    res[1] = end_point.Pos[1] - start_point.Pos[1];
    res[2] = end_point.Pos[2] - start_point.Pos[2];
  }
};
}

namespace vtkm
{
namespace filter
{

//-----------------------------------------------------------------------------
VTKM_CONT Lagrangian_deprecated::Lagrangian_deprecated()
  : initFlag(true)
  , extractFlows(false)
  , resetParticles(true)
  , stepSize(1.0f)
  , x_res(0)
  , y_res(0)
  , z_res(0)
  , cust_res(0)
  , SeedRes(vtkm::Id3(1, 1, 1))
  , writeFrequency(0)
{
}

VTKM_CONT
void Lagrangian_deprecated::SetCycle(vtkm::Id c)
{
  cycle = c;
}

VTKM_CONT
vtkm::Id Lagrangian_deprecated::GetCycle() const
{
  return cycle;
}


VTKM_CONT
void Lagrangian_deprecated::SetBasisParticles(
  const vtkm::cont::ArrayHandle<vtkm::Particle>& basisParticles)
{
  BasisParticles = basisParticles;
}

VTKM_CONT
vtkm::cont::ArrayHandle<vtkm::Particle> Lagrangian_deprecated::GetBasisParticles() const
{
  return BasisParticles;
}

VTKM_CONT
void Lagrangian_deprecated::SetBasisParticlesOriginal(
  const vtkm::cont::ArrayHandle<vtkm::Particle>& basisParticles)
{
  BasisParticlesOriginal = basisParticles;
}

VTKM_CONT
vtkm::cont::ArrayHandle<vtkm::Particle> Lagrangian_deprecated::GetBasisParticlesOriginal() const
{
  return BasisParticlesOriginal;
}

VTKM_CONT
void Lagrangian_deprecated::SetBasisParticleValidity(
  const vtkm::cont::ArrayHandle<vtkm::Id>& validity)
{
  BasisParticlesValidity = validity;
}

VTKM_CONT
vtkm::cont::ArrayHandle<vtkm::Id> Lagrangian_deprecated::GetBasisParticleValidity() const
{
  return BasisParticlesValidity;
}


//-----------------------------------------------------------------------------
void Lagrangian_deprecated::UpdateSeedResolution(const vtkm::cont::DataSet input)
{
  vtkm::cont::UnknownCellSet cell_set = input.GetCellSet();

  if (cell_set.CanConvert<vtkm::cont::CellSetStructured<1>>())
  {
    vtkm::cont::CellSetStructured<1> cell_set1 =
      cell_set.AsCellSet<vtkm::cont::CellSetStructured<1>>();
    vtkm::Id dims1 = cell_set1.GetPointDimensions();
    this->SeedRes[0] = dims1;
    if (this->cust_res)
    {
      this->SeedRes[0] = dims1 / this->x_res;
    }
  }
  else if (cell_set.CanConvert<vtkm::cont::CellSetStructured<2>>())
  {
    vtkm::cont::CellSetStructured<2> cell_set2 =
      cell_set.AsCellSet<vtkm::cont::CellSetStructured<2>>();
    vtkm::Id2 dims2 = cell_set2.GetPointDimensions();
    this->SeedRes[0] = dims2[0];
    this->SeedRes[1] = dims2[1];
    if (this->cust_res)
    {
      this->SeedRes[0] = dims2[0] / this->x_res;
      this->SeedRes[1] = dims2[1] / this->y_res;
    }
  }
  else if (cell_set.CanConvert<vtkm::cont::CellSetStructured<3>>())
  {
    vtkm::cont::CellSetStructured<3> cell_set3 =
      cell_set.AsCellSet<vtkm::cont::CellSetStructured<3>>();
    vtkm::Id3 dims3 = cell_set3.GetPointDimensions();
    this->SeedRes[0] = dims3[0];
    this->SeedRes[1] = dims3[1];
    this->SeedRes[2] = dims3[2];
    if (this->cust_res)
    {
      this->SeedRes[0] = dims3[0] / this->x_res;
      this->SeedRes[1] = dims3[1] / this->y_res;
      this->SeedRes[2] = dims3[2] / this->z_res;
    }
  }
}


//-----------------------------------------------------------------------------
void Lagrangian_deprecated::InitializeSeedPositions(const vtkm::cont::DataSet& input)
{
  vtkm::Bounds bounds = input.GetCoordinateSystem().GetBounds();

  Lagrangian_deprecated::UpdateSeedResolution(input);

  vtkm::Float64 x_spacing = 0.0, y_spacing = 0.0, z_spacing = 0.0;
  if (this->SeedRes[0] > 1)
    x_spacing = (double)(bounds.X.Max - bounds.X.Min) / (double)(this->SeedRes[0] - 1);
  if (this->SeedRes[1] > 1)
    y_spacing = (double)(bounds.Y.Max - bounds.Y.Min) / (double)(this->SeedRes[1] - 1);
  if (this->SeedRes[2] > 1)
    z_spacing = (double)(bounds.Z.Max - bounds.Z.Min) / (double)(this->SeedRes[2] - 1);
  // Divide by zero handling for 2D data set. How is this handled

  BasisParticles.Allocate(this->SeedRes[0] * this->SeedRes[1] * this->SeedRes[2]);
  BasisParticlesValidity.Allocate(this->SeedRes[0] * this->SeedRes[1] * this->SeedRes[2]);

  auto portal1 = BasisParticles.WritePortal();
  auto portal2 = BasisParticlesValidity.WritePortal();

  vtkm::Id id = 0;
  for (int z = 0; z < this->SeedRes[2]; z++)
  {
    vtkm::FloatDefault zi = static_cast<vtkm::FloatDefault>(z * z_spacing);
    for (int y = 0; y < this->SeedRes[1]; y++)
    {
      vtkm::FloatDefault yi = static_cast<vtkm::FloatDefault>(y * y_spacing);
      for (int x = 0; x < this->SeedRes[0]; x++)
      {
        vtkm::FloatDefault xi = static_cast<vtkm::FloatDefault>(x * x_spacing);
        portal1.Set(id,
                    vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(bounds.X.Min) + xi,
                                         static_cast<vtkm::FloatDefault>(bounds.Y.Min) + yi,
                                         static_cast<vtkm::FloatDefault>(bounds.Z.Min) + zi),
                                   id));
        portal2.Set(id, 1);
        id++;
      }
    }
  }
}

//-----------------------------------------------------------------------------
void Lagrangian_deprecated::InitializeCoordinates(const vtkm::cont::DataSet& input,
                                                  std::vector<Float64>& xC,
                                                  std::vector<Float64>& yC,
                                                  std::vector<Float64>& zC)
{
  vtkm::Bounds bounds = input.GetCoordinateSystem().GetBounds();

  vtkm::Float64 x_spacing = 0.0, y_spacing = 0.0, z_spacing = 0.0;
  if (this->SeedRes[0] > 1)
    x_spacing = (double)(bounds.X.Max - bounds.X.Min) / (double)(this->SeedRes[0] - 1);
  if (this->SeedRes[1] > 1)
    y_spacing = (double)(bounds.Y.Max - bounds.Y.Min) / (double)(this->SeedRes[1] - 1);
  if (this->SeedRes[2] > 1)
    z_spacing = (double)(bounds.Z.Max - bounds.Z.Min) / (double)(this->SeedRes[2] - 1);
  // Divide by zero handling for 2D data set. How is this handled

  for (int x = 0; x < this->SeedRes[0]; x++)
  {
    vtkm::FloatDefault xi = static_cast<vtkm::FloatDefault>(x * x_spacing);
    xC.push_back(bounds.X.Min + xi);
  }
  for (int y = 0; y < this->SeedRes[1]; y++)
  {
    vtkm::FloatDefault yi = static_cast<vtkm::FloatDefault>(y * y_spacing);
    yC.push_back(bounds.Y.Min + yi);
  }
  for (int z = 0; z < this->SeedRes[2]; z++)
  {
    vtkm::FloatDefault zi = static_cast<vtkm::FloatDefault>(z * z_spacing);
    zC.push_back(bounds.Z.Min + zi);
  }
}

//-----------------------------------------------------------------------------
VTKM_CONT vtkm::cont::DataSet Lagrangian_deprecated::DoExecute(const vtkm::cont::DataSet& input)
{
  if (cycle == 0)
  {
    InitializeSeedPositions(input);
    BasisParticlesOriginal.Allocate(this->SeedRes[0] * this->SeedRes[1] * this->SeedRes[2]);
    vtkm::cont::ArrayCopy(BasisParticles, BasisParticlesOriginal);
  }

  if (this->writeFrequency == 0)
  {
    throw vtkm::cont::ErrorFilterExecution(
      "Write frequency can not be 0. Use SetWriteFrequency().");
  }
  vtkm::cont::ArrayHandle<vtkm::Particle> basisParticleArray;
  vtkm::cont::ArrayCopy(BasisParticles, basisParticleArray);

  cycle += 1;
  const vtkm::cont::UnknownCellSet& cells = input.GetCellSet();
  const vtkm::cont::CoordinateSystem& coords =
    input.GetCoordinateSystem(this->GetActiveCoordinateSystemIndex());
  vtkm::Bounds bounds = input.GetCoordinateSystem().GetBounds();

  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::flow::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::flow::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::flow::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::flow::Stepper<RK4Type, GridEvalType>;

  vtkm::worklet::flow::ParticleAdvection particleadvection;
  vtkm::worklet::flow::ParticleAdvectionResult<vtkm::Particle> res;

  const auto field = input.GetField(this->GetActiveFieldName());
  FieldType velocities(field.GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::Vec3f>>(),
                       field.GetAssociation());

  GridEvalType gridEval(coords, cells, velocities);
  Stepper rk4(gridEval, static_cast<vtkm::Float32>(this->stepSize));

  res = particleadvection.Run(rk4, basisParticleArray, 1); // Taking a single step
  auto particles = res.Particles;

  vtkm::cont::DataSet outputData;
  vtkm::cont::DataSetBuilderRectilinear dataSetBuilder;

  if (cycle % this->writeFrequency == 0)
  {
    /* Steps to create a structured dataset */
    UpdateSeedResolution(input);
    vtkm::cont::ArrayHandle<vtkm::Vec3f> BasisParticlesDisplacement;
    BasisParticlesDisplacement.Allocate(this->SeedRes[0] * this->SeedRes[1] * this->SeedRes[2]);
    DisplacementCalculation displacement;
    this->Invoke(displacement, particles, BasisParticlesOriginal, BasisParticlesDisplacement);
    std::vector<Float64> xC, yC, zC;
    InitializeCoordinates(input, xC, yC, zC);
    outputData = dataSetBuilder.Create(xC, yC, zC);
    outputData.AddPointField("valid", BasisParticlesValidity);
    outputData.AddPointField("displacement", BasisParticlesDisplacement);

    if (this->resetParticles)
    {
      InitializeSeedPositions(input);
      BasisParticlesOriginal.Allocate(this->SeedRes[0] * this->SeedRes[1] * this->SeedRes[2]);
      vtkm::cont::ArrayCopy(BasisParticles, BasisParticlesOriginal);
    }
    else
    {
      vtkm::cont::ArrayCopy(particles, BasisParticles);
    }
  }
  else
  {
    ValidityCheck check(bounds);
    this->Invoke(check, particles, BasisParticlesValidity);
    vtkm::cont::ArrayCopy(particles, BasisParticles);
  }

  return outputData;
}

}
} // namespace vtkm::filter
