//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2017 Sandia Corporation.
//  Copyright 2017 UT-Battelle, LLC.
//  Copyright 2017 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//
//=============================================================================
#include <vtkm/worklet/ScatterPermutation.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/WorkletMapTopology.h>

#include <ctime>
#include <random>

namespace
{

class Worklet : public vtkm::worklet::WorkletMapCellToPoint
{
public:
  typedef void ControlSignature(CellSetIn cellset,
                                FieldOutPoint<IdType> outPointId,
                                FieldOutPoint<IdComponentType> outVisit);
  typedef void ExecutionSignature(InputIndex, VisitIndex, _2, _3);
  using InputDomain = _1;

  using ScatterType = vtkm::worklet::ScatterPermutation<>;

  Worklet(const vtkm::cont::ArrayHandle<vtkm::Id>& permutation)
    : Scatter(permutation)
  {
  }

  VTKM_EXEC void operator()(vtkm::Id pointId,
                            vtkm::IdComponent visit,
                            vtkm::Id& outPointId,
                            vtkm::IdComponent& outVisit) const
  {
    outPointId = pointId;
    outVisit = visit;
  }

  VTKM_CONT
  ScatterType GetScatter() const { return this->Scatter; }

private:
  ScatterType Scatter;
};

template <typename CellSetType>
void RunTest(const CellSetType& cellset, const vtkm::cont::ArrayHandle<vtkm::Id>& permutation)
{
  vtkm::cont::ArrayHandle<vtkm::Id> outPointId;
  vtkm::cont::ArrayHandle<vtkm::IdComponent> outVisit;

  Worklet worklet(permutation);
  vtkm::worklet::DispatcherMapTopology<Worklet>(worklet).Invoke(cellset, outPointId, outVisit);

  for (vtkm::Id i = 0; i < permutation.GetNumberOfValues(); ++i)
  {
    VTKM_TEST_ASSERT(outPointId.GetPortalConstControl().Get(i) ==
                       permutation.GetPortalConstControl().Get(i),
                     "output point ids do not match the permutation");
    VTKM_TEST_ASSERT(outVisit.GetPortalConstControl().Get(i) == 0, "incorrect visit index");
  }
}

void TestScatterPermutation()
{
  vtkm::cont::DataSet dataset = vtkm::cont::testing::MakeTestDataSet().Make2DUniformDataSet0();
  auto cellset = dataset.GetCellSet();
  vtkm::Id numberOfPoints = cellset.GetNumberOfPoints();

  vtkm::UInt32 seed = static_cast<vtkm::UInt32>(std::time(nullptr));

  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<vtkm::Id> countDistribution(1, 2 * numberOfPoints);
  std::uniform_int_distribution<vtkm::Id> ptidDistribution(0, numberOfPoints - 1);

  const int iterations = 5;
  std::cout << "Testing with random permutations " << iterations << " times\n";
  std::cout << "Seed: " << seed << std::endl;
  for (int i = 1; i <= iterations; ++i)
  {
    std::cout << "iteration: " << i << "\n";

    vtkm::Id count = countDistribution(generator);
    vtkm::cont::ArrayHandle<vtkm::Id> permutation;
    permutation.Allocate(count);

    auto portal = permutation.GetPortalControl();
    std::cout << "using permutation:";
    for (vtkm::Id j = 0; j < count; ++j)
    {
      auto val = ptidDistribution(generator);
      std::cout << " " << val;
      portal.Set(j, val);
    }
    std::cout << "\n";

    RunTest(cellset, permutation);
  }
}

} // anonymous namespace

int UnitTestScatterPermutation(int, char* [])
{
  return vtkm::cont::testing::Testing::Run(TestScatterPermutation);
}
