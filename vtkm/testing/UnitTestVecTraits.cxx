//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#include <vtkm/testing/VecTraitsTests.h>

#include <vtkm/testing/Testing.h>

namespace
{

static constexpr vtkm::Id MAX_VECTOR_SIZE = 5;
static constexpr vtkm::Id VecInit[MAX_VECTOR_SIZE] = { 42, 54, 67, 12, 78 };

struct TypeWithoutVecTraits
{
  vtkm::Id Value = -1;

  TypeWithoutVecTraits() = default;

  TypeWithoutVecTraits(vtkm::Id value)
    : Value(value)
  {
  }

  operator vtkm::Id() const { return this->Value; }
};

struct TestVecTypeFunctor
{
  template <typename T>
  void operator()(const T&) const
  {
    using Traits = vtkm::VecTraits<T>;
    using ComponentType = typename Traits::ComponentType;
    VTKM_TEST_ASSERT(Traits::NUM_COMPONENTS <= MAX_VECTOR_SIZE,
                     "Need to update test for larger vectors.");
    T inVector;
    for (vtkm::IdComponent index = 0; index < Traits::NUM_COMPONENTS; index++)
    {
      Traits::SetComponent(inVector, index, ComponentType(VecInit[index]));
    }
    T outVector;
    vtkm::testing::TestVecType<Traits::NUM_COMPONENTS>(inVector, outVector);
    vtkm::VecC<ComponentType> outVecC(outVector);
    vtkm::testing::TestVecType<Traits::NUM_COMPONENTS>(vtkm::VecC<ComponentType>(inVector),
                                                       outVecC);
    vtkm::VecCConst<ComponentType> outVecCConst(outVector);
    vtkm::testing::TestVecType<Traits::NUM_COMPONENTS>(vtkm::VecCConst<ComponentType>(inVector),
                                                       outVecCConst);
  }
};

void TestVecTraits()
{
  TestVecTypeFunctor test;
  vtkm::testing::Testing::TryTypes(test);
  std::cout << "vtkm::Vec<vtkm::FloatDefault, 5>" << std::endl;
  test(vtkm::Vec<vtkm::FloatDefault, 5>());
  std::cout << "TypeWithoutVecTraits" << std::endl;
  test(TypeWithoutVecTraits{});

  vtkm::testing::TestVecComponentsTag<vtkm::Id3>();
  vtkm::testing::TestVecComponentsTag<vtkm::Vec3f>();
  vtkm::testing::TestVecComponentsTag<vtkm::Vec4f>();
  vtkm::testing::TestVecComponentsTag<vtkm::VecC<vtkm::FloatDefault>>();
  vtkm::testing::TestVecComponentsTag<vtkm::VecCConst<vtkm::Id>>();
  vtkm::testing::TestScalarComponentsTag<vtkm::Id>();
  vtkm::testing::TestScalarComponentsTag<vtkm::FloatDefault>();
  vtkm::testing::TestScalarComponentsTag<TypeWithoutVecTraits>();
}

} // anonymous namespace

int UnitTestVecTraits(int argc, char* argv[])
{
  return vtkm::testing::Testing::Run(TestVecTraits, argc, argv);
}
