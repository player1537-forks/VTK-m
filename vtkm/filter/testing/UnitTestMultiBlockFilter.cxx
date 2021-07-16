//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/Math.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/cont/testing/Testing.h>

#include <vtkm/filter/CleanGrid.h>
#include <vtkm/filter/ClipWithField.h>
#include <vtkm/filter/Contour.h>
#include <vtkm/filter/Gradient.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/source/Tangle.h>

namespace
{
void ValidateResults(const vtkm::cont::PartitionedDataSet& truth,
                     const vtkm::cont::PartitionedDataSet& result,
                     const std::string& varName)
{
  VTKM_TEST_ASSERT(truth.GetNumberOfPartitions() == result.GetNumberOfPartitions());
  vtkm::Id numDS = truth.GetNumberOfPartitions();
  for (vtkm::Id i = 0; i < numDS; i++)
  {
    auto truthDS = truth.GetPartition(i);
    auto resultDS = result.GetPartition(i);
    VTKM_TEST_ASSERT(truthDS.GetNumberOfPoints() == resultDS.GetNumberOfPoints(),
                     "Wrong number of points");
    VTKM_TEST_ASSERT(truthDS.GetNumberOfCells() == resultDS.GetNumberOfCells(),
                     "Wrong number of cells");
    VTKM_TEST_ASSERT(resultDS.HasField(varName), "Missing field");

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> truthField, resultField;
    truthDS.GetField(varName).GetData().AsArrayHandle(truthField);
    resultDS.GetField(varName).GetData().AsArrayHandle(resultField);
    VTKM_TEST_ASSERT(truthField.GetNumberOfValues() == resultField.GetNumberOfValues(),
                     "Wrong number of field values");

    vtkm::Id numPts = truthField.GetNumberOfValues();
    const auto truthPortal = truthField.ReadPortal();
    const auto resultPortal = resultField.ReadPortal();
    for (vtkm::Id j = 0; j < numPts; j++)
      VTKM_TEST_ASSERT(truthPortal.Get(j) == resultPortal.Get(j), "Wrong value in field");
  }
}

void TestMultiBlockFilter()
{
  vtkm::cont::PartitionedDataSet pds;

  for (int i = 0; i < 50; i++)
  {
    vtkm::Id3 dims(9 + i, 9 + i, 9 + i);
    vtkm::source::Tangle tangle(dims);
    pds.AppendPartition(tangle.Execute());
  }

  std::cout << "ClipWithField" << std::endl;
  std::vector<vtkm::cont::PartitionedDataSet> results;
  std::vector<bool> flags = { false, true };
  for (const auto doThreading : flags)
  {
    vtkm::filter::ClipWithField clip;
    clip.SetRunMultiThreadedFilter(doThreading);
    clip.SetClipValue(0.0);
    clip.SetActiveField("nodevar");
    clip.SetFieldsToPass("nodevar", vtkm::cont::Field::Association::POINTS);
    auto result = clip.Execute(pds);
    VTKM_TEST_ASSERT(result.GetNumberOfPartitions() == pds.GetNumberOfPartitions());
    results.push_back(result);
  }
  ValidateResults(results[0], results[1], "nodevar");

  results.clear();

  /*

  std::cout << "Contour" << std::endl;
  vtkm::filter::Contour mc;
  mc.SetGenerateNormals(true);
  mc.SetIsoValue(0, 0.5);
  mc.SetActiveField("nodevar");
  mc.SetFieldsToPass(vtkm::filter::FieldSelection::MODE_NONE);
  result = mc.Execute(pds);
  VTKM_TEST_ASSERT(result.GetNumberOfPartitions() == pds.GetNumberOfPartitions());

  std::cout << "CleanGrid" << std::endl;
  vtkm::filter::CleanGrid clean;
  clean.SetCompactPointFields(true);
  clean.SetMergePoints(true);
  result = clean.Execute(pds);
  VTKM_TEST_ASSERT(result.GetNumberOfPartitions() == pds.GetNumberOfPartitions());


  std::cout << "Gradient" << std::endl;
  vtkm::filter::Gradient grad;
  grad.SetComputePointGradient(true);
  grad.SetActiveField("nodevar");
  grad.SetOutputFieldName("Gradient");
  result = grad.Execute(pds);
  VTKM_TEST_ASSERT(result.GetNumberOfPartitions() == pds.GetNumberOfPartitions());
  */
}

}; //namespace

int UnitTestMultiBlockFilter(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestMultiBlockFilter, argc, argv);
}
