//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_cont_testing_TestingCellLocatorXGCGrid_h
#define vtk_m_cont_testing_TestingCellLocatorXGCGrid_h

#include <random>
#include <string>

#include <vtkm/cont/ArrayHandleXGCCoordinates.h>
#include <vtkm/cont/CellLocatorXGCGrid.h>
#include <vtkm/cont/DataSet.h>

#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/cont/testing/Testing.h>

#include <vtkm/exec/CellLocatorXGCGrid.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/ScatterPermutation.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>


namespace
{
std::default_random_engine RandomGenerator;

class ParametricToWorldCoordinates : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn cellset,
                                FieldInPoint coords,
                                FieldInOutCell pcs,
                                FieldOutCell wcs);
  using ExecutionSignature = void(CellShape, _2, _3, _4);

  using ScatterType = vtkm::worklet::ScatterPermutation<>;

  VTKM_CONT
  static ScatterType MakeScatter(const vtkm::cont::ArrayHandle<vtkm::Id>& cellIds)
  {
    return ScatterType(cellIds);
  }

  template <typename CellShapeTagType, typename PointsVecType, typename PointType>
  VTKM_EXEC void operator()(CellShapeTagType cellShape,
                            PointsVecType points,
                            const PointType& pc,
                            PointType& wc) const
  {
    auto status = vtkm::exec::CellInterpolate(points, pc, cellShape, wc);
    std::cout << "Parametric to World: " << pc << " --> " << wc << std::endl;
    std::cout << " wedge points: " << std::endl;
    for (int i = 0; i < 6; i++)
      std::cout << i << ":  " << points[i] << std::endl;
    PointType wc2;
    vtkm::exec::ParametricCoordinatesToWorldCoordinates(points, pc, cellShape, wc2);
    std::cout << "  ---------------> convert back: " << wc2 << std::endl;

    if (status != vtkm::ErrorCode::Success)
      this->RaiseError(vtkm::ErrorString(status));
  }
};


void GenerateRandomInput(const vtkm::cont::DataSet& ds,
                         vtkm::Id count,
                         vtkm::cont::ArrayHandle<vtkm::Id>& cellIds,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& pcoords,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& wcoords)
{
  vtkm::Id numberOfCells = ds.GetNumberOfCells();

  std::uniform_int_distribution<vtkm::Id> cellIdGen(0, numberOfCells - 1);
  //std::uniform_real_distribution<vtkm::FloatDefault> pcoordGen(0.0f, 1.0f);
  std::uniform_real_distribution<vtkm::FloatDefault> pcoordGen(0.1f, 0.9f); //FIX ME

  cellIds.Allocate(count);
  pcoords.Allocate(count);
  wcoords.Allocate(count);

  auto cellIdPortal = cellIds.WritePortal();
  auto pcoordsPortal = pcoords.WritePortal();
  for (vtkm::Id i = 0; i < count; ++i)
  {
    vtkm::Id cid = cellIdGen(RandomGenerator);
    //cid = 0;
    cellIdPortal.Set(i, cid);
    auto p0 = pcoordGen(RandomGenerator);
    auto p1 = pcoordGen(RandomGenerator);
    auto p2 = pcoordGen(RandomGenerator);

    while (p0 + p1 > 1)
    {
      p0 = pcoordGen(RandomGenerator);
      p1 = pcoordGen(RandomGenerator);
    }

    vtkm::Vec3f pc{ p0, p1, p2 };
    std::cout << "cid= " << cellIdPortal.Get(i) << "  GEN RANDOMINPUT: " << pc << std::endl;
    pcoordsPortal.Set(i, pc);
  }

  vtkm::worklet::DispatcherMapTopology<ParametricToWorldCoordinates> dispatcher(
    ParametricToWorldCoordinates::MakeScatter(cellIds));
  dispatcher.Invoke(
    ds.GetCellSet(), ds.GetCoordinateSystem().GetDataAsMultiplexer(), pcoords, wcoords);

  std::cout << "****************************************" << std::endl;
  std::cout << "RandomInput" << std::endl << std::endl;
  std::ofstream ptF;
  ptF.open("points.txt");
  for (vtkm::Id i = 0; i < count; i++)
  {
    auto pt = wcoords.ReadPortal().Get(i);
    ptF << pt[0] << "," << pt[1] << "," << pt[2] << std::endl;
    std::cout << "   " << i << ": " << cellIdPortal.Get(i) << " " << pcoordsPortal.Get(i) << " --> "
              << wcoords.ReadPortal().Get(i) << std::endl;
  }
  std::cout << "****************************************" << std::endl;
  std::cout << "****************************************" << std::endl;
}


//-----------------------------------------------------------------------------
class FindCellWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn points,
                                ExecObject locator,
                                FieldOut cellIds,
                                FieldOut pcoords);
  using ExecutionSignature = void(_1, _2, _3, _4);

  template <typename LocatorType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const LocatorType& locator,
                            vtkm::Id& cellId,
                            vtkm::Vec3f& pcoords) const
  {
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      this->RaiseError(vtkm::ErrorString(status));
    }
  }
};


template <typename DeviceAdapter>
class TestingCellLocatorXGCGrid
{
public:
  using Algorithm = vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>;

  void TestLocator() const
  {
    vtkm::cont::DataSet dataset = vtkm::cont::testing::MakeTestDataSet().MakeXGCDataSet();
    vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();
    vtkm::cont::DynamicCellSet cellSet = dataset.GetCellSet();

    std::cout << "XGC Grid Test!!" << std::endl;
    vtkm::cont::CellLocatorXGCGrid locator;
    locator.SetCoordinates(coords);
    locator.SetCellSet(cellSet);

    locator.Update();

    //Generate some test points.
    const vtkm::Id numPoints = 15;
    vtkm::cont::ArrayHandle<vtkm::Id> expCellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> expPCoords;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> points;
    GenerateRandomInput(dataset, numPoints, expCellIds, expPCoords, points);

    vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
    vtkm::worklet::DispatcherMapField<FindCellWorklet> dispatcher;
    dispatcher.Invoke(points, locator, cellIds, pcoords);

    auto cellIdPortal = cellIds.ReadPortal();
    auto expCellIdsPortal = expCellIds.ReadPortal();
    auto pcoordsPortal = pcoords.ReadPortal();
    auto expPCoordsPortal = expPCoords.ReadPortal();
    for (vtkm::Id i = 0; i < numPoints; ++i)
    {
      std::cout << i << ":   TEST: " << pcoordsPortal.Get(i)
                << " should be: " << expPCoordsPortal.Get(i) << std::endl;
      std::cout << "    cellID: " << cellIdPortal.Get(i)
                << " should be: " << expCellIdsPortal.Get(i) << std::endl;
      /*
      VTKM_TEST_ASSERT(cellIdPortal.Get(i) == expCellIdsPortal.Get(i), "Incorrect cell ids");
      VTKM_TEST_ASSERT(test_equal(pcoordsPortal.Get(i), expPCoordsPortal.Get(i), 1e-3),
                       "Incorrect parameteric coordinates");
*/
    }


#if 0
    vtkm::Bounds bounds = coords.GetBounds();
    std::cout << "X bounds : " << bounds.X.Min << " to " << bounds.X.Max << std::endl;
    std::cout << "Y bounds : " << bounds.Y.Min << " to " << bounds.Y.Max << std::endl;
    std::cout << "Z bounds : " << bounds.Z.Min << " to " << bounds.Z.Max << std::endl;

    using StructuredType = vtkm::cont::CellSetStructured<3>;
    vtkm::Id3 cellDims =
      cellSet.Cast<StructuredType>().GetSchedulingRange(vtkm::TopologyElementTagCell());
    std::cout << "Dimensions of dataset : " << cellDims << std::endl;


    // Generate some sample points.
    using PointType = vtkm::Vec3f;
    std::vector<PointType> pointsVec;
    std::default_random_engine dre;
    std::uniform_real_distribution<vtkm::Float32> inBounds(0.0f, 4.0f);
    for (size_t i = 0; i < 10; i++)
    {
      PointType point = vtkm::make_Vec(inBounds(dre), inBounds(dre), inBounds(dre));
      pointsVec.push_back(point);
    }
    std::uniform_real_distribution<vtkm::Float32> outBounds(4.0f, 5.0f);
    for (size_t i = 0; i < 5; i++)
    {
      PointType point = vtkm::make_Vec(outBounds(dre), outBounds(dre), outBounds(dre));
      pointsVec.push_back(point);
    }
    std::uniform_real_distribution<vtkm::Float32> outBounds2(-1.0f, 0.0f);
    for (size_t i = 0; i < 5; i++)
    {
      PointType point = vtkm::make_Vec(outBounds2(dre), outBounds2(dre), outBounds2(dre));
      pointsVec.push_back(point);
    }

    // Add points right on the boundary.
    pointsVec.push_back(vtkm::make_Vec(0, 0, 0));
    pointsVec.push_back(vtkm::make_Vec(4, 4, 4));
    pointsVec.push_back(vtkm::make_Vec(4, 0, 0));
    pointsVec.push_back(vtkm::make_Vec(0, 4, 0));
    pointsVec.push_back(vtkm::make_Vec(0, 0, 4));
    pointsVec.push_back(vtkm::make_Vec(4, 4, 0));
    pointsVec.push_back(vtkm::make_Vec(0, 4, 4));
    pointsVec.push_back(vtkm::make_Vec(4, 0, 4));

    vtkm::cont::ArrayHandle<PointType> points =
      vtkm::cont::make_ArrayHandle(pointsVec, vtkm::CopyFlag::Off);
    // Query the points using the locators.
    vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
    vtkm::cont::ArrayHandle<PointType> parametric;
    vtkm::cont::ArrayHandle<bool> match;
    LocatorWorklet worklet(bounds, cellDims);
    vtkm::worklet::DispatcherMapField<LocatorWorklet> dispatcher(worklet);
    dispatcher.SetDevice(DeviceAdapter());
    dispatcher.Invoke(points, locator, cellIds, parametric, match);

    auto matchPortal = match.ReadPortal();
    for (vtkm::Id index = 0; index < match.GetNumberOfValues(); index++)
    {
      VTKM_TEST_ASSERT(matchPortal.Get(index), "Points do not match");
    }
    std::cout << "Test finished successfully." << std::endl;
#endif
  }

  void operator()() const
  {
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(DeviceAdapter());
    this->TestLocator();
  }
};

} // anonymous namespace

#endif
