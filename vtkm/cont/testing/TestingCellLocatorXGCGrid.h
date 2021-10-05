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
  ParametricToWorldCoordinates(bool isCyl)
    : IsCylindrical(isCyl)
  {
  }

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
    //wrap around....
    auto tmpPts = points;
    if (this->IsCylindrical && tmpPts[0][1] > tmpPts[3][1])
    {
      tmpPts[3][1] = static_cast<vtkm::FloatDefault>(vtkm::TwoPi());
      tmpPts[4][1] = static_cast<vtkm::FloatDefault>(vtkm::TwoPi());
      tmpPts[5][1] = static_cast<vtkm::FloatDefault>(vtkm::TwoPi());
    }

    auto status = vtkm::exec::CellInterpolate(tmpPts, pc, cellShape, wc);
    /*
    std::cout << "Parametric to World: " << pc << " --> " << wc << std::endl;
    std::cout << " wedge points: " << std::endl;

    for (int i = 0; i < 6; i++)
      std::cout << i << ":  " << tmpPts[i] << std::endl;
    PointType wc2;
    vtkm::exec::ParametricCoordinatesToWorldCoordinates(tmpPts, pc, cellShape, wc2);
    std::cout << "  ---------------> convert back: " << wc2 << std::endl;
*/
    if (status != vtkm::ErrorCode::Success)
      this->RaiseError(vtkm::ErrorString(status));
  }

private:
  bool IsCylindrical;
};

void GenerateRandomInput(const vtkm::cont::DataSet& ds,
                         vtkm::cont::ArrayHandle<vtkm::Id>& cellIds,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& pcoords,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& wcoords,
                         bool isCyl)
{
  vtkm::cont::DynamicCellSet cellSet = ds.GetCellSet();
  auto xgcCellSet = cellSet.Cast<vtkm::cont::CellSetExtrude>();

  vtkm::Id numCells = ds.GetNumberOfCells();
  if (!xgcCellSet.GetIsPeriodic())
    numCells -= xgcCellSet.GetNumberOfCellsPerPlane();

  auto cellIdGen = std::uniform_int_distribution<vtkm::Id>(0, numCells - 1);
  std::uniform_real_distribution<vtkm::FloatDefault> pcoordGen;

  //Cartesian is not as accurate. Points too close to a boundary will be outside.
  if (isCyl)
    pcoordGen = std::uniform_real_distribution<vtkm::FloatDefault>(0.0f, 1.0f);
  else
    pcoordGen = std::uniform_real_distribution<vtkm::FloatDefault>(0.1f, 0.9f);

  std::vector<vtkm::Id> cids;
  std::vector<vtkm::Vec3f> pc, wc;

  const vtkm::Id count = 20;
  for (vtkm::Id i = 0; i < count; i++)
  {
    vtkm::Id cid = cellIdGen(RandomGenerator);
    vtkm::FloatDefault p0 = pcoordGen(RandomGenerator);
    vtkm::FloatDefault p1 = pcoordGen(RandomGenerator);
    vtkm::FloatDefault p2 = pcoordGen(RandomGenerator);
    while (p0 + p1 > 1)
    {
      p0 = pcoordGen(RandomGenerator);
      p1 = pcoordGen(RandomGenerator);
    }

    vtkm::Vec3f p{ p0, p1, p2 };
    pc.push_back(p);
    cids.push_back(cid);
  }

  cellIds = vtkm::cont::make_ArrayHandle(cids, vtkm::CopyFlag::On);
  pcoords = vtkm::cont::make_ArrayHandle(pc, vtkm::CopyFlag::On);

  vtkm::cont::Invoker invoker;
  invoker(ParametricToWorldCoordinates{ isCyl },
          ParametricToWorldCoordinates::MakeScatter(cellIds),
          ds.GetCellSet(),
          ds.GetCoordinateSystem().GetDataAsMultiplexer(),
          pcoords,
          wcoords);
  /*
  std::cout << "****************************************" << std::endl;
  std::cout << "RandomInput" << std::endl << std::endl;
  std::ofstream ptF;
  ptF.open("points.txt");
  auto cidPortal = cellIds.ReadPortal();
  auto pcPortal = pcoords.ReadPortal();
  auto wcPortal = wcoords.ReadPortal();
  for (vtkm::Id i = 0; i < cidPortal.GetNumberOfValues(); i++)
  {
    auto pt = wcPortal.Get(i);
    ptF << pt[0] << "," << pt[1] << "," << pt[2] << std::endl;
    std::cout << "   " << i << ": cell= " << cidPortal.Get(i) << " param= " << pcPortal.Get(i) << " --> wc= "
              << wcPortal.Get(i) << std::endl;
  }
  std::cout << "****************************************" << std::endl;
  std::cout << "****************************************" << std::endl;
*/
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
      this->RaiseError(vtkm::ErrorString(status));
  }
};


template <typename DeviceAdapter>
class TestingCellLocatorXGCGrid
{
public:
  using Algorithm = vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>;

  void TestLocator(bool isCyl, vtkm::Id numPlanes) const
  {
    vtkm::cont::DataSet dataset =
      vtkm::cont::testing::MakeTestDataSet().MakeXGCDataSet(isCyl, numPlanes, true);
    vtkm::cont::CoordinateSystem coords = dataset.GetCoordinateSystem();
    vtkm::cont::DynamicCellSet cellSet = dataset.GetCellSet();

    std::cout << "XGC Grid Test. Cylindrical= " << isCyl << " numPlanes= " << numPlanes
              << std::endl;
    vtkm::cont::CellLocatorXGCGrid locator;
    locator.SetCoordinates(coords);
    locator.SetCellSet(cellSet);
    locator.Update();

    //Generate some test points.
    vtkm::cont::ArrayHandle<vtkm::Id> expCellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> expPCoords;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> points;
    GenerateRandomInput(dataset, expCellIds, expPCoords, points, isCyl);
    vtkm::Id numPoints = expCellIds.GetNumberOfValues();

    vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;

    vtkm::cont::Invoker invoker;
    invoker(FindCellWorklet{}, points, locator, cellIds, pcoords);

    auto cellIdPortal = cellIds.ReadPortal();
    auto expCellIdsPortal = expCellIds.ReadPortal();
    auto pcoordsPortal = pcoords.ReadPortal();
    auto expPCoordsPortal = expPCoords.ReadPortal();
    for (vtkm::Id i = 0; i < numPoints; ++i)
    {
      /*
      std::cout << i << ":   TEST: " << pcoordsPortal.Get(i)
                << " should be: " << expPCoordsPortal.Get(i) << std::endl;
      std::cout << ":       wcoords= : " << points.ReadPortal().Get(i) << std::endl;
      std::cout << "        cellID: " << cellIdPortal.Get(i)
                << " should be: " << expCellIdsPortal.Get(i) << std::endl;
      */
      VTKM_TEST_ASSERT(cellIdPortal.Get(i) == expCellIdsPortal.Get(i), "Incorrect cell ids");
      VTKM_TEST_ASSERT(test_equal(pcoordsPortal.Get(i), expPCoordsPortal.Get(i), 1e-3),
                       "Incorrect parameteric coordinates");
    }
  }

  void operator()() const
  {
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(DeviceAdapter());
    std::vector<bool> isCyl = { true, false };
    std::vector<vtkm::Id> numPlanes = { 8, 16, 32 };

    for (const auto& c : isCyl)
      for (const auto& np : numPlanes)
        this->TestLocator(c, np);
  }
};

} // anonymous namespace

#endif
