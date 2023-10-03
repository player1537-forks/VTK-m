#include <vtkm/Particle.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/CellLocatorUniformGrid.h>
#include <vtkm/cont/CellLocatorUniformBins.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ParticleArrayCopy.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/Timer.h>

#include "XGCParameters.h"
#include "RunPoincare2.h"
#include "Poincare2.h"
#include "Poincare2_multi_particle.h"

#include "ComputeAs.h"
//#include "ComputeAsCell.h"

#include <vtkm/exec/LocatorTimer.h>

vtkm::exec::LocatorTimer locTimer(30);

/*
extern void vtkm::exec::ResetTimers(int nThreads);
extern void UpdateTimer(std::size_t idx, const std::chrono::time_point<std::chrono::steady_clock>& start);
extern std::vector<std::vector<double>> vtkm::exec::timers;
extern std::vector<std::vector<vtkm::Id>> vtkm::exec::counters;
*/

void
RunSamplingTest(const vtkm::cont::DataSet& ds,
                const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& coeff_2D,
                const XGCParameters& xgcParams,
                std::map<std::string, std::vector<std::string>>& args)
{
  auto a = args["--AsGridSize"];
  vtkm::Id AsGridSizeX = 128, AsGridSizeY = 128;
  if (a.size() == 2 )
  {
    AsGridSizeX = std::atoi(a[0].c_str());
    AsGridSizeY = std::atoi(a[1].c_str());
  }
  else if (a.size() == 1)
  {
    AsGridSizeX = std::atoi(a[0].c_str());
    AsGridSizeY = AsGridSizeX;
  }

  //Compute the As values on a uniform grid.
  vtkm::Id3 dimsRZ(AsGridSizeX, AsGridSizeY, 1);

  //create a uniform RZ dataset.
  vtkm::FloatDefault dR = (xgcParams.eq_max_r-xgcParams.eq_min_r) / static_cast<vtkm::FloatDefault>(AsGridSizeX-1);
  vtkm::FloatDefault dZ = (xgcParams.eq_max_z-xgcParams.eq_min_z) / static_cast<vtkm::FloatDefault>(AsGridSizeY-1);
  vtkm::Vec3f originRZ(xgcParams.eq_min_r, xgcParams.eq_min_z, 0);
  vtkm::Vec3f spacingRZ(dR, dZ, 1);
  auto gridRZ = vtkm::cont::DataSetBuilderUniform::Create(dimsRZ, originRZ, spacingRZ);
  // std::cout<<"GridRZ**************************"<<std::endl;
  // gridRZ.PrintSummary(std::cout);

  //create a uniform psi-theta dataset.
  vtkm::FloatDefault psiNorm0 = xgcParams.psi_min / xgcParams.eq_x_psi;
  vtkm::FloatDefault psiNorm1 = xgcParams.psi_max / xgcParams.eq_x_psi;
  vtkm::FloatDefault dPsi = (psiNorm1 - psiNorm0) / static_cast<vtkm::FloatDefault>(AsGridSizeX-1);
  vtkm::FloatDefault dTheta = vtkm::TwoPi() / static_cast<vtkm::FloatDefault>(AsGridSizeY-1);
  vtkm::Vec3f originTP(0, psiNorm0, 0);
  vtkm::Vec3f spacingTP(dTheta, dPsi, 1);
  vtkm::Id3 dimsTP(AsGridSizeX, AsGridSizeY, 1); //xgcParams.numPlanes*2);
  auto gridTP = vtkm::cont::DataSetBuilderUniform::Create(dimsTP, originTP, spacingTP);
  // std::cout<<std::endl<<std::endl<<"GridTP**************************"<<std::endl;
  // gridTP.PrintSummary(std::cout);

  vtkm::cont::Invoker invoker;
  //vtkm::cont::CellLocatorTwoLevel locator2L;
  vtkm::cont::CellLocatorUniformBins locator2L;
  locator2L.SetDimensions({512, 512, 1});
  locator2L.SetCellSet(ds.GetCellSet());
  locator2L.SetCoordinates(ds.GetCoordinateSystem());
  locator2L.Update();
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> dAs_ff_rzp;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> dAsUniform;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> AsUniform;

  ComputeAsWorklet computeAs(xgcParams, AsGridSizeX, AsGridSizeY);

  AsUniform.AllocateAndFill(gridTP.GetNumberOfPoints(), 0.0);
  dAsUniform.AllocateAndFill(gridTP.GetNumberOfPoints(), 0.0);

  invoker(computeAs,
          gridRZ.GetCoordinateSystem(),
          locator2L,
          ds.GetCellSet(),
          coeff_2D,
          As_ff, dAs_ff_rzp,
          AsUniform, dAsUniform);

  std::string fname = "AsGrid.vtk";
  gridTP.AddField(vtkm::cont::make_FieldPoint("As_ff", AsUniform));
  gridTP.AddField(vtkm::cont::make_FieldPoint("dAs_ff", dAsUniform));
  gridRZ.AddField(vtkm::cont::make_FieldPoint("As_ff", AsUniform));

  // std::cout<<"Saving As file: "<<fname<<std::endl;
  vtkm::io::VTKDataSetWriter writer(fname.c_str());
  writer.WriteDataSet(gridTP);

/*
  {
    vtkm::io::VTKDataSetWriter writer("ds.vtk");
    writer.WriteDataSet(ds);
  }
*/


  return;


#if 0
  auto grid = vtkm::cont::DataSetBuilderUniform::Create(dims, origin, spacing);
  vtkm::cont::ArrayHandleUniformPointCoordinates uniform2DCoords({AsGridSizeX, AsGridSizeY, 1}, origin, spacing);
  vtkm::cont::CellSetStructured<2> uniform2DCells;
  uniform2DCells.SetPointDimensions(vtkm::Id2(AsGridSizeX, AsGridSizeY));

  // std::cout<<"********** NumCells= "<<grid.GetCellSet().GetNumberOfCells()<<std::endl;
  // grid.PrintSummary(std::cout);

  vtkm::cont::Invoker invoker;
  vtkm::cont::CellLocatorTwoLevel locator2L;
  locator2L.SetCellSet(ds.GetCellSet());
  locator2L.SetCoordinates(ds.GetCoordinateSystem());
  locator2L.Update();

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> dAs_ff_rzp;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> dAsUniform;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> AsUniform;

  // std::cout<<"********************* AS Point"<<std::endl;
  ComputeAsWorklet computeAs(xgcParams);
  vtkm::Id nPts = uniform2DCoords.GetNumberOfValues() * xgcParams.numPlanes * 2;
  AsUniform.Allocate(nPts);
  dAsUniform.Allocate(nPts);
  computeAs.NumSampsX = AsGridSizeX;
  computeAs.NumSampsY = AsGridSizeY;

  invoker(computeAs,
          uniform2DCoords,
          locator2L,
          ds.GetCellSet(),
          coeff_2D,
          As_ff, dAs_ff_rzp,
          AsUniform, dAsUniform);

  std::string fname = "AsGrid.vtk";
  grid.AddField(vtkm::cont::make_FieldPoint("As_ff", AsUniform));
  grid.AddField(vtkm::cont::make_FieldPoint("dAs_ff", dAsUniform));
  // grid.PrintSummary(std::cout);

  // std::cout<<"Saving As file: "<<fname<<std::endl;
  vtkm::io::VTKDataSetWriter writer(fname.c_str());
  writer.WriteDataSet(grid);



  //Check the errors.
  vtkm::cont::CellLocatorUniformGrid locatorU;
  locatorU.SetCoordinates(vtkm::cont::CoordinateSystem("coords", uniform2DCoords));
  locatorU.SetCellSet(uniform2DCells);
  locatorU.Update();

  AsErrorWorklet asErrorWorklet(xgcParams, uniform2DCoords.GetNumberOfValues());

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> AsError, dAsError, AsErrorPlane, dAsErrorMagPlane, dAsErrorRPlane, dAsErrorZPlane, dAsErrorPPlane;
  AsError.Allocate(ds.GetNumberOfPoints() * xgcParams.numPlanes * 2);
  dAsError.Allocate(ds.GetNumberOfPoints() * xgcParams.numPlanes * 2);
  invoker(asErrorWorklet,
          ds.GetCoordinateSystem(),
          As_ff,
          dAs_ff_rzp,

          uniform2DCells,
          locatorU,
          AsUniform,
          dAsUniform,

          AsError,
          dAsError,
          AsErrorPlane,
          dAsErrorMagPlane,
          dAsErrorRPlane,
          dAsErrorZPlane,
          dAsErrorPPlane);

  auto AsErrorTot = vtkm::cont::Algorithm::Reduce(AsError, 0.0f);
  auto AsErrorMin = vtkm::cont::Algorithm::Reduce(AsError, 1e10f, vtkm::Minimum());
  auto AsErrorMax = vtkm::cont::Algorithm::Reduce(AsError, -1e10f, vtkm::Maximum());
  // std::cout<<"AsError: m/M "<<AsErrorMin<<" "<<AsErrorMax<<" Avg= "<<AsErrorTot/(double)AsError.GetNumberOfValues()<<std::endl;

  auto dAsErrorTot = vtkm::cont::Algorithm::Reduce(dAsError, 0.0f);
  auto dAsErrorMin = vtkm::cont::Algorithm::Reduce(dAsError, 1e10f, vtkm::Minimum());
  auto dAsErrorMax = vtkm::cont::Algorithm::Reduce(dAsError, -1e10f, vtkm::Maximum());
  // std::cout<<"dAsError: m/M "<<dAsErrorMin<<" "<<dAsErrorMax<<" Avg= "<<dAsErrorTot/(double)dAsError.GetNumberOfValues()<<std::endl;

  auto dumpDS = ds;
  dumpDS.AddField(vtkm::cont::make_FieldPoint("AsError", AsErrorPlane));
  dumpDS.AddField(vtkm::cont::make_FieldPoint("dAsError", dAsErrorMagPlane));
  dumpDS.AddField(vtkm::cont::make_FieldPoint("dAsErrorR", dAsErrorRPlane));
  dumpDS.AddField(vtkm::cont::make_FieldPoint("dAsErrorZ", dAsErrorZPlane));
  dumpDS.AddField(vtkm::cont::make_FieldPoint("dAsErrorP", dAsErrorPPlane));

  /*
  // std::cout<<"Write file......"<<std::endl;
  vtkm::io::VTKDataSetWriter writer("AsError.vtk");
  writer.WriteDataSet(dumpDS);
  */

#endif
}


void
RunPoincare2(const vtkm::cont::DataSet& ds,
             std::vector<vtkm::Particle>& seeds,
             const XGCParameters& xgcParams,
             std::map<std::string, std::vector<std::string>>& args,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As_ff,
             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& dAs_ff_rzp,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& coeff_1D,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& coeff_2D,
             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& B_RZP,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& psi,
             vtkm::cont::ArrayHandle<vtkm::Vec3f>& tracesArr,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outR,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outZ,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outTheta,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outPsi,
             vtkm::cont::ArrayHandle<vtkm::Id>& outID)
{

  if (args.find("--AsGridSize") != args.end())
  {
    RunSamplingTest(ds, coeff_2D, xgcParams, args);
    return;
  }


  //Get all the arguments...
  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  if (args.find("--revDir") != args.end())
    stepSize = -stepSize;

  vtkm::Id numPunc = std::atoi(args["--numPunc"][0].c_str());

  bool useTraces = false;
  if (args.find("--traces") != args.end()) useTraces = std::atoi(args["--traces"][0].c_str());
  bool useBOnly = false;
  if (args.find("--useBOnly") != args.end()) useBOnly = true;
  bool useLinearB = false;
  if (args.find("--useLinearB") != args.end()) useLinearB = true;
  if (useLinearB)
  {
    useBOnly = true;
    // std::cout<<"Warning: Using linear B, forcing UseBOnly = true."<<std::endl;
  }
  bool usePrevCell = true;
  if (args.find("--UsePrevCell") != args.end()) usePrevCell = std::atoi(args["--UsePrevCell"][0].c_str());

  bool validateInterp = false;
  vtkm::Id validateInterpSkip = 1;
  if (args.find("--validateInterpolation") != args.end())
  {
    validateInterp = true;
    validateInterpSkip = static_cast<vtkm::Id>(std::stoi(args["--validateInterpolation"][0].c_str()));
  }

  bool useDeltaBScale = false, useBScale = false;
  vtkm::FloatDefault deltaBScale = 1.0, bScale = 1.0;
  if (args.find("--deltaBScale") != args.end())
  {
    useDeltaBScale = true;
    deltaBScale = std::atof(args["--deltaBScale"][0].c_str());
  }
  if (args.find("--BScale") != args.end())
  {
    useBScale = true;
    bScale = std::atof(args["--BScale"][0].c_str());
  }
  bool flipBPhi = false;
  if (args.find("--flipBPhi") != args.end())
    flipBPhi = true;

  auto cellSet = ds.GetCellSet().AsCellSet<vtkm::cont::CellSetSingleType<>>();

  vtkm::cont::CellLocatorTwoLevel locator2L;
  vtkm::cont::CellLocatorUniformBins locatorUB;

  bool useUB = false;
  if (args.find("--UniformBins") != args.end())
  {
    auto a = args["--UniformBins"];
    vtkm::Id nx = std::atoi(a[0].c_str());
    vtkm::Id ny = std::atoi(a[1].c_str());
    locatorUB.SetDimensions({nx, ny, 1});
    useUB = true;
  } else {
    locatorUB.SetDimensions({8000, 8000, 1});
    useUB = true;
  }

  bool useMultiPointWorklet = false;
  if (args.find("--MultiPointWorklet") != args.end())
    useMultiPointWorklet = true;

  if (args.find("--LocatorDensity") != args.end())
  {
    auto d = args["--LocatorDensity"];
    vtkm::FloatDefault density1 = std::atof(d[0].c_str());
    vtkm::FloatDefault density2 = std::atof(d[1].c_str());
    locator2L.SetDensityL1(density1);
    locator2L.SetDensityL2(density2);
    // std::cout<<"SetDensity: "<<density1<<" "<<density2<<std::endl;
  }

  auto startL = std::chrono::steady_clock::now();
  if (useUB)
  {
    locatorUB.SetCellSet(cellSet);
    locatorUB.SetCoordinates(ds.GetCoordinateSystem());
    locatorUB.Update();
  }
  else
  {
    locator2L.SetCellSet(cellSet);
    locator2L.SetCoordinates(ds.GetCoordinateSystem());
    locator2L.Update();
  }
  std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
  // std::cout<<"Locator build= "<<dt.count()<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
  auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);

  vtkm::cont::Invoker invoker;

  vtkm::Id maxItersPerPunc = 2500;
  if (args.find("--MaxItersPerPunc") != args.end())
    maxItersPerPunc = std::stoi(args["--MaxItersPerPunc"][0].c_str());

  bool bothDir = false;
  vtkm::Id2 fwdIdxRange(-1,-1);
  if (args.find("--bothDir") != args.end())
  {
    bothDir = true;
    fwdIdxRange = vtkm::Id2(0, static_cast<vtkm::Id>(seeds.size()/2));
  }

  PoincareWorklet2 worklet(numPunc, 0.0f, stepSize, maxItersPerPunc, useTraces, xgcParams, bothDir, fwdIdxRange);
  worklet.one_d_cub_dpsi_inv = 1.0/dPsi;
  worklet.UseLinearB = useLinearB;
  worklet.ValidateInterpolation = validateInterp;
  worklet.ValidateInterpolationSkip = validateInterpSkip;
  worklet.UseDeltaBScale = useDeltaBScale;
  worklet.DeltaBScale = deltaBScale;
  worklet.UseBScale = useBScale;
  worklet.BScale = bScale;
  worklet.UseBOnly = useBOnly;
  worklet.UsePrevCell = usePrevCell;
  worklet.FlipBPhi = flipBPhi;

  PoincareWorklet2_multi_particle workletMP(numPunc, 0.0f, stepSize, maxItersPerPunc, useTraces, xgcParams, bothDir, fwdIdxRange);
  workletMP.one_d_cub_dpsi_inv = 1.0/dPsi;
  workletMP.UseLinearB = useLinearB;
  workletMP.ValidateInterpolation = validateInterp;
  workletMP.ValidateInterpolationSkip = validateInterpSkip;
  workletMP.UseDeltaBScale = useDeltaBScale;
  workletMP.DeltaBScale = deltaBScale;
  workletMP.UseBScale = useBScale;
  workletMP.BScale = bScale;
  workletMP.UseBOnly = useBOnly;
  workletMP.UsePrevCell = usePrevCell;
  workletMP.FlipBPhi = flipBPhi;

  auto seedsArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

  if (useTraces)
  {
    std::vector<vtkm::Vec3f> t;
    t.resize(seedsArray.GetNumberOfValues()*worklet.MaxIter, {-100, -100, -100});
    // std::cout<<"Allocate TRACES: "<<t.size()<<std::endl;
    tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
  }

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayCopy(ds.GetCoordinateSystem().GetData(), coords);

  locTimer.Reset();
  if (useMultiPointWorklet)
  {
    //doesn't work....
    //auto groupArr = vtkm::cont::make_ArrayHandleGroupVec<10>(seeds);

    //using SeedArrType = vtkm::cont::ArrayHandle<vtkm::Particle>;
    vtkm::Id start = 0, end = seedsArray.GetNumberOfValues(), step=1;
    if (end % 10 != 0)
      throw vtkm::cont::ErrorBadType("Num seeds not divisible by 10");

    auto counter = vtkm::cont::make_ArrayHandleCounting(start, step, end);
    auto indexGroupArr = vtkm::cont::make_ArrayHandleGroupVec<10>(counter);
    // std::cout<<"PT0: "<<seedsArray.ReadPortal().Get(0)<<std::endl;

    if (useUB)
    {
      invoker(workletMP, seedsArray, indexGroupArr,
              locatorUB,
              cellSet, coords,
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              tracesArr, outR, outZ, outTheta, outPsi, outID);
    }
    else
    {
      invoker(workletMP, seedsArray, indexGroupArr,
              locator2L,
              cellSet, coords,
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              tracesArr, outR, outZ, outTheta, outPsi, outID);
    }
  }
  else
  {
    if (useUB)
    {
      invoker(worklet, seedsArray,
              locatorUB,
              cellSet, coords,
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              B_RZP, psi,
              tracesArr, outR, outZ, outTheta, outPsi, outID);
    }
    else
    {
      invoker(worklet, seedsArray,
              locator2L,
              cellSet, coords,
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              B_RZP, psi,
              tracesArr, outR, outZ, outTheta, outPsi, outID);
    }
  }

#if 0
  else if (false)
  {
    invoker(worklet, seedsArray,
            locatorUB, //locator2L,
            cellSet, coords,
            As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
            B_RZP, psi,
            tracesArr, outR, outZ, outTheta, outPsi, outID);
  }
  else
  {
  }
#endif

  outID.SyncControlArray();

  //Get cell timers...
  // std::cout<<std::endl<<"*********************************************************"<<std::endl;
  // std::cout<<"  Timers"<<std::endl;
  // std::cout<<std::setprecision(14);
  auto T = locTimer.GetTime();
  auto C = locTimer.GetCounters();
  std::vector<std::string> N = {"FindCell0     ",
                                "FindCell1     ",
                                "FindCell2     ",
                                "FindCellImpl0 ",
                                "FindCellImpl1 ",
                                "FindCellImpl2 ",
                                "PointInCell   "};
  for (std::size_t i = 0; i < vtkm::exec::LocatorTimer::NUM_IDX; i++)
  {
    // std::cout<<locTimer.Names[i]<<": "<<T[i]<<"  "<<C[i]<<" avg: "<<T[i]/C[i]<<std::endl;
  }
  // std::cout<<"*********************************************************"<<std::endl<<std::endl;


  //Look for min/max and avg num steps.
  vtkm::cont::ArrayHandle<vtkm::Vec3f> outP;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> outTime;
  vtkm::cont::ArrayHandle<vtkm::Id> _outID, outSteps;
  vtkm::cont::ArrayHandle<vtkm::ParticleStatus> outStatus;
  vtkm::cont::ParticleArrayCopy(seedsArray, outP, _outID, outSteps, outStatus, outTime);
  vtkm::Id minVal = vtkm::cont::Algorithm::Reduce(outSteps, 0, vtkm::Minimum());
  vtkm::Id maxVal = vtkm::cont::Algorithm::Reduce(outSteps, 0, vtkm::Maximum());
  vtkm::Id sum = vtkm::cont::Algorithm::Reduce(outSteps, 0);

  vtkm::FloatDefault avg = vtkm::FloatDefault(sum) / vtkm::FloatDefault(seedsArray.GetNumberOfValues());

  // std::cout<<"NumSteps: mM: "<<minVal<<" "<<maxVal<<" AVG= "<<avg<<std::endl;


}
