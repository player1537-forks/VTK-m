#include <vtkm/Particle.h>
#include <vtkm/cont/ArrayCopy.h>

#include "XGCParameters.h"
#include "RunPoincare2.h"
#include "Poincare2.h"


void
RunPoincare2(const vtkm::cont::DataSet& ds,
             std::vector<vtkm::Particle>& seeds,
             XGCParameters& xgcParams,
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
    std::cout<<"Warning: Using linear B, forcing UseBOnly = true."<<std::endl;
  }
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

  auto cellSet = ds.GetCellSet().AsCellSet<vtkm::cont::CellSetSingleType<>>();

  vtkm::cont::CellLocatorTwoLevel locator2L;
  if (args.find("--LocatorDensity") != args.end())
  {
    auto d = args["--LocatorDensity"];
    vtkm::FloatDefault density1 = std::atof(d[0].c_str());
    vtkm::FloatDefault density2 = std::atof(d[1].c_str());
    locator2L.SetDensityL1(density1);
    locator2L.SetDensityL2(density2);
    std::cout<<"SetDensity: "<<density1<<" "<<density2<<std::endl;
  }

  locator2L.SetCellSet(cellSet);
  locator2L.SetCoordinates(ds.GetCoordinateSystem());
  auto startL = std::chrono::steady_clock::now();
  locator2L.Update();
  std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
  std::cout<<"2L build= "<<dt.count()<<std::endl;

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

  auto seedsArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

  if (useTraces)
  {
    std::vector<vtkm::Vec3f> t;
    t.resize(seedsArray.GetNumberOfValues()*worklet.MaxIter, {-100, -100, -100});
    std::cout<<"Allocate TRACES: "<<t.size()<<std::endl;
    tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
  }

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayCopy(ds.GetCoordinateSystem().GetData(), coords);

  invoker(worklet, seedsArray,
          locator2L,
          cellSet, coords,
          As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
          B_RZP, psi,
          tracesArr, outR, outZ, outTheta, outPsi, outID);

  outID.SyncControlArray();
}
