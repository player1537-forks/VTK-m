#include <vtkm/Particle.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/CellSet.h>
#include <vtkm/cont/CellSetSingleType.h>
#include <vtkm/cont/CellLocatorUniformGrid.h>
#include <vtkm/cont/CellLocatorUniformBins.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ParticleArrayCopy.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/Timer.h>

#include "XGCParameters.h"
#include "SavePoincare.h"
#include "EvaluateAtNodes.h"

#include "ComputeAs.h"
//#include "ComputeAsCell.h"

#include <adios2.h>

static
void Save2Adios(const vtkm::cont::DataSet& ds,
                const vtkm::cont::CellSetSingleType<>& cells2D,
                vtkm::Id numNodes, vtkm::Id numPlanes,
                const std::string& fname)
{
  adios2::ADIOS adios;
  auto io = adios2::IO(adios.DeclareIO("output"));
  auto engine = io.Open(fname, adios2::Mode::Write);

  std::vector<std::size_t> shape2D, offset, size2D, shape3D, size3D;
  shape2D.push_back(numNodes);
  shape3D.push_back(numNodes*numPlanes);
  offset.push_back(0);
  size2D = shape2D;
  size3D = shape3D;

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  auto varR = io.DefineVariable<vtkm::FloatDefault>("R", shape2D, offset, size2D);
  auto varZ = io.DefineVariable<vtkm::FloatDefault>("Z", shape2D, offset, size2D);
  auto varB0_r = io.DefineVariable<vtkm::FloatDefault>("B0_r", shape2D, offset, size2D);
  auto varB0_z = io.DefineVariable<vtkm::FloatDefault>("B0_z", shape2D, offset, size2D);
  auto varB0_p = io.DefineVariable<vtkm::FloatDefault>("B0_p", shape2D, offset, size2D);
  /*
  auto var_dB_r = io.DefineVariable<vtkm::FloatDefault>("dB_r", shape2D, offset, size2D);
  auto var_dB_z = io.DefineVariable<vtkm::FloatDefault>("dB_z", shape2D, offset, size2D);
  auto var_dB_p = io.DefineVariable<vtkm::FloatDefault>("dB_p", shape2D, offset, size2D);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  auto varBr = io.DefineVariable<vtkm::FloatDefault>("Br", shape2D, offset, size2D);
  auto var_dBr_dr = io.DefineVariable<vtkm::FloatDefault>("dBr_dr", shape2D, offset, size2D);
  auto var_dBz_dz = io.DefineVariable<vtkm::FloatDefault>("dBz_dz", shape2D, offset, size2D);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  */

  auto varVec_r = io.DefineVariable<vtkm::FloatDefault>("vec_r", shape3D, offset, size3D);
  auto varVec_z = io.DefineVariable<vtkm::FloatDefault>("vec_z", shape3D, offset, size3D);
  auto varVec_p = io.DefineVariable<vtkm::FloatDefault>("vec_p", shape3D, offset, size3D);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0_rzp, vec_rzp, Br_rzp;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> Br, dBr_dr, dBz_dz;
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetPointField("B_rzp").GetData().AsArrayHandle(B0_rzp);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetPointField("B_dB_rzp").GetData().AsArrayHandle(vec_rzp);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.GetPointField("Br").GetData().AsArrayHandle(Br);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetPointField("dBr_dr").GetData().AsArrayHandle(dBr_dr);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetPointField("dBz_dz").GetData().AsArrayHandle(dBz_dz);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  //auto cellSet = ds.GetCellSet().AsCellSet<vtkm::cont::CellSetSingleType<>>();

  std::vector<vtkm::FloatDefault> R, Z, B0_r, B0_z, B0_p, dB_r, dB_z, dB_p, _dBr_dr, _dBz_dz, _Br;
  for (vtkm::Id i = 0; i < numNodes; i++)
  {
    R.push_back(coords.ReadPortal().Get(i)[0]);
    Z.push_back(coords.ReadPortal().Get(i)[1]);
    B0_r.push_back(B0_rzp.ReadPortal().Get(i)[0]);
    B0_z.push_back(B0_rzp.ReadPortal().Get(i)[1]);
    B0_p.push_back(B0_rzp.ReadPortal().Get(i)[2]);

    /*
    dB_r.push_back(Br_rzp.ReadPortal().Get(i)[0]);
    dB_z.push_back(Br_rzp.ReadPortal().Get(i)[1]);
    dB_p.push_back(Br_rzp.ReadPortal().Get(i)[2]);

    _dBr_dr.push_back(dBr_dr.ReadPortal().Get(i));
    _dBz_dz.push_back(dBz_dz.ReadPortal().Get(i));

    _Br.push_back(Br.ReadPortal().Get(i));
    */
  }
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  std::vector<vtkm::FloatDefault> vec_r, vec_z, vec_p;

  vtkm::Id idx=0;
  for (vtkm::Id p = 0; p < numPlanes; p++)
  {
    for (vtkm::Id i = 0; i < numNodes; i++)
    {
      vec_r.push_back(vec_rzp.ReadPortal().Get(idx)[0]);
      vec_z.push_back(vec_rzp.ReadPortal().Get(idx)[1]);
      vec_p.push_back(vec_rzp.ReadPortal().Get(idx)[2]);
      idx++;
    }
  }
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  std::vector<vtkm::Id> Conn;
  for (vtkm::Id i = 0; i < cells2D.GetNumberOfCells(); i++)
  {
    vtkm::Vec<vtkm::Id, 3> cellIds;
    cells2D.GetIndices(i, cellIds);
    for (int j = 0; j < 3; j++)
      Conn.push_back(cellIds[j]);
  }
  auto varConn = io.DefineVariable<vtkm::Id>("Conn", {Conn.size()}, {0}, {Conn.size()});


  engine.BeginStep();
  engine.Put<vtkm::FloatDefault>(varR, R.data());
  engine.Put<vtkm::FloatDefault>(varZ, Z.data());
  engine.Put<vtkm::Id>(varConn, Conn.data());

  engine.Put<vtkm::FloatDefault>(varB0_r, B0_r.data());
  engine.Put<vtkm::FloatDefault>(varB0_z, B0_z.data());
  engine.Put<vtkm::FloatDefault>(varB0_p, B0_p.data());
  engine.Put<vtkm::FloatDefault>(varVec_r, vec_r.data());
  engine.Put<vtkm::FloatDefault>(varVec_z, vec_z.data());
  engine.Put<vtkm::FloatDefault>(varVec_p, vec_p.data());

  /*
  engine.Put<vtkm::FloatDefault>(var_dB_r, dB_r.data());
  engine.Put<vtkm::FloatDefault>(var_dB_z, dB_z.data());
  engine.Put<vtkm::FloatDefault>(var_dB_p, dB_p.data());

  engine.Put<vtkm::FloatDefault>(var_dBr_dr, _dBr_dr.data());
  engine.Put<vtkm::FloatDefault>(var_dBz_dz, _dBz_dz.data());

  engine.Put<vtkm::FloatDefault>(varBr, _Br.data());
  */
  engine.EndStep();

  engine.Close();
}

void
SavePoincare(const vtkm::cont::DataSet& ds,
             const XGCParameters& xgcParams,
             std::map<std::string, std::vector<std::string>>& args)
{

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff, coeff_1D, coeff_2D, psi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_RZP, B_Norm_rzp, dAs_ff_rzp;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
  ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff_1D);
  ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);
  ds.GetField("psi2D").GetData().AsArrayHandle(psi);
  ds.GetField("B_RZP").GetData().AsArrayHandle(B_RZP);

  /*
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
  */

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
    std::cout<<"SetUniformBins: "<<nx<<" "<<ny<<std::endl;
  }

  if (args.find("--LocatorDensity") != args.end())
  {
    auto d = args["--LocatorDensity"];
    vtkm::FloatDefault density1 = std::atof(d[0].c_str());
    vtkm::FloatDefault density2 = std::atof(d[1].c_str());
    locator2L.SetDensityL1(density1);
    locator2L.SetDensityL2(density2);
    std::cout<<"SetDensity: "<<density1<<" "<<density2<<std::endl;
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
  std::cout<<"Locator build= "<<dt.count()<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
  auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);

  vtkm::cont::Invoker invoker;

  vtkm::Id maxItersPerPunc = 2500;
  if (args.find("--MaxItersPerPunc") != args.end())
    maxItersPerPunc = std::stoi(args["--MaxItersPerPunc"][0].c_str());

  //PoincareWorklet2 worklet(numPunc, 0.0f, stepSize, maxItersPerPunc, useTraces, xgcParams, bothDir, fwdIdxRange);
  //PoincareWorklet2_multi_particle worklet(numPunc, 0.0f, stepSize, maxItersPerPunc, useTraces, xgcParams, bothDir, fwdIdxRange);
  EvaluateAtNodes worklet(xgcParams);

  worklet.one_d_cub_dpsi_inv = 1.0/dPsi;
  worklet.UseDeltaBScale = false;
  worklet.DeltaBScale = 1.0;
  worklet.UseBScale = false;
  worklet.BScale = 1.0;
  worklet.UseBOnly = false;

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayCopy(ds.GetCoordinateSystem().GetData(), coords);

  if (useUB)
  {
    /*
    invoker(worklet, seedsArray,
            locatorUB,
            cellSet, coords,
            As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
            B_RZP, psi,
            tracesArr, outR, outZ, outTheta, outPsi, outID);
    */
  }
  else
  {
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> vecField_RZP, vec_B_RZP, vec_dB_RZP;
    vecField_RZP.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    vec_B_RZP.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    vec_dB_RZP.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> dBr_dr, Br, dBz_dz;
    dBr_dr.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    Br.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    dBz_dz.Allocate(coords.GetNumberOfValues() * xgcParams.numPlanes);
    //Evaluate at nodes...
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
    //DRP
    invoker(worklet, coords, vec_B_RZP, vec_dB_RZP, vecField_RZP,
            dBr_dr, Br, dBz_dz,
            locator2L, cellSet,
            As_ff, dAs_ff_rzp, coeff_1D, coeff_2D, B_RZP);
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
    //B_RZP, psi,
    //tracesArr, outR, outZ, outTheta, outPsi, outID);
    for (vtkm::Id i = 0; i < 25; i++)
      std::cout<<i<<":: "<<vecField_RZP.ReadPortal().Get(i)<<std::endl;

    if (1)
    {
      vtkm::cont::ArrayHandle<vtkm::FloatDefault> meow;
      vtkm::cont::ArrayHandle<vtkm::Vec3f> v;
      meow.Allocate(coords.GetNumberOfValues());
      v.Allocate(coords.GetNumberOfValues());
      for (vtkm::Id i = 0; i < coords.GetNumberOfValues(); i++)
      {
        meow.WritePortal().Set(i, (double)i);
        {
          vtkm::Vec3f vRZP = vecField_RZP.ReadPortal().Get(i);
          vtkm::Vec3f vRPZ(vRZP[0], vRZP[2], vRZP[1]);

          v.WritePortal().Set(i, vRPZ);
        }
      }

      std::cout<<"********************************** remove this...."<<std::endl;
      vtkm::Id NP = xgcParams.numPlanes;
      //NP = 2;

      //create a 3D dataset...
      vtkm::cont::ArrayHandle<vtkm::Vec3f> coords3D;
      coords3D.Allocate(coords.GetNumberOfValues() * NP);
      auto coords3DPortal = coords3D.WritePortal();
      auto coordsPortal = coords.ReadPortal();
      vtkm::Id idx = 0;

      vtkm::FloatDefault phi = 0;
      vtkm::FloatDefault dPhi = vtkm::TwoPi()/static_cast<double>(xgcParams.numPlanes);
      for (vtkm::Id i = 0; i < NP; i++)
      {
        for (vtkm::Id j = 0; j < coords.GetNumberOfValues(); j++)
        {
          vtkm::Vec3f ptRZ = coordsPortal.Get(j);
          auto R = ptRZ[0];
          auto Z = ptRZ[1];

          vtkm::Vec3f ptRPZ(ptRZ[0], phi, ptRZ[1]);
          //vtkm::Vec3f ptXYZ(R*vtkm::Cos(phi), R*vtkm::Sin(phi), Z);

          coords3DPortal.Set(idx, ptRPZ);
          idx++;
        }
        phi += dPhi;
      }

      //make the cellset.
      auto cells2D = cellSet;
      std::vector<vtkm::Id> wedgeIds;
      for (vtkm::Id p = 0; p < NP-1; p++)
      {
        for (vtkm::Id i = 0; i < cellSet.GetNumberOfCells(); i++)
        {
          vtkm::Vec<vtkm::Id, 3> cellIds;
          cells2D.GetIndices(i, cellIds);

          vtkm::Id offset = p*xgcParams.numNodes;
          for (vtkm::Id k = 0; k < 3; k++)
            wedgeIds.push_back(offset + cellIds[k]);

          offset = (p+1)*xgcParams.numNodes;
          //if (p == (NP-1)) offset = 0;

          for (vtkm::Id k = 0; k < 3; k++)
            wedgeIds.push_back(offset + cellIds[k]);
        }
      }
      auto wedgeConnArray = vtkm::cont::make_ArrayHandle(wedgeIds, vtkm::CopyFlag::On);
      vtkm::cont::CellSetSingleType<> cells3D;
      cells3D.Fill(coords3D.GetNumberOfValues(), vtkm::CellShapeTagWedge().Id, 6, wedgeConnArray);

      //DRP
      vtkm::cont::DataSet dsOut;
      dsOut.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords3D));
      dsOut.SetCellSet(cells3D);
      //dsOut.AddPointField("meow", meow);
      //dsOut.AddPointField("v", v);
      dsOut.AddPointField("B_dB_rzp", vecField_RZP);
      dsOut.AddPointField("B_rzp", vec_B_RZP);
      dsOut.AddPointField("dB", vec_dB_RZP);
      dsOut.AddPointField("dBr_dr", dBr_dr);
      dsOut.AddPointField("Br", Br);
      dsOut.AddPointField("dBz_dz", dBz_dz);

      Save2Adios(dsOut, cellSet, xgcParams.numNodes, xgcParams.numPlanes, "output.bp");
      dsOut.PrintSummary(std::cout);
      vtkm::io::VTKDataSetWriter writer("B_dB_field.vtk");
      writer.WriteDataSet(dsOut);
    }


#if 0
    //create a grouped vec..
    using SeedArrType = vtkm::cont::ArrayHandle<vtkm::Particle>;
    vtkm::Id start = 0, end = seedsArray.GetNumberOfValues(), step=1;
    auto counter = vtkm::cont::make_ArrayHandleCounting(start, step, end);
    auto groupArr = vtkm::cont::make_ArrayHandleGroupVec<10>(counter);
    //auto groupArr = vtkm::cont::ArrayHandleGroupVec<vtkm::cont::ArrayHandle<vtkm::Particle>, 10>(seedsArray);
    //auto groupArr = vtkm::cont::make_ArrayHandleGroupVec<10>(seedsArray);
    vtkm::cont::ArrayHandle<vtkm::Id> offsetArray;
    //auto groupArr = vtkm::cont::make_ArrayHandleGroupVecVariable(seedsArray, offsetArray);
    std::cout<<"********************************* Call worklet: "<<std::endl;
    for (int i = 0; i < groupArr.GetNumberOfValues(); i++)
    {
      auto x = groupArr.ReadPortal().Get(i);
      //for (int j = 0; j < x.GetNumberOfComponents(); j++)
        std::cout<<"  "<<i<<" : "<<x<<std::endl;
    }

    invoker(worklet, seedsArray, groupArr,
            locator2L,
            cellSet, coords,
            As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
            B_RZP, psi,
            tracesArr, outR, outZ, outTheta, outPsi, outID);
#endif
  }

  //outID.SyncControlArray();
}
