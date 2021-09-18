//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <typeinfo>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/filter/GhostCellClassify.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/EulerIntegrator.h>
#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/worklet/particleadvection/Stepper.h>
#include <vtkm/worklet/testing/GenerateTestDataSets.h>

#include <vtkm/filter/Gradient.h>

#include <vtkm/io/VTKDataSetWriter.h>

#include <fides/DataSetReader.h>

#include <adios2.h>
#include <random>
#include <chrono>
#include <mpi.h>

adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

int numPlanes = -1;
int numNodes = -1;
int numTri = -1;

const double psi_x = 0.2661956235889000;
const double xpoint_r = 1.557545038402000;
const double xpoint_z = -1.177067412978000;
const int nWedge = 1;

class adiosS
{
public:
    adiosS() {}
    adiosS(adios2::ADIOS *adiosPtr,
           const std::string &fn,
           const std::string &ioNm,
           const std::map<std::string, std::string> &args) : ioName(ioNm)
    {
        std::string pathNm = ".";

        auto x = args.find("--dir")->second;
        if (x.size() > 0) pathNm = x;
        this->fileName = pathNm + "/" + fn;
        std::cout<<"Open: "<<this->fileName<<std::endl;
        this->io = adios2::IO(adiosPtr->DeclareIO(this->ioName));
        this->engine = io.Open(fileName, adios2::Mode::Read);
    }
    ~adiosS() { engine.Close(); }
    adiosS& operator=(const adiosS &a)
    {
        ioName = a.ioName;
        fileName = a.fileName;
        io = a.io;
        engine = a.engine;
        return *this;
    }

    std::string ioName, fileName;
    adios2::IO io;
    adios2::Engine engine;
};

void
CalcBField(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> A_s;
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.GetField("B").GetData().AsArrayHandle(b);
  ds.GetField("apars").GetData().AsArrayHandle(A_s);

  auto bPortal = b.ReadPortal();
  auto aPortal = A_s.ReadPortal();

  vtkm::Id n = b.GetNumberOfValues();
  std::vector<vtkm::Vec3f> As_bHat(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f bn = vtkm::Normal(bPortal.Get(i));
    vtkm::FloatDefault As = aPortal.Get(i);
    As_bHat[i] = bn * As;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("As_bHat", vtkm::cont::make_ArrayHandle(As_bHat, vtkm::CopyFlag::On)));

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetComputeVorticity(true);
  gradient.SetActiveField("As_bHat");
  gradient.SetOutputFieldName("grad_As_bHat");
  std::cout<<"Compute Grad"<<std::endl;
  ds = gradient.Execute(ds);
  std::cout<<"Compute Grad DONE"<<std::endl;

#if 0
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curl;
  ds.GetField("Vorticity").GetData().AsArrayHandle(curl);
  auto cPortal = curl.ReadPortal();
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  //This is in cartesian coordintes.....
  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f c = cPortal.Get(i);
    V[i] = bPortal.Get(i) + c;
  }
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));
#endif

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f,3>> grad;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("grad_As_bHat").GetData().AsArrayHandle(grad);

  auto cPortal = coords.ReadPortal();
  auto gPortal = grad.ReadPortal();
  auto vec = grad.ReadPortal().Get(0);

  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f val = bPortal.Get(i);
    vtkm::FloatDefault R = cPortal.Get(i)[0];
    auto g = gPortal.Get(i);

    //R: (1/R * dAz/dT  - dAT/dZ)
    //T: dAr/dZ - dAz/dr
    //Z: 1/R [ d(rAt)/dr - dAr/dT]
    vtkm::FloatDefault rv, tv, zv;
    rv = 1/R * g[2][1] - g[1][2];
    tv = g[0][2] - g[2][1];
    zv = 1/R * (R*g[1][0] - g[0][1]);
    val[0] += rv;
    val[1] += tv;
    val[2] += zv;

    V[i] = val;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
}

void
ReadScalar(adiosS* stuff,
           vtkm::cont::DataSet& ds,
           const std::string& vname,
           bool isXYZ,
           bool cylAdd,
           std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  if (cylAdd)
  {
    for (int i = 0; i < numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname,
                                          vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadVec(adiosS* stuff,
        vtkm::cont::DataSet& ds,
        const std::string& vname,
        bool isXYZ,
        bool cylAdd,
        std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> vec;
  int NP = numPlanes;
  if (cylAdd)
    NP = NP+1;

  for (int p = 0; p < NP; p++)
  {
    //R,Z,T in file:
    for (int i = 0; i < numNodes; i++)
      vec.push_back(vtkm::Vec3f(tmp[i*3+0], tmp[i*3+2], tmp[i*3+1]));
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname,
                                          vtkm::cont::make_ArrayHandle(vec, vtkm::CopyFlag::On)));
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff, adiosS* dataStuff, bool isXYZ, bool cylAdd)
{
  std::vector<double> rz;
  std::vector<int> conn;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);

  double dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;
  double phi = 0.0;

  int NP = numPlanes;
  if (!isXYZ && cylAdd)
    NP++;

  //coords
  for (int p = 0; p < NP; p++)
  {
    if (p == NP-1)
      phi = vtkm::TwoPi();

    vtkm::Vec3f pt;
    for (int i = 0; i < numNodes; i++)
    {
      double R = rz[i*2 +0];
      double Z = rz[i*2 +1];

      if (isXYZ)
      {
        pt[0] = R*cos(phi);
        pt[1] = R*sin(phi);
        pt[2] = Z;
      }
      else
      {
        pt[0] = R;
        pt[1] = phi;
        pt[2] = Z;
      }

      coords.push_back(pt);
    }
    phi += dPhi;
  }

  //cells
  for (int p = 0; p < numPlanes; p++)
  {
    if (!isXYZ && p == (numPlanes-1) && !cylAdd)
      break;

    for (int i = 0; i < numTri*3; i+=3)
    {
      int off = p*numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);

      off = (p+1)*(numNodes);
      if (isXYZ && p == (numPlanes-1))
        off = 0;

      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  auto grid = dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");

//  vtkm::cont::DataSet ds;
//  ds.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", vtkm::cont::make_ArrayHandle(coords, vtkm::CopyFlag::On)));
//  ds.SetCellSet(cellSet);
  return grid;
}

void
RunPoincare(const vtkm::cont::DataSet& ds)
{
  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;

  FieldHandle BField;
  ds.GetField("V").GetData().AsArrayHandle(BField);
  const vtkm::FloatDefault stepSize = 0.05;
  FieldType velocities(BField);
  GridEvalType eval(ds, velocities);
  Stepper rk4(eval, stepSize);

  vtkm::worklet::Poincare p;
  vtkm::Plane<> plane({0,3,0}, {0,1,0});

  vtkm::Id numSeeds = 20;
  vtkm::Id maxPunctures = 500;
  vtkm::Id maxSteps = maxPunctures*1000;

  vtkm::FloatDefault x0 = 2.9, x1 = 3.5;
  vtkm::FloatDefault dx = (x1-x0) / (float)(numSeeds-1);
  vtkm::FloatDefault x = x0;
  std::vector<vtkm::Particle> seeds;
  for (vtkm::Id id = 0; id < numSeeds; id++, x+=dx)
    seeds.push_back({vtkm::Particle({x, 2.95, 0}, id)});
  auto seedsArr = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

  auto t1 = std::chrono::high_resolution_clock::now();
  auto res = p.Run(rk4, seedsArr, plane, maxSteps, maxPunctures, true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;


  if (1)
  {
    std::ofstream outPts;
    outPts.open("points.txt");
    int nPts = res.Positions.GetNumberOfValues();
    auto portal = res.Positions.ReadPortal();
    for (int i = 0; i < nPts; i++)
    {
      auto pt = portal.Get(i);
      outPts<<pt[0]<<", "<<pt[1]<<", "<<pt[2]<<std::endl;
    }
    outPts.close();
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  std::cout<<std::endl<<std::endl;

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;

  if (argc < 4)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  std::string dataDir = std::string(argv[2]);
  vtkm::Id numSeeds = std::stoi(argv[3]);
  vtkm::Id maxPunctures = std::stoi(argv[4]);
  std::map<std::string, std::string> args;
  args["--dir"] = dataDir;

  adios = new adios2::ADIOS;
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", args);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", args);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", args);

  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield"];
  meshStuff->engine.BeginStep();
  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);

  bool isXYZ = false;
  bool cylAdd = true;
  auto ds = ReadMesh(meshStuff, dataStuff, isXYZ, cylAdd);
  ReadScalar(dataStuff, ds, "dpot", isXYZ, cylAdd);
  ReadScalar(dataStuff, ds, "apars", isXYZ, cylAdd);
  ReadVec(bfieldStuff, ds, "B", isXYZ, cylAdd, "/node_data[0]/values");

  CalcBField(ds);
  ds.PrintSummary(std::cout);

  RunPoincare(ds);
  return 0;


  std::cout<<"**** Dump file...."<<std::endl;
  vtkm::io::VTKDataSetWriter writer("simple.vtk");
  writer.WriteDataSet(ds);

  return 0;

}
