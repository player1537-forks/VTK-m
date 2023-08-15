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
#include <iomanip>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayCopy.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/Particle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/density_estimate/ParticleDensityNearestGridPoint.h>

#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/Algorithm.h>
#ifdef VTKM_CUDA
#include <vtkm/cont/cuda/internal/ScopedCudaStackSize.h>
#include <vtkm/cont/cuda/internal/DeviceAdapterAlgorithmCuda.h>
#endif

#include <adios2.h>
#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include <variant>
#include <string>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "XGCParameters.h"
#include "FindMaxR.h"
#include "EvalField.h"
#include "RunPoincare2.h"
#include "SavePoincare.h"
#include "XGCHelpers.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <utility>


adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

//#define VALGRIND
#ifdef VALGRIND
#include <valgrind/callgrind.h>
#endif

#define BUILD_POINC2

#ifdef BUILD_POINC1
#include "Poincare.h"
#endif
#ifdef BUILD_POINC2
#include "Poincare2.h"
#endif
#ifdef BUILD_POINC3
#include "Poincare3.h"
#include "ComputeB.h"
#include "ComputeBCell.h"
#endif

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  out<<"[";
  for (const auto& x : v)
    out<<x<<" ";
  out<<"]";
  return out;
}


//using ArgumentType = std::variant<std::vector<int>, std::vector<float>, std::vector<std::string>>;
bool
ParseArgs(int argc, char **argv, std::map<std::string, std::vector<std::string>> &args)
{
  args.clear();
  /*
  args["--dir"] = {""};
  args["--numSeeds"] = {"10"};
  args["--numPunc"] = {"10"};
  args["--stepSize"] = {"0.01"};
  args["--varname"] = {"V"};
  args["--seed"] = {"0", "0", "0"};
  */

//  int i = 1;
  std::string a0;
  std::vector<std::string> a1;
  //while (i < argc)

  for (int i = 1; i < argc; i++)
  {
    std::string tmp(argv[i]);
    if (tmp.find("--") != std::string::npos)
    {
      if (!a0.empty())
      {
        args[a0] = a1;
        a1.clear();
      }

      a0 = tmp;
      continue;
    }
    else
      a1.push_back(tmp);
  }
  //last argument.
  if (!a0.empty())
    args[a0] = a1;

  std::cout<<"ARGS\n";
  for (const auto& it : args)
  {
    std::cout<<it.first<<" : {";
    for (const auto& jt : it.second)
      std::cout<<jt<<" ";
    std::cout<<"}\n";
  }

  bool retVal = true;

  auto requiredArgs = {"--dir", "--output"};

  for (const auto& a : requiredArgs)
    if (args.find(a) == args.end())
    {
      std::cerr<<"Error. Argument missing: "<<a<<std::endl;
      retVal = false;
    }

  return retVal;
}


class adiosS
{
public:
  adiosS() {}
  adiosS(adios2::ADIOS *adiosPtr,
         const std::string &fn,
         const std::string &ioNm,
         const std::map<std::string, std::string> &args)
    : ioName(ioNm)
  {
    std::string pathNm = ".";

    auto x = args.find("--dir")->second;
    if (x.size() > 0)
      pathNm = x;
    this->fileName = pathNm + "/" + fn;
    std::cout<<"Open: "<<this->fileName<<std::endl;
    this->io = adios2::IO(adiosPtr->DeclareIO(this->ioName));
    this->engine = io.Open(this->fileName, adios2::Mode::Read);
  }

  adiosS(adios2::ADIOS *adiosPtr,
         const std::string &fn,
         const std::string &ioNm,
         const adios2::Mode& mode)
    : ioName(ioNm)
    , fileName(fn)
  {
    std::cout<<"Open: "<<this->fileName<<std::endl;
    this->io = adios2::IO(adiosPtr->DeclareIO(this->ioName));
    this->engine = io.Open(this->fileName, mode);
  }

  ~adiosS()
  {
    std::cout<<"Closing "<<this->fileName<<std::endl;
    this->engine.Close();
  }

  adiosS& operator=(const adiosS &a)
  {
    this->ioName = a.ioName;
    this->fileName = a.fileName;
    this->io = a.io;
    this->engine = a.engine;
    return *this;
  }

  static bool Exists(const std::string& fname,
                     const std::map<std::string, std::string> &args)
  {
    std::string pathNm = ".";
    auto x = args.find("--dir")->second;
    if (x.size() > 0) pathNm = x;
    std::string fullName = pathNm + "/" + fname;

    std::ifstream ifile;
    ifile.open(fullName);
    if (ifile)
      return true;
    return false;
  }

  template <typename T>
  int GetVarInt(const std::string& varName)
  {
    auto var = this->io.InquireVariable<T>(varName);

    T val;
    this->engine.Get(var, &val, adios2::Mode::Sync);
    return static_cast<int>(val);
  }

  int GetVarInt(const std::string& varName)
  {
    if (this->io.InquireVariable<int>(varName))
      return this->GetVarInt<int>(varName);
    else if (this->io.InquireVariable<long>(varName))
      return this->GetVarInt<long>(varName);

    throw std::runtime_error("Unable read int var: " + varName);
  }

  std::string ioName, fileName;
  adios2::IO io;
  adios2::Engine engine;
};

void
ReadOther(adiosS* stuff,
          vtkm::cont::DataSet& ds,
          const std::string& vname,
          std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

#ifdef VTKM_USE_DOUBLE_PRECISION
  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WholeMesh,
                                     tmp,
                                     vtkm::CopyFlag::On));
#else
  //convert from double to vtkm::FloatDefault
  std::vector<vtkm::FloatDefault> tmpF(tmp.size());
  for (std::size_t i = 0; i < tmp.size(); i++)
    tmpF[i] = static_cast<vtkm::FloatDefault>(tmp[i]);

  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WholeMesh,
                                     tmpF,
                                     vtkm::CopyFlag::On));
#endif
}

void
ReadScalar(adiosS* stuff,
           XGCParameters& xgcParams,
           vtkm::cont::DataSet& ds,
           const std::string& vname,
           std::string fileName="",
           bool add3D=false,
           bool addExtra=false)
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  vtkm::Id numPts = ds.GetNumberOfPoints();
  std::vector<double> tmpPlane(numPts);
  for (int i = 0; i < numPts; i++)
    tmpPlane[i] = tmp[i];

  //Normalize psi by eq_x_psi
  /*
  if (vname == "psi")
  {
    for (int i = 0; i < numPts; i++)
      tmpPlane[i] /= eq_x_psi;
  }
  */

  if (addExtra && add3D)
  {
    for (int i = 0; i < xgcParams.numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D", vtkm::cont::make_ArrayHandle(tmpPlane, vtkm::CopyFlag::On)));
  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname, vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadB(adiosS* stuff,
      XGCParameters& xgcParams,
      vtkm::cont::DataSet& ds)
{
  std::string fileName = "/node_data[0]/values";

  auto Bvar = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(Bvar, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> b_rzp;
  for (int i = 0; i < xgcParams.numNodes; i++)
  {
    vtkm::Vec3f v(tmp[i*3+0], tmp[i*3+1], tmp[i*3+2]);
    b_rzp.push_back(v);
  }

  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP",vtkm::cont::make_ArrayHandle(b_rzp, vtkm::CopyFlag::On)));
}

void
printEntry(int zi, int ri, const std::vector<double>& v, int vtkmNotUsed(NZ), int NR)
{
  int idx = (zi*NR + ri)*16;

  std::cout<<"X["<<zi<<" "<<ri<<"] ="<<std::endl;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
      std::cout<<" "<<v[idx + i*4 + j]<<" ";
    std::cout<<std::endl;
  }
}

void
ReadPsiInterp(adiosS* eqStuff,
              adiosS* interpStuff,
              adiosS* unitsStuff,
              vtkm::cont::DataSet& ds,
              XGCParameters& xgcParams,
              std::map<std::string, std::vector<std::string>>& args)
{
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mpsi"), &xgcParams.eq_mpsi, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mr"), &xgcParams.eq_mr, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mz"), &xgcParams.eq_mz, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_r"), &xgcParams.eq_axis_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_z"), &xgcParams.eq_axis_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_r"), &xgcParams.eq_min_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_r"), &xgcParams.eq_max_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_z"), &xgcParams.eq_min_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_z"), &xgcParams.eq_max_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_psi"), &xgcParams.eq_x_psi, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_r"), &xgcParams.eq_x_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_z"), &xgcParams.eq_x_z, adios2::Mode::Sync);
  unitsStuff->engine.Get(unitsStuff->io.InquireVariable<int>("sml_wedge_n"), &xgcParams.sml_wedge_n, adios2::Mode::Sync);

  ReadOther(eqStuff, ds, "eq_I");
  ReadOther(eqStuff, ds, "eq_psi_grid");
  ReadOther(eqStuff, ds, "eq_psi_rz");
  ReadOther(interpStuff, ds, "coeff_1D", "one_d_cub_psi_acoef");
  ReadOther(interpStuff, ds, "coeff_2D", "psi_bicub_acoef");

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> psiGrid;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(psiGrid);
  xgcParams.itp_min_psi = psiGrid.ReadPortal().Get(0);
  xgcParams.itp_max_psi = psiGrid.ReadPortal().Get(xgcParams.eq_mpsi-1);


#if 0
  std::vector<double> tmp2D;
//  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("one_d_cub_psi_acoef"),
//                          tmp1D, adios2::Mode::Sync);
  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("psi_bicub_acoef"),
                          tmp2D, adios2::Mode::Sync);

  ds.AddField(vtkm::cont::make_Field("coeff_2D",
                                     vtkm::cont::Field::Association::WholeMesh,
                                     arr_coeff2D,
                                     vtkm::CopyFlag::On));

  /*
  vtkm::Id n = ds.GetField("eq_psi_grid").GetData().GetNumberOfValues();
  std::vector<std::vector<double>> coef1D;
  coef1D.resize(n);
  for (int i = 0; i < n; i++)
    coef1D[i].resize(4);

  int idx = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 4; j++)
    {
      coef1D[i][j] = tmp1D[idx];
      idx++;
    }
  */

  std::cout<<"****** READ: "<<tmp2D[0]<<" "<<tmp2D[16]<<" "<<tmp2D[32]<<std::endl;
  std::cout<<"          :: "<<tmp2D[1]<<" "<<tmp2D[17]<<std::endl;

  int idx = 0;
  int nr = xgcParams.eq_mr-1, nz = xgcParams.eq_mz-1;

  std::cout<<"   coeff_1D.size= "<<tmp2D.size()<<std::endl;

  printEntry(0, 250, tmp2D, nz, nr);
  printEntry(1023, 250, tmp2D, nz, nr);
  printEntry(1023, 0, tmp2D, nz, nr);


  printEntry(167, 213, tmp2D, nz, nr);
  printEntry(213, 167, tmp2D, nz, nr);


  std::vector<std::vector<std::vector<std::vector<double>>>> coef2D;
  coef2D.resize(nz);
  for (int i = 0; i < nz; i++)
  {
    coef2D[i].resize(nr);
    for (int j = 0; j < nr; j++)
    {
      coef2D[i][j].resize(4);
      for (int k = 0; k < 4; k++)
      {
        coef2D[i][j][k].resize(4);
        for (int m = 0; m < 4; m++)
        {
          coef2D[i][j][k][m] = tmp2D[idx];
          idx++;
        }
      }
    }
  }

  idx = 0;
  std::vector<vtkm::FloatDefault> arr_coeff2D(nz*nr*4*4);
  for (int i = 0; i < nz; i++)
    for (int j = 0; j < nr; j++)
      for (int k = 0; k < 4; k++)
        for (int m = 0; m < 4; m++)
        {
          //arr_coeff2D[idx] = coef2D[i][j][k][m];
          arr_coeff2D[idx] = coef2D[i][j][m][k];
          idx++;
        }

  ds.AddField(vtkm::cont::make_Field("coeff_2D",
                                     vtkm::cont::Field::Association::WholeMesh,
                                     arr_coeff2D,
                                     vtkm::CopyFlag::On));


  //Reorder coeff2d to see if it's wrong...
  idx = 0;
  std::vector<std::vector<std::vector<std::vector<double>>>> coeff_2D;
  coeff_2D.resize(nz);
  for (int i = 0; i < nz; i++)
  {
    coeff_2D[i].resize(nr);
    for (int j = 0; j < nr; j++)
    {
      coeff_2D[i][j].resize(4);
      for (int k = 0; k < 4; k++)
      {
        coeff_2D[i][j][k].resize(4);
        for (int m = 0; m < 4; m++)
        {
          coeff_2D[i][j][k][m] = tmp2D[idx];
          idx++;
        }
      }
    }
  }


  if (args.find("--dumpPsiGrid") != args.end())
  {
    //Put this on a 2D grid for debugging...
    vtkm::Vec2f origin2D(xgcParams.eq_min_r, xgcParams.eq_min_z);
    vtkm::Vec2f spacing2D((xgcParams.eq_max_r-xgcParams.eq_min_r)/double(xgcParams.eq_mr-1), (xgcParams.eq_max_z-xgcParams.eq_min_z)/double(xgcParams.eq_mz-1));
    auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(xgcParams.eq_mr, xgcParams.eq_mz),
                                                          origin2D, spacing2D);

    std::vector<std::vector<std::vector<vtkm::FloatDefault>>> cij(4);
    for (int k = 0; k < 4; k++) cij[k].resize(4);

    for (int i = 0; i < nz; i++)
      for (int j = 0; j < nr; j++)
      {
        for (int k = 0; k < 4; k++)
          for (int m = 0; m < 4; m++)
            cij[k][m].push_back(coeff_2D[i][j][k][m]);
      }

    for (int k = 0; k < 4; k++)
    {
      for (int m = 0; m < 4; m++)
      {
        char nm[32];
        sprintf(nm, "c%d%d", k,m);
        //std::cout<<"Add cij: "<<nm<<" "<<cij[k][m].size()<<std::endl;
        ds2D.AddCellField(nm, cij[k][m]);
      }
    }

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> b3d;
    ds.GetField("eq_psi_rz").GetData().AsArrayHandle(arr);
    ds2D.AddPointField("eq_psi_rz", arr);
    vtkm::io::VTKDataSetWriter writer("psiGrid.vtk");
    writer.WriteDataSet(ds2D);
  }
#endif


  if (args.find("--dumpPsiGrid") != args.end())
  {
    //Put this on a 2D grid for debugging...
    vtkm::Vec2f origin2D(xgcParams.eq_min_r, xgcParams.eq_min_z);
    vtkm::Vec2f spacing2D((xgcParams.eq_max_r-xgcParams.eq_min_r)/double(xgcParams.eq_mr-1), (xgcParams.eq_max_z-xgcParams.eq_min_z)/double(xgcParams.eq_mz-1));
    auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(xgcParams.eq_mr, xgcParams.eq_mz),
                                                          origin2D, spacing2D);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> carr;
    ds.GetField("coeff_2D").GetData().AsArrayHandle(carr);
    auto carrPortal = carr.ReadPortal();

    std::vector<std::vector<std::vector<vtkm::FloatDefault>>> cij(4);
    for (int k = 0; k < 4; k++) cij[k].resize(4);

    int nr = xgcParams.eq_mr-1, nz = xgcParams.eq_mz-1;
    for (int zi = 0; zi < nz; zi++)
      for (int ri = 0; ri < nr; ri++)
      {
        vtkm::Id idx = (zi*nr + ri)*16;
        for (int k = 0; k < 4; k++)
          for (int m = 0; m < 4; m++)
          {
            cij[k][m].push_back(carrPortal.Get(idx + k*4 + m));
          }
      }

    for (int k = 0; k < 4; k++)
    {
      for (int m = 0; m < 4; m++)
      {
        char nm[32];
        sprintf(nm, "c%d%d", k,m);
        //std::cout<<"Add cij: "<<nm<<" "<<cij[k][m].size()<<std::endl;
        ds2D.AddCellField(nm, cij[k][m]);
      }
    }

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> b3d;
    ds.GetField("eq_psi_rz").GetData().AsArrayHandle(arr);
    ds2D.AddPointField("eq_psi_rz", arr);
    vtkm::io::VTKDataSetWriter writer("psiGrid.vtk");
    writer.WriteDataSet(ds2D);
  }
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff, XGCParameters& xgcParams)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  //points.
  for (int i = 0; i < xgcParams.numNodes; i++)
  {
    double R = rz[i*2 +0];
    double Z = rz[i*2 +1];

    vtkm::Vec3f ptRZ(R,Z,0);
    coords.push_back(ptRZ);
  }

  //cells
  for (int i = 0; i < xgcParams.numTri*3; i+=3)
  {
    int p0 = conn[i+0];
    int p1 = conn[i+1];
    int p2 = conn[i+2];
    connIds.push_back(p0);
    connIds.push_back(p1);
    connIds.push_back(p2);
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagTriangle(), 3, connIds, "coords");
}

vtkm::cont::DataSet
ReadMesh3D(adiosS* meshStuff, bool addExtra, XGCParameters& xgcParams)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  double dPhi = vtkm::TwoPi()/static_cast<double>(xgcParams.numPlanes);

  //points.
  double phi = 0;
  int NP = xgcParams.numPlanes;
  if (addExtra) NP++;
  for (int p = 0; p < NP; p++)
  {
    std::cout<<"ReadMesh3D: phi= "<<phi<<std::endl;
    for (int i = 0; i < xgcParams.numNodes; i++)
    {
      double R = rz[i*2 +0];
      double Z = rz[i*2 +1];

      vtkm::Vec3f pt(R,Z,phi);
      coords.push_back(pt);
    }
    phi += dPhi;
  }

  //cells
  for (int p = 0; p < NP-1; p++)
  {
    for (int i = 0; i < xgcParams.numTri*3; i+=3)
    {
      int off = p*xgcParams.numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0+off);
      connIds.push_back(p1+off);
      connIds.push_back(p2+off);

      off = (p+1)*(xgcParams.numNodes);
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
}

static bool Exists(const std::string& /*fname*/)
{
  return false;
  /*
  std::ifstream ifile;
  ifile.open(fname);
  if (ifile)
    return true;
  return false;
  */
}

void WriteHeader(const std::string& fname, const std::string& header)
{
  std::ofstream f(fname, std::ofstream::out);
  f<<header<<std::endl;
  f.close();
}

void GetOutputDirFile(std::map<std::string, std::vector<std::string>>& args,
                      std::string& dirName,
                      std::string& fileName)
{
  dirName = "./";
  fileName = args["--output"][0];

  //create a directory name based on the args.
  if (args.find("--autoDir") != args.end())
  {
    std::stringstream ss;
    ss<<"h"<<(double)std::fabs(std::atof(args["--stepSize"][0].c_str()));
    ss<<"_p"<<std::atoi(args["--numPunc"][0].c_str());
    if (args.find("--psiRange") != args.end())
    {
      auto vals = args["--psiRange"];
      ss<<"_pR"<<(double)std::atof(vals[0].c_str())<<"_"<<(double)std::atof(vals[1].c_str())<<"_"<<std::atoi(vals[2].c_str());
    }
    else if (args.find("--psiVals") != args.end())
    {
      auto vals = args["--psiRange"];
      ss<<"_pV"<<(int)vals.size();
    }
    else
    {
      ss<<"_Sx";
    }

    if (args.find("--thetaRange") != args.end())
    {
      auto vals = args["--thetaRange"];
      if (vals.size() == 1)
      {
        ss<<"_tR"<<std::atoi(vals[0].c_str());
      }
      else
      {
        ss<<"_tR"<<std::atoi(vals[0].c_str())<<"_"<<(double)std::atof(vals[1].c_str())<<"_"<<std::atoi(vals[2].c_str());
      }
    }
    else if (args.find("--thetaVals") != args.end())
    {
      auto vals = args["--thetaVals"];
      ss<<"_tV"<<(int)vals.size();
    }

    if (args.find("--deltaBScale") != args.end())
    {
      auto deltaBScale = std::atof(args["--deltaBScale"][0].c_str());
      ss<<"_dB"<<deltaBScale;
    }

    dirName = ss.str();

    bool dirExists = false;
    struct stat info;
    if (stat(dirName.c_str(), &info) != 0)
      dirExists = false;
    else if (info.st_mode & S_IFDIR)
      dirExists = true;

    if (!dirExists)
    {
      std::cout<<"Dir is not there! Create it."<<std::endl;
      std::string cmd = "mkdir " + dirName;
      system(cmd.c_str());
    }
  }
  dirName = dirName + "/";
  std::cout<<"OUTPUT= "<<dirName + fileName<<std::endl;
}

void
DumpSeeds(const std::vector<vtkm::Particle>& seeds,
          const std::vector<vtkm::Vec2f>& thetaPsi,
          const XGCParameters& xgcParams)
{
  //Save to a text file.
  std::ofstream f("seeds.txt", std::ofstream::out);
  f<<"ID, R, Z, Phi"<<std::endl;
  f<<std::setprecision(12);
  for (std::size_t i = 0; i < seeds.size(); i++)
    f<<seeds[i].ID<<", "<<seeds[i].Pos[0]<<", "<<seeds[i].Pos[2]<<", "<<seeds[i].Pos[1]<<std::endl;
  f.close();

  auto N = seeds.size();
  //save to ADIOS
  auto io = adios2::IO(adios->DeclareIO("seeds"));
  auto engine = io.Open("poincare.seeds.bp", adios2::Mode::Write);

  std::vector<std::size_t> shape{N}, offset{0}, size=shape;
  auto varID = io.DefineVariable<vtkm::Id>("ID", shape, offset, size);
  auto varR = io.DefineVariable<vtkm::FloatDefault>("R", shape, offset, size);
  auto varP = io.DefineVariable<vtkm::FloatDefault>("Phi", shape, offset, size);
  auto varZ = io.DefineVariable<vtkm::FloatDefault>("Z", shape, offset, size);

  auto varTheta = io.DefineVariable<vtkm::FloatDefault>("Theta", shape, offset, size);
  auto varPsi = io.DefineVariable<vtkm::FloatDefault>("Psi", shape, offset, size);
  auto varPsiN = io.DefineVariable<vtkm::FloatDefault>("PsiN", shape, offset, size);

  std::vector<vtkm::Id> ID(N);
  std::vector<vtkm::FloatDefault> R(N), P(N), Z(N), Psi(N, 0.0), PsiN(N), Theta(N, 0.0);

  bool haveThetaPsi = !thetaPsi.empty();
  for (std::size_t i = 0; i < N; i++)
  {
    const auto& s = seeds[i];
    ID[i] = s.ID;
    R[i] = s.Pos[0];
    P[i] = s.Pos[1];
    Z[i] = s.Pos[2];
    if (haveThetaPsi)
    {
      Theta[i] = thetaPsi[i][0];
      Psi[i] = thetaPsi[i][1];
      PsiN[i] = Psi[i] / xgcParams.eq_x_psi;
    }
  }

  engine.BeginStep();
  engine.Put<vtkm::Id>(varID, ID.data());
  engine.Put<vtkm::FloatDefault>(varR, R.data());
  engine.Put<vtkm::FloatDefault>(varP, P.data());
  engine.Put<vtkm::FloatDefault>(varZ, Z.data());
  engine.Put<vtkm::FloatDefault>(varTheta, Theta.data());
  engine.Put<vtkm::FloatDefault>(varPsi, Psi.data());
  engine.Put<vtkm::FloatDefault>(varPsiN, PsiN.data());
  engine.EndStep();
}

void
SaveOutput(std::map<std::string, std::vector<std::string>>& args,
           const std::vector<std::vector<vtkm::Vec3f>>& traces,
           const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outR,
           const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outZ,
           const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outTheta,
           const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& outPsi,
           const vtkm::cont::ArrayHandle<vtkm::Id>& outID,
           int timeStep=0)
{

  std::string outFileName, outDir;
  GetOutputDirFile(args, outDir, outFileName);

  std::string tracesNm, puncNm, puncThetaPsiNm, adiosNm, csvNm;
  tracesNm = outDir + outFileName + ".traces.txt";
  puncNm = outDir + outFileName + ".punc.vtk";
  puncThetaPsiNm = outDir + outFileName + ".punc.theta_psi.vtk";
  adiosNm = outDir + outFileName + ".bp";
  csvNm = outDir + outFileName + ".csv";
  std::cout<<"ADIOS NAME= "<<adiosNm<<std::endl;

  bool tExists = Exists(tracesNm);
  bool pExists = Exists(puncNm);
  bool ptpExists = Exists(puncThetaPsiNm);
  bool csvExists = Exists(csvNm);
  std::cout<<"EXISTS: "<<tExists<<" "<<pExists<<" "<<ptpExists<<std::endl;

  //Write headers.
  if (!traces.empty())
  {
    if (!tExists)
      WriteHeader(tracesNm, "ID, R, Z, T");

    std::ofstream outTraces;
    outTraces.open(tracesNm, std::ofstream::app);
    outTraces<<std::setprecision(15);

    //write traces
    for (int i = 0; i < (int)traces.size(); i++)
    {
      for (const auto& pt : traces[i])
      {
        auto R = pt[0];
        auto Z = pt[2];
        auto PHI_N = pt[1];
        if (PHI_N > vtkm::TwoPi())
        {
          while (PHI_N > vtkm::TwoPi())
            PHI_N -= vtkm::TwoPi();
        }
        else
        {
          while (PHI_N < 0)
            PHI_N += vtkm::TwoPi();
        }
        outTraces<<i<<", "<<R<<", "<<Z<<", "<<PHI_N<<std::endl;
      }
    }
  }

  //save to adios
  bool firstTime = false;
  std::string argString = "";
  if (adiosStuff.find("output") == adiosStuff.end())
  {
    adiosStuff["output"] = new adiosS(adios, adiosNm, "output", adios2::Mode::Write);
    firstTime = true;

    argString = "";
    for (const auto& x : args)
    {
      std::string arg = x.first;
      for (const auto& y : x.second)
        arg = arg + " " + y;
      arg = arg + "\n";
      argString = argString + arg;
    }
  }

  adiosS* outputStuff = adiosStuff["output"];

  std::size_t nPunc = static_cast<std::size_t>(std::atoi(args["--numPunc"][0].c_str()));
  std::size_t nPts = static_cast<std::size_t>(outR.GetNumberOfValues());

  nPts /= nPunc;

  //create a (R,Z,ID) and (Theta, Psi, ID) arrays.
  auto RBuff = vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>(outR).GetReadPointer();
  auto ZBuff = vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>(outZ).GetReadPointer();
  auto TBuff = vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>(outTheta).GetReadPointer();
  auto PBuff = vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>(outPsi).GetReadPointer();
  auto IDBuff = vtkm::cont::ArrayHandleBasic<vtkm::Id>(outID).GetReadPointer();

  if (firstTime)
  {
    outputStuff->io.DefineAttribute<std::string>("Arguments", argString);
    std::vector<std::size_t> shape = {nPts, nPunc}, offset = {0,0}, size = {nPts, nPunc};
    outputStuff->io.DefineVariable<vtkm::FloatDefault>("R", shape, offset, size);
    outputStuff->io.DefineVariable<vtkm::FloatDefault>("Z", shape, offset, size);
    outputStuff->io.DefineVariable<vtkm::FloatDefault>("Theta", shape, offset, size);
    outputStuff->io.DefineVariable<vtkm::FloatDefault>("Psi", shape, offset, size);
    outputStuff->io.DefineVariable<int>("TimeStep", {1}, {0}, {1});

    outputStuff->io.DefineVariable<vtkm::Id>("ID", shape, offset, size);
  }
  auto vR = outputStuff->io.InquireVariable<vtkm::FloatDefault>("R");
  auto vZ = outputStuff->io.InquireVariable<vtkm::FloatDefault>("Z");
  auto vT = outputStuff->io.InquireVariable<vtkm::FloatDefault>("Theta");
  auto vP = outputStuff->io.InquireVariable<vtkm::FloatDefault>("Psi");
  auto vID = outputStuff->io.InquireVariable<vtkm::Id>("ID");
  auto vTS = outputStuff->io.InquireVariable<int>("TimeStep");

  outputStuff->engine.BeginStep();
  outputStuff->engine.Put<vtkm::FloatDefault>(vR, &(RBuff[0]));
  outputStuff->engine.Put<vtkm::FloatDefault>(vZ, &(ZBuff[0]));
  outputStuff->engine.Put<vtkm::FloatDefault>(vT, &(TBuff[0]));
  outputStuff->engine.Put<vtkm::FloatDefault>(vP, &(PBuff[0]));
  outputStuff->engine.Put<vtkm::Id>(vID, IDBuff);
  outputStuff->engine.Put<int>(vTS, &timeStep);
  outputStuff->engine.EndStep();

  {
    //save to csv
    std::ofstream outCSV;
    outCSV.open(csvNm, std::ofstream::app);
    outCSV << std::setprecision(15);

    if(!csvExists)
      outCSV << "ID,R,Z,T,P\n";

    auto rpID = outID.ReadPortal();
    auto rpR = outR.ReadPortal();
    auto rpZ = outZ.ReadPortal();
    auto rpT = outTheta.ReadPortal();
    auto rpP = outPsi.ReadPortal();

    for(int i = 0; i < outR.GetNumberOfValues(); ++i){
      outCSV << rpID.Get(i) << ",";
      outCSV <<  rpR.Get(i) << ",";
      outCSV <<  rpZ.Get(i) << ",";
      outCSV <<  rpT.Get(i) << ",";
      outCSV <<  rpP.Get(i) << "\n";
    }
    outCSV.close();
  }
}

void
Poincare(const vtkm::cont::DataSet& ds,
         XGCParameters& xgcParams,
         std::vector<vtkm::Particle>& seeds,
         std::map<std::string, std::vector<std::string>>& args,
         int timeStep=0)
{
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff, coeff_1D, coeff_2D, psi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_rzp, B_Norm_rzp, dAs_ff_rzp;
  //ds.GetField("As_phi_ff").GetData().AsArrayHandle(As_ff);
  std::cout<<"DRP: Fix me?? "<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
  ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff_1D);
  ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);
  ds.GetField("psi2D").GetData().AsArrayHandle(psi);
  ds.GetField("B_RZP").GetData().AsArrayHandle(B_rzp);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> tracesArr;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> outR, outZ, outTheta, outPsi;
  vtkm::cont::ArrayHandle<vtkm::Id> outID;
  vtkm::Id nSeeds = seeds.size();

  vtkm::Id numPunc = std::atoi(args["--numPunc"][0].c_str());
  vtkm::cont::ArrayHandleConstant<vtkm::Id> initIds(-1, numPunc*nSeeds);
  outR.Allocate(numPunc*nSeeds);
  outZ.Allocate(numPunc*nSeeds);
  outTheta.Allocate(numPunc*nSeeds);
  outPsi.Allocate(numPunc*nSeeds);
  outID.Allocate(numPunc*nSeeds);
  vtkm::cont::ArrayCopy(initIds, outID);

  auto start = std::chrono::steady_clock::now();
  RunPoincare2(ds, seeds, xgcParams, args,
               As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
               B_rzp, psi,
               tracesArr, outR, outZ, outTheta, outPsi, outID);
  auto end = std::chrono::steady_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout<<"PoincareTime= "<<elapsed_seconds.count()<<std::endl;
  //std::cout<<"vtkm::cont::Timer= "<<timer.GetElapsedTime()<<std::endl;
  std::cout<<"outputR.size()= "<<outR.GetNumberOfValues()<<std::endl;
  std::cout<<"outputTheta.size()= "<<outTheta.GetNumberOfValues()<<std::endl;
  std::cout<<"punctureID.size()= "<<outID.GetNumberOfValues()<<std::endl;

  //Save output
  std::vector<std::vector<vtkm::Vec3f>> traces;

  std::vector<std::vector<vtkm::Vec3f>> res;
  bool useTraces = false;
  if (args.find("--traces") != args.end()) useTraces = std::atoi(args["--traces"][0].c_str());
  if (useTraces)
  {
    traces.resize(1);
    auto portal = tracesArr.ReadPortal();
    vtkm::Id n = portal.GetNumberOfValues();
    for (vtkm::Id i = 0; i < n; i++)
    {
      auto v = portal.Get(i);
      if (v[2] > -1)
        traces[0].push_back(v);
    }
  }

  SaveOutput(args, traces, outR, outZ, outTheta, outPsi, outID, timeStep);

  if (args.find("--dumpDensity") != args.end())
  {
    std::string outFileName = args["--output"][0];
    auto vals = args["--dumpDensity"];
    vtkm::Id nx = std::atoi(vals[0].c_str());
    vtkm::Id ny = std::atoi(vals[1].c_str());

    //Create particle density
    vtkm::Id3 cellDims(nx, ny, 1);
    vtkm::FloatDefault dR = xgcParams.eq_max_r - xgcParams.eq_min_r;
    vtkm::FloatDefault dZ = xgcParams.eq_max_z - xgcParams.eq_min_z;
    vtkm::Vec3f origin(xgcParams.eq_min_r, xgcParams.eq_min_z, 0);
    vtkm::Vec3f spacing(dR/static_cast<vtkm::FloatDefault>(cellDims[0]-1),
                        dZ/static_cast<vtkm::FloatDefault>(cellDims[1]-1),
                        1);

    //Create particle dataset.
    vtkm::Id N = outR.GetNumberOfValues();
    vtkm::cont::ArrayHandle<vtkm::Vec3f> positions;
    positions.Allocate(N);
    auto pPortal = positions.WritePortal();
    auto rPortal = outR.ReadPortal();
    auto zPortal = outZ.ReadPortal();
    for (vtkm::Id i = 0; i < N; i++)
      pPortal.Set(i, vtkm::Vec3f(rPortal.Get(i), zPortal.Get(i), 0));
    //std::cout<<"Particle Pts: "<<positions.GetNumberOfValues()<<std::endl;

    vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleIndex(N), connectivity);
    auto particleDataSet = vtkm::cont::DataSetBuilderExplicit::Create(
      positions, vtkm::CellShapeTagVertex{}, 1, connectivity);

    //Save out the particles.
    //Add ID, psi0, punc
    std::string particleName = outFileName + "particles.vtk";
    vtkm::io::VTKDataSetWriter writerP(particleName);
    writerP.WriteDataSet(particleDataSet);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> mass;
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleConstant(1, N), mass);
    particleDataSet.AddCellField("mass", mass);

    vtkm::filter::density_estimate::ParticleDensityNearestGridPoint pDensityFilter;
    pDensityFilter.SetDimension(cellDims);
    pDensityFilter.SetOrigin(origin);
    pDensityFilter.SetSpacing(spacing);

    pDensityFilter.SetActiveField("mass");
    //pDensityFilter.SetComputeNumberDensity(true);
    //pDensityFilter.SetDivideByVolume(false);
    auto density3D = pDensityFilter.Execute(particleDataSet);

    //Density filter doesn't work for 2D, so create a 2D dataset.
    vtkm::Id2 dims2D(nx+1, ny+1);
    vtkm::Vec2f origin2D(origin[0], origin[1]);
    vtkm::Vec2f spacing2D(spacing[0], spacing[1]);
    auto density2D = vtkm::cont::DataSetBuilderUniform::Create(dims2D, origin2D, spacing2D);
    density2D.AddField(density3D.GetField("density"));

    std::string densityName = outFileName + ".vtk";
    vtkm::io::VTKDataSetWriter writer(densityName);
    writer.WriteDataSet(density2D);
  }
}

vtkm::cont::ArrayHandle<vtkm::FloatDefault>
CalcAs(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As, XGCParameters& xgcParams)
{
//  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

//  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);
  vtkm::Id numAs = As.GetNumberOfValues();
  auto asPortal = As.ReadPortal();

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_arr.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
      As_arr[p][i].resize(xgcParams.numNodes);
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      //vtkm::Vec3f cBhat = cbPortal.Get(n);
      for (int i = 0; i < 2; i++)
      {
        auto as = asPortal.Get(idx);
        //std::cout<<"As: "<<as<<" cBhat= "<<cBhat<<" :: "<<vtkm::Magnitude(cBhat)<<"  val= "<<val<<std::endl;
        //As_curlBhat[p][i][n] = val;
        As_arr[p][i][n] = as;
        idx++;
      }
    }
  }

  //flatten to 1d index.
  idx = 0;
  std::vector<vtkm::FloatDefault> arrAs(numAs);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < xgcParams.numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        idx++;
      }
    }
  }

  return vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);

  /*
  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WholeMesh,
                                     arrAs,
                                     vtkm::CopyFlag::On));
  */
}

vtkm::cont::ArrayHandle<vtkm::Vec3f>
Calc_dAs(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& dAsAH, XGCParameters& xgcParams)
{
  //Assumed that dAs_phi_ff is R,Z,Phi
  //vtkm::cont::ArrayHandle<vtkm::FloatDefault> dAsAH;
  //ds.GetField("dAs_phi_ff").GetData().AsArrayHandle(dAsAH);
  int nVals = dAsAH.GetNumberOfValues()/3;

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> dAs_ff;
  dAs_ff.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    dAs_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      dAs_ff[p][i].resize(xgcParams.numNodes);
  }


  int idx = 0;
  auto dAsPortal = dAsAH.ReadPortal();

  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      for (int c = 0; c < 3; c++) //vec coords
      {
        for (int i = 0; i < 2; i++) //idx
        {
          vtkm::FloatDefault val = dAsPortal.Get(idx);
          idx++;
          int cc = c;

          /*
          //Change to R,Phi,Z
          //swap phi and z
          if (c == 1) cc = 2;
          else if (c == 2) cc = 1;
          */

          dAs_ff[p][i][n][cc] = val;
        }
      }
    }
  }

  //
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> VEC_ff;
  VEC_ff.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    VEC_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC_ff[p][i].resize(xgcParams.numNodes);
  }

  idx = 0;
  std::vector<vtkm::Vec3f> arr_ff(nVals);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < xgcParams.numNodes; n++)
      {
        vtkm::Vec3f val = VEC_ff[p][i][n];
        val = dAs_ff[p][i][n];
        arr_ff[idx] = val;
        idx++;
      }
    }
  }

  //do some testing..
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    int off0 = p*xgcParams.numNodes*2;
    int off1 = p*xgcParams.numNodes*2 + xgcParams.numNodes;
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      auto V0 = dAs_ff[p][0][n];
      auto V1 = dAs_ff[p][1][n];

      auto X0 = arr_ff[off0 + n];
      auto X1 = arr_ff[off1 + n];
      auto diff = vtkm::Magnitude(V0-X0) + vtkm::Magnitude(V1-X1);
      if (diff > 0)
        std::cout<<"ERROR: "<<V0<<" : "<<V1<<"  diff0= "<<(V0-X0)<<" diff1= "<<(V1-X1)<<std::endl;

    }
  }

  return vtkm::cont::make_ArrayHandle(arr_ff, vtkm::CopyFlag::On);
  /*
  ds.AddField(vtkm::cont::make_Field("dAs_ff_rzp",
                                     vtkm::cont::Field::Association::WholeMesh,
                                     arr_ff,
                                     vtkm::CopyFlag::On));
  */
}

template <typename T>
void
UpdateField(vtkm::cont::DataSet& ds,
            const std::string& vName,
            vtkm::cont::ArrayHandle<T> &var)
{
  if (ds.HasField(vName))
  {
    vtkm::cont::ArrayHandle<T> oldVar;
    ds.GetField(vName).GetData().AsArrayHandle(oldVar);

    vtkm::cont::Algorithm::Copy(var, oldVar);
  }
  else
    ds.AddField(vtkm::cont::Field(vName,
                                  vtkm::cont::Field::Association::WholeMesh,
                                  var));
}

void
ReadTurbData(adiosS* turbStuff,
             std::map<std::string, std::vector<std::string>>& args,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As_arr,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& dAs_arr)
{
  std::string As_phiName = "As_phi_ff", dAs_phiName = "dAs_phi_ff";
  if (args.find("--AsVarName") != args.end())
    As_phiName = args["--AsVarName"][0];
  if (args.find("--dAsVarName") != args.end())
    dAs_phiName = args["--dAsVarName"][0];

  auto asV = turbStuff->io.InquireVariable<double>(As_phiName.c_str());
  auto dAsV = turbStuff->io.InquireVariable<double>(dAs_phiName.c_str());

  std::vector<double> arrAs, arrdAs;
  turbStuff->engine.Get(asV, arrAs, adios2::Mode::Sync);
  turbStuff->engine.Get(dAsV, arrdAs, adios2::Mode::Sync);

  bool useTurbulence = true;
  if (args.find("--turbulence") != args.end())
    useTurbulence = std::atoi(args["--turbulence"][0].c_str());
  if (useTurbulence == false)
  {
    for (auto& x : arrAs)
      x = 0.0;
    for (auto& x : arrdAs)
      x = 0.0;
  }

  As_arr = vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);
  dAs_arr = vtkm::cont::make_ArrayHandle(arrdAs, vtkm::CopyFlag::On);
  //auto AsPhi = vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);
  //auto dAsPhi = vtkm::cont::make_ArrayHandle(arrdAs, vtkm::CopyFlag::On);
  //UpdateField(ds, "As_phi_ff", AsPhi);
  //UpdateField(ds, "dAs_phi_ff", dAsPhi);
}

vtkm::cont::DataSet
ReadStaticData(std::map<std::string, std::vector<std::string>>& args,
               XGCParameters& xgcParams,
               const std::string& coeffFile)
{
  if (adios != nullptr)
  {
    std::cerr<<"Re-reading static data!!!"<<std::endl;
  }

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  if (adios == nullptr)
  {
    if (args.find("--xml") == args.end())
      adios = new adios2::ADIOS;
    else
      adios = new adios2::ADIOS(args["--xml"][0]);
  }

  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "bfield", adiosArgs);
  adiosStuff["coeff"] = new adiosS(adios, coeffFile, "coeff", adiosArgs);
  adiosStuff["units"] = new adiosS(adios, "xgc.units.bp", "units", adiosArgs);

  auto meshStuff = adiosStuff["mesh"];
  auto equilStuff = adiosStuff["equil"];
  auto coeffStuff = adiosStuff["coeff"];
  auto bfieldStuff = adiosStuff["bfield"];
  auto unitsStuff = adiosStuff["units"];

  meshStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  coeffStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();
  unitsStuff->engine.BeginStep();

  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &xgcParams.numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &xgcParams.numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  xgcParams.psi_min = xgcParams.psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < xgcParams.psi_min) xgcParams.psi_min = p;
    if (p > xgcParams.psi_max) xgcParams.psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<xgcParams.psi_min<<" "<<xgcParams.psi_max<<std::endl;
  std::cout<<"********                  PSIN m/M = "<<xgcParams.psi_min*xgcParams.eq_x_psi<<" "<<xgcParams.psi_max*xgcParams.eq_x_psi<<std::endl;

  auto ds = ReadMesh(meshStuff, xgcParams);
  ReadPsiInterp(equilStuff, coeffStuff, unitsStuff, ds, xgcParams, args);
  ReadScalar(meshStuff, xgcParams, ds, "psi");
  ReadB(bfieldStuff, xgcParams, ds);

  //Add a normalized psi
  for (auto& psi : psiVals)
    psi /= xgcParams.eq_x_psi;
  ds.AddField(vtkm::cont::make_FieldPoint("psiNorm",
                                          vtkm::cont::make_ArrayHandle(psiVals, vtkm::CopyFlag::On)));

  meshStuff->engine.EndStep();
  equilStuff->engine.EndStep();
  coeffStuff->engine.EndStep();
  bfieldStuff->engine.EndStep();

  return ds;
}

vtkm::cont::DataSet
ReadDataSet_ORIG(std::map<std::string, std::vector<std::string>>& args, XGCParameters& xgcParams)
{
  auto ds = ReadStaticData(args, xgcParams, "xgc.bfield-all.bp");

  //Get the data...
  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);

  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield-all"];

  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  xgcParams.numPlanes = dataStuff->GetVarInt("nphi");

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
  ReadTurbData(bfieldStuff, args, As_arr, dAs_arr);

  dataStuff->engine.EndStep();
  bfieldStuff->engine.EndStep();

  auto As = CalcAs(As_arr, xgcParams);
  auto dAs = Calc_dAs(dAs_arr, xgcParams);
  UpdateField(ds, "As_ff", As);
  UpdateField(ds, "dAs_ff_rzp", dAs);

  return ds;
}

vtkm::cont::DataSet
ReadDataSet(std::map<std::string, std::vector<std::string>>& args, XGCParameters& xgcParams)
{
  vtkm::cont::DataSet ds;

  if (args.find("--test") != args.end())
  {
    ds = ReadDataSet_ORIG(args, xgcParams);
  }
  else if (args.find("--InputVTK") != args.end())
  {
    std::string fname = args["--InputVTK"][0];
    vtkm::io::VTKDataSetReader reader(fname);
    std::cout<<"Reading input data from: "<<fname<<std::endl;

    ds = reader.ReadDataSet();
    ds.PrintSummary(std::cout);
  }
  else
  {
    std::string coeffFile = "xgc.poincare_init.bp";
    if (args.find("--CoeffFile") != args.end())
      coeffFile = args["--CoeffFile"][0];

    ds = ReadStaticData(args, xgcParams, coeffFile);

    if (args.find("--AsFilePath") != args.end())
    {
      auto xgc3DFilePath = args["--AsFilePath"][0];
      adiosStuff["data"] = new adiosS(adios, xgc3DFilePath, "3d", adios2::Mode::Read);
    }
    else
    {
      std::string xgc3DFile = "xgc.3d.bp";
      if (args.find("--3DFile") != args.end())
        xgc3DFile = args["--3DFile"][0];

      std::map<std::string, std::string> adiosArgs;
      adiosArgs["--dir"] = args["--dir"][0];
      adiosStuff["data"] = new adiosS(adios, xgc3DFile, "3d", adiosArgs);
    }

    auto dataStuff = adiosStuff["data"];
    dataStuff->engine.BeginStep();
    xgcParams.numPlanes = dataStuff->GetVarInt("nphi");

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
    ReadTurbData(dataStuff, args, As_arr, dAs_arr);
    dataStuff->engine.EndStep();

    auto As = CalcAs(As_arr, xgcParams);
    auto dAs = Calc_dAs(dAs_arr, xgcParams);
    UpdateField(ds, "As_ff", As);
    UpdateField(ds, "dAs_ff_rzp", dAs);
  }

  if (args.find("--dumpGrid") != args.end())
  {
    ds.PrintSummary(std::cout);
    vtkm::io::VTKDataSetWriter writer("grid.vtk");
    writer.WriteDataSet(ds);
  }

  return ds;
}

int oneDBlocks = 16;
int threadsPerBlock = 16;
#ifdef VTKM_CUDA
vtkm::cont::cuda::ScheduleParameters
mySchedParams(char const* name,
              int major,
              int minor,
              int multiProcessorCount,
              int maxThreadsPerMultiProcessor,
              int maxThreadsPerBlock)
{
  vtkm::cont::cuda::ScheduleParameters p;
  p.one_d_blocks = oneDBlocks;
  p.one_d_threads_per_block = threadsPerBlock;

  return p;
}
#endif

/*
template <typename Coeff_1DType>
VTKM_EXEC
vtkm::FloatDefault I_interpol(const vtkm::FloatDefault& psi,
                              const int& ideriv,
                              const vtkm::FloatDefault& one_d_cub_dpsi_inv,
                              const vtkm::Id& ncoeff,
                              const Coeff_1DType& coeff_1D)
{
  vtkm::FloatDefault pn = psi * one_d_cub_dpsi_inv;
  int ip=floor(pn);
  ip=std::min(std::max(ip,0), int(ncoeff-1));
  vtkm::FloatDefault wp=pn-(vtkm::FloatDefault)(ip);

  int idx = ip*4;

  //vtkm::FloatDefault acoef[4];
  //acoef[0] = one_d_cub_acoef(ip).coeff[0];
  //acoef[1] = one_d_cub_acoef(ip).coeff[1];
  //acoef[2] = one_d_cub_acoef(ip).coeff[2];
  //acoef[3] = one_d_cub_acoef(ip).coeff[3];

  const vtkm::FloatDefault acoef[4] = {coeff_1D.Get(idx+0),
                                       coeff_1D.Get(idx+1),
                                       coeff_1D.Get(idx+2),
                                       coeff_1D.Get(idx+3)};

  vtkm::FloatDefault iVal = 0.0;
  if (ideriv==0)
    iVal = acoef[0]+(acoef[1]+(acoef[2]+acoef[3]*wp)*wp)*wp;
  else if (ideriv==1)
    iVal = (acoef[1]+(2.0*acoef[2]+3.0*acoef[3]*wp)*wp)*one_d_cub_dpsi_inv;

  return iVal; // * this->sml_bp_sign; = -1 !!
}
*/

template <typename CoeffType>
vtkm::FloatDefault
InterpolatePsi(const vtkm::Vec2f& ptRZ,
               const CoeffType& coeff,
               const vtkm::Id& ncoeff,
               const vtkm::Id2& nrz,
               const vtkm::Vec2f& rzmin,
               const vtkm::Vec2f& drz,
               const XGCParameters& /*xgcParams*/,
               vtkm::Vec3f& B0)
{
  vtkm::FloatDefault psi = 0;
  /*
  std::vector<double> eq_rgrid, eq_zgrid, rc_cub, zc_cub;
  double dr = (xgcParams.eq_max_r - xgcParams.eq_min_r) / double(xgcParams.eq_mr-1);
  double dz = (xgcParams.eq_max_z - xgcParams.eq_min_z) / double(xgcParams.eq_mz-1);
  for (int i = 0; i < xgcParams.eq_mr; i++)
  {
    double v = xgcParams.eq_min_r + dr * double(i);
    eq_rgrid.push_back(v);
  }
  for (int i = 0; i < xgcParams.eq_mz; i++)
  {
    double v = xgcParams.eq_min_z + dz * double(i);
    eq_zgrid.push_back(v);
  }

  for (int i = 0; i < xgcParams.eq_mr-1; i++)
    rc_cub.push_back((eq_rgrid[i] + eq_rgrid[i+1])/2.0);
  for (int i = 0; i < xgcParams.eq_mz-1; i++)
    zc_cub.push_back((eq_zgrid[i] + eq_zgrid[i+1])/2.0);
  */

  int r_i = XGCHelper::GetIndex(ptRZ[0], nrz[0], rzmin[0], 1.0/drz[0]);
  int z_i = XGCHelper::GetIndex(ptRZ[1], nrz[1], rzmin[1], 1.0/drz[1]);
  vtkm::FloatDefault Rc = rzmin[0] + (vtkm::FloatDefault)(r_i)*drz[0];
  vtkm::FloatDefault Zc = rzmin[1] + (vtkm::FloatDefault)(z_i)*drz[1];
  //auto rc2 = rc_cub[r_i];
  //auto zc2 = zc_cub[z_i];

  auto Rc_1 = Rc + drz[0];
  auto Zc_1 = Zc + drz[1];
  Rc = (Rc + Rc_1) * 0.5;
  Zc = (Zc + Zc_1) * 0.5;

  //vtkm::Id offset = (r_i * ncoeff + z_i) * 16; //DRP
  vtkm::Id offset = (z_i * ncoeff + r_i) * 16;

  /*
  std::cout<<"InterpolatePsi: "<<ptRZ<<std::endl;
  std::cout<<"  i/j= "<<r_i<<" "<<z_i<<std::endl;
  std::cout<<"  ncoeff= "<<ncoeff<<std::endl;
  std::cout<<"  offset= "<<offset<<std::endl;
  std::cout<<"  Rc/Zc= "<<Rc<<" "<<Zc<<std::endl;
  */

  vtkm::FloatDefault dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
  XGCHelper::EvalBicub2(ptRZ[0], ptRZ[1], Rc, Zc, offset, coeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);

  auto R = ptRZ[0];
  B0[0] = -dpsi_dz / R;
  B0[1] = dpsi_dr / R;
  B0[2] = 0;

  /*
  std::cout<<" psi= "<<psi<<std::endl;
  std::cout<<" dPsi_d= "<<dpsi_dr<<" "<<dpsi_dz<<std::endl;
  std::cout<<" d2Psi_d= "<<d2psi_drdz<<" "<<d2psi_d2r<<" "<<d2psi_d2z<<std::endl;
  */

  if (psi < 0) //DRP
    psi = 0;

  return psi;
}

template <typename CoeffType>
vtkm::FloatDefault
InterpolatePsi(const vtkm::Vec2f& ptRZ,
               const CoeffType& coeff,
               const vtkm::Id& ncoeff,
               const vtkm::Id2& nrz,
               const vtkm::Vec2f& rzmin,
               const vtkm::Vec2f& drz,
               const XGCParameters& xgcParams)
{
  vtkm::Vec3f B0;
  return InterpolatePsi(ptRZ, coeff, ncoeff, nrz, rzmin, drz, xgcParams, B0);
}

std::vector<vtkm::Particle>
GenerateNormalizedFromThetaPsiSeedsPairs(std::map<std::string, std::vector<std::string>>& args,
                      const vtkm::cont::DataSet& ds,
                      XGCParameters& xgcParams,
                      std::vector<vtkm::Vec2f>& seedsThetaPsi)

{
  std::vector<vtkm::Particle> seeds;
  std::vector<std::string> vals;
  std::vector<vtkm::FloatDefault> thetaVals, psiVals;
  const vtkm::FloatDefault degToRad = vtkm::Pi()/180.0;

  vals = args["--psiVals"];
  for (const auto& v: vals)
    psiVals.push_back(std::atof(v.c_str()) * xgcParams.eq_x_psi);

  vals = args["--thetaVals"];
  for (const auto& v: vals)
    thetaVals.push_back(std::atof(v.c_str()) * degToRad);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxR;
  FindMaxR(ds, xgcParams, thetaVals, maxR);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  ds.GetField("coeff_2D").GetData().AsArrayHandle(arr);
  auto coeff = arr.ReadPortal();

  auto maxRPortal = maxR.ReadPortal();
  vtkm::Id ncoeff = xgcParams.eq_mr-1;
  vtkm::Id2 nrz(xgcParams.eq_mr, xgcParams.eq_mz);
  vtkm::Vec2f rzmin(xgcParams.eq_min_r, xgcParams.eq_min_z);
  vtkm::Vec2f drz((xgcParams.eq_max_r-xgcParams.eq_min_r)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mr-1),
                  (xgcParams.eq_max_z-xgcParams.eq_min_z)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mz-1));

  const vtkm::FloatDefault dR = 0.000001;

  vtkm::Id ID = 0;
  for (std::size_t i = 0; i < thetaVals.size(); i++)
  {
    auto psiTarget = psiVals[i];
    auto theta = thetaVals[i];
    auto cost = vtkm::Cos(theta), sint = vtkm::Sin(theta);

    vtkm::FloatDefault r0 = 0;
    vtkm::FloatDefault r1 = r0;
    vtkm::FloatDefault psi = 0;
    vtkm::FloatDefault R, Z;

    r1 = r0+dR;
    bool found = false;
    while (!found && r1 < maxRPortal.Get(i)) {
      R = xgcParams.eq_axis_r + r1*cost;
      Z = xgcParams.eq_axis_z + r1*sint;
      psi = InterpolatePsi({R,Z}, coeff, ncoeff, nrz, rzmin, drz, xgcParams);
      if (psi >= psiTarget) {
        found = true;
        break;
      }

      r0 += dR;
      r1 += dR;
    }
    if (!found) continue;
    auto diffPsi = psi - psiTarget;

    //Now, do a binary search to find psi between (r0, r1)
    int cnt = 0;
    while (diffPsi > 1e-6 && cnt < 100) {
      auto rMid = (r0+r1)/2.0;
      R = xgcParams.eq_axis_r + rMid*cost;
      Z = xgcParams.eq_axis_z + rMid*sint;
      psi = InterpolatePsi({R,Z}, coeff, ncoeff, nrz, rzmin, drz, xgcParams);

      if (psi < psiTarget) r0 = rMid;// mid is inside, range = (rmid, r1)
      else r1 = rMid; //mid is outside, range = (r0, rmid)
          
      diffPsi = vtkm::Abs(psi - psiTarget);
      cnt++;
    }

    vtkm::Vec3f pt_rpz(R, 0, Z);
    vtkm::Particle p({pt_rpz, ID++});
    seeds.push_back(p);
    seedsThetaPsi.push_back({theta, psiTarget});
  }

  return seeds;
}


std::vector<vtkm::Particle>
GenerateThetaPsiSeeds(std::map<std::string, std::vector<std::string>>& args,
                      const vtkm::cont::DataSet& ds,
                      XGCParameters& xgcParams,
                      std::vector<vtkm::Vec2f>& seedsThetaPsi)

{
  std::vector<vtkm::Particle> seeds;

  std::vector<std::string> vals;
  std::vector<vtkm::FloatDefault> thetaVals, psiVals;
  const vtkm::FloatDefault degToRad = vtkm::Pi()/180.0;

  if (args.find("--psiRange") != args.end())
  {
    vals = args["--psiRange"];
    vtkm::FloatDefault psiNorm0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault psiNorm1 = std::atof(vals[1].c_str());
    int numPts = std::atoi(vals[2].c_str());

    auto psiMin = psiNorm0 * xgcParams.eq_x_psi, psiMax = psiNorm1 * xgcParams.eq_x_psi;
    auto dPsi = (psiMax-psiMin) / static_cast<vtkm::FloatDefault>(numPts-1);

    auto psi = psiMin;
    for (int i = 0; i < numPts; i++, psi += dPsi)
      psiVals.push_back(psi);
  }
  else if (args.find("--psiVals") != args.end())
  {
    vals = args["--psiVals"];
    for (const auto& v: vals)
      psiVals.push_back(std::atof(v.c_str()) * xgcParams.eq_x_psi);
  }

  if (args.find("--thetaRange") != args.end())
  {
    vals = args["--thetaRange"];
    vtkm::FloatDefault theta0 = 0, theta1 = 360, dTheta;
    int numThetas = -1;
    if (vals.size() == 1)
    {
      numThetas = std::atoi(vals[0].c_str());
      dTheta = (theta1-theta0) / static_cast<vtkm::FloatDefault>(numThetas);
    }
    else
    {
      theta0 = std::atof(vals[0].c_str());
      theta1 = std::atof(vals[1].c_str());
      numThetas = std::atoi(vals[2].c_str());
      dTheta = (theta1-theta0) / static_cast<vtkm::FloatDefault>(numThetas-1);
    }


    auto theta = theta0;
    for (int i = 0; i < numThetas; i++, theta += dTheta)
    {
      thetaVals.push_back(theta*degToRad);
      //std::cout<<"Theta_i= "<<theta<<std::endl;
    }
  }
  else if (args.find("--thetaVals") != args.end())
  {
    vals = args["--thetaVals"];
    for (const auto& v: vals)
      thetaVals.push_back(std::atof(v.c_str()) * degToRad);
  }

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxR;
  FindMaxR(ds, xgcParams, thetaVals, maxR);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  ds.GetField("coeff_2D").GetData().AsArrayHandle(arr);
  auto coeff = arr.ReadPortal();

  auto maxRPortal = maxR.ReadPortal();
  vtkm::Id ncoeff = xgcParams.eq_mr-1;
  vtkm::Id2 nrz(xgcParams.eq_mr, xgcParams.eq_mz);
  vtkm::Vec2f rzmin(xgcParams.eq_min_r, xgcParams.eq_min_z);
  vtkm::Vec2f drz((xgcParams.eq_max_r-xgcParams.eq_min_r)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mr-1),
                  (xgcParams.eq_max_z-xgcParams.eq_min_z)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mz-1));

  const vtkm::FloatDefault dR = 0.000001;

  vtkm::Id ID = 0;
  for (std::size_t i = 0; i < thetaVals.size(); i++)
  {
    auto theta = thetaVals[i];
    auto cost = vtkm::Cos(theta), sint = vtkm::Sin(theta);

    //std::cout<<"Theta_"<<i<<" = "<<theta<<" psiN: "<<psiNorm0<<" "<<psiNorm1<<std::endl;

    vtkm::FloatDefault r0 = 0;
    for (const auto& psiTarget : psiVals)
    {
      vtkm::FloatDefault r1 = r0;
      vtkm::FloatDefault psi = 0;

      //std::cout<<" PsiTarget= "<<psiTarget<<std::endl;
      vtkm::FloatDefault R, Z;

      r1 = r0+dR;
      bool found = false;
      while (!found && r1 < maxRPortal.Get(i))
      {
        R = xgcParams.eq_axis_r + r1*cost;
        Z = xgcParams.eq_axis_z + r1*sint;
        psi = InterpolatePsi({R,Z}, coeff, ncoeff, nrz, rzmin, drz, xgcParams);
        //std::cout<<"R: "<<r0<<" "<<r1<<" --> "<<vtkm::Vec2f(R,Z)<<" psi(r1)= "<<std::setprecision(15)<<psi<<std::endl;
        if (psi >= psiTarget)
        {
          found = true;
          break;
        }

        r0 += dR;
        r1 += dR;
      }
      if (!found)
        continue;
      auto diffPsi = psi - psiTarget;
      //std::cout<<"   Pt_"<<j<<" RZ= "<<R<<" "<<Z<<"  psiN= "<<psi /xgcParams.eq_x_psi;
      //std::cout<<std::setprecision(15)<<" diff= "<<diffPsi/xgcParams.eq_x_psi<<std::endl;
      //std::cout<<std::setprecision(15)<<"  psiTarget= "<<psiTarget<<"  psi= "<<psi<<std::endl;

      //Now, do a binary search to find psi between (r0, r1)
      int cnt = 0;
      while (diffPsi > 1e-10 && cnt < 100)
      {
        auto rMid = (r0+r1)/2.0;
        R = xgcParams.eq_axis_r + rMid*cost;
        Z = xgcParams.eq_axis_z + rMid*sint;
        psi = InterpolatePsi({R,Z}, coeff, ncoeff, nrz, rzmin, drz, xgcParams);

        if (psi < psiTarget) // mid is inside, range = (rmid, r1)
          r0 = rMid;
        else
          r1 = rMid;  //mid is outside, range = (r0, rmid)

        diffPsi = vtkm::Abs(psi - psiTarget);
        //std::cout<<std::setprecision(15)<<"        "<<cnt<<":  R=("<<r0<<" "<<r1<<")  psi= "<<std::setprecision(15)<<psi<<" "<<dPsi<<std::endl;
        cnt++;
      }

      //std::cout<<"   **** Pt_"<<j<<" RZ= "<<R<<" "<<Z<<"  psiN= "<<psi /xgcParams.eq_x_psi;
      //std::cout<<std::setprecision(15)<<" diff= "<<diffPsi/xgcParams.eq_x_psi<<std::endl;
      //std::cout<<"CNT= "<<cnt<<std::setprecision(15)<<" dPsi= "<<dPsi<<std::endl;

      //std::cout<<i<<" "<<j<<" ID: "<<ID<<" PT= "<<(psi/xgcParams.eq_x_psi)<<" "<<theta<<std::endl;
      //std::cout<<"  RZ= "<<R<<" "<<Z<<" cnt= "<<cnt<<" "<<"psiN= "<<(psi/xgcParams.eq_x_psi)<<" diffPsi= "<<diffPsi<<std::endl;

      vtkm::Vec3f pt_rpz(R, 0, Z);
      vtkm::Particle p({pt_rpz, ID++});
      seeds.push_back(p);
      seedsThetaPsi.push_back({theta, psiTarget});
    }
  }

  return seeds;
}

std::vector<vtkm::Particle>
GenerateSeeds(const vtkm::cont::DataSet& ds,
              XGCParameters& xgcParams,
              std::map<std::string, std::vector<std::string>>& args)
{
  std::vector<vtkm::Particle> seeds;
  std::vector<vtkm::Vec2f> seedsThetaPsi;

  if (args.find("--range") != args.end())
  {
    vtkm::Id numSeeds = std::stoi(args["--numSeeds"][0]);

    auto vals = args["--range"];
    vtkm::FloatDefault r0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault r1 = std::atof(vals[1].c_str());

    vtkm::FloatDefault dr = (r1-r0) / (float)(numSeeds-1);
    vtkm::FloatDefault r = r0;

    for (vtkm::Id id = 0; id < numSeeds; id++, r+=dr)
      seeds.push_back({{r, .1, 0}, id});
    //std::cout<<"SEEDS= "<<seeds<<std::endl;
  }
  else if (args.find("--psiRange") != args.end() ||
           args.find("--thetaRange") != args.end() ||
           args.find("--psiVals") != args.end() ||
           args.find("--thetaVals") != args.end())
  {
    seeds = GenerateThetaPsiSeeds(args, ds, xgcParams, seedsThetaPsi);
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({{r, t, z}, 0});
  }
  else if (args.find("--jong1") != args.end())
  {
    //first point in traces.v2
    //zone = 3132, incident nodes: 1521, 1612, 1613
    seeds = {{{3.029365, 6.183185, 0.020600}, 0}};
  }
  else if (args.find("--jong6") != args.end())
  {
    auto vals = args["--jong6"];
    std::cout<<"VALS= "<<vals<<std::endl;

    std::vector<vtkm::Vec3f> allSeeds;
    allSeeds = {
      {3.351443028564415449, 0.0, -0.451648806402756176}, //pt_0, ID= 10670   blue 3 islands
      {3.187329423521033878, 0.0, -0.665017624967372267}, //pt_1, ID= 12000
      {1.992020349316277139, 0.0, -0.126203396421661285}, //pt_2, ID= 13100
      {3.018666196722858963, 0.0, 0.073864239629065770}, //pt_3, ID= 0
      {3.176582679765305173, 0.0, -0.220557108925872658}, //pt_4, ID= 4000    stochastic problem area
      {2.129928300491922499, 0.0, -0.176508860570331577}, //pt_5, ID= 10153  semi stoch. xgc spread out, vtkm thin. good with 0.001
      {2.568671712782164995, 0.0, 0.050249128799423198}, //pt_6, ID= 100
      {2.934624677179501262, 0.0, 0.220686132855778427}, //pt_7, ID= 500
      {2.959288366738244580, 0.0, 0.448869653975662142}, //pt_8, ID= 5000
    };

    if (vals.size() == 0) //all the points.
    {
      std::size_t n = allSeeds.size();
      for (std::size_t i = 0; i < n; i++)
        seeds.push_back({allSeeds[i], static_cast<vtkm::Id>(i)});
    }
    else
    {
      int n = std::stoi(vals[0]);
      if (n >= (int)allSeeds.size())
        std::cout<<"Bad seed!!! #allseeds= "<<allSeeds.size()<<std::endl;

      seeds = {{allSeeds[n], 0}};
    }
  }
  else if (args.find("--jongrz") != args.end())
  {
    //Seed from the blue island.
    //seeds = {{3.351443,             0, -0.451649}};
    seeds = {{{3.351443028564415449, 0, -0.45164880640275617552}, 0}};
  }
  else if (args.find("--afterN") != args.end())
  {
    seeds = {
      {{3.321888620239255019, 0.0, 0.478933972623594384}, 0},
      {{2.568934684085571352, 0.0, 0.731290913908178353}, 1},
      {{3.493628202658771720, 0.0, 0.433951677589735296}, 2},
      {{2.862485694515508605, 0.0, 0.208737305948038576}, 3},
      {{2.905837753215041008, 0.0, -0.397811882628356817}, 4},
      {{3.391834939600261389, 0.0, -0.350011953142094434}, 5}
    };
  }
  else if (args.find("--parseAdios") != args.end())
  {
    auto dir = args["--dir"][0];
    auto vals = args["--parseAdios"];
    auto fname = dir + "/" + vals[0];
    int skip = 1;
    if (vals.size() > 1)
      skip = std::atoi(vals[1].c_str());

    std::cout<<"Reading seeds from "<<fname<<" skip= "<<skip<<std::endl;
    if (adios == nullptr)
    {
      if (args.find("--xml") == args.end())
        adios = new adios2::ADIOS;
      else
        adios = new adios2::ADIOS(args["--xml"][0]);
    }

    auto io = adios2::IO(adios->DeclareIO("seedsIO"));
    auto engine = io.Open(fname, adios2::Mode::Read);

    auto v = io.InquireVariable<float>("ephase");
    std::vector<float> tmp;
    engine.Get(v, tmp, adios2::Mode::Sync);
    engine.Close();

    int n = tmp.size();
    int id = 0;
    for (int i = 0; i < n; i += (9*skip))
    {
      auto R = static_cast<vtkm::FloatDefault>(tmp[i + 0]);
      auto Z = static_cast<vtkm::FloatDefault>(tmp[i + 1]);

      vtkm::Particle p(vtkm::Vec3f(R, 0.0, Z), id++);
      seeds.push_back(p);
    }
  }
  else if (args.find("--parse") != args.end())
  {
    //Generaate the seed list by running jongAll.py
    std::cout<<"READING: "<<args["--parse"][0]<<std::endl;
    std::ifstream seedFile;
    seedFile.open(args["--parse"][0]);
    std::string line;
    vtkm::Id pid = 0;
    while (std::getline(seedFile, line))
    {
      vtkm::FloatDefault r, z, p;
#ifdef VTKM_USE_DOUBLE_PRECISION
      sscanf(line.c_str(), "%lf, %lf, %lf", &r, &z, &p);
#else
      sscanf(line.c_str(), "%f, %f, %f", &r, &z, &p);
#endif
      seeds.push_back({{r,p,z}, pid});
      pid++;
    }
    std::cout<<"Seeds= "<<seeds<<std::endl;
  }

  if (args.find("--bothDir") != args.end())
  {
    seeds.reserve(2*seeds.size());
    seeds.insert(seeds.end(), seeds.begin(), seeds.end());

    seedsThetaPsi.reserve(2*seedsThetaPsi.size());
    seedsThetaPsi.insert(seedsThetaPsi.end(), seedsThetaPsi.begin(), seedsThetaPsi.end());
  }
  std::cout<<" ****** Num Seeds= "<<seeds.size()<<std::endl;

  if (args.find("--dumpSeeds") != args.end())
    DumpSeeds(seeds, seedsThetaPsi, xgcParams);

  return seeds;
}

void
StreamingPoincare(std::map<std::string, std::vector<std::string>>& args)
{
  XGCParameters xgcParams;

  std::string coeffFile = "xgc.poincare_init.bp";
  if (args.find("--CoeffFile") != args.end())
    coeffFile = args["--CoeffFile"][0];

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];
  auto ds = ReadStaticData(args, xgcParams, coeffFile);

  std::string fName = args["--streaming"][0];
  std::string outputFile = args["--output"][0];

  adiosStuff["data"] = new adiosS(adios, fName, "3d", adios2::Mode::Read);
  auto dataStuff = adiosStuff["data"];

  std::vector<vtkm::Particle> seeds;
  adios2::StepStatus status;
  int step = 0;
  while (true)
  {
    status = dataStuff->engine.BeginStep(adios2::StepMode::Read, 1000.0);
    if (status != adios2::StepStatus::OK)
    {
      std::cerr<<"Failed to read "<<fName<<" step "<<step<<". Exiting now"<<std::endl;
      break;
    }

    //Initialize the num planes variable.
    if (step == 0)
      xgcParams.numPlanes = dataStuff->GetVarInt("nphi");

    int timeStep = step;
    auto tsVar = dataStuff->io.InquireVariable<int>("tindex");
    if (tsVar)
      dataStuff->engine.Get(tsVar, &timeStep, adios2::Mode::Sync);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
    ReadTurbData(dataStuff, args, As_arr, dAs_arr);
    auto As = CalcAs(As_arr, xgcParams);
    auto dAs = Calc_dAs(dAs_arr, xgcParams);
    UpdateField(ds, "As_ff", As);
    UpdateField(ds, "dAs_ff_rzp", dAs);
    if (step == 0 && args.find("--dumpGrid") != args.end())
    {
      ds.PrintSummary(std::cout);
      vtkm::io::VTKDataSetWriter writer("grid.vtk");
      writer.WriteDataSet(ds);
    }

    std::cout<<step<<": Read data timeStep= "<<timeStep<<std::endl;
    dataStuff->engine.EndStep();

    if (step == 0) //create the seeds...
      seeds = GenerateSeeds(ds, xgcParams, args);

    auto seedsCopy = seeds;

    std::cout<<"Dump to "<<outputFile<<std::endl;
    Poincare(ds, xgcParams, seedsCopy, args, timeStep);

    step++;
  }
}

void InteractivePoincare(std::map<std::string, std::vector<std::string>>& args)
{
  FILE *fdread;
  FILE *fdwrite;
  int n;
  const char *fifoName = args["--mkfifoName"][0].c_str();
  float psiVal, thetaVal;
  int numPunc;
  float stepSize;
  int bytes_read = 0; 
  
  mkfifo(fifoName, 0666);
  XGCParameters xgcParams;
  auto ds = ReadDataSet(args, xgcParams);

  //printf("POINCARE INTERACTIVE READY\n");
  //fflush(stdout);

  

  fdwrite = fopen(fifoName, "w");
  fputs("ready", fdwrite);
  fclose(fdwrite);

  for(;;){
    if ((fdread = fopen(fifoName, "r")) == NULL){
      exit(EXIT_FAILURE);
    }


    //delete what might already be there so insert succeeds
    args.erase("--psiVals");
    args.erase("--thetaVals");
    args.erase("--numPunc");
    args.erase("--stepSize");

    std::vector<std::string> strPsiVals;
    std::vector<std::string> strThetaVals;
    std::vector<std::string> strNumPunc;
    std::vector<std::string> strStepSize;

    //read in psiVals and thetaVals
    while(fscanf(fdread, "%f%f", &psiVal, &thetaVal) == 2) {
      strPsiVals.push_back(std::to_string(psiVal));
      strThetaVals.push_back(std::to_string(thetaVal));
    }

    //read in semicolon followed by numPunc and stepSizes
    if(fscanf(fdread, ";%d%f", &numPunc, &stepSize) == 2){
      strNumPunc.push_back(std::to_string(numPunc));
      strStepSize.push_back(std::to_string(stepSize));
    } else {
      strNumPunc.push_back("1000");
      strStepSize.push_back("0.01");
    }        
    
    args.insert(std::pair<std::string, std::vector<std::string>>("--psiVals", strPsiVals));
    args.insert(std::pair<std::string, std::vector<std::string>>("--thetaVals", strThetaVals));
    args.insert(std::pair<std::string, std::vector<std::string>>("--numPunc", strNumPunc));
    args.insert(std::pair<std::string, std::vector<std::string>>("--stepSize", strStepSize));
    
    std::vector<vtkm::Particle> seeds;
    std::vector<vtkm::Vec2f> seedsThetaPsi;
    seeds = GenerateNormalizedFromThetaPsiSeedsPairs(args, ds, xgcParams, seedsThetaPsi);
    return;
    Poincare(ds, xgcParams, seeds, args, 0);


    fclose(fdread);
    if((fdwrite = fopen(fifoName, "w")) == NULL){ 
      exit(EXIT_FAILURE);
    }
    fputs("done", fdwrite);
    fclose(fdwrite);
  }
}

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank, numRanks;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::map<std::string, std::vector<std::string>> args;
  if (!ParseArgs(argc, argv, args))
    return 0;

  if (args.find("--gpuParams") != args.end())
  {
    oneDBlocks = std::atoi(args["--gpuParams"][0].c_str());
    threadsPerBlock = std::atoi(args["--gpuParams"][1].c_str());
  }
  std::cout<<"GPUParams: "<<oneDBlocks<<" "<<threadsPerBlock<<std::endl;

#ifdef VTKM_CUDA
  if (args.find("--gpu") != args.end())
    vtkm::cont::cuda::InitScheduleParameters(mySchedParams);
#endif

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  if (argc < 7)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures stepSize poincVar [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  if (args.find("--gpu") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});
  if (args.find("--kokkos") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagKokkos{});
  else if (args.find("--openmp") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP{});
  else if (args.find("--serial") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});
  else
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});


  if (args.find("--debug") != args.end())
  {
    /*
    XGCParameters xgcParams;
    auto ds = ReadDataSet(args, xgcParams);
    auto vals = args["--debug"];

    vtkm::FloatDefault R, Z;
    if (vals.size() == 1) //point index.
    {
      vtkm::Id idx = std::atoi(vals[0].c_str());
      vtkm::cont::ArrayHandle<vtkm::Vec3f> c;
      ds.GetCoordinateSystem().GetData().AsArrayHandle(c);
      R = c.ReadPortal().Get(idx)[0];
      Z = c.ReadPortal().Get(idx)[1];
    }
    else
    {
      R = std::atof(vals[0].c_str());
      Z = std::atof(vals[1].c_str());
    }
    */
    /*
    vtkm::Vec3f pt(R,Z,0);

    std::cout<<std::setprecision(15)<<"Debug pt= "<<pt<<std::endl;
    vtkm::Vec2f ptRZ(R,Z);
    vtkm::Vec3f B0;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
    ds.GetField("coeff_2D").GetData().AsArrayHandle(arr);
    auto coeff = arr.ReadPortal();
    vtkm::Id ncoeff = xgcParams.eq_mr-1;
    vtkm::Id2 nrz(xgcParams.eq_mr-1, xgcParams.eq_mz-1); //DRP
    vtkm::Vec2f rzmin(xgcParams.eq_min_r, xgcParams.eq_min_z);
    vtkm::Vec2f drz((xgcParams.eq_max_r-xgcParams.eq_min_r)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mr-1),
                    (xgcParams.eq_max_z-xgcParams.eq_min_z)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mz-1));

    auto psiHigh = InterpolatePsi(ptRZ,  coeff, ncoeff, nrz, rzmin, drz, xgcParams, B0);
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> coeff1D;
    ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff1D);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
    ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
    auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);
    vtkm::FloatDefault one_d_cub_dpsi_inv = 1.0/dPsi;
    auto Ipsi = I_interpol(psiHigh, 0, one_d_cub_dpsi_inv, ncoeff, coeff1D.ReadPortal());  //rgn=1
    B0[2] = Ipsi / R;


    vtkm::cont::CellLocatorTwoLevel locator2L;
    locator2L.SetCellSet(ds.GetCellSet());
    locator2L.SetCoordinates(ds.GetCoordinateSystem());
    locator2L.Update();
    auto psiLinear = EvalScalar(pt, locator2L, ds, "psi2D");
    auto B0Linear = EvalVec(pt, locator2L, ds, "B_RZP");

    std::cout<<std::setprecision(15);
    std::cout<<"PsiLinear = "<<psiLinear<<" norm: "<<(psiLinear/xgcParams.eq_x_psi)<<std::endl;
    std::cout<<"PsiHigh   = "<<psiHigh<<" norm: "<<(psiHigh/xgcParams.eq_x_psi)<<std::endl;
    std::cout<<"B0Linear  = "<<B0Linear<<std::endl;
    std::cout<<"B0High    = "<<B0<<std::endl;

    return 0;
    */
  }

  if (args.find("--streaming") != args.end())
  {
    StreamingPoincare(args);
  }
  else if (args.find("--saveVecField") != args.end())
  {
    XGCParameters xgcParams;
    auto ds = ReadDataSet(args, xgcParams);
    //ds.PrintSummary(std::cout);
    //vtkm::io::VTKDataSetWriter writer2("grid.vtk");
    //writer2.WriteDataSet(ds);
    SavePoincare(ds, xgcParams, args);
  }
  else if (args.find("--interactive") != args.end())
  {
    InteractivePoincare(args);
  }
  else
  {
    XGCParameters xgcParams;
    auto ds = ReadDataSet(args, xgcParams);
    //ds.PrintSummary(std::cout);
    //vtkm::io::VTKDataSetWriter writer2("grid.vtk");
    //writer2.WriteDataSet(ds);
    auto seeds = GenerateSeeds(ds, xgcParams, args);
    Poincare(ds, xgcParams, seeds, args);
  }

  //close all the adios stuff.
  for (auto& i : adiosStuff)
    delete i.second;

  return 0;
}

/*
//cyclone case

./examples/poincare/Poincare  --vField B --dir ../data/run_1920 --traces 0 --useHighOrder --turbulence 1 --jong6 --openmp  --output bumm --numPunc 303 --gpuParams 256 128 --stepSize 0.1 --test


//iter case
./examples/poincare/Poincare  --vField B --dir ../data/XGC_GB/test_GB_small_su455 --traces 0 --useHighOrder --turbulence 1 --psiRange .1 .9 4 2 --openmp  --output bumm --numPunc 303 --gpuParams 256 128 --stepSize 0.1

//streaming case
./examples/poincare/Poincare  --vField B --dir ../data/XGC_GB/test_GB_small_su455 --traces 0 --useHighOrder --turbulence 1 --psiRange .1 .9 4 2 --openmp  --output bumm --numPunc 500 --gpuParams 256 128 --stepSize 0.05 --useLinearB --streaming xgc.3d.panout.1.bp



./examples/poincare/Poincare  --vField B --dir ../data/XGC_GB/test_GB_small_su455 --traces 0 --useHighOrder --turbulence 1 --psiRange .1 .7 10 4 --openmp  --output bumm --numPunc 100 --gpuParams 256 128 --stepSize 0.05 --streaming xgc.3d.panout.5.bp --userLinearB



*/


/*

ITER:
0 point?      : 6.362548633399328, 0.5648666773221587
interior point:  6.234733194519893, 0.5953562419309504

i=6500 : 5.197881014778669 -1.9994832294333764

i=7300 : 4.866152226365543 -2.042357129245992

*/


/*
psi theta range:

./examples/poincare/Poincare --dir ../data/summitRunJong-v3 --openmp --output OUT --numPunc 500 --stepSize 0.01 --dumpSeeds --psiThetaRange .95 .99 50 260 262 265 267 270 272 275 277 280 285 290 300 310 320 330 340 350  240 245 250 45 90 135 180


Need to save out seed-psi0 for each particle.

 */
