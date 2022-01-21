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
#include <vtkm/Geometry.h>
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

#include <vtkm/cont/CellLocatorGeneral.h>

//#include <vtkm/filter/Gradient.h>

#include <vtkm/io/VTKDataSetWriter.h>

//#include <fides/DataSetReader.h>

#include <adios2.h>

#include <random>
#include <chrono>
#include <variant>
#include <filesystem>
#include <mpi.h>

/*
radius values: max range: 2.5 3.7
for interesting psi: 3.5  3.7
*/

/*
TODO:
Make sure that wrap around works.
Get the Bs in 3D working.

*/


adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

bool useTurb=true;
int numPlanes = -1;
int numNodes = -1;
int numTri = -1;
float XScale = 1;
vtkm::FloatDefault eq_axis_r = 2.8, eq_axis_z = 0.0, eq_x_psi = 0.0697345;
vtkm::FloatDefault eq_min_r = 1.60014, eq_max_r = 3.99986;
vtkm::FloatDefault eq_min_z = -1.19986, eq_max_z = 1.19986;
vtkm::FloatDefault psi_min = -1, psi_max = -1;
int eq_mr = -1, eq_mz = -1;
//  vtkm::FloatDefault eq_x_psi = 0.0697345, eq_x_r = 2.8, eq_x_z = -0.99988;


using Ray3f = vtkm::Ray<vtkm::FloatDefault, 3, true>;

#include "Poincare.h"

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  out<<"[";
  for (const auto& x : v)
    out<<x<<" ";
  out<<"]";
  return out;
}

std::vector<vtkm::Vec3f>
ConvertToThetaPsi(const std::vector<vtkm::Vec3f>& pts)
{
  std::vector<vtkm::Vec3f> output;
  for (const auto& p : pts)
  {
    auto R = p[0];
    auto Z = p[2];
    auto psi = vtkm::Sqrt(((R-eq_axis_r)*(R-eq_axis_r) + Z*Z));
    auto theta = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
    if (theta < 0) theta += vtkm::TwoPi();

    output.push_back({theta, psi, 0});
  }

  return output;
}

void
GetPlaneIdx(const vtkm::FloatDefault& phi,
            const vtkm::Id& nPlanes,
            vtkm::FloatDefault& phiN,
            vtkm::Id& plane0,
            vtkm::Id& plane1,
            vtkm::FloatDefault& phi0,
            vtkm::FloatDefault& phi1,
            vtkm::Id& numRevs,
            vtkm::FloatDefault& T)
{
  vtkm::FloatDefault dPhi = vtkm::TwoPi()/static_cast<double>(nPlanes);
  std::vector<vtkm::FloatDefault> PHIs;

  /*
  vtkm::FloatDefault p = 0;
  for (int i = 0; i < nPlanes+1; i++, p+= dPhi)
    PHIs.push_back(p);
  std::cout<<"PHIs= "<<PHIs<<std::endl;
  */

  numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
  //rem = std::fmod(vtkm::Abs(phi), vtkm::TwoPi());
  phiN = phi;
  if (phi < 0)
  {
    //rem = -rem;
    phiN += ((1+numRevs) * vtkm::TwoPi());
  }

  plane0 = vtkm::Floor(phiN / dPhi);
  plane1 = plane0 + 1;
  phi0 = static_cast<vtkm::FloatDefault>(plane0)*dPhi;
  phi1 = static_cast<vtkm::FloatDefault>(plane1)*dPhi;
  if (plane1 == nPlanes)
    plane1 = 0;
  T = (phiN-phi0) / (phi1-phi0);

//  std::cout<<phi<<"  : ("<<phi0<<" "<<phi1<<") :: ["<<plane0<<" "<<plane1<<"]"<<std::endl;

  return;

  /*
  plane0 = -1;
  plane1 = -1;
  vtkm::FloatDefault phi0 = -1.0f, phi1 = -1.0f;
  for (int i = 0; i < PHIs.size()-1; i++)
  {
    if (phiN > PHIs[i] && phiN <= PHIs[i+1])
    {
      plane0 = i; //(i % nPlanes);
      plane1 = i+1; //((i+1) % nPlanes);
      break;
    }
  }
  auto S = phiN / dPhi;
  auto S0 = vtkm::Floor(S);
  auto S1 = S0 + 1;
  if (S1 == nPlanes)
    S1 = 0;
  std::cout<<" **** S= "<<S<<" ("<<S0<<" "<<S1<<")"<<std::endl;

  phi0 = (plane0*dPhi);
  phi1 = (plane1*dPhi);
  if (plane1 == nPlanes)
    plane1 = 0;
  T = (phiN-phi0) / (phi1-phi0);
  */


  std::cout<<"phiN= "<<phiN<<" ";
  std::cout<<"plane0,plane1= "<<plane0<<" "<<plane1<<std::endl;
  std::cout<<"Phi0, Phi1= "<<phi0<<" "<<phi1<<std::endl;
  std::cout<<"      T= "<<T<<std::endl;

/*
  vtkm::Id idx = vtkm::Floor(phiN / dPhi);
  std::cout<<"    idx= "<<idx<<" :: "<<rem<<" "<<dPhi<<" phiN= "<<phiN<<std::endl;
  if (idx < 0)
    idx += nPlanes;

  if (idx == nPlanes-1)
  {
    plane0 = 0;
    plane1 = idx;
  }
  else
  {
    plane0 = idx;
    plane1 = plane0+1;
  }
*/
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

  return true;
}


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
ComputeV(vtkm::cont::DataSet& ds)
{
#if 0
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

  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f val = bPortal.Get(i);
    vtkm::FloatDefault R = cPortal.Get(i)[0];
    auto g = gPortal.Get(i);

    //From: https://www.therightgate.com/deriving-curl-in-cylindrical-and-spherical/
    //R: (1/R * dAz/dT  - dAT/dZ)
    //T: dAr/dZ - dAz/dr
    //Z: Az/R + dAt/dr - 1/R dAr/dT]
    vtkm::FloatDefault rv, tv, zv;
    rv = 1/R * g[2][1] - g[1][2];
    tv = g[0][2] - g[2][1];
    zv = As_bHat[i][1]/R + g[1][0] - 1/R*g[0][1];
    val[0] += rv;
    val[1] += tv;
    val[2] += zv;

    V[i] = val;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
#endif
}

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

  if ((vname == "As_phi_ff" || vname == "dAs_phi_ff")&& useTurb == false)
  {
    for (auto& x : tmp)
      x = 0.0;
  }

  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     tmp,
                                     vtkm::CopyFlag::On));
}

void
ReadScalar(adiosS* stuff,
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

  if (addExtra && add3D)
  {
    for (int i = 0; i < numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D", vtkm::cont::make_ArrayHandle(tmpPlane, vtkm::CopyFlag::On)));
  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname, vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadPsiInterp(adiosS* eqStuff,
              adiosS* interpStuff,
              vtkm::cont::DataSet& ds)
{
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mr"), &eq_mr, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mz"), &eq_mz, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_r"), &eq_axis_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_z"), &eq_axis_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_r"), &eq_min_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_r"), &eq_max_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_z"), &eq_min_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_z"), &eq_max_z, adios2::Mode::Sync);

  ReadOther(eqStuff, ds, "eq_I");
  ReadOther(eqStuff, ds, "eq_psi_grid");
  ReadOther(eqStuff, ds, "eq_psi_rz");
  ReadOther(interpStuff, ds, "coeff_1D", "one_d_cub_psi_acoef");

  std::vector<double> tmp2D;
//  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("one_d_cub_psi_acoef"),
//                          tmp1D, adios2::Mode::Sync);
  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("psi_bicub_acoef"),
                          tmp2D, adios2::Mode::Sync);

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
  int nr = eq_mr-1, nz = eq_mz-1;
  //int ni = 150, nj = 150;
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
          arr_coeff2D[idx] = coef2D[i][j][k][m];
          idx++;
        }

  ds.AddField(vtkm::cont::make_Field("coeff_2D",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
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


  vtkm::Vec2f origin2D(eq_min_r, eq_min_z);
  vtkm::Vec2f spacing2D((eq_max_r-eq_min_r)/150., (eq_max_z-eq_min_z)/150.);
  auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(151, 151),
                                                        origin2D, spacing2D);
  std::vector<vtkm::FloatDefault> c00;
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> cij(4);
  for (int k = 0; k < 4; k++) cij[k].resize(4);

  for (int i = 0; i < nz; i++)
    for (int j = 0; j < nr; j++)
    {
      c00.push_back(coeff_2D[i][j][0][0]);
      for (int k = 0; k < 4; k++)
        for (int m = 0; m < 4; m++)
          cij[k][m].push_back(coeff_2D[i][j][k][m]);
    }

  //ds2D.AddPointField("c00", c00);
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

  //ds2D.AddCellField("bum", cij[0][0]);
  //ds.PrintSummary(std::cout);
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b3d;
  ds.GetField("eq_psi_rz").GetData().AsArrayHandle(arr);
  ds2D.AddPointField("eq_psi_rz", arr);

//  ds.GetField("B_RZP").GetData().AsArrayHandle(b3d);
//  std::vector<vtkm::Vec3f> B2D(

  vtkm::io::VTKDataSetWriter writer("psiGrid.vtk");
  writer.WriteDataSet(ds2D);
  vtkm::io::VTKDataSetWriter writer2("grid.vtk");
  writer2.WriteDataSet(ds);


  /*
  auto v = stuff->io.InquireVariable<double>("eq_psi_rz");
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  std::cout<<"eq_psi_rz: "<<tmp.size()<<std::endl;
  for (int i = eq_mr; i < 2*eq_mr; i++)
    std::cout<<"eq_psi_rz["<<i<<"] = "<<tmp[i]<<std::endl;
  */


  //Let's evaluate the b field.
  vtkm::FloatDefault R = 2, Z = 0;

}

//-----------------------------------------------------------------------------
class NormalizeWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn input,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2);
  using InputDomain = _1;

  template <typename VectorType>
  VTKM_EXEC void operator()(const VectorType& input,
                            VectorType& output) const
  {
    output = vtkm::Normal(input);
  }
};

//Sets it as RZP
vtkm::cont::ArrayHandle<vtkm::Vec3f>
ComputeCurl(const vtkm::cont::DataSet& inDS)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  inDS.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);
  return curlBNorm;

#if 0
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  inDS.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B_RZP").GetData().AsArrayHandle(B);

  //Get the gradients.
  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("B_RZP");
  auto ds = gradient.Execute(inDS);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  ds.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    auto B0 = bPortal.Get(i);
    auto ptRPZ = cPortal.Get(i);
    auto R = ptRPZ[0];
    //auto Z = ptRPZ[2];

    auto inv_r = 1.0/R;
    auto Bmag = vtkm::Magnitude(B0);
    auto over_Bmag = 1.0/Bmag;
    auto over_Bmag2 = over_Bmag * over_Bmag;

    auto br = B0[0];
    auto bz = B0[1];
    auto bphi = B0[2];

    auto GRAD = gPortal.Get(i);

    auto dbrdr = GRAD[0][0];
    auto dbzdr = GRAD[0][1];
    auto dbpdr = GRAD[0][2];

    auto dbrdz = GRAD[1][0];
    auto dbzdz = GRAD[1][1];
    auto dbpdz = GRAD[1][2];

    auto dbrdp = GRAD[2][0];
    auto dbzdp = GRAD[2][1];
    //auto dbpdp = GRAD[2][2];

    //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    //dbdphi=0D0  ! no B perturbation
    auto dbdr = (br*dbrdr + bphi*dbpdr + bz*dbzdr) * over_Bmag;
    auto dbdz = (br*dbrdz + bphi*dbpdz + bz*dbzdz) * over_Bmag;
    //auto dbdphi = 0;


    vtkm::Vec3f curl_B;
    //R curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
    curl_B[0] =          dbzdp*inv_r - dbpdz;
    //Z curl_B(2)  = fld%bphi*inv_r + fld%dbpdr - fld%dbrdp*inv_r
    curl_B[1] =          bphi*inv_r +     dbpdr -     dbrdp*inv_r;
    //phi curl_B(3)  = fld%dbrdz - fld%dbzdr
    curl_B[2] =            dbrdz -     dbzdr;

    vtkm::Vec3f curl_nb;

    //R,Z,Phi
    //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2

    curl_nb[0] = curl_B[0]*over_Bmag + (bphi * dbdz) * over_Bmag2;
    curl_nb[1] = curl_B[1]*over_Bmag + (-bphi * dbdr) * over_Bmag2;
    curl_nb[2] = curl_B[2]*over_Bmag + (bz * dbdr - br * dbdz) * over_Bmag2;

    portal.Set(i, curl_nb);
  }

  return curlBNorm;


  //old way...
#if 0
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(B);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;

  //wrong way...
  /*
    vtkm::filter::Gradient gradient;
    gradient.SetComputePointGradient(true);
    gradient.SetComputeVorticity(true);
    gradient.SetActiveField("B2D_Norm");
    auto tmpDS = gradient.Execute(ds);
    tmpDS.GetField("Vorticity").GetData().AsArrayHandle(curlBNorm);
   */

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("BRZP_Norm");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    vtkm::Vec3f ptRZ = cPortal.Get(i);
    auto R = ptRZ[0];
    auto invR = 1.0f / R;
    auto BPhi = bPortal.Get(i)[2]; //This needs to be the BPhi at the position ?
    //std::cout<<__LINE__<<" fix bPhi!"<<std::endl;

    vtkm::Vec<vtkm::Vec3f, 3> grad = gPortal.Get(i);
    // mesh is in R,Z space.
    auto _dBr = grad[0];
    auto _dBz = grad[1];
    auto _dBphi = grad[2];

//    if (R > 3.0 && R < 3.1)
//      std::cout<<"B= "<<bPortal.Get(i)<<" GRAD: "<<_dBr<<" "<<_dBz<<" "<<_dBphi<<std::endl;

    //((dBr/dR, dBz/dR, dBphi/dR) (dBr/dz, dBz/dz, dBphi/dz) (dBr/dPhi, dBz/dPhi, dBphi/dPhi))
    auto dBr_dR     = grad[0][0];
    auto dBz_dR     = grad[0][1];
    auto dBphi_dR   = grad[0][2];

    auto dBr_dZ     = grad[1][0];
    auto dBz_dZ     = grad[1][1];
    auto dBphi_dZ   = grad[1][2];

    auto dBr_dPhi   = grad[2][0];
    auto dBz_dPhi   = grad[2][1];
    auto dBphi_dPhi = grad[2][2];

    vtkm::Vec3f curl;
    //curl_R = 1/R * dBz/dPhi - dBphi/dZ
    curl[0] = invR* dBz_dPhi - dBphi_dZ;

    //curl_Phi = dBr/dZ - dBz/dR
    curl[1] = dBr_dZ - dBz_dR;

    //curl_Z = BPhi / R + dPhi_dR - dBr_dPhi *1/R
    //curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
    //curl_B.Z = BPhi * invR + dBphi_dR - dBr_dphi*invR
    curl[2] = BPhi * invR  + dBphi_dR  - dBr_dPhi * invR;

    //std::cout<<"dBr/dPhi= "<<dBr_dPhi<<std::endl;
    //curl[2] = invR*(BPhi + R*dBphi_dR - dBr_dPhi);
    //curl(2) = 1/r*(B_phi + dB_phi_dR * R - dB_r/dPhi);

    //std::cout<<"********CURL: "<<curl<<" "<<vtkm::Magnitude(curl)<<std::endl;
    portal.Set(i, curl);
  }

  return curlBNorm;
#endif
#endif
}

//This uses the curl_b
vtkm::cont::ArrayHandle<vtkm::Vec3f>
ComputeCurl2(const vtkm::cont::DataSet& ds)
{
#if 0
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, Brzp;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(Brzp);

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("B_RZP_2D");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = Brzp.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    vtkm::Vec3f ptRZ = cPortal.Get(i);
    auto R = ptRZ[0];
    auto invR = 1.0f / R;
    auto B = bPortal.Get(i);
    auto BR = B[0];
    auto BZ = B[1];
    auto BPhi = B[2];
    auto BMag = vtkm::Magnitude(B);
    auto over_B = 1.0/BMag;
    auto over_B2 = over_B * over_B;

    vtkm::Vec<vtkm::Vec3f, 3> grad = gPortal.Get(i);
    if (R > 3.0 && R < 3.1)
      std::cout<<"B= "<<bPortal.Get(i)<<" GRAD: "<<grad[0]<<" "<<grad[1]<<" "<<grad[2]<<std::endl;

    //((dBr/dR, dBz/dR, dBphi/dR) (dBr/dz, dBz/dz, dBphi/dz) (dBr/dPhi, dBz/dPhi, dBphi/dPhi))
    auto dBr_dR   = grad[0][0];
    auto dBz_dR   = grad[0][1];
    auto dBphi_dR = grad[0][2];

    auto dBr_dZ   = grad[1][0];
    auto dBz_dZ   = grad[1][1];
    auto dBphi_dZ = grad[1][2];

    auto dBr_dPhi = grad[2][0];
    auto dBz_dPhi = grad[2][1];
    //auto dBphi_dPhi = grad[2][2];


    //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    auto db_dr = (BR*dBr_dR + BPhi*dBphi_dR + BZ*dBz_dR) * over_B;

    //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    auto db_dz = (BR*dBr_dZ + BPhi*dBphi_dZ + BZ*dBz_dZ) * over_B;

    //dbdphi=0D0  ! no B perturbation
    //auto db_dphi = 0;

    vtkm::Vec3f curl_B;

    //Setting curl in R,Phi,Z

    //curl_B.R = 1/R * dBz/dPhi - dBphi/dZ
    curl_B[0] = invR* dBz_dPhi - dBphi_dZ;

    //curl_B.Phi = dBr/dZ - dBz/dR
    curl_B[1] = dBr_dZ - dBz_dR;

    //curl_Z = BPhi / R + dPhi_dR - dBr_dPhi *1/R
    //curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
    //curl_B.Z = BPhi * invR + dBphi_dR - dBr_dphi*invR
    curl_B[2] = BPhi * invR  + dBphi_dR  - dBr_dPhi * invR;


    vtkm::Vec3f curl_nb;

    //curl_nb.R
    //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb[0] = curl_B[0]*over_B + (BPhi * db_dz) * over_B2;

    //curl_nb.Phi
    //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    curl_nb[1] = curl_B[1]*over_B + (BZ * db_dr - BR * db_dz) * over_B2;


    //curl_nb.Z
    //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb[2] = curl_B[1]*over_B + (BPhi * db_dr) * over_B2;

    portal.Set(i, curl_nb);
  }

  return curlBNorm;
#endif
}

void
ReadB(adiosS* stuff,
      vtkm::cont::DataSet& ds)
{
  std::string fileName = "/node_data[0]/values";

  auto Bvar = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(Bvar, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> b_rzp;
  for (int i = 0; i < numNodes; i++)
  {
    vtkm::Vec3f v(tmp[i*3+0], tmp[i*3+1], tmp[i*3+2]);
    b_rzp.push_back(v);
  }

  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP",vtkm::cont::make_ArrayHandle(b_rzp, vtkm::CopyFlag::On)));

  vtkm::cont::ArrayHandle<vtkm::Vec3f> bhat_rzp, b;
  vtkm::cont::Invoker invoker;

  ds.GetField("B_RZP").GetData().AsArrayHandle(b);
  invoker(NormalizeWorklet{}, b, bhat_rzp);
  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP_Norm",bhat_rzp));

  /*
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm_rzp;
  curlBNorm_rzp = ComputeCurl(ds);
  ds.AddField(vtkm::cont::make_FieldPoint("curl_nb_rzp", curlBNorm_rzp));
  */
}

void
ReadVec(adiosS* stuff,
        vtkm::cont::DataSet& ds,
        const std::string& vname,
        std::string fileName="",
        bool add3D=false,
        bool addExtra=false)
{
  if (fileName.empty())
    fileName = vname;

  bool isB = (vname == "B");

  auto var = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(var, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> vec, vecRZP, vec2d;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < numNodes; i++)
    {
      //B is saved as: R,Z,Phi
      int vidx = (isB ? i : (p*numNodes+i));
      //Swap to R,Phi,Z
      vtkm::Vec3f v(tmp[vidx*3+0], tmp[vidx*3+2], tmp[vidx*3+1]);
      //vtkm::Vec3f v(tmp[vidx*3+0], tmp[vidx*3+1], tmp[vidx*3+2]);
      vec.push_back(v);
      if (p == 0)
      {
        vec2d.push_back(v);
        vecRZP.push_back(vtkm::Vec3f(tmp[vidx*3+0], tmp[vidx*3+1], tmp[vidx*3+2]));
      }
    }
  }

  if (addExtra && add3D)
  {
    for (int i = 0; i < numNodes; i++)
      vec.push_back(vec[i]);
  }

  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname,vtkm::cont::make_ArrayHandle(vec, vtkm::CopyFlag::On)));
  else
  {
    ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D",vtkm::cont::make_ArrayHandle(vec2d, vtkm::CopyFlag::On)));
    ds.AddField(vtkm::cont::make_FieldPoint("B_RZP_2D",vtkm::cont::make_ArrayHandle(vecRZP, vtkm::CopyFlag::On)));

    if (isB)
    {
      vtkm::cont::ArrayHandle<vtkm::Vec3f> b, brzp, bhat, brzphat;
      vtkm::cont::Invoker invoker;
      ds.GetField("B2D").GetData().AsArrayHandle(b);
      invoker(NormalizeWorklet{}, b, bhat);
      ds.AddField(vtkm::cont::make_FieldPoint("B2D_Norm",bhat));

      ds.GetField("B_RZP_2D").GetData().AsArrayHandle(brzp);
      invoker(NormalizeWorklet{}, brzp, brzphat);
      ds.AddField(vtkm::cont::make_FieldPoint("BRZP_Norm", brzphat));


      vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
      curlBNorm = ComputeCurl(ds);


      ds.AddField(vtkm::cont::make_FieldPoint("Curl_B2D_Norm", curlBNorm));
      ds.AddField(vtkm::cont::make_FieldPoint("curl_nb_rzp", curlBNorm));

      //Add B3D.
      std::vector<vtkm::Vec3f> B3D;
      auto portal = b.ReadPortal();
      for (vtkm::Id p = 0; p < numPlanes; p++)
        for (vtkm::Id n = 0; n < numNodes; n++)
          B3D.push_back(portal.Get(n));
      ds.AddField(vtkm::cont::make_Field("B3D",
                                         vtkm::cont::Field::Association::WHOLE_MESH,
                                         B3D,
                                         vtkm::CopyFlag::On));
    }
  }
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff)
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
  for (int i = 0; i < numNodes; i++)
  {
    double R = rz[i*2 +0];
    double Z = rz[i*2 +1];

    vtkm::Vec3f ptRZ(R,Z,0);
    coords.push_back(ptRZ);
  }

  //cells
  for (int i = 0; i < numTri*3; i+=3)
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
ReadMesh3D(adiosS* meshStuff, bool addExtra)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  double dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);

  //points.
  double phi = 0;
  int NP = numPlanes;
  if (addExtra) NP++;
  for (int p = 0; p < NP; p++)
  {
    std::cout<<"ReadMesh3D: phi= "<<phi<<std::endl;
    for (int i = 0; i < numNodes; i++)
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
    for (int i = 0; i < numTri*3; i+=3)
    {
      int off = p*numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0+off);
      connIds.push_back(p1+off);
      connIds.push_back(p2+off);

      off = (p+1)*(numNodes);
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
}

class CalculateABhat : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn B,
                                FieldIn apars,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2, _3);

  VTKM_EXEC void operator()(const vtkm::Vec3f& bvec,
                            const vtkm::FloatDefault& apars,
                            vtkm::Vec3f& output) const
  {
    output = vtkm::Normal(bvec) * apars;
  }
};

class CalculateVecField : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn coords,
                                FieldIn B,
                                FieldIn apars,
                                FieldIn gradient,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);

  template <typename GradientType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const vtkm::Vec3f& bvec,
                            const vtkm::Vec3f& a_bHat,
                            const GradientType& grad,
                            vtkm::Vec3f& output) const
  {
    const auto& R = point[0];
    output = bvec;

    //From: https://www.therightgate.com/deriving-curl-in-cylindrical-and-spherical/
    //R: (1/R * dAz/dT  - dAT/dZ)
    //T: dAr/dZ - dAz/dr
    //Z: Az/R + dAt/dr - 1/R dAr/dT]
    auto rv = 1/R * grad[2][1] - grad[1][2];
    auto tv = grad[0][2] - grad[2][1];
    auto zv = a_bHat[1]/R + grad[1][0] - 1/R*grad[0][1];

    output[0] += rv;
    output[1] += tv;
    output[2] += zv;
  }
};

//-----------------------------------------------------------------------------
class FindCellWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn points,
                                WholeCellSetIn<> cellSet,
                                ExecObject locator,
                                FieldOut cellIds,
                                FieldOut pcoords);
//                                FieldOut ptIndices);
  using ExecutionSignature = void(_1, _2, _3, _4, _5); //, _6);
  using InputDomain = _1;

  template <typename LocatorType, typename CellSetType> //, typename PtIndexType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const CellSetType& vtkmNotUsed(cellSet),
                            const LocatorType& locator,
                            vtkm::Id& cellId,
                            vtkm::Vec3f& pcoords) const
  //                            PtIndexType& ptIndices ) const
  {
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
//      std::cout<<"Cell not found! "<<point<<std::endl;
//      std::cout<<"    ***** Try reducing the step size."<<std::endl;
      this->RaiseError(vtkm::ErrorString(status));
    }
    //ptIndices = cellSet.GetIndices(cellId);
    //auto x = cellSet.GetIndices(cellId);
    //ptIndices = x;
  }
};


std::vector<std::vector<vtkm::Vec3f>>
ConvertPuncturesToThetaPsi(const std::vector<std::vector<vtkm::Vec3f>>& puncturesRZ,
                           const vtkm::cont::DataSet& ds)
{
  std::cout<<"Convert Punctures: "<<puncturesRZ.size()<<std::endl;

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<std::vector<vtkm::Vec3f>> puncturesTP(puncturesRZ.size());
  for (int p = 0; p < (int)puncturesRZ.size(); p++)
  {
    auto pts = puncturesRZ[p];
    for (auto& pt : pts)
    {
      auto R = pt[0];
      auto Z = pt[2];
      pt = vtkm::Vec3f(R,Z,0);
    }
    std::cout<<"RZ_pts.size()= "<<pts.size()<<std::endl;
    auto RZpoints = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

    vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
    vtkm::cont::Invoker invoker;

    invoker(FindCellWorklet{}, RZpoints, ds.GetCellSet(), locator, cellIds, pcoords);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> psi, theta;
    ds.GetField("psi2D").GetData().AsArrayHandle(psi);
    ds.GetField("theta2D").GetData().AsArrayHandle(theta);

    auto cPortal = cellIds.ReadPortal();
    auto pPortal = pcoords.ReadPortal();
    auto psiPortal = psi.ReadPortal();
    auto thetaPortal = theta.ReadPortal();
    auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

    auto pRZpoints = RZpoints.ReadPortal();
    vtkm::Id numPts = RZpoints.GetNumberOfValues();
    puncturesTP[p].resize(numPts);
    for (vtkm::Id i = 0; i < numPts; i++)
    {
      auto R = pRZpoints.Get(i)[0];
      auto Z = pRZpoints.Get(i)[1];
      vtkm::Id vIds[3];
      vtkm::Id cid = cPortal.Get(i);
      cs.GetCellPointIds(cid, vIds);

      vtkm::VecVariable<vtkm::FloatDefault, 3> pVals, tVals;
      for (vtkm::Id j = 0; j < 3; j++)
      {
        pVals.Append(psiPortal.Get(vIds[j]));
        tVals.Append(thetaPortal.Get(vIds[j]));
      }

      vtkm::FloatDefault psiI, thetaI;
      vtkm::exec::CellInterpolate(pVals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), psiI);
      vtkm::exec::CellInterpolate(tVals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), thetaI);
      thetaI += vtkm::Pi();

      auto thetaVal = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
      if (thetaVal < 0)
        thetaVal += vtkm::TwoPi();

      //puncturesTP[p][i][0] = thetaI;
      puncturesTP[p][i][0] = thetaVal;
      puncturesTP[p][i][1] = psiI / eq_x_psi;
    }
  }

  return puncturesTP;
}

void
InterpVector(const vtkm::cont::DataSet& ds,
             const vtkm::cont::CellLocatorGeneral& locator,
             const std::vector<vtkm::Vec3f>& pts,
             const std::string& vName,
             std::vector<vtkm::Vec3f>& out,
             bool is3D = false)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> points;
  std::vector<vtkm::Id> offset(pts.size(), 0);
  if (is3D)
  {
    vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);

    auto pts2 = pts;
    for (std::size_t i = 0; i < pts.size(); i++)
    {
      //std::cout<<"Wrap around check: line= "<<__LINE__<<std::endl;
      vtkm::Id off = static_cast<vtkm::Id>(vtkm::Floor(pts2[i][1] / phiSpacing));
      pts2[i][1] = pts2[i][2];
      pts2[i][2] = 0;
      offset[i] = off * numNodes;
//      std::cout<<"******* Offset:  "<<pts2[i]<<" "<<off<<" --> "<<offset[i]<<std::endl;
    }
    points = vtkm::cont::make_ArrayHandle(pts2, vtkm::CopyFlag::On);
  }
  else
    points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  //Find the cell on the RZ plane.
  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
//    std::cout<<"CID: "<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
//    std::cout<<"  vals: "<<vPortal.Get(vIds[0]+offset[i])<<" "<<vPortal.Get(vIds[1]+offset[i])<<" "<<vPortal.Get(vIds[2]+offset[i])<<std::endl;
    for (vtkm::Id j = 0; j < 3; j++)
      vals.Append(vPortal.Get(vIds[j]+offset[i]));

    vtkm::Vec3f v;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
    out.push_back(v);
//    std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
//    std::cout<<" -----> "<<v<<std::endl;
  }
}

std::vector<std::vector<vtkm::Vec3f>>
EvalTensor(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName,
           const std::vector<int>& offset)
{
  for (std::size_t i = 0; i < pts.size(); i++)
    if (pts[i][2] != 0)
      std::cout<<"********************************************************** FIX ME: "<<__LINE__<<std::endl;

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::Off);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  //std::cout<<"EvalVector("<<vName<<"): "<<pts<<" offset= "<<offset<<std::endl;
  //Find the cell on the RZ plane.
  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f,3>> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  std::vector<std::vector<vtkm::Vec3f>> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);

    auto v0 = vPortal.Get(vIds[0]+offset[i]);
    auto v1 = vPortal.Get(vIds[1]+offset[i]);
    auto v2 = vPortal.Get(vIds[2]+offset[i]);

    std::vector<vtkm::Vec3f> grads(3);
    for (int j = 0; j < 3; j++)
    {
      vtkm::VecVariable<vtkm::Vec3f, 3> vals;
      vals.Append(v0[j]);
      vals.Append(v1[j]);
      vals.Append(v2[j]);

      vtkm::Vec3f v;
      vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
      grads[j] = v;
    }
    out.push_back(grads);
  }

  return out;
}

bool
PtLoc(const vtkm::cont::DataSet& ds,
      const vtkm::cont::CellLocatorGeneral& locator,
      const vtkm::Vec3f& pt,
      vtkm::Vec3f& param,
      vtkm::Vec<vtkm::Id, 3>& vIds)
{
  if (pt[2] != 0)
    std::cout<<"********************************************************** FIX ME: "<<__LINE__<<std::endl;

  std::vector<vtkm::Vec3f> pts = {pt};

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::Off);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  vtkm::Id ids[3];
  vtkm::Id cid = cPortal.Get(0);
  cs.GetCellPointIds(cid, ids);
  vIds[0] = ids[0];
  vIds[1] = ids[1];
  vIds[2] = ids[2];

  param = pPortal.Get(0);

  return true;
}

template <typename PortalType>
vtkm::Vec3f
EvalV(const PortalType& vPortal,
      const vtkm::Id& offset,
      const vtkm::Vec<vtkm::Id, 3>& vId,
      const vtkm::Vec3f& param)
{
  vtkm::VecVariable<vtkm::Vec3f, 3> vals;
  vals.Append(vPortal.Get(vId[0]+offset));
  vals.Append(vPortal.Get(vId[1]+offset));
  vals.Append(vPortal.Get(vId[2]+offset));

  vtkm::Vec3f v;
  vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), v);
  return v;
}

template <typename PortalType>
vtkm::FloatDefault
EvalS(const PortalType& sPortal,
      const vtkm::Id& offset,
      const vtkm::Vec<vtkm::Id, 3>& vId,
      const vtkm::Vec3f& param)
{
  vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
  vals.Append(sPortal.Get(vId[0]+offset));
  vals.Append(sPortal.Get(vId[1]+offset));
  vals.Append(sPortal.Get(vId[2]+offset));

  vtkm::FloatDefault s;
  vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), s);
  return s;
}

std::vector<vtkm::Vec3f>
EvalVector(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName,
           const std::vector<int>& offset)
{
  for (std::size_t i = 0; i < pts.size(); i++)
    if (pts[i][2] != 0)
      std::cout<<"********************************************************** FIX ME: "<<__LINE__<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::Vec3f> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);
  auto vPortal = V.ReadPortal();

  std::vector<vtkm::Vec3f> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Vec3f param;
    vtkm::Vec<vtkm::Id, 3> vIds;
    PtLoc(ds, locator, pts[i], param, vIds);

    vtkm::Vec3f v = EvalV(vPortal, (vtkm::Id)offset[i], vIds, param);
    out.push_back(v);
  }
  return out;
}

std::vector<vtkm::Vec3f>
EvalVector(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName)
{
  std::vector<int> offset(pts.size(), 0);
  return EvalVector(ds, locator, pts, vName, offset);
}


std::vector<vtkm::FloatDefault>
InterpScalar(const vtkm::cont::DataSet& ds,
             const vtkm::cont::CellLocatorGeneral& locator,
             const std::vector<vtkm::Vec3f>& pts,
             const std::string& vName,
             const std::vector<int>& offsets)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> points;
  points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);
  auto vPortal = V.ReadPortal();


  std::vector<vtkm::FloatDefault> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Vec3f param;
    vtkm::Vec<vtkm::Id, 3> vIds;
    PtLoc(ds, locator, pts[i], param, vIds);

    vtkm::FloatDefault v = EvalS(vPortal, (vtkm::Id)offsets[i], vIds, param);
    out.push_back(v);
  }
  return out;
}

//Input points are R,Phi,Z
void
Evaluate(const vtkm::cont::DataSet& ds,
         const vtkm::cont::CellLocatorGeneral& locator,
         const std::vector<vtkm::Vec3f>& pts,
         const std::string& vField,
         std::vector<vtkm::Vec3f>& output)
{
  bool isB = vField == "B";

  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_rzp, B_Norm_rzp, dAs_ff_rzp, Curl_nb_rzp;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff;
  ds.GetField("B_RZP").GetData().AsArrayHandle(B_rzp);
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("B_RZP_Norm").GetData().AsArrayHandle(B_Norm_rzp);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
  //ds.GetField("AsCurlBHat_RZP").GetData().AsArrayHandle(AsCurlBHat_rzp);
  ds.GetField("curl_nb_rzp").GetData().AsArrayHandle(Curl_nb_rzp);

  auto pB_rzp = B_rzp.ReadPortal();
  auto pAs_ff = As_ff.ReadPortal();
  auto pB_Norm_rzp = B_Norm_rzp.ReadPortal();
  auto pdAs_ff_rzp = dAs_ff_rzp.ReadPortal();
  //auto pAsCurlBHat_rzp = AsCurlBHat_rzp.ReadPortal();
  auto pCurl_nb_rzp = Curl_nb_rzp.ReadPortal();

  for (const auto& x : pts)
  {
    auto ptRPZ = x;
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    GetPlaneIdx(Phi, numPlanes, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);

    vtkm::Vec3f ptRZ(R, Z, 0);

    vtkm::Vec3f particlePos_param;
    vtkm::Vec<vtkm::Id,3> particlePos_vids;
    PtLoc(ds, locator, ptRZ, particlePos_param, particlePos_vids);
    auto B0_rzp = EvalV(pB_rzp, 0, particlePos_vids, particlePos_param);

    std::cout<<"Meow: "<<ptRPZ<<std::endl;
    std::cout<<"    B0= "<<B0_rzp<<std::endl;

    if (isB)
    {
      vtkm::Vec3f vec_rpz(B0_rzp[0], B0_rzp[2]/R, B0_rzp[1]);
      output.push_back(vec_rpz);
      continue;
    }

    vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

    vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray_rpz({R, phiN, Z}, B0_rpz);

    //Get point on mid plane.  Use the R,Z for this point for triangle finds.
    vtkm::FloatDefault RP_T;
    vtkm::Vec3f ptOnMidPlane_rpz;
    bool tmp;
    midPlane.Intersect(ray_rpz, RP_T, ptOnMidPlane_rpz, tmp);

    //Now, interpolate between Phi_i and Phi_i+1
    vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
    vtkm::FloatDefault T10 = 1.0f - T01;

    //Get vec at Phi0 and Phi1.
    //x_ff is in rzp
    vtkm::Vec3f x_ff_rzp(ptOnMidPlane_rpz[0], ptOnMidPlane_rpz[2], 0);
    std::vector<int> offsets(2);
    offsets[0] = planeIdx0*numNodes*2;
    offsets[1] = planeIdx0*numNodes*2 + numNodes;


    const vtkm::FloatDefault basis = 0.0f;
    auto B0_R = B0_rzp[0];
    auto B0_Z = B0_rzp[1];
    auto x_ff_R = x_ff_rzp[0];
    //auto x_ff_Z = x_ff_rzp[1];

    //gradPsi: pt on mid plane?  (question)
    //dPsi/dR = B0_Z * R
    //dPsi/dZ = -B0_R * R;
    vtkm::Vec3f gradPsi_rzp(B0_Z * x_ff_R, -B0_R * x_ff_R, 0);
    vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi_rzp);

    vtkm::Vec2f rvec(0,0), zvec(0,0);
    rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];
    rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi_rzp[1];
    zvec[0] =         (1.0-basis) * gammaPsi * (-gradPsi_rzp[1]);
    zvec[1] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];

    //Get the vectors in the ff coordinates.
    //auto dAs_ff_rzp = EvalVector(ds, locator, {x_ff_rzp, x_ff_rzp}, "dAs_ff_rzp", offsets);
    //auto dAs_ff0_rzp = dAs_ff_rzp[0];
    //auto dAs_ff1_rzp = dAs_ff_rzp[1];

    vtkm::Vec3f x_ff_param;
    vtkm::Vec<vtkm::Id,3> x_ff_vids;
    PtLoc(ds, locator, x_ff_rzp, x_ff_param, x_ff_vids);
    auto dAs_ff0_rzp = EvalV(pdAs_ff_rzp, offsets[0], x_ff_vids, x_ff_param);
    auto dAs_ff1_rzp = EvalV(pdAs_ff_rzp, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault wphi[2] = {T10, T01}; //{T01, T10};
    vtkm::Vec3f gradAs_rpz;

    //vec.r = wphi[0]*( rvec[0]*V.r[0] + zvec[0]*V.z[0]) +
    //        wphi[1]*( rvec[0]*V.r[1] + zvec[0]*v.z[1]);
    //vec.p = wphi[0]*V.phi[0] +
    //        whpi[1]*V.phi[1];
    //vec.z = wphi[0]*( rvec[1]*V.r[0] + zvec[1]*V.z[0]) +
    //        wphi[1]*( rvec[1]*V.r[1] + zvec[1]*V.Z[1]);
    gradAs_rpz[0] = wphi[0]*(rvec[0]*dAs_ff0_rzp[0] + zvec[0]*dAs_ff0_rzp[1]) +
                    wphi[1]*(rvec[0]*dAs_ff1_rzp[0] + zvec[0]*dAs_ff1_rzp[1]);
    gradAs_rpz[1] = wphi[0] * dAs_ff0_rzp[2] +
                    wphi[1] * dAs_ff1_rzp[2];
    gradAs_rpz[2] = wphi[0]*(rvec[1]*dAs_ff0_rzp[0] + zvec[1]*dAs_ff0_rzp[1]) +
                    wphi[1]*(rvec[1]*dAs_ff1_rzp[0] + zvec[1]*dAs_ff1_rzp[1]);

    vtkm::FloatDefault BMag = vtkm::Magnitude(B0_rzp);
    //project using bfield.
    //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
    gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*B0_rzp[0] - gradAs_rpz[2]*B0_rzp[1]) / B0_rzp[2];

    //deltaB = AsCurl(bhat) + gradAs x bhat.
    //std::vector<int> off = {planeIdx0*numNodes};
    //vtkm::Vec3f AsCurl_bhat_rzp = EvalVector(ds, locator, {x_ff_rzp}, "AsCurlBHat_RZP", off)[0];
    //auto AsCurl_bhat_rzp = EvalV(pAsCurlBHat_rzp, 0, x_ff_vids, x_ff_param);

    //vtkm::Vec3f curl_nb_rzp = EvalVector(ds, locator, {ptRZ}, "curl_nb_rzp")[0];
    std::cout<<"    pos_ids= "<<particlePos_vids<<std::endl;
    std::cout<<"    pos_parms= "<<particlePos_param<<std::endl;
    auto curl_nb_rzp = EvalV(pCurl_nb_rzp, 0, particlePos_vids, particlePos_param);

    //auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
    //vtkm::FloatDefault As_ff0 = As_ff[0];
    //vtkm::FloatDefault As_ff1 = As_ff[1];
    auto As_ff0 = EvalS(pAs_ff, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = EvalS(pAs_ff, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * curl_nb_rzp;
    std::cout<<"    As= "<<As<<std::endl;
    std::cout<<"    curl_nb_rzp= "<<curl_nb_rzp<<std::endl;
    std::cout<<"    curl_nb_rzp.size()= "<<pCurl_nb_rzp.GetNumberOfValues()<<std::endl;
    std::cout<<"    curl_nb_rzp.Get(3819)= "<<pCurl_nb_rzp.Get(3819)<<std::endl;
    std::cout<<"    curl_nb_rzp.Get(v0)= "<<pCurl_nb_rzp.Get(particlePos_vids[0])<<std::endl;
    std::cout<<"    curl_nb_rzp.Get(v1)= "<<pCurl_nb_rzp.Get(particlePos_vids[1])<<std::endl;
    std::cout<<"    curl_nb_rzp.Get(v2)= "<<pCurl_nb_rzp.Get(particlePos_vids[2])<<std::endl;
    std::cout<<"    AsCurl_bhat_rzp= "<<AsCurl_bhat_rzp<<std::endl;

    //vtkm::Vec3f bhat_rzp = EvalVector(ds, locator, {ptRZ}, "B_RZP_Norm")[0];
    auto bhat_rzp = EvalV(pB_Norm_rzp, 0, particlePos_vids, particlePos_param);
    std::cout<<"    bhat_rzp= "<<bhat_rzp<<std::endl;

    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    std::cout<<"    gradAs_rzp= "<<gradAs_rzp<<std::endl;
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);

    std::cout<<"    deltaB= "<<deltaB_rzp<<std::endl<<std::endl;

    deltaB_rzp[2] /= R;
    B0_rzp[2] /= R;

    vtkm::Vec3f vec_rzp = B0_rzp + deltaB_rzp;

    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    std::cout<<"    vec_rpz= "<<vec_rpz<<std::endl<<std::endl;
    output.push_back(vec_rpz);
  }


#if 0
  /*
  for (std::size_t i = 0; i < pts.size(); i++)
    output.push_back({0, -.1, 0});
  return;
  */

  /*
  vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);
  std::cout<<"\n\n********************************************************"<<std::endl;
  std::cout<<"phiSpacing= "<<phiSpacing<<std::endl;
  std::cout<<"NumPlanes= "<<numPlanes<<" spacing= "<<phiSpacing<<std::endl;
  std::cout<<"Plane, Phi"<<std::endl;
  for (int i = 0; i < numPlanes; i++)
    std::cout<<i<<", "<<((float)i * phiSpacing)<<"  deg= "<<(i*phiSpacing*57.92958)<<std::endl;
  std::cout<<"** pts= "<<pts<<std::endl;
  */

  bool isB = vField == "B3D";
  bool isV = vField == "V";
  bool isV2 = vField == "V2";
  bool isX = vField == "X";

  for (const auto& x : pts)
  {
    auto pt = x;
    vtkm::FloatDefault phi = pt[1];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    GetPlaneIdx(phi, numPlanes, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);

    //std::cout<<pt<<" : "<<phiN<<" Pln: "<<planeIdx0<<" "<<planeIdx1<<" Phi: "<<Phi0<<" "<<Phi1<<std::endl;
    if (isX || isB)
    {
      vtkm::Vec3f ptRZ(pt[0], pt[2], 0);
      std::vector<vtkm::Vec3f> P = {ptRZ};
      auto B0_rzp = EvalVector(ds, locator, P, "B_RZP_2D")[0];
      auto B = B0_rzp[0];
      B0_rzp[2] = B0_rzp[2] / pt[0];
      if (isB)
      {
        output.push_back(B0_rzp);
        continue;
      }

      //Calculate X
      vtkm::Vec3f rayPt(pt[0], phiN, pt[2]);
      Ray3f ray0(rayPt, -B), ray1(rayPt, B);

      //std::cout<<"Ray: "<<rayPt<<" "<<B<<std::endl;
      vtkm::Plane<> Plane0({0,Phi0,0}, {0,-1,0}), Plane1({0,Phi1,0}, {0,-1,0});

      vtkm::Vec3f ptOnPlane0, ptOnPlane1;
      vtkm::FloatDefault T0, T1;
      bool tmp;
      Plane0.Intersect(ray0, T0, ptOnPlane0, tmp);
      Plane1.Intersect(ray1, T1, ptOnPlane1, tmp);

      auto dist01 = vtkm::Magnitude(ptOnPlane1-ptOnPlane0);
      auto dist0i = vtkm::Magnitude(pt-ptOnPlane0) / dist01;
      //auto disti1 = vtkm::Magnitude(pt-ptOnPlane1) / dist01;

      //Eval X(p0_rz, p1_rz)
      std::vector<vtkm::Vec3f> P2 = { {ptOnPlane0[0], ptOnPlane0[2], 0}, {ptOnPlane1[0], ptOnPlane1[2], 0} };
      std::vector<int> offsets(2);
      offsets[0] = (int)(planeIdx0 * numNodes);
      offsets[1] = (int)(planeIdx1 * numNodes);

      auto X = EvalVector(ds, locator, P2, vField, offsets);
      //std::cout<<"     Eval X @ "<<P2<<" ---> "<<X<<std::endl;

      auto res = vtkm::Lerp(X[0], X[1], dist0i);
      res[1] /= pt[0];
      res = res+B;
      output.push_back(res);
    }
    //For V
    else
    {
      //Wrap around case....
      if (planeIdx0 == numPlanes-1 && planeIdx1 == 0)
      {
        //std::cout<<"Wrap around: "<<planeIdx0<<" phi= "<<phiN<<std::endl;
      }

      vtkm::Vec3f particleRZ(pt[0], pt[2], 0);
      //This is RZP
      vtkm::Vec3f B0_rzp = EvalVector(ds, locator, {particleRZ}, "B_RZP_2D")[0];
      vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

      vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
      vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
      Ray3f ray_rpz({pt[0], phiN, pt[2]}, B0_rpz);

      //Get point on mid plane.  Use the R,Z for this point for triangle finds.
      vtkm::FloatDefault RP_T;
      vtkm::Vec3f ptOnMidPlane_rpz;
      bool tmp;
      midPlane.Intersect(ray_rpz, RP_T, ptOnMidPlane_rpz, tmp);

      //Now, interpolate between Phi_i and Phi_i+1
      vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
      vtkm::FloatDefault T10 = 1.0f - T01;

      //Get vec at Phi0 and Phi1.
      //x_ff is in rzp
      vtkm::Vec3f x_ff_rzp(ptOnMidPlane_rpz[0], ptOnMidPlane_rpz[2], 0);
      std::vector<int> offsets(2);
      offsets[0] = planeIdx0*numNodes*2;
      offsets[1] = planeIdx0*numNodes*2 + numNodes;

      vtkm::Vec3f vec_phi0, vec_phi1;
      vtkm::Vec3f vec;
      if (isV)
      {
        //DRP: Are these vecs in RZP ??
        auto vecs = EvalVector(ds, locator, {x_ff_rzp, x_ff_rzp}, vField, offsets);
        vec_phi0 = vecs[0];
        vec_phi1 = vecs[1];
        vec = vec_phi0 * T01 + vec_phi1 * T10;
      }
      else if (isV2)
      {
        const vtkm::FloatDefault basis = 0.0f;

        auto R = pt[0];
        auto Phi = pt[1];
        auto x_ff_R = x_ff_rzp[0];
        auto x_ff_Z = x_ff_rzp[1];
        auto Z = pt[2];
        auto B0_R = B0_rzp[0];
        auto B0_Z = B0_rzp[1];

        //gradPsi: pt on mid plane?  (question)
        //dPsi/dR = B0_Z * R
        //dPsi/dZ = -B0_R * R;
        vtkm::Vec3f gradPsi_rzp(B0_Z * x_ff_R, -B0_R * x_ff_R, 0);
        vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi_rzp);

        vtkm::Vec2f rvec(0,0), zvec(0,0);
        rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];
        rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi_rzp[1];
        zvec[0] =         (1.0-basis) * gammaPsi * (-gradPsi_rzp[1]);
        zvec[1] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];

        //Get the vectors in the ff coordinates.
        auto dAs_ff_rpz = EvalVector(ds, locator, {x_ff_rzp, x_ff_rzp}, "dAs_ff_rpz", offsets);
        auto dAs_ff0_rpz = dAs_ff_rpz[0];
        auto dAs_ff1_rpz = dAs_ff_rpz[1];
        //std::cout<<x_ff<<":  dAs_ff= ::: "<<dAs_ff0<<" "<<dAs_ff1<<" mag: "<<vtkm::Magnitude(dAs_ff0)<<" "<<vtkm::Magnitude(dAs_ff1)<<std::endl;

        vtkm::FloatDefault wphi[2] = {T10, T01}; //{T01, T10};
        vtkm::Vec3f gradAs_rpz;
        //std::cout<<"T::: "<<phiN<<" ("<<Phi0<<" "<<Phi1<<") wphi: "<<wphi[0]<<" "<<wphi[1]<<std::endl;

        //vec.r = wphi[0]*( rvec[0]*V.r[0] + zvec[0]*V.z[0]) +
        //        wphi[1]*( rvec[0]*V.r[1] + zvec[0]*v.z[1]);
        //vec.p = wphi[0]*V.phi[0] +
        //        whpi[1]*V.phi[1];
        //vec.z = wphi[0]*( rvec[1]*V.r[0] + zvec[1]*V.z[0]) +
        //        wphi[1]*( rvec[1]*V.r[1] + zvec[1]*V.Z[1]);
        gradAs_rpz[0] = wphi[0]*(rvec[0]*dAs_ff0_rpz[0] + zvec[0]*dAs_ff0_rpz[2]) +
                        wphi[1]*(rvec[0]*dAs_ff1_rpz[0] + zvec[0]*dAs_ff1_rpz[2]);
        gradAs_rpz[1] = wphi[0] * dAs_ff0_rpz[1] +
                        wphi[1] * dAs_ff1_rpz[1];
        gradAs_rpz[2] = wphi[0]*(rvec[1]*dAs_ff0_rpz[0] + zvec[1]*dAs_ff0_rpz[2]) +
                        wphi[1]*(rvec[1]*dAs_ff1_rpz[0] + zvec[1]*dAs_ff1_rpz[2]);

        vtkm::FloatDefault BMag = vtkm::Magnitude(B0_rzp);
        //project using bfield.
        //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
        gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*B0_rzp[0] - gradAs_rpz[2]*B0_rzp[1]) / B0_rzp[2];

        //deltaB = AsCurl(bhat) + gradAs x bhat.
        std::vector<int> off = {planeIdx0*numNodes};
        vtkm::Vec3f AsCurl_bhat = EvalVector(ds, locator, {x_ff_rzp}, "AsCurlBHat", off)[0];

        vtkm::Vec3f curl_bhat_rzp = EvalVector(ds, locator, {particleRZ}, "curl_nb_rzp")[0];
        auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
        vtkm::FloatDefault As_ff0 = As_ff[0];
        vtkm::FloatDefault As_ff1 = As_ff[1];

        vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
        AsCurl_bhat_rzp = As * curl_bhat_rzp;

        //std::cout<<"As*Curl(Bhat):: "<<AsCurl_bhat<<" "<<vtkm::Magnitude(AsCurl_bhat)<<std::endl;
        vtkm::Vec3f bhat_rzp = EvalVector(ds, locator, {particleRZ}, "BRZP_Norm")[0];
        //std::cout<<"bhat_rzp= "<<bhat_rzp<<" "<<vtkm::Magnitude(bhat_rzp)<<std::endl;

        vtkm::Vec3f deltaB = AsCurl_bhat_rzp + vtkm::Cross(gradAs, bhat_rzp);

        vec = deltaB;
        std::cout<<"Meow: deltaB= "<<deltaB_rzp<<std::endl;

        //std::cout<<"deltaB: "<<deltaB<<" :: "<<vtkm::Magnitude(deltaB)<<std::endl;

        std::cout<<"Meow: pt= "<<pt<<" "<<particleRZ<<" x_rzp= "<<x_ff_rzp<<" "<<phiN<<" "<<Phi0<<" "<<Phi1<<" wphi=("<<wphi[0]<<", "<<wphi[1]<<")"<<std::endl;
        std::cout<<"           B0_rzp= "<<B0_rzp<<"  /R= "<<vtkm::Vec3f(B0_rzp[0], B0_rzp[1]/particleRZ[0], B0_rzp[2])<<std::endl;
        std::cout<<"    curl_bhat_rzp= "<<curl_bhat_rzp<<std::endl;
        std::cout<<"           As= "<<As<<std::endl;
        std::cout<<"Bsa (As*curl_nb)= "<<AsCurl_bhat<<std::endl;
        std::cout<<"          dAs_rpz= "<<dAs_ff0_rpz<<" "<<dAs_ff1_rpz<<" --> "<<gradAs<<std::endl;
        std::cout<<"\n\n"<<std::endl;
      }

      B0_rzp[1] = B0_rzp[1] / particleRZ[0];
      vec[1] = vec[1] / particleRZ[0];
      auto result = vec + B0_rzp;
      output.push_back(result);
    }
  }
#endif
}

static int rk4_counter = 0;

std::vector<vtkm::Vec3f>
RK4(const vtkm::cont::DataSet& ds,
    const vtkm::cont::CellLocatorGeneral& locator,
    const std::vector<vtkm::Vec3f>& pts,
    const std::string& vField,
    const std::vector<bool>& /*pointMask*/,
    vtkm::FloatDefault h)
{
  std::vector<vtkm::Vec3f> k1;
  Evaluate(ds, locator, pts, vField, k1);

  std::vector<vtkm::Vec3f> tmp(pts.size());
  vtkm::FloatDefault h_2 = h/2.0;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k1[i]*h_2;
  std::vector<vtkm::Vec3f> k2;
  Evaluate(ds, locator, tmp, vField, k2);

  std::vector<vtkm::Vec3f> k3;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k2[i]*h_2;
  Evaluate(ds, locator, tmp, vField, k3);

  std::vector<vtkm::Vec3f> k4;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k3[i]*h;
  Evaluate(ds, locator, tmp, vField, k4);

  vtkm::FloatDefault h_6 = h/6.0;
  std::vector<vtkm::Vec3f> newPts(pts.size());
  for (std::size_t i = 0; i < pts.size(); i++)
  {
    newPts[i] = pts[i] + h_6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    //auto rkv = h_6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    //std::cout<<"*******************  RK4_"<<rk4_counter<<" : "<<pts[i]<<" === "<<rkv<<" norm: "<<vtkm::Normal(rkv)<<" =======> "<<newPts[i]<<std::endl;

    /*
    //Wrap around.
    if (newPts[i][1] < 0)
      newPts[i][1] += vtkm::TwoPi();
    else if (newPts[i][1] > vtkm::TwoPi())
      newPts[i][1] -= vtkm::TwoPi();
    */
  }

  rk4_counter++;
  return newPts;
}

std::vector<std::vector<vtkm::Vec3f>>
Poincare(const vtkm::cont::DataSet& ds,
         std::vector<vtkm::Vec3f>& pts,
         const std::string& vField,
         vtkm::FloatDefault h,
         int numPunc,
         bool useWorklet,
         bool useBOnly,
         bool useHighOrder,
         std::vector<std::vector<vtkm::Vec3f>>* traces=nullptr)

{
  if (useWorklet)
  {
    vtkm::cont::CellLocatorGeneral locator;
    locator.SetCellSet(ds.GetCellSet());
    locator.SetCoordinates(ds.GetCoordinateSystem());
    locator.Update();

    PoincareWorklet worklet(numPunc, 0.0f, h, (traces!=nullptr));
    worklet.UseBOnly = useBOnly;
    worklet.UseHighOrder = useHighOrder;

    vtkm::cont::Invoker invoker;
    std::vector<vtkm::Particle> s;
    for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
      s.push_back(vtkm::Particle(pts[i], i));
    auto seeds = vtkm::cont::make_ArrayHandle(s, vtkm::CopyFlag::On);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff, coeff_1D, coeff_2D;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> B_rzp, B_Norm_rzp, dAs_ff_rzp;//, AsCurlBHat_rzp, curl_nb_rzp;
    ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
    ds.GetField("B_RZP").GetData().AsArrayHandle(B_rzp);
    ds.GetField("B_RZP_Norm").GetData().AsArrayHandle(B_Norm_rzp);
    ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
    //ds.GetField("AsCurlBHat_RZP").GetData().AsArrayHandle(AsCurlBHat_rzp);
    //ds.GetField("curl_nb_rzp").GetData().AsArrayHandle(curl_nb_rzp);
    ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff_1D);
    ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);

    vtkm::cont::ArrayHandle<vtkm::Vec3f> tracesArr;
    std::vector<vtkm::Vec3f> o, t;
    o.resize(numPunc*pts.size(), {-100, -100, -100});
    auto output = vtkm::cont::make_ArrayHandle(o, vtkm::CopyFlag::On);
    if (traces != nullptr)
    {
      t.resize(pts.size()*worklet.MaxIter, {-100, -100, -100});
      std::cout<<"TRACES: "<<t.size()<<std::endl;
      tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
    }

    std::vector<vtkm::Id> puncID(o.size(), -1);
    vtkm::cont::ArrayHandle<vtkm::Id> punctureID = vtkm::cont::make_ArrayHandle(puncID, vtkm::CopyFlag::On);

    auto start = std::chrono::steady_clock::now();
    invoker(worklet, seeds, locator, ds.GetCellSet(),
            B_rzp, B_Norm_rzp, /*curl_nb_rzp,*/ As_ff, dAs_ff_rzp,
            coeff_1D, coeff_2D,
            tracesArr, output, punctureID);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout<<"PoincareTime= "<<elapsed_seconds.count()<<std::endl;
    std::cout<<"output.size()= "<<output.GetNumberOfValues()<<std::endl;
    std::cout<<"punctureID.size()= "<<punctureID.GetNumberOfValues()<<std::endl;
    //vtkm::cont::printSummary_ArrayHandle(output, std::cout);
    //vtkm::cont::printSummary_ArrayHandle(punctureID, std::cout);

    std::vector<std::vector<vtkm::Vec3f>> res;
    if (traces)
    {
      auto portal = tracesArr.ReadPortal();
      vtkm::Id n = portal.GetNumberOfValues();
      for (vtkm::Id i = 0; i < n; i++)
      {
        auto v = portal.Get(i);
        if (v[2] > -1)
          (*traces)[0].push_back(v);
      }
    }

    auto portal = output.ReadPortal();
    vtkm::Id n = portal.GetNumberOfValues();
    res.resize(pts.size());
    for (vtkm::Id i = 0; i < n; i++)
    {
      vtkm::Id id = punctureID.ReadPortal().Get(i);
      if (id >= 0)
      {
        auto p = portal.Get(i);
        res[id].push_back(p);
      }
    }
    return res;
  }

//  const vtkm::FloatDefault planeVal = 2.0f;
//  const vtkm::FloatDefault planeVal = vtkm::Pi();
  const vtkm::FloatDefault planeVal = 0; //vtkm::TwoPi();

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<std::vector<vtkm::Vec3f>> punctures(pts.size());
  std::vector<int> puncCount(pts.size(), 0);
  std::vector<bool> pointMask(pts.size(), true);

  if (traces)
    for (int i = 0; i < (int)pts.size(); i++)
      (*traces)[i].push_back(pts[i]);
/*
  auto thetaPsi = ConvertToThetaPsi(pts);
  std::cout<<"Poincare: "<<pts<<"  theta,psi= "<<thetaPsi<<std::endl;

  vtkm::Vec3f p(pts[0][0], pts[0][2], 0);
  //B0: R,phi,Z
  auto B0 = (EvalVector(ds, locator, {p}, "B_RZP"))[0];
  std::cout<<"B0= "<<B0<<std::endl;
  auto B0_pol = vtkm::Sqrt(B0[0]*B0[0] + B0[2]*B0[2]);
  std::cout<<"B0_pol = "<<B0_pol<<std::endl;

  std::cout<<"B0_phi/R = qB0_pol/r_minor"<<std::endl;
  std::cout<<B0[1]<<" / "<<pts[0][0]<<" = q "<<B0_pol<<" / "<<thetaPsi[0][1]<<std::endl;
  auto x = B0[1] / pts[0][0];
  auto y = B0_pol / thetaPsi[0][1];
  std::cout<<"----> "<<x<<" "<<y<<std::endl;
  std::cout<<std::endl<<std::endl<<std::endl;
*/


  int maxIter = numPunc*1000000;
  //maxIter = 15000;
  for (int i = 0; i < maxIter; i++)
  {
    auto newPts = RK4(ds, locator, pts, vField, pointMask, h);

    for (std::size_t j = 0; j < pts.size(); j++)
    {
      //We puncture the plane if the points are on opposite sides of the plane.
      int nRevs0 = vtkm::Floor(vtkm::Abs(pts[j][1] / vtkm::TwoPi()));
      int nRevs1 = vtkm::Floor(vtkm::Abs(newPts[j][1] / vtkm::TwoPi()));
      //std::cout<<" PCHECK: "<<pts[j][1]<<" "<<newPts[j][1]<<" planeVal= "<<planeVal<<" nREVS0= "<<nRevs0<<" "<<nRevs1<<std::endl;

      if (i % 500 == 0) std::cout<<"Poinc iter: "<<i<<" "<<pts[j]<<" nRevs: "<<nRevs0<<" "<<nRevs1<<std::endl;
      //if (newPts[j][1] < 0) pointMask[j] = false;

      if (nRevs1 > nRevs0)
      {
        punctures[j].push_back(newPts[j]);
        puncCount[j]++;
        std::cout<<"PUNC: "<<pts[j]<<" --> "<<newPts[j]<<"  planeVal= "<<planeVal<<" punc= "<<puncCount[j]<<std::endl;
        std::cout<<j<<":    "<<newPts[j]<<" "<<puncCount[j]<<std::endl;
        if (puncCount[j] == numPunc)
          pointMask[j] = false;
      }

      if (traces)
        (*traces)[j].push_back(newPts[j]);
    }

    //All points are done.
    if (std::accumulate(pointMask.begin(), pointMask.end(), 0) == 0)
      break;

    pts = std::move(newPts);
  }

  return punctures;
}

void TestLocator(const vtkm::cont::DataSet& ds)
{
  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<vtkm::Vec3f> pts = {{2, 0, 0}, {3,.1,0}, {2.4, -.5, 0}};

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id, 3>> vertIds;
  vtkm::cont::Invoker invoker;

  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);
  std::cout<<"Points:  "; vtkm::cont::printSummary_ArrayHandle(points, std::cout, true);
  std::cout<<"CellIds: "; vtkm::cont::printSummary_ArrayHandle(cellIds, std::cout, true);
  std::cout<<"Param:   "; vtkm::cont::printSummary_ArrayHandle(pcoords, std::cout, true);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> B;
  ds.GetField("B2D").GetData().AsArrayHandle(B);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    std::cout<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
    for (vtkm::Id j = 0; j < 3; j++)
      vals.Append(bPortal.Get(vIds[0]));

    vtkm::Vec3f interp;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), interp);
    std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
    std::cout<<" -----> "<<interp<<std::endl;
  }
}

void
CalcV(vtkm::cont::DataSet& ds)
{
  std::cout<<__FILE__<<" "<<__LINE__<<" fix me. only works for 1 index. (see below)"<<std::endl;
  //DeltaB = As * curl(Bhat) + grad_As x bhat
  //Advect: DeltaB + B0

  //ds.PrintSummary(std::cout);
  vtkm::cont::ArrayHandle<vtkm::Vec3f> AsCurlBhat, gradAs, B0;
  ds.GetField("AsCurlBHat").GetData().AsArrayHandle(AsCurlBhat);
  ds.GetField("gradAs").GetData().AsArrayHandle(gradAs);
  ds.GetField("B2D").GetData().AsArrayHandle(B0);

  std::vector<vtkm::Vec3f> vField((numPlanes * numNodes * 2), vtkm::Vec3f(1,0,0));
  auto B = vField;
  auto Bn = vField;
  auto acb = vField;
  auto gas = vField;

  auto acbPortal = AsCurlBhat.ReadPortal();
  auto gasPortal = gradAs.ReadPortal();
  auto b0Portal = B0.ReadPortal();

  vtkm::Id idx = 0, idx2 = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        auto b = b0Portal.Get(n);
        auto bn = vtkm::Normal(b);

        auto v1 = acbPortal.Get(idx);
        auto v2 = gasPortal.Get(idx);
        auto val = v1 + vtkm::Cross(v2, bn);
        vField[idx] = val;

#if 0
        if (i == 0) //<-------------------- Only do this for ONE plane.
        {
          auto v1 = acbPortal.Get(idx);
          auto v2 = gasPortal.Get(idx);
          //auto cross = vtkm::Cross(v2, bn);
          auto val  = v1 + vtkm::Cross(v2, bn); // + b;

          //R,Phi,Z.
          //vtkm::FloatDefault dPhi = (float)(p+1) * -0.01;
          //val = {0, dPhi, 0};
          //Need to store as R,Z,Phi
          //val = b;
          vField[idx2] = val;
          B[idx2] = b;
          Bn[idx2] = bn;
          acb[idx2] = v1;
          gas[idx2] = v2;

//          if (vtkm::Magnitude(v2) > 100)
//            std::cout<<"************************ "<<val<<" :: "<<v1<<" "<<v2<<" "<<bn<<" "<<b<<std::endl;

          idx2++;
        }
#endif
        idx++;
      }
    }
  }
  std::cout<<"Calc V: "<<idx2<<" "<<vField.size()<<std::endl;

#if 0
  //Now duplicate plane 0 to plane N.
  for (int n = 0; n < numNodes; n++)
  {
    vField.push_back(vField[n]);
    B.push_back(B[n]);
    Bn.push_back(Bn[n]);
    acb.push_back(acb[n]);
    gas.push_back(gas[n]);
  }


  ds3D.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(vField, vtkm::CopyFlag::On)));

  if (0)
  {
    ds3D.AddField(vtkm::cont::make_FieldPoint("B", vtkm::cont::make_ArrayHandle(B, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("Bn", vtkm::cont::make_ArrayHandle(Bn, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("ACB", vtkm::cont::make_ArrayHandle(acb, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("GAS", vtkm::cont::make_ArrayHandle(gas, vtkm::CopyFlag::On)));

    std::cout<<"Dumping file"<<std::endl;
    vtkm::io::VTKDataSetWriter writer("dsWithV.vtk");
    writer.WriteDataSet(ds3D);
  }
#endif

  ds.AddField(vtkm::cont::make_Field("V",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     vField, vtkm::CopyFlag::On));
}

void
CalcX(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B, Bs, X;

  bool isB2D = ds.HasField("B2D");

  if (isB2D)
    ds.GetField("B2D").GetData().AsArrayHandle(B);
  else
    ds.GetField("B").GetData().AsArrayHandle(B);

  ds.GetField("Bs").GetData().AsArrayHandle(Bs);

  X.Allocate(Bs.GetNumberOfValues());
  auto BP = B.ReadPortal();
  auto BsP = Bs.ReadPortal();
  auto XP = X.WritePortal();

  vtkm::Id idx = 0;
  for (vtkm::Id p = 0; p < numPlanes; p++)
    for (vtkm::Id i = 0; i < numNodes; i++, idx++)
    {
      vtkm::Id bidx = (isB2D ? i : idx);
      XP.Set(idx, BsP.Get(idx)-BP.Get(bidx));
    }

  ds.AddField(vtkm::cont::make_FieldPoint("X", X));
}

void
CalcAs(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);
  vtkm::Id numAs = As.GetNumberOfValues();
  auto asPortal = As.ReadPortal();

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_arr.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
      As_arr[p][i].resize(numNodes);
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
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
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arrAs,
                                     vtkm::CopyFlag::On));
}

void
Calc_dAs(vtkm::cont::DataSet& ds)
{
  //Assumed that dAs_phi_ff is R,Z,Phi
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> dAsAH;
  ds.GetField("dAs_phi_ff").GetData().AsArrayHandle(dAsAH);
  int nVals = dAsAH.GetNumberOfValues()/3;

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> dAs_ff;
  dAs_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    dAs_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      dAs_ff[p][i].resize(numNodes);
  }


  int idx = 0;
  auto dAsPortal = dAsAH.ReadPortal();

  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
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
  VEC_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    VEC_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC_ff[p][i].resize(numNodes);
  }

  idx = 0;
  std::vector<vtkm::Vec3f> arr_ff(nVals);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f val = VEC_ff[p][i][n];
        val = dAs_ff[p][i][n];
        arr_ff[idx] = val;
        idx++;
      }
    }
  }

  //do some testing..
  for (int p = 0; p < numPlanes; p++)
  {
    int off0 = p*numNodes*2;
    int off1 = p*numNodes*2 + numNodes;
    for (int n = 0; n < numNodes; n++)
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

  ds.AddField(vtkm::cont::make_Field("dAs_ff_rzp",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_ff,
                                     vtkm::CopyFlag::On));
}

void
CalcAsCurlBHat(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBhat;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);
  vtkm::Id numAs = As.GetNumberOfValues();
  auto asPortal = As.ReadPortal();


  ds.GetField("curl_nb_rzp").GetData().AsArrayHandle(curlBhat);
  auto cbPortal = curlBhat.ReadPortal();


  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> As_curlBhat;
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_curlBhat.resize(numPlanes);
  As_arr.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    As_curlBhat[p].resize(2);
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
    {
      As_curlBhat[p][i].resize(numNodes);
      As_arr[p][i].resize(numNodes);
    }
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
    {
      vtkm::Vec3f cBhat = cbPortal.Get(n);
      for (int i = 0; i < 2; i++)
      {
        auto as = asPortal.Get(idx);
        auto val = as * cBhat;
        //std::cout<<"As: "<<as<<" cBhat= "<<cBhat<<" :: "<<vtkm::Magnitude(cBhat)<<"  val= "<<val<<std::endl;
        As_curlBhat[p][i][n] = val;
        As_arr[p][i][n] = as;
        idx++;
      }
    }
  }

  //flatten to 1d index.
  idx = 0;
  std::vector<vtkm::Vec3f> arr(numAs);
  std::vector<vtkm::FloatDefault> arrAs(numAs);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        vtkm::Vec3f val = As_curlBhat[p][i][n];
        arr[idx] = val;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("AsCurlBHat_RZP",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr,
                                     vtkm::CopyFlag::On));
  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arrAs,
                                     vtkm::CopyFlag::On));
}

void
CalcGradAs(vtkm::cont::DataSet& ds)
{
  //Assumed that dAs_phi_ff is R,Z,Phi
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> dAsAH;
  ds.GetField("dAs_phi_ff").GetData().AsArrayHandle(dAsAH);
  int nVals = dAsAH.GetNumberOfValues()/3;

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> dAs_ff;
  dAs_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    dAs_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      dAs_ff[p][i].resize(numNodes);
  }


  int idx = 0;
  auto dAsPortal = dAsAH.ReadPortal();

  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
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
  VEC_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    VEC_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC_ff[p][i].resize(numNodes);
  }

  idx = 0;
  std::vector<vtkm::Vec3f> arr_ff(nVals);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f val = VEC_ff[p][i][n];
        val = dAs_ff[p][i][n];
        arr_ff[idx] = val;
        idx++;
      }
    }
  }

  //do some testing..
  for (int p = 0; p < numPlanes; p++)
  {
    int off0 = p*numNodes*2;
    int off1 = p*numNodes*2 + numNodes;
    for (int n = 0; n < numNodes; n++)
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

  ds.AddField(vtkm::cont::make_Field("dAs_ff_rzp",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_ff,
                                     vtkm::CopyFlag::On));


  /*
  std::cout<<"Print Out dAs_phi_ff"<<std::endl;
  for (int p = 0; p < numPlanes; p++)
    for (int n = 0; n < numNodes; n++)
    {
      auto v0 = dAs_ff[p][0][n];
      auto v1 = dAs_ff[p][1][n];
      if (vtkm::Magnitude(v0) > 1e-7 && vtkm::Magnitude(v1) > 1e-7)
        std::cout<<p<<" :: dAs_phi_ff["<<n<<"]= "<<dAs_ff[p][0][n]<<" "<<dAs_ff[p][1][n]<<std::endl;
    }
  */


#if 0
  //B is R,Z,Phi
  std::cout<<"****************************************************************   Get R,Z,Phi thing right....."<<std::endl;
  //calculate gradPsi;
  std::vector<vtkm::Vec3f> gradPsi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B_RZP").GetData().AsArrayHandle(B0);
  gradPsi.resize(numNodes);
  auto cPortal = coords.ReadPortal();
  auto b0Portal = B0.ReadPortal();
  for (int n = 0; n < numNodes; n++)
  {
    auto b = b0Portal.Get(n);
    auto pt = cPortal.Get(n);
    vtkm::Vec3f gv(b[2] * pt[0], -b[0] * pt[0], 0);
    gradPsi[n] = gv;
  }

  std::vector<std::vector<std::vector<vtkm::Vec3f>>> VEC;
  VEC.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    VEC[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC[p][i].resize(numNodes);
  }

  vtkm::FloatDefault basis = 0.0f;
  std::vector<std::vector<vtkm::FloatDefault>> wphi = {{1.0f, 0.0f}, {0.0f, 1.0f}};

  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f gPsi = gradPsi[n];
        vtkm::Vec3f dAs = dAs_ff[p][i][n];
        vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gPsi);

        vtkm::Vec2f rvec(0,0), zvec(0,0);
        rvec[0] = basis + (1.0-basis) * gammaPsi *   gPsi[0];
        rvec[1] =         (1.0-basis) * gammaPsi *   gPsi[1];
        zvec[0] =         (1.0-basis) * gammaPsi * (-gPsi[1]);
        zvec[1] = basis + (1.0-basis) * gammaPsi *   gPsi[0];

        //R,Phi,Z
        //vec.r = vec[0] = rvec[0]*field.R[0] + zvec[0] * field.Z[0]  //plane 0
        //vec.r +=         rvec[0]*field.R[1] + zvec[0] * field.Z[1]  //plane 1
        //vec.z = vec[2] = rvec[1]*field.R[0] + zvec[1] * field.Z[0]  //plane 0
        //vec.z +=         rvec[1]*field.R[1] + zvec[1] * field.Z[1]  //plane 1
        //vec.phi = vec[1] = field.Phi[0]                             //plane 0
        //vec.phi += vec[1] = field.Phi[1]                            //plane 1
        vtkm::Vec3f vec;
        vec[0] = rvec[0] * dAs[0] + zvec[0] * dAs[2];
        vec[2] = rvec[1] * dAs[0] + zvec[1] * dAs[2];
        vec[1] = dAs[1];

        auto B = b0Portal.Get(n);
        auto Bmag = vtkm::Magnitude(B);

        //parallel derivative to phi derivative
        //vec.phi = (vec.phi*Bmag - vec.R*B.r - vec.z*b.z) / B.phi
        vec[1] = (vec[1]*Bmag - vec[0]*B[0] - vec[2]*B[2]) / B[1];
        VEC[p][i][n] = vec;
      }
    }
  }


  //flatten to 1d index
  idx = 0;
  std::vector<vtkm::Vec3f> arr(nVals);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f val = VEC[p][i][n];
        arr[idx] = val;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("gradAs",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr,
                                     vtkm::CopyFlag::On));
#endif
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

void
SaveOutput(const std::vector<std::vector<vtkm::Vec3f>>& traces,
           const std::vector<std::vector<vtkm::Vec3f>>& punctures,
           const std::vector<std::vector<vtkm::Vec3f>>& puncturesTP,
           const std::string& outFileName = "")
{
  std::string tracesNm, puncNm, puncThetaPsiNm;
  if (outFileName.empty())
  {
    tracesNm = "./traces.txt";
    puncNm = "./punctures.txt";
    puncThetaPsiNm = "./punctures.theta_psi.txt";
  }
  else
  {
    tracesNm = outFileName + ".traces.txt";
    puncNm = outFileName + ".punc.txt";
    puncThetaPsiNm = outFileName + ".punc.theta_psi.txt";
  }

  bool tExists = Exists(tracesNm);
  bool pExists = Exists(puncNm);
  bool ptpExists = Exists(puncThetaPsiNm);
  std::cout<<"EXISTS: "<<tExists<<" "<<pExists<<" "<<ptpExists<<std::endl;

  //Write headers.
  if (!tExists)
    WriteHeader(tracesNm, "ID, R, Z, T");
  if (!pExists)
    WriteHeader(puncNm, "ID, R, Z");
  if (!ptpExists)
    WriteHeader(puncThetaPsiNm, "ID, THETA, PSI");

  std::ofstream outTraces, outPunc, outPuncThetaPsi;

  outTraces.open(tracesNm, std::ofstream::app);
  outPunc.open(puncNm, std::ofstream::app);
  outPuncThetaPsi.open(puncThetaPsiNm, std::ofstream::app);

  //write traces
  for (int i = 0; i < (int)traces.size(); i++)
  {
    for (const auto& pt : traces[i])
    {
      auto R = pt[0];
      auto Z = pt[2];
      auto PHI_N = pt[1];
      while (PHI_N < 0)
        PHI_N += vtkm::TwoPi();

      outTraces<<i<<", "<<R<<", "<<Z<<", "<<PHI_N<<std::endl;
    }
  }

  //write punctures
  for (int i = 0; i < (int)punctures.size(); i++)
  {
    for (const auto& p : punctures[i])
      outPunc<<std::setprecision(12)<<i<<", "<<p[0]<<", "<<p[2]<<std::endl;

    for (const auto& p : puncturesTP[i])
      outPuncThetaPsi<<std::setprecision(12)<<i<<", "<<p[0]<<", "<<p[1]<<std::endl;
  }
}

vtkm::cont::DataSet
ReadData(std::map<std::string, std::vector<std::string>>& args)
{
  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  adios = new adios2::ADIOS;
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);
  //adiosStuff["psi_bicub_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "psi_bicub_acoef", adiosArgs);
  //adiosStuff["one_d_cub_psi_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "one_d_cub_acoef", adiosArgs);
  adiosStuff["interp_coeff"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "interp_coeff", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);

  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield"];
  auto bfield_allStuff = adiosStuff["bfield-all"];
  auto interp_coeffStuff = adiosStuff["interp_coeff"];
  auto equilStuff = adiosStuff["equil"];
  //auto one_dcub_acoefStuff = adiosStuff["one_d_cub_psi_acoef"];

  meshStuff->engine.BeginStep();
  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();
  bfield_allStuff->engine.BeginStep();
  interp_coeffStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  //one_dcub_acoefStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  psi_min = psiVals[0];
  psi_max = psiVals[psiVals.size()-1];

  //Try and do everything in cylindrical coords and worklets.
  auto ds = ReadMesh(meshStuff);
  ReadScalar(dataStuff, ds, "dpot");
  ReadScalar(dataStuff, ds, "apars", "apars", true);
  ReadScalar(meshStuff, ds, "psi");
  ReadScalar(meshStuff, ds, "theta");
  ReadB(bfieldStuff, ds);

  /*
  ReadVec(bfieldStuff, ds, "B", "/node_data[0]/values");
  ReadVec(bfieldStuff, ds, "B3D", "/node_data[0]/values", true, false);
  ReadVec(bfield_allStuff, ds, "Bs", "Bs", true);
  */
  ReadOther(bfieldStuff, ds, "As_phi_ff");
  ReadOther(bfieldStuff, ds, "dAs_phi_ff");
  ReadPsiInterp(equilStuff, interp_coeffStuff, ds);

  //CalcX(ds);

  //ds.PrintSummary(std::cout);
  CalcAs(ds);
  Calc_dAs(ds);
  //CalcV(ds);

  if (0)
  {
    auto ds3d = ReadMesh3D(meshStuff, true);
    ReadScalar(dataStuff, ds3d, "apars", "apars", true, true);
    ReadVec(bfieldStuff, ds3d, "B", "/node_data[0]/values", true, true);
    ReadVec(bfield_allStuff, ds3d, "Bs", "Bs", true, true);
    //CalcX(ds3d);
    vtkm::io::VTKDataSetWriter writer("debug.vtk");
    writer.WriteDataSet(ds3d);
  }

  //ds.PrintSummary(std::cout);
//  vtkm::io::VTKDataSetWriter writer("debug.vtk");
//  writer.WriteDataSet(ds);

  return ds;
}

void
GradientTest()
{
#if 0
  std::vector<vtkm::Vec3f> coords;

  coords.push_back(vtkm::Vec3f(0,0,0));
  coords.push_back(vtkm::Vec3f(1,0,0));
  coords.push_back(vtkm::Vec3f(0,1,0));
  coords.push_back(vtkm::Vec3f(1,1,0));

  std::vector<vtkm::Id> conn;

  //tri 0
  conn.push_back(0);
  conn.push_back(2);
  conn.push_back(1);

  //tri1
  conn.push_back(1);
  conn.push_back(2);
  conn.push_back(3);

  vtkm::cont::DataSetBuilderExplicit dsb;
  auto ds = dsb.Create(coords, vtkm::CellShapeTagTriangle(), 3, conn);

  std::vector<vtkm::Vec3f> vecs;
  vecs.push_back(vtkm::Vec3f(0,0,2));
  vecs.push_back(vtkm::Vec3f(1,0,1));
  vecs.push_back(vtkm::Vec3f(0,1,1));
  vecs.push_back(vtkm::Vec3f(1,1,1));

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(vecs, vtkm::CopyFlag::On)));
  //ds.PrintSummary(std::cout);

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("V");
  gradient.SetOutputFieldName("gradV");
  auto out = gradient.Execute(ds);
  //out.PrintSummary(std::cout);
#endif
}

void
Debug(const vtkm::cont::DataSet& inDS)
{
#if 0
  auto ds = inDS;

  /*
  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetComputeVorticity(true);
  gradient.SetActiveField("B_RZP_2D");
  auto ds = gradient.Execute(inDS);
  */

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  //auto cPortal = coords.ReadPortal();

  vtkm::Vec3f pt(3.029365, 0.020600, 0);

  //Use node position.
  //pt = cPortal.Get(1521);


  vtkm::Vec3f B0 = EvalVector(ds, locator, {pt}, "B_RZP")[0];
  vtkm::Vec3f curlBh = EvalVector(ds, locator, {pt}, "curl_nb_rzp")[0];
  std::cout<<"B0_rzp("<<pt<<")= "<<B0<<std::endl;
  std::cout<<"curlBhat("<<pt<<")= "<<curlBh<<std::endl;
  std::cout<<"**************** DO Z,Phi need to swap????? ***************"<<std::endl;

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  //gradient.SetActiveField("BRZP_Norm");
  gradient.SetActiveField("B_RZP");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);
  auto GRAD = EvalTensor(tmpDS, locator, {pt}, "Gradients", {0})[0];
/*
  std::cout<<"dBr/dr= "<<GRAD[0][0]<<std::endl;
  std::cout<<"dBz/dr= "<<GRAD[0][1]<<std::endl;
  std::cout<<"dBp/dr= "<<GRAD[0][2]<<std::endl;

  std::cout<<"dBr/dz= "<<GRAD[1][0]<<std::endl;
  std::cout<<"dBz/dz= "<<GRAD[1][1]<<std::endl;
  std::cout<<"dBp/dz= "<<GRAD[1][2]<<std::endl;

  std::cout<<"dBr/dp= "<<GRAD[2][0]<<std::endl;
  std::cout<<"dBz/dp= "<<GRAD[2][1]<<std::endl;
  std::cout<<"dBp/dp= "<<GRAD[2][2]<<std::endl;
*/

  auto R = pt[0];
  //auto Z = pt[1];
  auto inv_r = 1.0/R;
  auto Bmag = vtkm::Magnitude(B0);
  auto over_Bmag = 1.0/Bmag;
  auto over_Bmag2 = over_Bmag * over_Bmag;

  auto br = B0[0];
  auto bz = B0[1];
  auto bphi = B0[2];

  auto dbrdr = GRAD[0][0];
  auto dbzdr = GRAD[0][1];
  auto dbpdr = GRAD[0][2];

  auto dbrdz = GRAD[1][0];
  auto dbzdz = GRAD[1][1];
  auto dbpdz = GRAD[1][2];

  auto dbrdp = GRAD[2][0];
  auto dbzdp = GRAD[2][1];
  auto dbpdp = GRAD[2][2];

  std::cout<<"dbrdr= "<<dbrdr<<std::endl;
  std::cout<<"dbrdz= "<<dbrdz<<std::endl;
  std::cout<<"dbrdp= "<<dbrdp<<std::endl;

  std::cout<<"dbzdr= "<<dbzdr<<std::endl;
  std::cout<<"dbzdz= "<<dbzdz<<std::endl;
  std::cout<<"dbzdp= "<<dbzdp<<std::endl;

  std::cout<<"dbpdr= "<<dbpdr<<std::endl;
  std::cout<<"dbpdz= "<<dbpdz<<std::endl;
  std::cout<<"dbpdp= "<<dbpdp<<std::endl;


  //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  //dbdphi=0D0  ! no B perturbation
  auto dbdr = (br*dbrdr + bphi*dbpdr + bz*dbzdr) * over_Bmag;
  auto dbdz = (br*dbrdz + bphi*dbpdz + bz*dbzdz) * over_Bmag;
  auto dbdphi = 0;

  auto div = dbrdr + br/R + dbzdz;
  std::cout<<"Check divervgence: "<<div<<std::endl;

  vtkm::Vec3f curl_B;
  //R curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
  curl_B[0] =          dbzdp*inv_r - dbpdz;
  //Z curl_B(2)  = fld%bphi*inv_r + fld%dbpdr - fld%dbrdp*inv_r
  curl_B[1] =          bphi*inv_r +     dbpdr -     dbrdp*inv_r;
  std::cout<<"    curl_b.z: "<<(bphi*inv_r)<<" + "<<dbpdr<<" - "<<(dbrdp*inv_r)<<std::endl;
  //phi curl_B(3)  = fld%dbrdz - fld%dbzdr
  curl_B[2] =            dbrdz -     dbzdr;

  std::cout<<"curl_B= "<<curl_B<<std::endl;
  std::cout<<"  dbdr= "<<dbdr<<std::endl;
  std::cout<<"  dbdz= "<<dbdz<<std::endl;
  std::cout<<"  dbdphi= "<<dbdphi<<std::endl;


  vtkm::Vec3f curl_nb;

  //R,Z,Phi
  //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
  //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
  //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2

  curl_nb[0] = curl_B[0]*over_Bmag + (bphi * dbdz) * over_Bmag2;
  curl_nb[1] = curl_B[1]*over_Bmag + (-bphi * dbdr) * over_Bmag2;
  curl_nb[2] = curl_B[2]*over_Bmag + (bz * dbdr - br * dbdz) * over_Bmag2;

  std::cout<<"curl_nb= "<<curl_nb<<std::endl;

#if 0
  vtkm::Id offset = 2673266;
  std::vector<vtkm::Id> vids = {1521, 1612, 1613};

  vtkm::cont::ArrayHandle<vtkm::Vec3f> b, coords, curl_bhat;
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(b);
  ds.GetField("Curl_B2D_Norm").GetData().AsArrayHandle(curl_bhat);
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  std::vector<vtkm::Vec3f> B_vals, curl_nb;
  for (const auto& id : vids)
  {
    B_vals.push_back(b.ReadPortal().Get(id));
    curl_nb.push_back(curl_bhat.ReadPortal().Get(id));
  }

  std::vector<vtkm::Vec<vtkm::Vec3f, 3>> grad;
  std::vector<vtkm::Vec3f> pts;
  for (const auto& id : vids)
  {
    grad.push_back(gradients.ReadPortal().Get(id));
    pts.push_back(coords.ReadPortal().Get(id));
  }

  for (int i = 0; i < 3; i++)
  {
    auto c = pts[i];
    auto g = grad[i];

    auto dBr_dR     = g[0][0];
    auto dBz_dR     = g[0][1];
    auto dBphi_dR   = g[0][2];

    auto dBr_dZ     = g[1][0];
    auto dBz_dZ     = g[1][1];
    auto dBphi_dZ   = g[1][2];

    auto dBr_dPhi   = g[2][0];
    auto dBz_dPhi   = g[2][1];
    auto dBphi_dPhi = g[2][2];

    std::cout<<i<<": "<<c<<std::endl;
    std::cout<<"    B= "<<B_vals[i]<<std::endl;
    std::cout<<"    Curl_bhat= "<<curl_nb[i]<<std::endl;
    std::cout<<"    dB*_dR= "<<g[0]<<std::endl;
    std::cout<<"    dB*_dZ= "<<g[1]<<std::endl;
    std::cout<<"    dB*_dP= "<<g[2]<<std::endl;

    std::cout<<"   dBp_dr= "<<dBphi_dR<<std::endl;
    std::cout<<"   dBr_dp= "<<dBr_dPhi<<std::endl;
    std::cout<<"   dBp_dz= "<<dBphi_dZ<<std::endl;
    std::cout<<"   dBz_dr= "<<dBz_dR<<std::endl;
    std::cout<<"   dBr_dz= "<<dBr_dZ<<std::endl;
    std::cout<<"\n\n";
  }
#endif
#endif
}

void
SaveStuff(const vtkm::cont::DataSet& inDS)
{
  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(inDS.GetCoordinateSystem());
  ds.SetCellSet(inDS.GetCellSet());

  ds.AddField(inDS.GetField("B_RZP_Norm"));
  ds.AddField(inDS.GetField("curl_nb_rzp"));
  ds.AddField(inDS.GetField("B_RZP"));

  //Add gradPsi
  vtkm::Id nPts = ds.GetNumberOfPoints();
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B_RZP").GetData().AsArrayHandle(B0);
  auto cPortal = coords.ReadPortal();
  auto b0Portal = B0.ReadPortal();

  std::vector<vtkm::Vec3f> gradPsi(nPts);
  for (vtkm::Id i = 0; i < nPts; i++)
  {
    gradPsi[i][0] = b0Portal.Get(i)[2] * cPortal.Get(i)[0];
    gradPsi[i][1] = -b0Portal.Get(i)[0] * cPortal.Get(i)[0];
    gradPsi[i][2] = 0;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("gradPsi", vtkm::cont::make_ArrayHandle(gradPsi, vtkm::CopyFlag::On)));

  vtkm::io::VTKDataSetWriter writer("stuff.vtk");
  writer.WriteDataSet(ds);
}

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  std::map<std::string, std::vector<std::string>> args;
  ParseArgs(argc, argv, args);

  if (argc < 7)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures stepSize poincVar [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  if (args.find("--gpu") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});
  else if (args.find("--openmp") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP{});
  else if (args.find("--serial") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});

  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  int numPunc = std::atoi(args["--numPunc"][0].c_str());
  std::string vField = args["--vField"][0];

  std::vector<vtkm::Vec3f> seeds;
  bool useWorklet = std::atoi(args["--worklet"][0].c_str());
  bool useTraces = std::atoi(args["--traces"][0].c_str());
  std::string outFileName = args["--output"][0];
  useTurb = true;
  if (args.find("--turbulence") != args.end())
    useTurb = std::atoi(args["--turbulence"][0].c_str());
  auto ds = ReadData(args);
  //ds.PrintSummary(std::cout);
  //return 0;
  /*
  SaveStuff(ds);
  Debug(ds);
  return 0;
  */


  bool useBOnly = false, useHighOrder = false;
  if (args.find("--useBOnly") != args.end()) useBOnly = true;
  if (args.find("--useHighOrder") != args.end()) useHighOrder = true;

  if (args.find("--range") != args.end())
  {
    vtkm::Id numSeeds = std::stoi(args["--numSeeds"][0]);

    auto vals = args["--range"];
    vtkm::FloatDefault r0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault r1 = std::atof(vals[1].c_str());

    vtkm::FloatDefault dr = (r1-r0) / (float)(numSeeds-1);
    vtkm::FloatDefault r = r0;

    for (vtkm::Id id = 0; id < numSeeds; id++, r+=dr)
      seeds.push_back({r, .1, 0});
    std::cout<<"SEEDS= "<<seeds<<std::endl;
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({r, t, z});
  }
  else if (args.find("--jong1") != args.end())
  {
    //first point in traces.v2
    //zone = 3132, incident nodes: 1521, 1612, 1613
    seeds = {{3.029365, 6.183185, 0.020600}};

    //take position at vid=1512
    //seeds = {{2.98872, 0, -0.113239}};

    //first point in log file w/ dpsi_dr
    //seeds = {{3.030292040947820, 0, 0}};

    //testing seed...
    //seeds = {{2,0,0}};
    //seeds = {{2,0,.5}};
    //seeds = {{2,0,-.5}};
    //seeds = {{2.8, 0, -.99}};

    /*
      Jong's changes in diagnosis.F90
      R,Z = 3.0, 0.0: psi,dpsi_dr,dpsi_dz: 1.0227020266024015E-002   4.6889284669612431E-002  -4.6889284669613417E-002
     */
//    seeds = {{3, 0, 0}};

    //first point in MEOW MEOW MEOW
    //Particle: RZP= 2.986310165629829 0.1353622587347799 6.183185307179587
    //seeds = {{2.986310165629829, 6.183185307179587, 0.1353622587347799}};
//    seeds = {{3.019020736951258, 6.183185307179587, 7.1164151319539654E-002}};

    //seed for runs on summit.
    //seeds = {{2.728835848680459, 6.183185307179587, 0.2190207369512611}};
/*
DRP:  derivs_sp MEOW MEOW MEOW ***********************************************
 Particle: RZP=     3.019020736951258        7.1164151319539654E-002     6.183185307179587
  B=   -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860      mag=   0.2200670521901169   4.544069591735594
 ****************************************************************
 ****************************************************************
 Calc curl_B
   over_B    4.544069591735594
   inv_r   0.3312332332668388
   fld.br/bz/phi  -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860
   fld.dbRdx:    7.5098966316331705E-003  -8.0453390321500493E-002     0.000000000000000
   fld.dbZdx:    5.9158099337405949E-002  -5.5705426429988473E-003     0.000000000000000
   fld.dbPdx:    7.2622853854721864E-002   -0.000000000000000     0.000000000000000
   dbdr =   -6.7708980323630388E-002
   dbdz =    1.6843565038781637E-003
   dbdp =     0.000000000000000
   As  =     2.0422619291435113E-006
   Bsa2 =   -1.5573090758542901E-008  -6.2601835972401922E-007   -1.3466549373597934E-006
   Bsa =    9.8643536496446089E-006  -7.1673502824034881E-005   -7.4497222472069995E-006
   dAs =    7.1300391181276443E-005   9.9533406652970541E-006   -4.4521462973275313E-007
fld stuff....
   fld.rzp=     3.019020736951258        7.1164151319539654E-002     6.183185307179587
   fld.psi=    6.9734500328945031E-003
   fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860
   fld.As=    2.0422619291435113E-006
   fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006    4.4521462973275313E-007
   dpsi_dr=    5.4401681238839920E-002
   dpsi_dz=    1.7676215185990881E-002
   deltaB?=     0.000000000000000         0.000000000000000     0.000000000000000
 ****************************************************************
 ****************************************************************


 * bummy ***************************************************************
 * bummy ***************************************************************
 DRP: here in efidld_common.F90
   fld.rzp=     3.019020736951258        7.1164151319539654E-002
    6.183185307179587
   fld.psi=    6.9734500328945031E-003
   fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002
  -0.2192499017639860
   fld.As=    2.0422619291435113E-006            ******** VTKm = As= 2.01509e-06
   fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006  4.4521462973275313E-007
vTKm=  gradAs_rpz= [-7.13004e-05,4.45215e-07,-9.95334e-06]

   dpsi_dr=    5.4401681238839920E-002
   dpsi_dz=    1.7676215185990881E-002
   gamma_psi=    17.48211272627623
   wp=    2.9502307908158221E-002
   iphi=            47
   node=          1617
   irhoX=            47           47
   irho0p1=            47
   wphi=    0.7639437268410916        0.2360562731589084
   rvec=    0.9510563239163462        0.3090175864554083
   zvec=   -0.3090175864554083        0.9510563239163462
   deltaB?=     0.000000000000000         0.000000000000000     0.000000000000000
   As calculation: wrho0=     1.000000000000000         0.000000000000000
                   wrho1=     1.000000000000000         0.000000000000000
                   irho0=            47
                   irho1=            47
                   irho0p1=            47
                   irho1p1=            47
                   node=          1617
                   VECPOTS(0,irho0, node)    2.4447478377585865E-006
                   VECPOTS(1,irho1, node)    2.8654162153181949E-006
                   VECPOTS(0,irho0p1,node)    2.4447478377585865E-006
                   VECPOTS(1,irho1p1,node)    2.8654162153181949E-006

 * zoommy ***************************************************************
 * zoommy ***************************************************************


 */





    /*
      B=  -1.1258751808180879E-002   1.5496320905027214E-002 -0.2216514572458676
      dAs   7.3166248873948126E-005   9.5118595155724718E-005   4.0429622421441507E-006
      PTL1(R,Z,P,rho)    2.986310165629829  0.1353622587347799         6.183185307179587       -4.5917006108385793E-004
      Bsa   9.5002789220787541E-005 -7.3646160159566971E-005  -1.1920069817798782E-005
      Bvec  -1.1258751808180879E-002 1.5496320905027214E-002  -0.2216514572458676
      As   3.0485639983994535E-006
      curl_nb  -1.4504793857055040E-002 -0.3136527384215291       -0.6593914097083449
      curl_B    0.000000000000000 -4.7905885599546799E-018  -0.1419851267560749

      Calc curl_B
        over_B    4.494835378368263
        inv_r   0.3348613990298943
        BVALUE: fld.br/bz/phi  -1.1258751808180879E-002   1.5496320905027214E-002 -0.2216514572458676
        fld.dbRdx:    1.2882152089190432E-002  -7.6545946509187848E-002  0.000000000000000
        fld.dbZdx:    6.5439180246887108E-002  -9.1120307073726311E-003  0.000000000000000
        fld.dbPdx:    7.4222517070366034E-002   -0.000000000000000  0.000000000000000

        dbdr =   -7.0040769970213385E-002
        dbdz =    3.2390182056756529E-003
        dbdp =     0.000000000000000

     //stuff from efield_common.F90
     DRP: here in efidld_common.F90
     fld.rzp=     3.019020736951258 7.1164151319539654E-002  6.183185307179587
     fld.psi=    6.9734500328945031E-003
     fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002 -0.2192499017639860
     fld.As=    2.0422619291435113E-006
     fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006 4.4521462973275313E-007
     dpsi_dr=    5.4401681238839920E-002
     dpsi_dz=    1.7676215185990881E-002
     gamma_psi=    17.48211272627623

     wp=    2.9502307908158221E-002
     iphi=            47
     node=          1617
     irhoX=            47           47
     irho0p1=            47
     wphi=    0.7639437268410916        0.2360562731589084
     rvec=    0.9510563239163462        0.3090175864554083
     zvec=   -0.3090175864554083        0.9510563239163462
     deltaB?=     0.000000000000000         0.000000000000000    0.000000000000000
     */

    /*
      x_ff debugging:
      DRP: field_following_pos2()
      x=    3.019020736951258 7.1164151319539654E-002 -->
            3.021638837619697 6.2535188082394388E-002
      dphi=   1.7275076525105959E-002
      phi   6.183185307179587  -->   6.217735460229799

////use this one.
pIn     2.728835848680459        0.2190207369512611        6.183185307179587
 *************************************************
  i=             1
  x=    2.728835848680459        0.2190207369512611      phi=   6.183185307179587
     bvec_interpol:     2.728835848680459        0.2190207369512611       -->   -1.9935856993796582E-002  -6.4775663199158279E-003  -0.2425649752146412
 b/bphi=   8.2187698269940740E-002   2.6704458523675777E-002 x(1)=   2.728835848680459      dx=    0.2242767373595473      7.2872083739006915E-002
    bvec_interpol: dx=   0.2242767373595473        7.2872083739006915E-002
  dx1=   0.2242767373595473        7.2872083739006915E-002
 k1_rpz= 0.2242767373595448,     0.07287208373900525]                        GOOD

   x_tmp=    2.730773047580803      0.2196501723628287
    tmp= [2.730773047580803,6.181090142925015,0.2196501723628287]



     bvec_interpol:     2.730773047580803        0.2196501723628287       -->
  -1.9978785301185673E-002  -6.2967293967192660E-003  -0.2423929006426939
 b/bphi=   8.2423145431292824E-002   2.5977367241465283E-002 x(1)=
    2.730773047580803      dx=    0.2250789040406072
   7.0938294310101874E-002
    bvec_interpol: dx=   0.2250789040406072        7.0938294310101874E-002
  dx2=   0.2250789040406072        7.0938294310101874E-002
   x_tmp=    2.730779976326204        0.2196334691826448
     bvec_interpol:     2.730779976326204        0.2196334691826448       -->
  -1.9977575372099224E-002  -6.2961967126177248E-003  -0.2423922856247466
 b/bphi=   8.2418362946694565E-002   2.5975235541798633E-002 x(1)=
    2.730779976326204      dx=    0.2250664152164190
   7.0932653097900436E-002
    bvec_interpol: dx=   0.2250664152164190        7.0932653097900436E-002
  dx3=   0.2250664152164190        7.0932653097900436E-002
   x_tmp=    2.732723888226554        0.2202461039616561
     bvec_interpol:     2.732723888226554        0.2202461039616561       -->
  -2.0018867953906675E-002  -6.1149800853820989E-003  -0.2422198608691359
 b/bphi=   8.2647508268210379E-002   2.5245576739414608E-002 x(1)=
    2.732723888226554      dx=    0.2258528201469401
   6.8989190627854940E-002
    bvec_interpol: dx=   0.2258528201469401        6.8989190627854940E-002
  dx4=   0.2258528201469401        6.8989190627854940E-002
   ** x=    2.732723950718343        0.2202461248374213
 *************************************************


     */
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


      /*
      {3.351443028564415449, 0.0, -0.451648806402756176}, //blue 3 islands.
      {3.187329423521033878, 0.0, -0.665017624967372267},
      {1.992020349316277139, 0.0, -0.126203396421661285},
      {3.018666196722858963, 0.0, 0.073864239629065770},
      {3.176582679765305173, 0.0, -0.220557108925872658},  //stochastic region 4000
      {2.179226604128697176, 0.0, 0.291539359807166554},

      {2.552260904008052389, 0, -0.003112355795000767}, //stochastic: 100
      {2.904147293825020348, 0, 0.233676309266207888},  //stochastic: 500
      {2.834116132387578091, 0, 0.301602359825420996},  //stochastic: 1000
      {3.221313019897226848, 0, -0.169787551332608172}, //stochastic: 5000
      {2.1299283004919225 0, -0.17650886057033158},     //semi stoch: pt 10153: xgc spread out, vtkm thin.
      */

    if (vals.size() == 0) //all the points.
      seeds = allSeeds;
    else
    {
      int n = std::stoi(vals[0]);
      if (n >= allSeeds.size())
        std::cout<<"Bad seed!!! #allseeds= "<<allSeeds.size()<<std::endl;

      seeds = {allSeeds[n]};
    }

    /*
allSeeds[0] 50 punctures
--stepSize 0.01   Total error=  0.0120200778585 maxErr= 0.000484904284968
--stepSize 0.005  Total error=  0.00764654588018 maxErr= 0.000306355820869
--stepSize 0.001  Total error=  0.00458809221615 maxErr= 0.000220362884491

-stepSize  0.00050 Total error=  0.00410665950475 maxErr= 0.000212761006711
-stepSize  0.00025 Total error=  0.00395742797169 maxErr= 0.000212790194159
-stepSize  0.00010 Total error=  0.00384771565615 maxErr= 0.000209736086882

allSeeds[1] 50 punctures
--stepSize 0.01   Total error=  0.0246758872498 maxErr= 0.00102710677529
--stepSize 0.005  Total error=  0.0204393203537 maxErr= 0.000897389164912
--stepSize 0.001  Total error=  0.0168112860277 maxErr= 0.00089728753628

allSeeds[2] 50 punctures
--stepSize 0.01  Total error=  0.020815736223 maxErr= 0.000694933702884
--stepSize 0.005 Total error=  0.0157798325015 maxErr= 0.000548627689596
--stepSize 0.001 Total error=  0.012512912312 maxErr= 0.000501483340821


allSeeds[3] 50 punctures
--stepSize 0.01  Total error=  0.00930387190718 maxErr= 0.000376516678999
--stepSize 0.005 Total error=  0.00660720169504 maxErr= 0.000267672744868
--stepSize 0.001 Total error=  0.00501206967247 maxErr= 0.000244691060404


allSeeds[4] 50 punctures ##Stochastic region
--stepSize 0.01    Total error=  3.11058237394 maxErr= 0.417674988212
--stepSize 0.005   Total error=  3.10614570903 maxErr= 0.422411249578
--stepSize 0.001   Total error=  3.0891389475 maxErr= 0.424464351612
--stepSize 0.0005  Total error=  3.10299880864 maxErr= 0.423057887001
--stepSize 0.00025 Total error=  3.10476342536 maxErr= 0.421858627004


10 punctures:
--stepSize 0.000001


allSeeds[5] 50 punctures
--stepSize 0.01  Total error=  0.0173742171359 maxErr= 0.000733842147223
--stepSize 0.005 Total error=  0.0130991456975 maxErr= 0.000589881209904
--stepSize 0.001 Total error=  0.00999171349076 maxErr= 0.000566749317839

allSeeds[5] PID=10153 50 punctures
--stepSize 0.001  Total error=  0.00630002849371 maxErr= 0.000329808548949
--stepSize 0.0005 Total error=  0.00625729353866 maxErr= 0.000329798904507



     */

    //traces.v2 pt near begining.
    //seeds = {{3.024768, 6.070249, 0.049700}};

    //seeds from data/sku_8000/jong.py
    //pts in: xgc_theta_psi.txt, xgc_punctures.txt
//    seeds = {
//      {3.351443028564415449, 0.0, -0.451648806402756176}, //blue 3 islands.
//      {3.187329423521033878, 0.0, -0.665017624967372267},
//      {1.992020349316277139, 0.0, -0.126203396421661285},
//      {3.018666196722858963, 0.0, 0.073864239629065770},
//      {3.176582679765305173, 0.0, -0.220557108925872658},  //stochastic region
//      {2.179226604128697176, 0.0, 0.291539359807166554},

      //stochastic region, i=4000
      /*
s         3.176582679765305173, -0.220557108925872658
p1        2.673415694283511446, -0.404035922651817592
p2        2.366230557551912916, -0.054420066925619182
p3        2.602919798501353910, 0.404861298966848748
p4        3.164266550057113658, 0.258207175353824703
p5        3.087268470563441003, -0.311153651480581217
p6        2.556488556701177028, -0.370963671520662341
p7        2.363954167247493743, 0.104848146584189700
p8        2.744892558803513793, 0.421606828485864227
p9        3.244749979529466088, 0.034975199218512609
p10       2.780834806912544810, -0.455027987386564192
p11       2.329460849125147615, -0.073678279004152566
       */
//    };
  }
  else if (args.find("--jongrz") != args.end())
  {
    //Seed from the blue island.
    //seeds = {{3.351443,             0, -0.451649}};
    seeds = {{3.351443028564415449, 0, -0.45164880640275617552}};
  }
  else if (args.find("--afterN") != args.end())
  {
    seeds = {
      {3.321888620239255019, 0.0, 0.478933972623594384},
      {2.568934684085571352, 0.0, 0.731290913908178353},
      {3.493628202658771720, 0.0, 0.433951677589735296},
      {2.862485694515508605, 0.0, 0.208737305948038576},
      {2.905837753215041008, 0.0, -0.397811882628356817},
      {3.391834939600261389, 0.0, -0.350011953142094434},
    };
  }
  else if (args.find("--parse") != args.end())
  {
//    ./examples/poincare/Simple2.3 --vField B --dir ../data/sku_8000/POINC --worklet 1 --traces 0 --useHighOrder --parse ../data/sku_8000/seeds.txt   --output bumm --numPunc 1000 --stepSize 0.001

    //Generaate the seed list by running jongAll.py
    std::cout<<"READING: "<<args["--parse"][0]<<std::endl;
    std::ifstream seedFile;
    seedFile.open(args["--parse"][0]);
    std::string line;
    while (std::getline(seedFile, line))
    {
      vtkm::FloatDefault r, p, z;
      sscanf(line.c_str(), "%lf, %lf, %lf", &r, &p, &z);
      seeds.push_back({r,p,z});
    }
    std::cout<<"NumSeeds= "<<seeds.size()<<std::endl;
  }

  std::vector<std::vector<vtkm::Vec3f>> traces(seeds.size());
  auto punctures = Poincare(ds, seeds, vField, stepSize, numPunc, useWorklet, useBOnly, useHighOrder, (useTraces ? &traces : nullptr));

  auto puncturesTP = ConvertPuncturesToThetaPsi(punctures, ds);

  std::cout<<"TRACES: "<<traces.size()<<std::endl;

  SaveOutput(traces, punctures, puncturesTP, outFileName);

#if 0
  std::ofstream outPts, outPtsPsiTheta;

  outPts.open("punctures.txt");
  outPts<<"ID,R,Z,T"<<std::endl;
  for (int i = 0; i < (int)punctures.size(); i++)
    for (const auto& p : punctures[i])
    {
      outPts<<i<<", "<<p[0]<<","<<p[2]<<","<<p[1]<<std::endl;
    }

  std::ofstream outTraces, RZ, thetaPsi;
  outTraces.open("traces.txt"), RZ.open("rz.txt"), thetaPsi.open("thetaPsi.txt");
  outTraces<<"ID,R,Z,T,THETA,PSI"<<std::endl;
  RZ<<"ID,R,Z,T"<<std::endl;
  thetaPsi<<"ID,theta,psi,Z"<<std::endl;
  for (int i = 0; i < (int)traces.size(); i++)
  {
    int idx = 0;
    for (const auto& p : traces[i])
    {
      auto R = p[0];
      auto Z = p[2];
      auto PHI = p[1]; //std::fmod(p[1], vtkm::TwoPi());
      int numRevs = 0;
      auto PHI_N = PHI;
      while (PHI_N < 0)
      {
        PHI_N += vtkm::TwoPi();
        numRevs++;
      }
      //auto PHI_N = PHI + (numRevs*vtkm::TwoPi());
//      PHI = p[1];
//      if (PHI < 0) PHI = -PHI;
      auto theta = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
      if (theta < 0) theta += vtkm::TwoPi();
      auto psi = vtkm::Sqrt(((R-eq_axis_r)*(R-eq_axis_r) + Z*Z));

      outTraces<<idx<<", "<<R<<", "<<Z<<", "<<PHI_N<<", "<<theta<<", "<<psi<<std::endl;
//      outTraces<<idx<<", "<<PHI<<" "<<PHI_N<<" nr= "<<numRevs<<std::endl;
      RZ<<idx<<", "<<p[0]<<", "<<p[2]<<", 0"<<std::endl;
      thetaPsi<<idx<<", "<<theta<<", "<<psi<<", 0"<<std::endl;
      idx++;
    }
  }
#endif

  return 0;
}


/*
XGC:
3.029365, 0.020600, 6.183185
3.029293, 0.021400, 6.180291



mine:
3.02936, 0.0206, 6.18318
3.02934, 0.0209196, 6.17947

delta= -0.000047, .0004804, .000821


 */


#if 0
  float phi = -0.01;
  int cnt = 0;
  numPlanes = 48;
  float dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);
  int idx = 0;
  std::cout<<"dPhi= "<<dPhi<<std::endl;
  int rev0 = 0;

  vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);
  std::cout<<"\n\n********************************************************"<<std::endl;
  std::cout<<"phiSpacing= "<<phiSpacing<<std::endl;
  std::cout<<"NumPlanes= "<<numPlanes<<" spacing= "<<phiSpacing<<std::endl;
  std::cout<<"Plane, Phi"<<std::endl;
  for (int i = 0; i < numPlanes; i++)
    std::cout<<i<<", "<<((float)i * phiSpacing)<<"  deg= "<<(i*phiSpacing*57.2958)<<std::endl;
  std::cout<<std::endl<<std::endl;

  while (cnt < 10)
  {

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault T, phi0, phi1;
    GetPlaneIdx(phi, numPlanes, planeIdx0, planeIdx1, phi0, phi1, numRevs, T);
/*
    int  numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
    vtkm::FloatDefault rem = std::fmod(phi, vtkm::TwoPi());

    vtkm::Id planeIdx = static_cast<vtkm::Id>(vtkm::Floor(rem / dPhi));
    vtkm::Id planeIdx0, planeIdx1;
    if (planeIdx < 0)
      planeIdx += numPlanes;
    if (planeIdx == numPlanes-1)
    {
      planeIdx0 = 0;
      planeIdx1 = planeIdx;
    }
    else
    {
      planeIdx0 = planeIdx;
      planeIdx1 = planeIdx0+1; //no //B is going in the NEGATIVE phi direction.
    }
*/
    std::cout<<idx<<":  phi= "<<phi<<" #Rev= "<<numRevs<<" planes: ("<<planeIdx0<<" "<<planeIdx1<<") T= "<<T<<std::endl;
    std::cout<<"*******************************************\n\n\n"<<std::endl;

    if (numRevs > rev0)
    {
      rev0 = numRevs;
      cnt++;
    }

    phi -= 0.05;
    idx++;
  }
  return 0;
#endif
