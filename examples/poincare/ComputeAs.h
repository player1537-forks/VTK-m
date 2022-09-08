#include <vtkm/worklet/WorkletMapTopology.h>

#ifndef ComputeAS_h
#define ComputeAS_h

//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
class ComputeAsWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn CoordsRZ,
                                ExecObject locator2L,
                                WholeCellSetIn<> cellSetExplicit,
                                WholeArrayIn Coeff_2D,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayInOut AsOut,
                                WholeArrayInOut dAsOut);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  using InputDomain = _1;

  ~ComputeAsWorklet()
  {
  }

  ComputeAsWorklet(const XGCParameters& xgcParams, int nX, int nY)
    : NumSampsX(nX)
    , NumSampsY(nY)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;

    this->nr = xgcParams.eq_mr-1;
    this->nz = xgcParams.eq_mz-1;
    this->rmin = xgcParams.eq_min_r;
    this->rmax = xgcParams.eq_max_r;
    this->zmin = xgcParams.eq_min_z;
    this->zmax = xgcParams.eq_max_z;
    this->eq_axis_r = xgcParams.eq_axis_r;
    this->eq_axis_z = xgcParams.eq_axis_z;
    this->eq_x_r = xgcParams.eq_x_r;
    this->eq_x_z = xgcParams.eq_x_z;
    this->eq_x_psi = xgcParams.eq_x_psi;
    this->eq_x2_r = xgcParams.eq_axis_r;
    this->eq_x2_z = xgcParams.eq_max_z;
    this->eq_x2_psi = xgcParams.eq_x_psi;
    this->eq_x_slope = -(this->eq_x_r - this->eq_axis_r) / (this->eq_x_z - this->eq_axis_z);
    this->eq_x2_slope = -(this->eq_x2_r - this->eq_axis_r) / (this->eq_x2_z - this->eq_axis_z);

    this->dr = (xgcParams.eq_max_r - xgcParams.eq_min_r) / vtkm::FloatDefault(this->nr);
    this->dz = (xgcParams.eq_max_z - xgcParams.eq_min_z) / vtkm::FloatDefault(this->nz);
    this->dr_inv = 1.0/this->dr;
    this->dz_inv = 1.0/this->dz;

    this->ncoeff_r = xgcParams.eq_mr-1;
    this->ncoeff_z = xgcParams.eq_mz-1;
    this->ncoeff_psi = xgcParams.eq_mpsi-1;
    this->min_psi = xgcParams.psi_min / this->eq_x_psi;
    this->max_psi = xgcParams.psi_max / this->eq_x_psi;
    this->itp_min_psi = xgcParams.itp_min_psi;
    this->itp_max_psi = xgcParams.itp_max_psi;
    this->one_d_cub_dpsi_inv = 1.0 / ((this->max_psi-this->min_psi)/vtkm::FloatDefault(this->ncoeff_psi));

    this->dpsi = (this->max_psi - this->min_psi) / vtkm::FloatDefault(this->NumSampsX);
    this->dtheta = vtkm::TwoPi() / vtkm::FloatDefault(this->NumSampsY);
    this->dpsi_inv = 1.0 / this->dpsi;
    this->dtheta_inv = 1.0 / this->dtheta;

    this->SpacingRZ = vtkm::Vec3f((this->rmax-this->rmin)/vtkm::FloatDefault(this->NumSampsX-1),
                                  (this->zmax-this->zmin)/vtkm::FloatDefault(this->NumSampsY-1),
                                  1);
    this->SpacingPT = vtkm::Vec3f((this->max_psi-this->min_psi)/vtkm::FloatDefault(this->NumSampsX-1),
                                  vtkm::TwoPi()/vtkm::FloatDefault(this->NumSampsY-1),
                                  1);

    this->InvSpacingRZ = vtkm::Vec3f(1.0f/this->SpacingRZ[0],
                                     1.0f/this->SpacingRZ[1],
                                     1.0f/this->SpacingRZ[2]);
    this->InvSpacingPT = vtkm::Vec3f(1.0f/this->SpacingPT[0],
                                     1.0f/this->SpacingPT[1],
                                     1.0f/this->SpacingPT[2]);
  }

  template <typename CoordsType, typename CoeffType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType, typename AsOutType, typename dAsOutType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& uniformRZ,          //uniform2D  gridRZ
                            const LocatorType& locator2L,         //locator2L
                            const CellSetType& cellSetExplicit,   //explicit cellset
                            const CoeffType& Coeff_2D,            // bicubic interpolation coeff
                            const AsPhiType& As_phi_ff,
                            const dAsPhiType& dAs_phi_ff_RZP,
                            AsOutType& AsOut,
                            dAsOutType& dAsOut_RZP) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f param;

    auto R = uniformRZ[0];
    auto Z = uniformRZ[1];

    //R = 7.90; //7.80503;
    //Z = 0.00; //0.831567;

    //R = 6.20;
    //Z = 0.0;
    vtkm::Vec3f CoordRZ(R,Z,0);

    /*
    if (idx == 5049)
    {
      std::cout<<"Weird land...."<<std::endl;
      vtkm::Vec3f invSp(1.0/(double)this->NumSampsX, 1.0/(double)this->NumSampsY, 1.0);
      vtkm::Vec3f origin(this->rmin, this->zmin, 0);
      std::cout<<"invSp: "<<invSp<<" O= "<<origin<<std::endl;

      auto tmp = CoordRZ - origin;
      std::cout<<"   "<<CoordRZ<<" --> "<<tmp<<std::endl;
      tmp = tmp * this->InvSpacingRZ;
      std::cout<<"       ----> "<<tmp<<std::endl;
      vtkm::Id3 logicalCell(tmp);
      vtkm::Id idx2 = (logicalCell[2] * this->NumSampsY + logicalCell[1]) * this->NumSampsX + logicalCell[0];

      std::cout<<"      "<<idx2<<std::endl;
    }
    */


    //Find the triangle.
    vtkm::ErrorCode status = locator2L.FindCell(CoordRZ, cellId, param);

    auto dR = (this->rmax-this->rmin), dZ = (this->zmax-this->zmin);
    auto drInv = dR / vtkm::FloatDefault(this->NumSampsX);
    auto dzInv = dZ / vtkm::FloatDefault(this->NumSampsY);
    auto R_i = this->GetIndex(CoordRZ[0], this->NumSampsX, this->rmin, this->InvSpacingRZ[0]);
    auto Z_i = this->GetIndex(CoordRZ[1], this->NumSampsY, this->zmin, this->InvSpacingRZ[1]);
    /*
    std::cout<<"**********************************************************"<<std::endl;
    std::cout<<"   "<<idx<<" RZ= "<<CoordRZ<<std::endl;
    std::cout<<"        RZ_i= "<<R_i<<" "<<Z_i<<std::endl;
    */
    if (status != vtkm::ErrorCode::Success)
    {
      return;





      /*
      vtkm::Vec3f tmp(-10000,-1, -1);
      for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
      {
        vtkm::Id outIdx = n*(this->NumSampsX * this->NumSampsY) + idx;
        if (outIdx >= AsOut.GetNumberOfValues())
        {
          std::cout<<"Crashing coming at: "<<idx<<" outIdx= "<<outIdx<<std::endl;
        }
        //vtkm::Id outIdx = n * this->Num2DPts + idx;
        AsOut.Set(outIdx, -1);
        dAsOut_RZP.Set(outIdx, tmp);
      }
      return;
      */
    }

    auto tmp =  cellSetExplicit.GetIndices(cellId);
    vtkm::Vec<vtkm::Id,3> vIds;
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    R_i = this->GetIndex(CoordRZ[0], this->NumSampsX, this->rmin, this->InvSpacingRZ[0]);
    Z_i = this->GetIndex(CoordRZ[1], this->NumSampsY, this->zmin, this->InvSpacingRZ[1]);

    vtkm::Vec3f ptPT;
    vtkm::Id PTidx;
    this->CalcPsiTheta(CoordRZ, Coeff_2D, ptPT, PTidx);
    //ptPT[0] = 0.75;

    vtkm::Id idx_x = idx / this->NumSampsX;
    vtkm::Id idx_y = idx % this->NumSampsY;

    auto P_i = this->GetIndex(ptPT[0], this->NumSampsX, this->min_psi, this->InvSpacingPT[0]);
    auto T_i = this->GetIndex(ptPT[1], this->NumSampsY, 0, this->InvSpacingPT[1]);

    //vtkm::Id idxOut = T_i*idx_x + P_i;
    //vtkm::Id idxOut = (this->NumSampsY+T_i) * this->NumSampsX + P_i;
    //vtkm::Id idxOut = (this->NumSampsY+T_i) * this->NumSampsX + P_i;
    vtkm::Id idxOut = P_i*this->NumSampsX + T_i;
    idxOut = T_i*this->NumSampsX + P_i;
    idxOut = P_i*this->NumSampsY + T_i;

    if (idxOut == 5049)
    {
//      std::cout<<"****************************** bum bum"<<std::endl;
//      std::cout<<CoordRZ<<" "<<R_i<<" "<<Z_i<<" idx= "<<idx<<std::endl;
//      std::cout<<"  "<<ptPT<<" "<<P_i<<" "<<T_i<<std::endl;
    }

    AsOut.Set(idxOut, AsOut.Get(idxOut) + ptPT[0]); //ptPT[1]);
    return;



//    if (idx == 10)
//      std::cout<<"*********************** pPT= "<<ptPT<<std::endl;
    /*
    //P_i = 27;
    if (idx % 3 == 0)
      P_i = 16;
    else
      P_i = 29;
    */

//    std::cout<<"        PT_i= "<<P_i<<" "<<T_i<<std::endl;

    //vtkm::Id idx2 =   (this->NumSampsY + logicalCell[1]) * this->NumSampsX + logicalCell[0];
    vtkm::Id outIdx = (this->NumSampsY+T_i) * this->NumSampsX + P_i;

    vtkm::FloatDefault as;
    as = 10;
    AsOut.Set(outIdx, as);

    /*
    for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
    {
      vtkm::Id rzIdx = (Z_i*this->NumSampsX) + R_i;

      vtkm::Id outIdx = (n*this->NumSampsX*this->NumSampsY) + idx;      //this works. Make the red outline shape.
      //vtkm::Id outIdx = (n*this->NumSampsX*this->NumSampsY) + rzIdx;
      //outIdx = (n*this->NumSampsX*this->NumSampsY) + rzIdx;

      auto ptIdx = (T_i * this->NumSampsX) + P_i;
      //ptIdx = (P_i*256) + T_i;
      outIdx = (n*this->NumSampsX*this->NumSampsY) + ptIdx;


      vtkm::FloatDefault as = ptPT[0];
      AsOut.Set(outIdx, as);
    }
    */
    return;

//    vtkm::Vec3f ptPT;
//    vtkm::Id PTidx;
    this->CalcPsiTheta(CoordRZ, Coeff_2D, ptPT, PTidx);

    //Use the param values to set the point for each plane.
    for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
    {
      vtkm::FloatDefault as = 0;
      vtkm::Vec3f das_rzp;

      vtkm::Id offset = n*this->NumNodes;
      vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
      vals.Append(As_phi_ff.Get(vIds[0]+offset));
      vals.Append(As_phi_ff.Get(vIds[1]+offset));
      vals.Append(As_phi_ff.Get(vIds[2]+offset));
      vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), as);

      vtkm::VecVariable<vtkm::Vec3f, 3> valv;
      valv.Append(dAs_phi_ff_RZP.Get(vIds[0]+offset));
      valv.Append(dAs_phi_ff_RZP.Get(vIds[1]+offset));
      valv.Append(dAs_phi_ff_RZP.Get(vIds[2]+offset));
      vtkm::exec::CellInterpolate(valv, param, vtkm::CellShapeTagTriangle(), das_rzp);

      vtkm::Id outIdx = n*(this->NumSampsX * this->NumSampsY) + idx; //works
      //vtkm::Id outIdx = n*(this->NumSampsX * this->NumSampsY) + PTidx;
      as = idx;
      as = ptPT[0];
      as = 10;
      AsOut.Set(outIdx, as);
      dAsOut_RZP.Set(outIdx, das_rzp);
    }
  }

  template <typename CoeffType>
  VTKM_EXEC
  void CalcPsiTheta(const vtkm::Vec3f& ptRZ,
                    const CoeffType& Coeff_2D,
                    vtkm::Vec3f& ptPT,
                    vtkm::Id& idxPT) const
  {
    auto R = ptRZ[0], Z = ptRZ[1];
    auto r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    auto z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = this->rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = this->zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;
    vtkm::Id offset = (z_i * this->ncoeff_r + r_i) * 16;
    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    this->EvalBicub2(R, Z, Rc, Zc, offset, Coeff_2D, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);

    auto theta = vtkm::ATan2(Z-this->eq_axis_z, R-this->eq_axis_r);
    if (theta < 0)
      theta += vtkm::TwoPi();

    //calc idx.
    auto P_i = this->GetIndex(psi, this->NumSampsX, this->min_psi, this->InvSpacingPT[0]);
    auto T_i = this->GetIndex(theta, this->NumSampsY, 0, this->InvSpacingPT[1]);
    /*
    if (P_i < 0 || P_i >= this->NumSampsX)
      std::cout<<" issue with P_i: "<<P_i<<std::endl;
    if (T_i < 0 || T_i >= this->NumSampsY)
      std::cout<<" issue with T_i: "<<T_i<<std::endl;
    */

    psi = psi / this->eq_x_psi;

    ptPT[0] = psi;
    ptPT[1] = theta;
    ptPT[2] = 0;

    idxPT = T_i * this->NumSampsX + P_i;
    //idxPT = P_i * this->NumSampsY + T_i;
  }

  VTKM_EXEC
  int GetIndex(const vtkm::FloatDefault& x,
               const int& nx,
               const vtkm::FloatDefault& xmin,
               const vtkm::FloatDefault& dx_inv) const
  {
    int A = nx-1;
    auto B = (x-xmin)*dx_inv;
    vtkm::FloatDefault B0 = (x-xmin);
    vtkm::FloatDefault B1 = B0*dx_inv;
    int C = int(B1);
    int D = std::min(nx-1, C);
    int idx = std::max(0, std::min(nx-1,
                                   int((x-xmin)*dx_inv)) );
    return idx;

  }

  template <typename CoeffType>
  VTKM_EXEC
  void EvalBicub2(const vtkm::FloatDefault& x,
                  const vtkm::FloatDefault& y,
                  const vtkm::FloatDefault& xc,
                  const vtkm::FloatDefault& yc,
                  const vtkm::Id& offset,
                  const CoeffType& Coeff_2D,
                  vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                  vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02) const
  {
    vtkm::FloatDefault dx = x - xc;
    vtkm::FloatDefault dy = y - yc;

    //fortran code.

    f00 = f01 = f10 = f11 = f20 = f02 = 0.0f;
    vtkm::FloatDefault xv[4] = {1, dx, dx*dx, dx*dx*dx};
    vtkm::FloatDefault yv[4] = {1, dy, dy*dy, dy*dy*dy};
    vtkm::FloatDefault fx[4] = {0,0,0,0};
    vtkm::FloatDefault dfx[4] = {0,0,0,0};
    vtkm::FloatDefault dfy[4] = {0,0,0,0};
    vtkm::FloatDefault dfx2[4] = {0,0,0,0};
    vtkm::FloatDefault dfy2[4] = {0,0,0,0};

    /*
    for (int j = 0; j < 4; j++)
    {
      std::cout<<"acoeff_"<<j<<": ";
      for (int i = 0; i < 4; i++)
        std::cout<<Coeff_2D.Get(offset + j*4 + i)<<" ";
      std::cout<<std::endl;
    }
    */

    for (int j=0; j<4; j++)
    {
      for (int i=0; i<4; i++)
        fx[j] = fx[j] + xv[i]*Coeff_2D.Get(offset + j*4 + i); //acoeff[i][j];
      for (int i=1; i<4; i++)
        dfx[j] = dfx[j] + vtkm::FloatDefault(i)*xv[i-1]*Coeff_2D.Get(offset + j*4 + i); //acoeff[i][j];
      for (int i=2; i<4; i++)
        dfx2[j] = dfx2[j] + vtkm::FloatDefault(i*(i-1))*xv[i-2]*Coeff_2D.Get(offset + j*4 + i); //acoeff[i][j];
    }

    for (int j = 0; j < 4; j++)
    {
      f00 = f00 + fx[j]*yv[j];
      f10 = f10 + dfx[j]*yv[j];
      f20 = f20 + dfx2[j]*yv[j];
    }

    for (int j = 1; j < 4; j++)
    {
      dfy[j] = vtkm::FloatDefault(j)*yv[j-1];
      f01 = f01 + fx[j]*dfy[j];
      f11 = f11 + dfx[j]*dfy[j];
    }

    for (int j = 2; j < 4; j++)
    {
      dfy2[j] = vtkm::FloatDefault(j*(j-1))*yv[j-2];
      f02 = f02 + fx[j]*dfy2[j];
    }
  }

  vtkm::Vec3f SpacingRZ, SpacingPT;
  vtkm::Vec3f InvSpacingRZ, InvSpacingPT;
  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;

  vtkm::Id NumSampsX = 256;
  vtkm::Id NumSampsY = 256;

  int nr, nz;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;
  vtkm::FloatDefault eq_axis_z, eq_axis_r, eq_axis_psi;
  vtkm::FloatDefault eq_x_r, eq_x_z, eq_x_psi, eq_x2_z, eq_x2_r, eq_x2_psi;
  vtkm::FloatDefault eq_x_slope, eq_x2_slope;

  int ncoeff_r, ncoeff_z, ncoeff_psi;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault itp_min_psi, itp_max_psi;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault sml_bp_sign = -1.0f;
  vtkm::FloatDefault dpsi, dtheta, dpsi_inv, dtheta_inv;
};


class AsErrorWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn CoordsRZ,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,

                                WholeCellSetIn<> cellSetUniform,
                                ExecObject locatorUniform,
                                WholeArrayIn As_phi_ff_uniform,
                                WholeArrayIn dAs_phi_ff_RZP_uniform,
                                WholeArrayOut AsError,
                                WholeArrayOut dAsError,
                                FieldOut AsErrorPlane,
                                FieldOut dAsErrorMagPlane,
                                FieldOut dAsErrorRPlane,
                                FieldOut dAsErrorZPlane,
                                FieldOut dAsErrorPPlane);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);
  using InputDomain = _1;

  AsErrorWorklet(const XGCParameters& xgcParams, vtkm::Id num2dPts)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;
    this->Num2DPts = num2dPts;
  }

  template <typename CoordsType, typename CellSetUniformType, typename AsType, typename dAsType, typename LocatorType, typename AsErrorType, typename dAsErrorType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& coordRZ,
                            const AsType& As_phi_ff,
                            const dAsType& dAs_phi_ff,

                            const CellSetUniformType& cellSetUniform,
                            const LocatorType& locatorUniform,
                            const AsType& As_phi_ff_uniform,
                            const dAsType& dAs_phi_ff_uniform,
                            AsErrorType& AsError,
                            dAsErrorType& dAsError,
                            vtkm::FloatDefault& AsErrorPlane,
                            vtkm::FloatDefault& dAsErrorMagPlane,
                            vtkm::FloatDefault& dAsErrorRPlane,
                            vtkm::FloatDefault& dAsErrorZPlane,
                            vtkm::FloatDefault& dAsErrorPPlane) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f param;
    auto statusU = locatorUniform.FindCell(coordRZ, cellId, param);
    /*
    if (statusU != vtkm::ErrorCode::Success)
    {
      std::cout<<"Locator Error at "<<idx<<" "<<coordRZ<<std::endl;
    }
    */

    auto vIds = cellSetUniform.GetIndices(cellId);

    AsErrorPlane = 0;
    dAsErrorMagPlane = 0;
    dAsErrorRPlane = dAsErrorZPlane = dAsErrorPPlane = 0;
    for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
    {
      vtkm::Id offset = n*this->NumNodes;
      vtkm::Id offsetU = n * this->Num2DPts;

      auto As = As_phi_ff.Get(idx + offset);
      auto dAs = dAs_phi_ff.Get(idx + offset);

      vtkm::VecVariable<vtkm::FloatDefault, 4> vals;
      vtkm::VecVariable<vtkm::Vec3f, 4> valsVec;
      for (vtkm::Id i = 0; i < 4; i++)
      {
        vals.Append(As_phi_ff_uniform.Get(vIds[i]+offsetU));
        valsVec.Append(dAs_phi_ff_uniform.Get(vIds[i]+offsetU));
      }

      vtkm::FloatDefault AsU;
      vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagQuad(), AsU);

      vtkm::Vec3f dAsU;
      vtkm::exec::CellInterpolate(valsVec, param, vtkm::CellShapeTagQuad(), dAsU);

      auto errAs = vtkm::Abs(As - AsU);
      AsErrorPlane += errAs;

      AsError.Set(idx + offset, errAs);

      vtkm::FloatDefault errdAs;
      /*
      auto dAsMag = vtkm::Magnitude(dAsU);
      if (dAsMag > 0)
        errdAs = vtkm::Magnitude(dAs - dAsU) / dAsMag;
      else
        errdAs = vtkm::Magnitude(dAs - dAsU);
      */
      auto d_dAs = dAs - dAsU;
      errdAs = vtkm::Magnitude(d_dAs);

      dAsError.Set(idx+offset, errdAs);
      dAsErrorMagPlane += errdAs;
      dAsErrorRPlane += vtkm::Abs(d_dAs[0]);
      dAsErrorZPlane += vtkm::Abs(d_dAs[1]);
      dAsErrorPPlane += vtkm::Abs(d_dAs[2]);
    }
  }

  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;
  vtkm::Id Num2DPts;
};




#endif // ComputeAS_h
