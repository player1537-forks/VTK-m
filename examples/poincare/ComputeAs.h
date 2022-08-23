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
                                ExecObject locatorRZ,
                                WholeCellSetIn<> cellSetRZ,
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

  ComputeAsWorklet(const XGCParameters& xgcParams)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;
  }

  template <typename CoordsType, typename CoeffType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType, typename AsOutType, typename dAsOutType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& CoordRZ,            //uniform2D
                            const LocatorType& locatorRZ,         //uniformLocator
                            const CellSetType& cellSetRZ,         //explicit cellset
                            const CoeffType& Coeff_2D,            // bicubic interpolation coeff
                            const AsPhiType& As_phi_ff,
                            const dAsPhiType& dAs_phi_ff_RZP,
                            AsOutType& AsOut,
                            dAsOutType& dAsOut_RZP) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f param;

    //Find the triangle.
    vtkm::ErrorCode status = locatorRZ.FindCell(CoordRZ, cellId, param);

    if (status != vtkm::ErrorCode::Success)
    {
      vtkm::Vec3f tmp(0,0,0);
      for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
      {
        vtkm::Id outIdx = n * this->Num2DPts + idx;
        AsOut.Set(outIdx, 0);
        dAsOut_RZP.Set(outIdx, tmp);
      }

      return;
    }

    auto tmp =  cellSetRZ.GetIndices(cellId);
    vtkm::Vec<vtkm::Id,3> vIds;
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

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

      vtkm::Id outIdx = n * this->Num2DPts + idx;
      AsOut.Set(outIdx, as);
      dAsOut_RZP.Set(outIdx, das_rzp);
    }
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

  vtkm::Id Num2DPts;
  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;
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
    if (statusU != vtkm::ErrorCode::Success)
    {
      std::cout<<"Locator Error at "<<idx<<" "<<coordRZ<<std::endl;
    }

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
