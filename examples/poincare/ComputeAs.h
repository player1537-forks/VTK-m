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
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayInOut AsOut,
                                WholeArrayInOut dAsOut);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7);
  using InputDomain = _1;

  ~ComputeAsWorklet()
  {
  }

  ComputeAsWorklet(const XGCParameters& xgcParams)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;
  }

  template <typename CoordsType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType, typename AsOutType, typename dAsOutType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& CoordRZ,
                            const LocatorType& locatorRZ,
                            const CellSetType& cellSetRZ,
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
