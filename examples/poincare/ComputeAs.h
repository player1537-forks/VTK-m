
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
  using ControlSignature = void(FieldIn Coords,
                                ExecObject locator,
                                WholeCellSetIn<> cellSet,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                FieldOut AsOut,
                                FieldOut dAsOut);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7);
  using InputDomain = _1;

  ~ComputeAsWorklet()
  {
  }

  ComputeAsWorklet()
  {
    this->NumNodes = numNodes;
    this->NumPlanes = numPlanes;
  }

  template <typename CoordsType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& CoordsRZN,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const AsPhiType& As_phi_ff,
                            const dAsPhiType& dAs_phi_ff_RZP,
                            vtkm::FloatDefault& AsOut,
                            vtkm::Vec3f& dAsOut_RZP) const
  {
    vtkm::FloatDefault R = CoordsRZN[0];
    vtkm::FloatDefault Z = CoordsRZN[1];
    vtkm::Id N = static_cast<vtkm::Id>(CoordsRZN[2]);

    if (idx == 0)
    {
      R = 3.0;
      Z = 0.1;
      N = 12;
    }

    vtkm::Id cellId;
    vtkm::Vec3f ptRZ(R,Z,0), param;
    vtkm::Vec<vtkm::Id,3> vIds;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param);

    AsOut = 0;
    dAsOut_RZP = {0,0,0};
    if (status != vtkm::ErrorCode::Success)
      return;

    auto tmp =  cellSet.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    vtkm::Id offset = N*this->NumNodes;

    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    vals.Append(As_phi_ff.Get(vIds[0]+offset));
    vals.Append(As_phi_ff.Get(vIds[1]+offset));
    vals.Append(As_phi_ff.Get(vIds[2]+offset));
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), AsOut);

    vtkm::VecVariable<vtkm::Vec3f, 3> valv;
    valv.Append(dAs_phi_ff_RZP.Get(vIds[0]+offset));
    valv.Append(dAs_phi_ff_RZP.Get(vIds[1]+offset));
    valv.Append(dAs_phi_ff_RZP.Get(vIds[2]+offset));
    vtkm::exec::CellInterpolate(valv, param, vtkm::CellShapeTagTriangle(), dAsOut_RZP);

    if (idx == 0)
    {
      std::cout<<idx<<" Pt: "<<ptRZ<<" N= "<<N<<" param= "<<param<<std::endl;
      std::cout<<"    cid: "<<cellId<<" vids= "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
      std::cout<<"    As:  "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<" --> "<<AsOut<<std::endl;
      std::cout<<"    dAs: "<<valv[0]<<" "<<valv[1]<<" "<<valv[2]<<" --> "<<dAsOut_RZP<<std::endl;
    }
  }

  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;
};
