

//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif


//#define DO_TRACES 0

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class PoincareWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldInOut particles,
                                ExecObject locator,
                                WholeCellSetIn<> cellSet,
                                WholeArrayIn Coords,
                                WholeArrayIn B_RZP,
                                WholeArrayIn B_Norm_RZP,
                                WholeArrayIn Curl_NB_RZP,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayIn AsCurlBHat_RZP,
                                WholeArrayInOut output,
                                WholeArrayInOut traces);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);
  using InputDomain = _1;

  PoincareWorklet(vtkm::Id maxPunc, vtkm::FloatDefault planeVal, vtkm::FloatDefault stepSize)
    : MaxIter(maxPunc * 1000000)
    , MaxPunc(maxPunc)
    , PlaneVal(planeVal)
    , StepSize(stepSize)
  {
    this->NumPlanes = numPlanes;
    this->NumNodes = numNodes;
    this->dPhi = vtkm::TwoPi()/static_cast<vtkm::FloatDefault>(this->NumPlanes);
    this->StepSize_2 = this->StepSize / 2.0;
    this->StepSize_6 = this->StepSize / 6.0;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename OutputType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            vtkm::Particle& particle,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const CoordsType& /*coords*/,
                            const BFieldType& B_RZP,
                            const BFieldType& B_Norm_RZP,
                            const BFieldType& Curl_NB_RZP,
                            const AsFieldType& AsPhiFF_RZP,
                            const DAsFieldType& DAsPhiFF_RZP,
                            OutputType& output,
                            OutputType& traces) const
  {
    DBG("Begin: "<<particle<<std::endl);
    while (true)
    {
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particle.Pos<<" #s= "<<particle.NumSteps<<std::endl);
      if (!this->TakeRK4Step(particle, locator, cellSet,
                             B_RZP, B_Norm_RZP, Curl_NB_RZP,
                             AsPhiFF_RZP, DAsPhiFF_RZP, newPos))
      {
        break;
      }

      DBG("     *** Step--> "<<newPos<<std::endl);
      vtkm::Id numRevs0 = vtkm::Floor(vtkm::Abs(particle.Pos[1] / vtkm::TwoPi()));
      vtkm::Id numRevs1 = vtkm::Floor(vtkm::Abs(newPos[1] / vtkm::TwoPi()));

      particle.Pos = newPos;
      particle.NumSteps++;

#ifdef DO_TRACES
      if (particle.NumSteps < this->MaxIter)
        traces.Set(idx*this->MaxIter + particle.NumSteps, particle.Pos);
#endif

      if (numRevs1 > numRevs0)
      {
        vtkm::Id i = (idx * this->MaxPunc) + particle.NumPunctures;
        output.Set(i, particle.Pos);
        particle.NumPunctures++;
        if (idx == 0 && particle.NumPunctures%10 == 0 ) std::cout<<" ***** PUNCTURE n= "<<particle.NumPunctures<<std::endl;
        DBG("************* PUNCTURE n= "<<particle.NumPunctures<<std::endl);
      }

      if (particle.NumSteps >= this->MaxIter || particle.NumPunctures >= this->MaxPunc)
        break;
    }
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType>
  bool TakeRK4Step(vtkm::Particle& particle,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const BFieldType& B_RZP,
                   const BFieldType& B_Norm_RZP,
                   const BFieldType& Curl_NB_RZP,
                   const AsFieldType& AsPhiFF_RZP,
                   const DAsFieldType& DAsPhiFF_RZP,
                   vtkm::Vec3f& res) const
  {
    vtkm::Vec3f tmp, k1, k2, k3, k4, p0;

    p0 = particle.Pos;
    DBG("    ****** K1"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF_RZP, DAsPhiFF_RZP, k1))
      return false;
    tmp = p0 + k1*this->StepSize_2;

    DBG("    ****** K2"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF_RZP, DAsPhiFF_RZP, k2))
      return false;
    tmp = p0 + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF_RZP, DAsPhiFF_RZP, k3))
      return false;
    tmp = p0 + k3*this->StepSize_2;

    DBG("    ****** K4"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF_RZP, DAsPhiFF_RZP, k4))
      return false;
    tmp = p0 + k4*this->StepSize;

    res = p0 + this->StepSize_6*(k1 + 2*k2 + 2*k3 + k4);

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType>
  bool Evaluate(vtkm::Vec3f& ptRPZ,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const BFieldType& B_RZP,
                const BFieldType& B_Norm_RZP,
                const BFieldType& Curl_NB_RZP,
                const AsFieldType& AsPhiFF_RZP,
                const DAsFieldType& DAsPhiFF_RZP,
                vtkm::Vec3f& res) const
  {

    vtkm::Vec3f ptRZ(ptRPZ[0], ptRPZ[2], 0);

    vtkm::Vec3f B0_rzp;
    if (!this->EvaluateVec(ptRZ, locator, cellSet, B_RZP, 0, B0_rzp))
      return false;

    B0_rzp[2] /= ptRZ[0];
    res = vtkm::Vec3f(B0_rzp[0], B0_rzp[1], B0_rzp[2]);

    return true;

#if 0
    std::cout<<" ####################### offset0 > offset1.    how to  handle (47,0)"<<std::endl;
    //vtkm::FloatDefault xxx = AsPhiFF.Get(0);

    //auto x = AsPhiFF.Get(0);
    vtkm::Id offset = 0;
    vtkm::Vec3f B0; //R, Phi, Z
    if (!this->EvaluateVec(pos,
                           locator,
                           cellSet,
                           B0Field,
                           offset,
                           B0))
    {
      return false;
    }
    std::cout<<"  ********$$$$$$$$$$$$$$$$ B0= "<<B0<<std::endl;

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    auto phi = pos[1];
    this->GetPlaneIdx(phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);
    auto PhiMid = Phi0 + (Phi1-Phi0)/2.0;

    DBG("   phiN="<<phiN<<" planes: ("<<planeIdx0<<" "<<planeIdx1<<") :: "<<Phi0<<" "<<Phi1<<" T= "<<T<<" midPhi= "<<PhiMid<<std::endl);

    //Point on midplane.
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray({pos[0], phiN, pos[2]}, B0);

    vtkm::FloatDefault RP_T;
    vtkm::Vec3f ptOnMidPlane, ptOnPlane0, ptOnPlane1;
    bool tmp;
    if (!midPlane.Intersect(ray, RP_T, ptOnMidPlane, tmp) || tmp)
      return false;

    vtkm::Plane<> plane0({0,Phi0, 0}, {0, -1, 0}), plane1({0,Phi1, 0}, {0, -1, 0});
    Ray3f ray0(ptOnMidPlane, -B0), ray1(ptOnMidPlane, B0);
    vtkm::FloatDefault T0, T1;
    if ((!plane0.Intersect(ray0, T0, ptOnPlane0, tmp) || tmp) ||
        (!plane1.Intersect(ray1, T1, ptOnPlane1, tmp) || tmp))
    {
      return false;
    }

    vtkm::Vec3f ptN(pos[0], phiN, pos[2]);
    auto dist01 = vtkm::Magnitude(ptOnPlane1-ptOnPlane0);
    auto dist0i = vtkm::Magnitude(ptN-ptOnPlane0) / dist01;
    auto T01 = dist0i;
    auto T10 = 1-T01;

    DBG("         midP_int= "<<ptOnMidPlane<<" T01: "<<T01<<"  T10: "<<T10<<std::endl);

    vtkm::Id cellId;
    vtkm::Vec3f pcoords;
    vtkm::Vec3f posRZ(ptOnMidPlane[0], ptOnMidPlane[2], 0);
    auto status = locator.FindCell(posRZ, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      DBG("FindCell FAILED: "<<posRZ<<std::endl);
      return false;
    }

    //Get nodeIndices for the triangle.
    //********************************************* assumes offset0 < offset1 !!!!
    auto indices = cellSet.GetIndices(cellId);
    // As_phi_ff: (nPlanes, nNodes, 0:1)
    //Offset for phi_i
    vtkm::Id offset0 = planeIdx0 * this->NumNodes * 2;
    vtkm::Id offset1 = offset0 + this->NumNodes;
    DBG("          indices= "<<indices[0]<<", "<<indices[1]<<", "<<indices[2]<<std::endl);
    DBG("          offset_0= "<<offset<<std::endl);

    //Evaluate As_phi_ff at the P0 and P1.
    vtkm::FloatDefault vals0[3], vals1[3];
    for(int k = 0; k < 3; k++) vals0[k] = AsPhiFF.Get(indices[k]+offset0);
    for(int k = 0; k < 3; k++) vals1[k] = AsPhiFF.Get(indices[k]+offset1);

    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    for (int k = 0; k < 3; k++) vals.Append(vals0[k]*T01 + vals1[k]*T10);

    vtkm::FloatDefault As_ff;
    vtkm::exec::CellInterpolate(vals, pcoords, vtkm::CellShapeTagTriangle(), As_ff);

    DBG("         As_ff at pt: "<<As_ff<<std::endl);
    DBG("          vals= "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<" pcoords= "<<pcoords<<std::endl);
    //Calculate gradpsi
    //We evaluated B0 up above.
    // gradPsi = (dPsi/dR, dPsi/dZ) = (B0.z * R, -B0.r * R)
    vtkm::Vec3f gradPsi(B0[2] * pos[0], -B0[0]*pos[0], 0);
    std::cout<<"              gradPsi= "<<gradPsi<<std::endl;
//    DBG("           gradPsi= "<<gradPsi<<std::end);

    vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi);
    vtkm::FloatDefault basis = 0.0f;

    vtkm::Vec3f rvec(0,0,0), zvec(0,0,0);
    rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi[0];
    rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi[1];
    zvec[0] = basis + (1.0-basis) * gammaPsi * (-gradPsi[1]);
    zvec[1] =         (1.0-basis) * gammaPsi *   gradPsi[0];

    std::cout<<"        rvec= "<<rvec<<" zvec= "<<zvec<<std::endl;

    vtkm::Id voffset0 = planeIdx0 * this->NumNodes * 3;
    vtkm::Id voffset1 = planeIdx1 * this->NumNodes * 3 + this->NumNodes * 3;

    vtkm::Vec3f dAs0[3], dAs1[3];
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
      {
        dAs0[l][k] = dAsPhiFF.Get(voffset0 + indices[l] + k);
        dAs1[l][k] = dAsPhiFF.Get(voffset1 + indices[l] + k);
      }


    //vtkm::Vec3f vec;
    //vec[0] =

    res = {0,-.2,0};
    return true;
#endif
  }

  template <typename LocatorType, typename CellSetType, typename VecFieldType>
  bool EvaluateVec(const vtkm::Vec3f& pos,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const VecFieldType& VecField,
                   const vtkm::Id& offset,
                   vtkm::Vec3f& out) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f pcoords;
    vtkm::Vec3f posRZ(pos[0], pos[2], 0);
    auto status = locator.FindCell(posRZ, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      DBG("FindCell FAILED: "<<posRZ<<std::endl);
      return false;
    }

    auto indices = cellSet.GetIndices(cellId);
    //auto pts = vtkm::make_VecFromPortalPermute(&indices, coords);

    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(VecField.Get(indices[0]+offset));
    vals.Append(VecField.Get(indices[1]+offset));
    vals.Append(VecField.Get(indices[2]+offset));
    vtkm::Vec3f vec;
    vtkm::exec::CellInterpolate(vals, pcoords, vtkm::CellShapeTagTriangle(), vec);

    out[0] = vec[0];
    out[1] = vec[2];
    /*
    if (this->UseCylindrical)
      out[1] /= posRZ[0];
    */
    out[2] = vec[1];

    DBG("          EvaluateVec("<<pos<<") off= "<<offset<<" = "<<vec<<"  divByR= "<<out<<std::endl);

    return true;

  }


  template <typename LocatorType, typename CellSetType, typename CoordsType, typename BFieldType>
  bool EvaluateB(const vtkm::Vec3f& pos,
                 const LocatorType& locator,
                 const CellSetType& cellSet,
                 const CoordsType& coords,
                 const BFieldType& B0Field,
                 vtkm::Vec3f& B0) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f pcoords;
    vtkm::Vec3f posRZ(pos[0], pos[2], 0);
    auto status = locator.FindCell(posRZ, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      DBG("FindCell FAILED: "<<posRZ<<std::endl);
      return false;
    }

    auto indices = cellSet.GetIndices(cellId);
    //auto pts = vtkm::make_VecFromPortalPermute(&indices, coords);

    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(B0Field.Get(indices[0]));
    vals.Append(B0Field.Get(indices[1]));
    vals.Append(B0Field.Get(indices[2]));
    vtkm::Vec3f vec;
    vtkm::exec::CellInterpolate(vals, pcoords, vtkm::CellShapeTagTriangle(), vec);

    B0[0] = vec[0];
    B0[1] = vec[2] / posRZ[0];
    B0[2] = vec[1];

    return true;
  }

  void
  GetPlaneIdx(const vtkm::FloatDefault& phi,
              vtkm::FloatDefault& phiN,
              vtkm::Id& plane0,
              vtkm::Id& plane1,
              vtkm::FloatDefault& phi0,
              vtkm::FloatDefault& phi1,
              vtkm::Id& numRevs,
              vtkm::FloatDefault& T) const
  {
    numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
    //rem = std::fmod(vtkm::Abs(phi), vtkm::TwoPi());
    phiN = phi;
    if (phi < 0)
    {
      //rem = -rem;
      phiN += ((1+numRevs) * vtkm::TwoPi());
    }

    plane0 = vtkm::Floor(phiN / this->dPhi);
    plane1 = plane0 + 1;
    phi0 = static_cast<vtkm::FloatDefault>(plane0)*this->dPhi;
    phi1 = static_cast<vtkm::FloatDefault>(plane1)*this->dPhi;

    DBG("             ***** phi= "<<phi<<" phiN= "<<phiN<<" "<<plane0<<" "<<plane1<<" "<<phi0<<" "<<phi1<<std::endl);

    if (plane1 == this->NumPlanes)
      plane1 = 0;
    T = (phiN-phi0) / (phi1-phi0);
  }

  vtkm::Id MaxIter = 0;
  vtkm::Id MaxPunc = 0;
  vtkm::FloatDefault PlaneVal = 0.0f;
  vtkm::FloatDefault StepSize;
  vtkm::FloatDefault StepSize_2;
  vtkm::FloatDefault StepSize_6;

  vtkm::Id NumPlanes;
  vtkm::Id NumNodes;
  vtkm::FloatDefault dPhi;

  bool UseBOnly = false;
};
