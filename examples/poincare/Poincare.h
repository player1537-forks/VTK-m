

//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

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
                                WholeArrayInOut traces,
                                WholeArrayInOut output,
                                WholeArrayInOut punctureID);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12);
  using InputDomain = _1;

PoincareWorklet(vtkm::Id maxPunc, vtkm::FloatDefault planeVal, vtkm::FloatDefault stepSize, bool saveTraces)
    : MaxIter(maxPunc * 1000000)
    , MaxPunc(maxPunc)
    , PlaneVal(planeVal)
    , StepSize(stepSize)
    , SaveTraces(saveTraces)
  {
    this->NumPlanes = numPlanes;
    this->NumNodes = numNodes;
    this->dPhi = vtkm::TwoPi()/static_cast<vtkm::FloatDefault>(this->NumPlanes);
    this->StepSize_2 = this->StepSize / 2.0;
    this->StepSize_6 = this->StepSize / 6.0;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename OutputType, typename IdType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            vtkm::Particle& particle,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const CoordsType& /*coords*/,
                            const BFieldType& B_RZP,
                            const BFieldType& B_Norm_RZP,
                            const BFieldType& Curl_NB_RZP,
                            const AsFieldType& AsPhiFF,
                            const DAsFieldType& DAsPhiFF_RZP,
                            OutputType& traces,
                            OutputType& output,
                            IdType punctureID) const
  {
    DBG("Begin: "<<particle<<std::endl);
    while (true)
    {
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particle.Pos<<" #s= "<<particle.NumSteps<<std::endl);
      if (!this->TakeRK4Step(particle, locator, cellSet,
                             B_RZP, B_Norm_RZP, Curl_NB_RZP,
                             AsPhiFF, DAsPhiFF_RZP, newPos))
      {
        break;
      }

      DBG("     *** Step--> "<<newPos<<std::endl);
      vtkm::Id numRevs0 = vtkm::Floor(vtkm::Abs(particle.Pos[1] / vtkm::TwoPi()));
      vtkm::Id numRevs1 = vtkm::Floor(vtkm::Abs(newPos[1] / vtkm::TwoPi()));

      particle.Pos = newPos;
      particle.NumSteps++;

      if (this->SaveTraces)
        traces.Set(idx*this->MaxIter + particle.NumSteps, particle.Pos);

      if (numRevs1 > numRevs0)
      {
        vtkm::Id i = (idx * this->MaxPunc) + particle.NumPunctures;
        output.Set(i, particle.Pos);
        punctureID.Set(i, idx);
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
                   const AsFieldType& AsPhiFF,
                   const DAsFieldType& DAsPhiFF_RZP,
                   vtkm::Vec3f& res) const
  {
    vtkm::Vec3f tmp, k1, k2, k3, k4, p0;

    p0 = particle.Pos;
    DBG("    ****** K1"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, k1))
      return false;
    tmp = p0 + k1*this->StepSize_2;

    DBG("    ****** K2"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, k2))
      return false;
    tmp = p0 + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, k3))
      return false;
    tmp = p0 + k3*this->StepSize_2;

    DBG("    ****** K4"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, k4))
      return false;

    res = p0 + this->StepSize_6*(k1 + 2*k2 + 2*k3 + k4);

    return true;
  }

  template <typename PortalType>
    vtkm::FloatDefault
    EvalS(const PortalType& sPortal,
          const vtkm::Id& offset,
          const vtkm::Vec<vtkm::Id, 3>& vId,
          const vtkm::Vec3f& param) const
  {
    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    vals.Append(sPortal.Get(vId[0]+offset));
    vals.Append(sPortal.Get(vId[1]+offset));
    vals.Append(sPortal.Get(vId[2]+offset));

    vtkm::FloatDefault s;
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), s);
    return s;
  }

  template <typename PortalType>
    vtkm::Vec3f
    EvalV(const PortalType& vPortal,
          const vtkm::Id& offset,
          const vtkm::Vec3f& param,
          const vtkm::Vec<vtkm::Id, 3>& vId) const

  {
    //std::cout<<"    ******** vid= "<<vId[0]<<" "<<vId[1]<<" "<<vId[2]<<std::endl;
    //std::cout<<"    ******** vec= "<<vPortal.Get(vId[0])<<" "<<vPortal.Get(vId[1])<<" "<<vPortal.Get(vId[2])<<std::endl;
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(vPortal.Get(vId[0]+offset));
    vals.Append(vPortal.Get(vId[1]+offset));
    vals.Append(vPortal.Get(vId[2]+offset));

    vtkm::Vec3f v;
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), v);
    return v;
  }

  template <typename LocatorType, typename CellSetType>
  bool PtLoc(const vtkm::Vec3f& ptRZ,
             const LocatorType& locator,
             const CellSetType& cs,
             vtkm::Vec3f& param,
             vtkm::Vec<vtkm::Id, 3>& vIds) const
  {
    vtkm::Id cellId;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param);
    if (status != vtkm::ErrorCode::Success)
    {
      std::cout<<"Point not found: "<<ptRZ<<std::endl;
      return false;
    }

    //vtkm::VecVariable<vtkm::Id, 3> tmp;
    auto tmp =  cs.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType>
  bool Evaluate(vtkm::Vec3f& ptRPZ,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const BFieldType& B_RZP,
                const BFieldType& B_Norm_RZP,
                const BFieldType& Curl_NB_RZP,
                const AsFieldType& AsPhiFF,
                const DAsFieldType& DAsPhiFF_RZP,
                vtkm::Vec3f& res) const
  {
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    this->GetPlaneIdx(Phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);

    vtkm::Vec3f ptRZ(R, Z, 0);

    vtkm::Vec3f particlePos_param;
    vtkm::Vec<vtkm::Id,3> particlePos_vids;
    if (!this->PtLoc(ptRZ, locator, cellSet, particlePos_param, particlePos_vids))
      return false;

    auto B0_rzp = this->EvalV(B_RZP, 0, particlePos_param, particlePos_vids);
    //std::cout<<"Meow: "<<ptRPZ<<std::endl;
    //std::cout<<"    B0= "<<B0_rzp<<std::endl;

    if (this->UseBOnly)
    {
      B0_rzp[2] /= ptRZ[0];
      res = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);
      return true;
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
    this->PtLoc(x_ff_rzp, locator, cellSet, x_ff_param, x_ff_vids);
    auto dAs_ff0_rzp = this->EvalV(DAsPhiFF_RZP, offsets[0], x_ff_param, x_ff_vids);
    auto dAs_ff1_rzp = this->EvalV(DAsPhiFF_RZP, offsets[1], x_ff_param, x_ff_vids);

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
    //auto AsCurl_bhat_rzp = this->EvalV(AsCurlBHat_RZP, 0, x_ff_vids, x_ff_param);

    //vtkm::Vec3f curl_nb_rzp = EvalVector(ds, locator, {ptRZ}, "curl_nb_rzp")[0];
    //std::cout<<"    pos_ids= "<<particlePos_vids<<std::endl;
    //std::cout<<"    pos_parms= "<<particlePos_param<<std::endl;
    auto curl_nb_rzp = this->EvalV(Curl_NB_RZP, 0, particlePos_param, particlePos_vids);

    //auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
    //vtkm::FloatDefault As_ff0 = As_ff[0];
    //vtkm::FloatDefault As_ff1 = As_ff[1];
    auto As_ff0 = this->EvalS(AsPhiFF, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = this->EvalS(AsPhiFF, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * curl_nb_rzp;
    //std::cout<<"    As= "<<As<<std::endl;
    //std::cout<<"    curl_nb_rzp= "<<curl_nb_rzp<<std::endl;
    //std::cout<<"    curl_nb_rzp.size()= "<<Curl_NB_RZP.GetNumberOfValues()<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v0)= "<<Curl_NB_RZP.Get(particlePos_vids[0])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v1)= "<<Curl_NB_RZP.Get(particlePos_vids[1])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v2)= "<<Curl_NB_RZP.Get(particlePos_vids[2])<<std::endl;
    //std::cout<<"    AsCurl_bhat_rzp= "<<AsCurl_bhat_rzp<<std::endl;

    //vtkm::Vec3f bhat_rzp = EvalVector(ds, locator, {ptRZ}, "B_RZP_Norm")[0];
    auto bhat_rzp = this->EvalV(B_Norm_RZP, 0, particlePos_param, particlePos_vids);
    //std::cout<<"    bhat_rzp= "<<bhat_rzp<<std::endl;

    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    //std::cout<<"    gradAs_rzp= "<<gradAs_rzp<<std::endl;
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);

    //std::cout<<"    deltaB= "<<deltaB_rzp<<std::endl;

    deltaB_rzp[2] /= R;
    B0_rzp[2] /= R;

    vtkm::Vec3f vec_rzp = B0_rzp + deltaB_rzp;
    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    res = vec_rpz;
    //std::cout<<"    vec_rpz= "<<vec_rpz<<std::endl<<std::endl;

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename VecFieldType>
  bool EvaluateVec(const vtkm::Vec3f& pos_rpz,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const VecFieldType& VecField,
                   const vtkm::Id& offset,
                   vtkm::Vec3f& out) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f pcoords;
    vtkm::Vec3f posRZ(pos_rpz[0], pos_rpz[2], 0);
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
    std::cout<<"   EvaluateVec: "<<vec<<std::endl;

    out[0] = vec[0];
    out[1] = vec[2];
    /*
    if (this->UseCylindrical)
      out[1] /= posRZ[0];
    */
    out[2] = vec[1];

    DBG("          EvaluateVec("<<pos_rpz<<") off= "<<offset<<" = "<<vec<<"  divByR= "<<out<<std::endl);

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
  bool SaveTraces = false;
};
