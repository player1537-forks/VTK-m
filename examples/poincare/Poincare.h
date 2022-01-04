//
//./examples/poincare/Simple2.3 --vField B --dir ../data/sku_8000/POINC --worklet 1 --traces 1 --numPunc 2 --stepSize 0.01 --useHighOrder --jong1 --output bumm



//XGC code:  /gpfs/alpine/proj-shared/csc143/pugmire/XGC-Devel/XGC-Devel-poincare-debug
// make xgc-eem
// make xgc-eem-cpp
//run: bsub job-summit.sh
/*

Jong's branch w/ coeff and debug stuff.
/gpfs/alpine/proj-shared/csc143/jyc/summit/exp-xgc-poincare/exp-poincare-8000-dave

//run: bsub job-summit.sh


//old code (that runs...)
/gpfs/alpine/proj-shared/csc143/pugmire/XGC-Devel/XGC-Devel-poincare-debug
launched job: 1714079
wrong number of procs.
fixed script and ran again. Looks like it works.

run with DDT

bsub -q debug -nnodes 6  -P PHY122 -W 0:30 -Is $SHELL -l

[ -z $JOBSIZE ] && JOBSIZE=$(((LSB_DJOB_NUMPROC-1)/42))
[ -z $ISTEP ] && ISTEP=3000
WDIR=exp-poincare-$ISTEP
SDIR=poincare_plot_for_su419.org
export OMP_NUM_THREADS=2
export COLUMNS=512
export OMPI_MCA_coll_ibm_collselect_mode_barrier=failsafe

date



ddt --connect jsrun -n $((JOBSIZE*32)) -a1 -c1 -g0 -r32 -brs /usr/bin/stdbuf -oL -eL ./xgc-eem 2>&1 | tee run-$JOBID.log



*/


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
                                WholeArrayIn coeff_1D,
                                WholeArrayIn coeff_2D,
                                WholeArrayInOut traces,
                                WholeArrayInOut output,
                                WholeArrayInOut punctureID);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);
  using InputDomain = _1;

PoincareWorklet(vtkm::Id maxPunc, vtkm::FloatDefault planeVal, vtkm::FloatDefault stepSize, bool saveTraces)
    : MaxIter(maxPunc * 100000000)
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


    std::cout<<"These are hard coded!!! Fix me!!"<<std::endl;

    this->nr = 150; //eq_mr;
    this->nz = 150; //eq_mz;
    this->rmin = eq_min_r;
    this->rmax = eq_max_r;
    this->zmin = eq_min_z;
    this->zmax = eq_max_z;
    this->dr = (eq_max_r - eq_min_r) / vtkm::FloatDefault(this->nr);
    this->dz = (eq_max_z - eq_min_z) / vtkm::FloatDefault(this->nz);
    this->dr_inv = 1.0/this->dr;
    this->dz_inv = 1.0/this->dz;


    this->min_psi = 0.0;
    this->max_psi = 0.0697345;
    this->one_d_cub_dpsi_inv = 1.0/.0004649;
    this->ncoeff = 150;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC bool HighOrderB(const vtkm::Vec3f& ptRPZ,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const BFieldType& B_RZP,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            vtkm::Vec3f& res) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2];

//    std::cout<<"***************************************"<<std::endl;
//    std::cout<<"HighOrderB"<<std::endl;
//    std::cout<<" ptRPZ= "<<ptRPZ<<std::endl;

    vtkm::Vec3f ptRZ(R,Z,0);

    if (0)
    {
      vtkm::Vec3f particlePos_param;
      vtkm::Vec<vtkm::Id,3> particlePos_vids;
      if (!this->PtLoc(ptRZ, locator, cellSet, particlePos_param, particlePos_vids))
        return false;
      std::cout<<"B: ids= "<<particlePos_vids<<std::endl;
      std::cout<<"B: par= "<<particlePos_param<<std::endl;
      auto B0_rzp = this->EvalV(B_RZP, 0, particlePos_param, particlePos_vids);
      B0_rzp[2] /= R;
      auto res_rpz = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

      std::cout<<"B_rpz("<<R<<" "<<Z<<") = "<<res_rpz<<std::endl;
    }

    int r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);
    //std::cout<<"r_i, z_i = "<<r_i<<" "<<z_i<<std::endl;
    /*
    r_i = 76-1;
    z_i = 76-1;
    r_i = 88-1;
    z_i = 76-1;
    */

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;
    /*
    Rc = this->rmin + (this->rmax-this->rmin)/(this->nr-1) * vtkm::FloatDefault(r_i+1);
    Zc = this->zmin + (this->zmax-this->zmin)/(this->nr-1) * vtkm::FloatDefault(z_i+1);
    auto Rc_1 = this->rmin + (this->rmax-this->rmin)/(this->nr-1) * vtkm::FloatDefault(r_i+2);
    auto Zc_1 = this->zmin + (this->zmax-this->zmin)/(this->nr-1) * vtkm::FloatDefault(z_i+2);
    Rc = (Rc+Rc_1)/2.;
    Zc = (Zc+Zc_1)/2.;

    std::cout<<"Fix me: "<<__LINE__<<std::endl;
    std::cout<<"Rc, Zc= "<<Rc<<" "<<Zc<<std::endl;

    Rc = 2.999976000000000;
    Zc = 7.9990400000000666E-003;

    std::cout<<"SHOULD BE: Rc, Zc= "<<Rc<<" "<<Zc<<std::endl;
    std::cout<<"(r_i,z_i)= "<<r_i<<" "<<z_i<<"  Rc,Zc= "<<Rc<<" "<<Zc<<std::endl;
    */

    //Get the coeffcients (z,r,4,4)
    vtkm::Matrix<vtkm::FloatDefault, 4, 4> acoeff;
    //offset = ri * nz + zi
    vtkm::Id offset = (r_i * this->ncoeff + z_i) * 16;
    //offset = (z_i * this->ncoeff + r_i)*16;
    vtkm::Id idx = 0;
    //std::cout<<"Offset= "<<(offset/16)<<" 16: "<<offset<<std::endl;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        acoeff[ii][jj] = Coeff_2D.Get(offset+idx);
        //acoeff[jj][ii] = Coeff_2D.Get(offset+idx); //z,r
        //std::cout<<"c_"<<ii<<jj<<"= "<<Coeff_2D.Get(offset+idx)<<std::endl;
        idx++;
      }

    double psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);
    /*
    std::cout<<" psi= "<<psi<<std::endl;
    std::cout<<" dpsi_dr = "<<dpsi_dr<<std::endl;
    std::cout<<" dpsi_dz = "<<dpsi_dz<<std::endl;
    */

    vtkm::FloatDefault fld_I = this->I_interpol(psi, 0, Coeff_1D);
    //vtkm::FloatDefault fld_dIdpsi = this->I_interpol(psi, 1, Coeff_1D);

    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault Br = -dpsi_dz * over_r;
    vtkm::FloatDefault Bz = dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r * this->sml_bp_sign;

    res = vtkm::Vec3f(Br, Bp, Bz);
    /*
    std::cout<<"BUM_B_rpz= "<<res<<std::endl;
    std::cout<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    */

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType, typename OutputType, typename IdType>
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
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            OutputType& traces,
                            OutputType& output,
                            IdType punctureID) const
  {
    DBG("Begin: "<<particle<<std::endl);
    /*
    //values for dopri
    vtkm::FloatDefault facold = 1e-4;
    vtkm::FloatDefault hlamb = 0.0;
    vtkm::FloatDefault h = this->StepSize;
    */

    while (true)
    {
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particle.Pos<<" #s= "<<particle.NumSteps<<std::endl);
      if (!this->TakeRK4Step(particle, locator, cellSet,
                             B_RZP, B_Norm_RZP, Curl_NB_RZP,
                             AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, newPos))
      {
        std::cout<<"*****************       All done: RK step failed."<<std::endl;
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
      {
        std::cout<<"************************************* All done: "<<particle<<std::endl;
        break;
      }
    }
    std::cout<<"Particle done: "<<idx<<std::endl;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType>
  bool TakeDoPriStep(vtkm::Particle& particle,
                     const LocatorType& locator,
                     const CellSetType& cellSet,
                     const BFieldType& B_RZP,
                     const BFieldType& B_Norm_RZP,
                     const BFieldType& Curl_NB_RZP,
                     const AsFieldType& AsPhiFF,
                     const DAsFieldType& DAsPhiFF_RZP,
                     vtkm::FloatDefault& h,
                     vtkm::FloatDefault& facold,
                     vtkm::FloatDefault& hlamb,
                     vtkm::Vec3f& res) const
  {
    int n_accepted = 0, n_rejected = 0, n_steps = 0, n_eval = 0;
    static const double safe = 0.9;
    static const double epsilon = std::numeric_limits<double>::epsilon();
    static const double facl = 0.2;
    static const double facr = 10.0;
    static const double beta = 0.04;
    static const unsigned int nstiff = 100;
    static const double a21=0.2, a31=3.0/40.0, a32=9.0/40.0, a41=44.0/45.0,
                    a42=-56.0/15.0, a43=32.0/9.0, a51=19372.0/6561.0,
                    a52=-25360.0/2187.0,
                    a53=64448.0/6561.0, a54=-212.0/729.0, a61=9017.0/3168.0,
                    a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0,
                    a65=-5103.0/18656.0,
                    a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0,
                    a75=-2187.0/6784.0, a76=11.0/84.0;

    static const double e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0,
      e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;

    auto t_local = particle.Time;

    // stepsize underflow??
    if (0.1*h <= t_local*epsilon)
      return false;

    auto p0 = particle.Pos;
    auto t = particle.Time;
    vtkm::Vec3f yCur = p0;
    vtkm::Vec3f vCur;
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, vCur))
      return false;

    bool reject = false;
    int iasti = 0, nonsti = 0;
    vtkm::Vec3f y_stiff;
    vtkm::FloatDefault direction = 1.0;

    vtkm::Vec3f k1 = vCur;

    vtkm::Vec3f k2 = a21*k1;
    auto y_new = yCur + h*k2;

    vtkm::Vec3f k3 = a31*k1 + a32*k2;
    y_new = yCur + h*k3;

    vtkm::Vec3f k4 = a41*k1 + a42*k2 + a43*k3;
    y_new = yCur + h*k4;

    vtkm::Vec3f k5 = a51*k1 + a52*k2 + a53*k3 + a54*k4;
    y_new = yCur + h*k5;

    vtkm::Vec3f k6 = a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5;
    y_new = yCur + h * k6;

    vtkm::Vec3f k7 = a71*k1 + a73*k3 + a74*k4 + a75*k5 + a76*k6;
    y_new = yCur + h * k7;


    vtkm::Vec3f ee = h * ( e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7 );
    double sk, sqr;

    vtkm::FloatDefault err = 0.0, h_new, fac11;
    for( size_t i=0; i<3; i++ )
    {
      sk = abstol + reltol * std::max(std::abs(yCur[i]), std::abs(y_new[i]));
      sqr = ee[i]/sk;
      err += sqr*sqr;
    }

    err = vtkm::Sqrt(err / 3.0);

    // compute next potential stepsize
    fac11 = pow( err, 0.2 - beta*0.75 );

    // Lund-stabilization
    double fac = fac11 / pow( facold, beta );

    // we require facl <= h_new/h <= facr
    fac = std::max( 1.0/facr, std::min( 1.0/facl, fac/safe ) );

    h_new = h / fac;
    std::cout<<"h_new= "<<h_new<<" was "<<h<<" fac= "<<fac<<std::endl;

    if( h_new > std::numeric_limits<double>::max() )
      h_new = std::numeric_limits<double>::max();

    if( h_new < -std::numeric_limits<double>::max() )
      h_new = std::numeric_limits<double>::max();

    if( err <= 1.0 )
    {
      // step accepted
      facold = std::max( err, 1.0e-4 );
      n_accepted++;

      // stiffness detection
      if( !(n_accepted % nstiff) || (iasti > 0) )
      {
        double stnum = 0.0, stden = 0.0, sqr_;

        for( size_t i=0; i < 3; i++ )
        {
          sqr_ = k7[i] - k6[i];
          stnum += sqr_ * sqr_;
          sqr_ = y_new[i] - y_stiff[i];
          stden += sqr_ * sqr_;
        }

        if( stden > 0.0 )
          hlamb = h * sqrt( stnum/stden );

        if( fabs(hlamb) > 3.25 )
        {
          nonsti = 0;
          iasti++;

          if( iasti == 15 )
          {
            return false;
            /*
            if (DebugStream::Level5())
            {
              debug5 << "\tavtIVPDopri5::Step(): exiting at t = "
                     << t << ", problem seems stiff (y = " << yCur
                     << ")\n";
            }
            return avtIVPSolver::STIFFNESS_DETECTED;
            */
          }
        }
        else
        {
          nonsti++;
          if( nonsti == 6 )
            iasti = 0;
        }

        // --- step looks ok - prepare for return
        if( reject )
          h_new = direction * std::min( std::abs(h_new), std::abs(h) );

        yCur = y_new;
        vCur = k7;
        t = t+h;

//        if( period && last )
//          t += epsilon;

        // Set the step size on sucessful step.
        h = h_new;
        return true;
      }
      else
      {
        // step rejected
        h_new = h / std::min( 1.0/facl, fac11/safe );
        reject = true;

        if( n_accepted >= 1 )
          n_rejected++;

        // Update the step size to the new step size.
        h = h_new;
      }
    }

    res = y_new;
    return true;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  bool TakeRK4Step(vtkm::Particle& particle,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const BFieldType& B_RZP,
                   const BFieldType& B_Norm_RZP,
                   const BFieldType& Curl_NB_RZP,
                   const AsFieldType& AsPhiFF,
                   const DAsFieldType& DAsPhiFF_RZP,
                   const Coeff_1DType& Coeff_1D,
                   const Coeff_2DType& Coeff_2D,
                   vtkm::Vec3f& res) const
  {
    vtkm::Vec3f tmp, k1, k2, k3, k4, p0;

    //k1 = F(p)
    //k2 = F(p+hk1/2)
    //k3 = F(p+hk2/2)
    //k4 = F(p+hk3)
    //Yn+1 = Yn + 1/6 h (k1+2k2+2k3+k4)
    p0 = particle.Pos;
    DBG("    ****** K1"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k1))
      return false;
    tmp = p0 + k1*this->StepSize_2;

    DBG("    ****** K2"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k2))
      return false;
    tmp = p0 + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k3))
      return false;
    tmp = p0 + k3*this->StepSize;

    DBG("    ****** K4"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, Curl_NB_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k4))
      return false;

    vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4)/6.0;
    res = p0 + this->StepSize * vec;

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

  int GetIndex(const vtkm::FloatDefault& x,
               const int& nx,
               const vtkm::FloatDefault& xmin,
               const vtkm::FloatDefault& dx_inv) const
  {
    //std::cout<<"GetIndex: "<<x<<" "<<nx<<" min= "<<xmin<<" dx_inv "<<dx_inv<<std::endl;
    //return std::max(0, std::min(nx-1, (int)((x-xmin)*dx_inv)) );
    //return std::max(0, std::min(nx-1,    (int)((x-xmin)*dx_inv)) );
    int idx = std::max(1, std::min(nx  , 1 + int ((x-xmin)*dx_inv)) );
    return idx-1;
  }

  bool eval_bicub_2(const vtkm::FloatDefault& x,
                    const vtkm::FloatDefault& y,
                    const vtkm::FloatDefault& xc,
                    const vtkm::FloatDefault& yc,
                    const vtkm::Matrix<vtkm::FloatDefault, 4, 4>& acoeff,
                    vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                    vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02) const
  {
    double dx = x - xc;
    double dy = y - yc;

    //fortran code.

    f00 = f01 = f10 = f11 = f20 = f02 = 0.0f;
    double xv[4] = {1, dx, dx*dx, dx*dx*dx};
    double yv[4] = {1, dy, dy*dy, dy*dy*dy};
    double fx[4] = {0,0,0,0};
    double dfx[4] = {0,0,0,0};
    double dfy[4] = {0,0,0,0};
    double dfx2[4] = {0,0,0,0};
    double dfy2[4] = {0,0,0,0};

    for (int j=0; j<4; j++)
    {
      for (int i=0; i<4; i++)
        fx[j] = fx[j] + xv[i]*acoeff[i][j];
      for (int i=1; i<4; i++)
        dfx[j] = dfx[j] + double(i)*xv[i-1]*acoeff[i][j];
      for (int i=2; i<4; i++)
        dfx2[j] = dfx2[j] + double(i*(i-1))*xv[i-2]*acoeff[i][j];
    }

    for (int j = 0; j < 4; j++)
    {
      f00 = f00 + fx[j]*yv[j];
      f10 = f10 + dfx[j]*yv[j];
      f20 = f20 + dfx2[j]*yv[j];
    }

    for (int j = 1; j < 4; j++)
    {
      dfy[j] = double(j)*yv[j-1];
      f01 = f01 + fx[j]*dfy[j];
      f11 = f11 + dfx[j]*dfy[j];
    }

    for (int j = 2; j < 4; j++)
    {
      dfy2[j] = double(j*(j-1))*yv[j-2];
      f02 = f02 + fx[j]*dfy2[j];
    }






    //c++ code.
    /*
    double fx_i, dfx_i, dfx2_i;

    f00 = 0;
    f01 = 0;
    f02 = 0;

    fx_i = ((acoeff[0][3]*dx + acoeff[0][2])*dx + acoeff[0][1])*dx + acoeff[0][0];
    f00 = f00 + fx_i;

    fx_i = ((acoeff[1][3]*dx + acoeff[1][2])*dx + acoeff[1][1])*dx + acoeff[1][0];
    f00 = f00 + dy*fx_i;
    f01 = f01 +    fx_i;

    fx_i = ((acoeff[2][3]*dx + acoeff[2][2])*dx + acoeff[2][1])*dx + acoeff[2][0];
    f00 = f00 + dy*(dy*fx_i);
    f01 = f01 + 2.0*(dy*fx_i);
    f02 = f02 + 2.0*fx_i;

    fx_i = ((acoeff[3][3]*dx + acoeff[3][2])*dx + acoeff[3][1])*dx + acoeff[3][0];
    f00 = f00 + dy*(dy*(dy*fx_i));
    f01 = f01 + 3.0*(dy*(dy*fx_i));
    f02 = f02 + 6.0*(dy*fx_i);

    f10 = 0;
    f11 = 0;

    dfx_i = (acoeff[0][3]*3.0*dx + acoeff[0][2]*2.0)*dx + acoeff[0][1];
    f10 = f10 + dfx_i;

    dfx_i = (acoeff[1][3]*3.0*dx + acoeff[1][2]*2.0)*dx + acoeff[1][1];
    f10 = f10 + dy*dfx_i;
    f11 = f11 +    dfx_i;

    dfx_i = (acoeff[2][3]*3.0*dx + acoeff[2][2]*2.0)*dx + acoeff[2][1];
    f10 = f10 + dy*(dy*dfx_i);
    f11 = f11 + 2.0*(dy*dfx_i);

    dfx_i = (acoeff[3][3]*3.0*dx + acoeff[3][2]*2.0)*dx + acoeff[3][1];
    f10 = f10 + dy*(dy*dy*dfx_i);
    f11 = f11 + 3.0*(dy*dy*dfx_i);

    f20 = 0;

    dfx2_i = (3.*2.*acoeff[0][3]*dx + 2.*1.*acoeff[0][2]);
    f20 = f20 + dfx2_i;

    dfx2_i = (3.*2.*acoeff[1][3]*dx + 2.*1.*acoeff[1][2]);
    f20 = f20 + dy*dfx2_i;

    dfx2_i = (3.*2.*acoeff[2][3]*dx + 2.*1.*acoeff[2][2]);
    f20 = f20 + (dy*dy)*dfx2_i;

    dfx2_i = (3.*2.*acoeff[3][3]*dx + 2.*1.*acoeff[3][2]);
    f20 = f20 + (dy*dy)*dy*dfx2_i;
    */
    return true;
  }

  template <typename Coeff_1DType>
  vtkm::FloatDefault I_interpol(const vtkm::FloatDefault& psi,
                                const int& ideriv,
                                const Coeff_1DType& coeff_1D) const
  {
    double pn = psi * this->one_d_cub_dpsi_inv;
    int ip=floor(pn);
    ip=std::min(std::max(ip,0),this->ncoeff-1);
    double wp=pn-(double)(ip);

    int idx = ip*4;

    //double acoef[4];
    //acoef[0] = one_d_cub_acoef(ip).coeff[0];
    //acoef[1] = one_d_cub_acoef(ip).coeff[1];
    //acoef[2] = one_d_cub_acoef(ip).coeff[2];
    //acoef[3] = one_d_cub_acoef(ip).coeff[3];

    double acoef[4] = {coeff_1D.Get(idx+0),
                       coeff_1D.Get(idx+1),
                       coeff_1D.Get(idx+2),
                       coeff_1D.Get(idx+3)};

    vtkm::FloatDefault iVal = 0.0;
    if (ideriv==0)
      iVal = acoef[0]+(acoef[1]+(acoef[2]+acoef[3]*wp)*wp)*wp;
    else if (ideriv==1)
      iVal = (acoef[1]+(2.0*acoef[2]+3.0*acoef[3]*wp)*wp)*one_d_cub_dpsi_inv;

    return iVal;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename Coeff_1DType, typename Coeff_2DType>
  bool HighOrderEval(vtkm::Vec3f& ptRPZ,
                     const LocatorType& locator,
                     const CellSetType& cellSet,
                     const BFieldType& B_RZP,
                     const Coeff_1DType& coeff_1D,
                     const Coeff_2DType& coeff_2D,
                     vtkm::Vec3f& res) const
  {

    //res is R,P,Z
    bool val = this->HighOrderB(ptRPZ, locator, cellSet, B_RZP, coeff_1D, coeff_2D, res);
    res[1] /= ptRPZ[0];
    return val;
/*
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];


    int i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int j = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = rmin + (vtkm::FloatDefault)(i)*dr;
    vtkm::FloatDefault Zc = zmin + (vtkm::FloatDefault)(j)*dz;

    //vtkm::FloatDefault f00, f10, f01, f11, f20, f02;
    //eval_bicub_2(x,y, rc(i), zc(j), acoeff_all(j,i).coeff, f00,f10,f01,f11,f20,f02 );
    vtkm::Matrix<vtkm::FloatDefault, 4, 4> acoeff;
    vtkm::Id offset = i*this->nr + j, idx = 0;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        //acoeff[ii][jj] = coeff_2D.Get(offset+idx);
        acoeff[jj][ii] = coeff_2D.Get(offset+idx);
        idx++;
      }

    double psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;

    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);
    std::cout<<" psi= "<<psi<<std::endl;
    std::cout<<" dpsi_dr = "<<dpsi_dr<<std::endl;
    std::cout<<" dpsi_dz = "<<dpsi_dz<<std::endl;

    //double fld_I =      I_interp.interpol(psi,0, is_in_rgn3, equil.xpt_psi, bt_sign);
    //double fld_dIdpsi = I_interp.interpol(psi,1, is_in_rgn3, equil.xpt_psi, bt_sign);

    //All in region 1
    //double psi = min(max(min_psi,psi_in),max_psi);
    //double r8value = interpolation(psi, ideriv );
    //return bt_sign*r8value;

    vtkm::FloatDefault fld_I = this->I_interpol(psi, 0, coeff_1D);
    vtkm::FloatDefault fld_dIdpsi = this->I_interpol(psi, 1, coeff_1D);


    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault Br = -dpsi_dz * over_r;
    vtkm::FloatDefault Bz = dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r;

    res = vtkm::Vec3f(Br, Bz, Bp);

    if (1)
    {
      vtkm::Vec3f ptRZ(R, Z, 0);
      vtkm::Vec3f particlePos_param;
      vtkm::Vec<vtkm::Id,3> particlePos_vids;
      if (!this->PtLoc(ptRZ, locator, cellSet, particlePos_param, particlePos_vids))
        return false;
      auto B0_rzp = this->EvalV(B_RZP, 0, particlePos_param, particlePos_vids);
      B0_rzp[2] /= R;
      auto res0 = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

      std::cout<<"B_rzp("<<R<<" "<<Z<<" "<<Phi<<") = "<<res<<" orig= "<<res0<<std::endl;
    }

    return true;
*/
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  bool Evaluate(vtkm::Vec3f& ptRPZ,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const BFieldType& B_RZP,
                const BFieldType& B_Norm_RZP,
                const BFieldType& Curl_NB_RZP,
                const AsFieldType& AsPhiFF,
                const DAsFieldType& DAsPhiFF_RZP,
                const Coeff_1DType& Coeff_1D,
                const Coeff_2DType& Coeff_2D,
                vtkm::Vec3f& res) const
  {
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    if (this->UseHighOrder)
      return this->HighOrderEval(ptRPZ, locator, cellSet, B_RZP, Coeff_1D, Coeff_2D, res);

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

  //dopri
  vtkm::FloatDefault h_init = 0.0;
  vtkm::FloatDefault reltol = 1e-4;
  vtkm::FloatDefault abstol = 1e-6;

  vtkm::Id NumPlanes;
  vtkm::Id NumNodes;
  vtkm::FloatDefault dPhi;

  bool UseBOnly = false;
  bool UseHighOrder = false;
  bool SaveTraces = false;

  int nr, nz;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;

  int ncoeff;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault sml_bp_sign = -1.0f;
};
