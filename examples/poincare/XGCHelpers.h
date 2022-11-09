#ifndef XGCHelpers_h
#define XGCHelpers_h

namespace XGCHelper
{

VTKM_EXEC
inline int GetIndex(const vtkm::FloatDefault& x,
             const int& nx,
             const vtkm::FloatDefault& xmin,
             const vtkm::FloatDefault& dx_inv)
{
  int idx = std::max(0, std::min(nx-1,
                                 int((x-xmin)*dx_inv)) );
  return idx;

}

template <typename Coeff_2DType>
VTKM_EXEC
void EvalBicub2(const vtkm::FloatDefault& x,
                const vtkm::FloatDefault& y,
                const vtkm::FloatDefault& xc,
                const vtkm::FloatDefault& yc,
                const vtkm::Id& offset,
                const Coeff_2DType& Coeff_2D,
                vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02)
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

}

#endif //XGCHelpers_h
