#include <TMB.hpp>
using namespace Eigen;

template <class Type>
Type objective_function <Type>:: operator () ()
{
  // data ( input from R)
  DATA_VECTOR (bs2_G);
  DATA_VECTOR (bs2_GT);
  DATA_VECTOR (os22_G);
  DATA_VECTOR (os22_GT);
  DATA_VECTOR (gamma_j);

  // parameters ( input from R)
  PARAMETER (tau2);
  PARAMETER (rho);

  int n = gamma_j.size();
  Type nll = 0.0;
  matrix<Type> V_inv(2,2);
  matrix<Type> V(2,2);

  for(int j=0; j<n; j++) {
    //Construct the Vj_inv matrix
    V_inv(0,0) = os22_G(j) + Type(1.0)/(tau2*(Type(1.0)-rho*rho));
    V_inv(0,1) = - rho/(tau2*(Type(1.0)-rho*rho));
    V_inv(1,0) = - rho/(tau2*(Type(1.0)-rho*rho));
    V_inv(1,1) = os22_GT(j) + Type(1.0)/(tau2*(Type(1.0)-rho*rho));

    Type det = V_inv(0,0)*V_inv(1,1)-V_inv(0,1)*V_inv(1,0);
    // try{
    V(0,0) = (Type(1.0)/det)*V_inv(1,1);
    V(1,0) = -(Type(1.0)/det)*V_inv(0,1);
    V(0,1) = -(Type(1.0)/det)*V_inv(1,0);
    V(1,1) = (Type(1.0)/det)*V_inv(0,0);
    // }
    // catch (const char* msg)
    // {
    //   printf ("Pseudo Inverse is used");
    //   V = V_inv.completeOrthogonalDecomposition().pseudoInverse();
    // };
    nll -= gamma_j(j)*(-Type(2.0)*log(tau2)-log(Type(1.0)-rho*rho)-log(det)+
      bs2_G(j)*bs2_G(j)*V(0,0)+2*bs2_G(j)*bs2_GT(j)*V(0,1)+bs2_GT(j)*bs2_GT(j)*V(1,1));
  }
  return nll;
}
