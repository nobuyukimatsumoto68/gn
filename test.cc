#include "typedefs_and_globalfuncs.h"
#include "linkconfig.h"

int main(){

  const int Nc=2;
  LinkConfig ell(Nc);
  // ell.update( MC::Random(Nc, Nc) );
  std::cout << ell.W << std::endl;

  // std::cout << std::exp(I*ell.theta) * ell.Phi*ell.U << std::endl;
  // std::cout << ell.U.determinant() << std::endl;
  
  std::cout << ell.J() << std::endl;
  std::cout << ell.J().determinant() << std::endl;

  return 0;
}
