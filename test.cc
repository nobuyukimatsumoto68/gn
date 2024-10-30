#include "global.h"
#include "force.h"
#include "gauge.h"
#include "action.h"
#include "hmc.h"

#include <unsupported/Eigen/MatrixFunctions>

int main(){

  using Force = ForceSingleLink;
  using Gauge = LinkConfig;
  using Action = GaussianAction;

  const int Nc=2;
  const double beta = 2.0;
  Gauge W(Nc);
  Action S(beta);

  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );
  std::cout << W << std::endl;

  {
    Force pi(Nc);
    pi.rand();

    HMC<Force, Gauge, Action> hmc(S);

    std::cout << pi << std::endl;
    std::cout << W << std::endl;

    std::cout << S.d( W ) << std::endl;
    std::cout << 0.5 * S.d( W ) << std::endl;
    // pi += 0.5 * S.d( W );

    hmc.leapfrog_explicit(pi, W);

    // std::cout << pi << std::endl;
    // std::cout << W << std::endl;
  }


  return 0;
}
