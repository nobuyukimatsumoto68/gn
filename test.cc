#include "global.h"
#include "force.h"
#include "gauge.h"
#include "action.h"
#include "hmc.h"

#include <unsupported/Eigen/MatrixFunctions>

int main(){

  const int Nc=2;

  using Force = ForceSingleLink;

  using Gauge = LinkConfig;
  Gauge W(Nc);

  using Action = GaussianAction;
  const double beta = 2.0;
  Action S(beta);

  // srand((unsigned int) time(0));
  W.update_W_from( MC::Random(Nc, Nc) );
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
