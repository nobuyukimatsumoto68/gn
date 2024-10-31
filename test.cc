#include "global.h"
#include "force.h"
#include "gauge.h"
#include "action.h"
#include "hmc.h"

#include <unsupported/Eigen/MatrixFunctions>

int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

  int nsteps = 10;
  if (argc>1){ nsteps = atoi(argv[1]); }

  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  const int Nc=2;
  Gauge W(Nc);

  // using Action = GaussianAction;
  // const double beta = 2.0;
  // Action S(beta);

  // using Action = GaussianPhiAction;
  // const double lambda = 2.0;
  // Action S(lambda);

  using Action = WilsonGaussianAction;
  const double beta = 3.0;
  const double lambda = 2.0;
  Action S(beta, lambda);

  // ------------------


  // std::cout << S(W) << std::endl;
  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );
  // W.theta -= M_PI;
  // W.update_others();
  // std::cout << S(W) << std::endl;

  {
    Force pi(Nc);
    pi.rand();

    const double stot = 1.0;
    HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

    const double Hinit = hmc.H(pi,W);
    hmc.leapfrog_explicit(pi, W);
    const double Hfin = hmc.H(pi,W);

    const double diff = Hfin-Hinit;
    std::cout << hmc.tau << "\t" << diff << std::endl;
  }


  return 0;
}
