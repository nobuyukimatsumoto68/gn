#include <cmath>
#include <fstream>
#include <filesystem>

#include "global.h"
#include "rnd.h"
#include "force.h"
#include "gauge.h"
#include "kernel.h"
#include "action.h"
#include "hmc.h"
#include "obs.h"

#include <unsupported/Eigen/MatrixFunctions>



int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

  int seed = 0;
  if (argc>1){ seed = atoi(argv[1]); }

  // ------------------

  using Force = ForceSingleLink;
  using Gauge = LinkConfig;
  using Action = WilsonGaussianAction;
  using Kernel = TrivialKernel;

  const int Nc=2;
  Gauge W(Nc);

  std::vector<Obs<double, Gauge>*> obslist;
  std::string data_path="./obs/";
  std::filesystem::create_directory(data_path);

  // ------------------

  double beta = 2.3;
  if (argc>2){ beta = atof(argv[2]); }
  const double lambda = 1.0;
  Action S(beta, lambda);

  // ------------------

  TrivialKernel K;

  // ------------------

  const double lambda_0 = 2.0 * std::cyl_bessel_i( 1, beta ) / beta;
  const double lambda_F = 2.0 * std::cyl_bessel_i( 2, beta ) / beta;
  double exact = 0.0;
  if( std::abs(lambda_0)>1.0e-14 ) exact = 2.0*lambda_F/lambda_0;

  Obs<double, Gauge> retrU( "retrU", beta, [](const Gauge& W ){ return W.U.trace().real(); }, exact );
  obslist.push_back(&retrU);

  // ------------------

  mt.seed( seed );
  W.randomize( gaussian );

  // ------------------


  {
    const double stot = 1.0;
    const int nsteps = 2;

    int ntherm=1000;
    int niter=10000;
    if (argc>3){ ntherm = atoi(argv[3]); }
    if (argc>4){ niter = atoi(argv[4]); }
    const int interval=5;

    HMC<Force, Gauge, Action, Kernel> hmc(S, K, stot, nsteps);
    double dH, r;
    bool is_accept;

    for(int n=0; n<ntherm; n++){
      hmc.run(W, r, dH, is_accept, &gaussian, &uniform);
    }

    for(int n=0; n<niter; n++){
      hmc.run(W, r, dH, is_accept, &gaussian, &uniform);
      std::clog << "n = " << n
		<< ", r = " << r
		<< ", dH = " << dH
		<< ", is_accept = " << is_accept
		<< std::endl;

      if(n%interval==0){
	for(auto pt : obslist) pt->meas( W );
      }
    }

    std::cout << "# beta \t mean \t err \t exact" << std::endl;
    for(auto pt : obslist){
      double mn, er;
      const int binsize = 20;
      pt->jackknife( mn, er, binsize );
      std::cout << pt->description << "\t\t"
		<< pt->param << "\t"
		<< mn << "\t" << er << "\t" << pt->exact <<  std::endl;

      // ------------------------
      
      // std::ofstream of;
      // of.open( data_path+pt->description+".dat", std::ios::out | std::ios::app);
      // if(!of) assert(false);
      // of << std::scientific << std::setprecision(15);
      // of << pt->param << "\t"
      // 	 << mn << "\t"
      // 	 << er << std::endl;
    }
  }


  return 0;
}
