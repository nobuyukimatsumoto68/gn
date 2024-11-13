#include <cmath>
#include <fstream>
#include <filesystem>

#include "global.h"
#include "rnd.h"
#include "force.h"
#include "gauge.h"
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

  const int Nc=3;
  Gauge W(Nc);

  std::vector<Obs<double, Gauge>*> obslist;
  std::string data_path="./obs_exc/";
  std::filesystem::create_directory(data_path);

  // ------------------

  // using Action = GaussianAction;
  // const double beta = 2.0;
  // Action S(beta);
  // Obs<double, Gauge> trW( "trW", [](const Gauge& W ){ return W().trace().real(); }, 0.0 );
  // obslist.push_back(&trW);
  // Obs<double, Gauge> trWsq( "trWsq", [](const Gauge& W ){ return W().squaredNorm(); }, 2.0*Nc*Nc*Nc / beta );
  // obslist.push_back(&trWsq);

  // ------------------
  
  // using Action = GaussianPhiAction;
  // const double lambda = 2.0;
  // Action S(lambda);

  // ------------------

  using Action = WilsonGaussianAction;
  const double beta = 2.5;
  double lambda = 1.0;
  if (argc>2){ lambda = atof(argv[2]); }
  Action S(beta, lambda);

  Obs<double, Gauge> phi_tr_re( "retrPhi", lambda, [](const Gauge& W ){
    return ( W.Phi - W.id() ).trace().real();
  }, 0.0 );
  obslist.push_back(&phi_tr_re);
  Obs<double, Gauge> phi_tr_im( "imtrPhi", lambda, [](const Gauge& W ){
    return ( W.Phi - W.id() ).trace().imag();
  }, 0.0 );
  obslist.push_back(&phi_tr_im);
  Obs<double, Gauge> phi_norm( "trPhisq", lambda, [](const Gauge& W ){
    return ( W.Phi - W.id() ).squaredNorm();
  }, 0.0 );
  obslist.push_back(&phi_norm);

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

    HMC<Force, Gauge, Action> hmc(S, stot, nsteps);
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
      
      std::ofstream of;
      of.open( data_path+pt->description+std::to_string(pt->param)+".dat", std::ios::out | std::ios::app);
      if(!of) assert(false);
      of << std::scientific << std::setprecision(15);
      for(auto iv=pt->v.begin(); iv!=pt->v.end(); ++iv) of << pt->param << "\t" << *iv << std::endl;
    }
  }


  return 0;
}
