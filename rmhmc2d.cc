#include <iostream>
#include <iomanip>
// #include <cmath>
// #include <fstream>
#include <filesystem>

// #include "global.h"
#include "lattice.h"
#include "rng.h"
#include "force.h"
#include "gauge.h"
#include "kernel.h"
#include "action.h"
#include "integrator.h"
#include "hmc.h"
#include "obs.h"

#include <unsupported/Eigen/MatrixFunctions>



int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);

  int seed = 0;
  if (argc>1){ seed = atoi(argv[1]); }

  // ------------------

  using Force = Force2D;
  using Gauge = Dim2Gauge;
  using Action = WilsonGaussianAndDet2D;

  // using Kernel = IdpWHW;
  // using Integrator = ImplicitLeapfrog<Force,Gauge,Action,Kernel>;
  using Rng = ParallelRng;
  // using HMC = HMC<Force,Gauge,Integrator,Rng>;

  // ---------------

  const Lattice lat( Lattice::Coord{{ 4,4 }} );

  // ---------------

  const int Nc=2;
  Gauge W(lat, Nc);

  std::vector<Obs<double, Gauge>*> obslist;
  std::string data_path="./obs/";
  std::filesystem::create_directory(data_path);

  // ------------------

  double beta = 3.3;
  if (argc>2){ beta = atof(argv[2]); }
  const double lambda = 2.0;
  const double kappa = 4.0;
  Action S(lat, beta, lambda, kappa);

  // ------------------
  
  // const double alpha = 0.001;
  // Kernel K(Nc, alpha);

  // // ------------------

  Rng rng(lat);
  W.randomize( rng );

  std::cout << "S = " << S(W) << std::endl;
  
  



  // // ------------------

  // const double lambda_0 = 2.0 * std::cyl_bessel_i( 1, beta ) / beta;
  // const double lambda_F = 2.0 * std::cyl_bessel_i( 2, beta ) / beta;
  // double exact = 0.0;
  // if( std::abs(lambda_0)>1.0e-14 ) exact = 2.0*lambda_F/lambda_0;
  // std::cout << "exact = " << exact << std::endl;

  // Obs<double, Gauge> retrU( "retrU", beta, [](const Gauge& W ){ return W.U.trace().real(); }, exact );
  // obslist.push_back(&retrU);
  // Obs<double, Gauge> phi_norm( "trPhisq", beta, [](const Gauge& W ){
  //   return ( W.Phi - W.id() ).squaredNorm();
  // }, 0.0 );
  // obslist.push_back(&phi_norm);
  // Obs<double, Gauge> phi_det_abs( "absdetPhi", beta, [](const Gauge& W ){
  //   return std::abs(W.Phi.determinant());
  // }, 0.0 );
  // obslist.push_back(&phi_det_abs);
  // Obs<double, Gauge> detK( "absdetK", beta, [&](const Gauge& W ){
  //   return std::abs(K.det(W));
  // }, 0.0 );
  // obslist.push_back(&detK);


  // // ------------------

  // rng.seed( seed );
  // W.randomize( [&](){ return rng.gaussian(); },
  // 	       [&](){ return rng.gaussian(); }
  // 	       );


  // // ------------------


  // const double stot = 1.0;
  // const int nsteps = 10;
  // Integrator md(S, K, stot, nsteps);
  // HMC hmc(md, rng, stot, nsteps);

  // {
  //   int ntherm=1000;
  //   int niter=10000;
  //   if (argc>3){ ntherm = atoi(argv[3]); }
  //   if (argc>4){ niter = atoi(argv[4]); }
  //   const int interval = 10;
  //   const int binsize = 10;

  //   double dH, r;
  //   bool is_accept;

  //   for(int n=0; n<ntherm; n++) hmc.run(W, r, dH, is_accept);

  //   for(int n=0; n<niter; n++){
  //     hmc.run(W, r, dH, is_accept);
  //     std::clog << "n = " << n
  //       	<< ", r = " << r
  //       	<< ", dH = " << dH
  //       	<< ", is_accept = " << is_accept
  //       	<< std::endl;

  //     if(n%interval==0){
  //       for(auto pt : obslist) {
  //         pt->meas( W );
  //         std::cout << pt->description << "\t"
  //                   << *(pt->v.end()-1) << std::endl;
  //       }
  //     }
  //   }

  //   std::cout << "# beta \t mean \t err \t exact" << std::endl;
  //   for(auto pt : obslist){
  //     double mn, er;
  //     pt->jackknife( mn, er, binsize );
  //     std::cout << pt->description << "\t\t"
  //       	<< pt->param << "\t"
  //       	<< mn << "\t" << er << "\t" << pt->exact <<  std::endl;

  //     // ------------------------
      
  //     {
  // 	std::ofstream of;
  // 	of.open( data_path+pt->description+".dat", std::ios::out | std::ios::app);
  // 	if(!of) assert(false);
  // 	of << std::scientific << std::setprecision(15);
  // 	of << pt->param << "\t"
  // 	   << mn << "\t"
  // 	   << er << "\t"
  // 	   << pt->exact << std::endl;
  //     }

  //     // ------------------------

  //     {
  // 	std::ofstream of;
  // 	of.open( data_path+pt->description+std::to_string(pt->param)+"_.dat", std::ios::out | std::ios::app);
  // 	if(!of) assert(false);
  // 	of << std::scientific << std::setprecision(15);
  // 	for(auto iv=pt->v.begin(); iv!=pt->v.end(); ++iv) of << pt->param << "\t" << *iv << std::endl;
  //     }
  //   }
  // }


  return 0;
}
