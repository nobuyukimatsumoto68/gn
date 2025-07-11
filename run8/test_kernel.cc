#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <omp.h>
constexpr int nparallel=3;
#define Nc 2
#define NA (Nc*Nc-1)
#define NG (2*Nc*Nc)
#define NH (Nc*Nc)
#define DIM 2

#include "lattice.h"
#include "rng.h"
#include "generators.h"
#include "force.h"
#include "gauge.h"
#include "kernel.h"
#include "action.h"
#include "integrator.h"
#include "hmc.h"
// #include "obs.h"

#include <unsupported/Eigen/MatrixFunctions>


int main( int argc, char *argv[] ){
  std::cout << std::scientific << std::setprecision(15);
  omp_set_dynamic(0);
  omp_set_num_threads(nparallel);
  Eigen::setNbThreads(1);
  std::cout << "# openmp threads = " << omp_get_num_threads() << std::endl;
  std::cout << "# eigen threads = " << Eigen::nbThreads() << std::endl;

  int seed = 0;
  if (argc>1){ seed = atoi(argv[1]); }

  // ------------------

  using Force = LinkForce;
  using Gauge = LinkConf;
  using Rng = SingleRng;

  using Action = WilsonGaussianAndDet2;

  using Kernel = TrivialKernel2;
  using Integrator = ImplicitLeapfrog<Force,Gauge,Action,Kernel>;
  using HMC = HMC<Force,Gauge,Integrator,Rng>;

  // ---------------

  // const Lattice lat( Lattice::Coord{{ 4, 4 }} );

  // ---------------

  Gauge W;

  Rng rng(1);
  W.randomize( rng, 1.0 );

  Kernel K;
  Force f = K.gen( W, rng );
  Force f0 = f;

  // ------------------

  double beta = 2.0; // 4.0
  const double lambda = 1.0;
  const double kappa = 5.0;
  Action S(beta, lambda, kappa);

  // ------------------



  // ------------------

  // const double alpha = 0.001;
  // Kernel K(lat);
  // Force f(lat);

  // Force Wf = K.act( W, f );

  // Laplacian::MN mat = K.matrix_free();
  // std::cout << "# herm: " << (mat-mat.adjoint()).norm() << std::endl;

  // Eigen::SelfAdjointEigenSolver<Laplacian::MN> ev( mat );
  // std::cout << "# ev = " << std::endl
  //           << ev.eigenvalues() << std::endl;

  // ------------------


  // std::cout << "f.pi = " << std::endl
  //           << f.pi.transpose() << std::endl;
  // std::cout << "f.pi0 = " << std::endl
  //           << f.pi0 << std::endl;
  // std::cout << "f.rho = " << std::endl
  //           << f.rho.transpose() << std::endl;
  // std::cout << "f.rho0 = " << std::endl
  //           << f.rho0 << std::endl;

  // Force::VG wbasis = f.wbasis( W.J() );
  // f.update_from( wbasis, W.J() );

  // std::cout << "check. diff = " << ( f-f0 ).norm() << std::endl;

  // Force::MC piM = f.get_pi();
  // f.update_pi( piM );
  // std::cout << "check. diff = " << ( f-f0 ).norm() << std::endl;

  using Complex = std::complex<double>;
  Complex I = Complex(0.0, 1.0);


  {

    // std::cout << "S = " << S(ell) << std::endl;
    Force dS = S.d( W );
    Force::VG check = dS.wbasis( W.J() );

    std::cout << "dS (check) = " << std::endl;
    for(int i=0; i<Nc; i++){
      for(int j=0; j<Nc; j++){
        const double eps = 1.0e-5;
        Gauge WP = W;
        Gauge WM = W;

        WP.W(i,j) += eps;
        WP.update_others();
        WM.W(i,j) -= eps;
        WM.update_others();

        std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
        std::cout << check(Nc*i + j) << std::endl;
      }
    }
    for(int i=0; i<Nc; i++){
      for(int j=0; j<Nc; j++){
        const double eps = 1.0e-5;
        Gauge WP = W;
        Gauge WM = W;

        WP.W(i,j) += I*eps;
        WP.update_others();
        WM.W(i,j) -= I*eps;
        WM.update_others();

        std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
        std::cout << check(Nc*Nc + Nc*i + j) << std::endl;
      }
    }
    std::cout << std::endl;
  }

  const double stot = 1.0;
  int nsteps=10;
  Integrator md(S, K, stot, nsteps);
  // HMC hmc(md, rng, stot, nsteps);

  { // dHdW
    Force p = K.gen( W, rng );
    const Force dHdW = md.dHdW( p, W );
    Force::VG check = dHdW.wbasis( W.J() );
    const double eps = 1.0e-5;

    // for(int qij=0; qij<2*Nc*Nc; qij++){
    //   std::cout << tmp[qij] << std::endl;

    //   Gauge Wp = W;
    //   Gauge Wm = W;
    //   Wp[qij] += eps;
    //   Wm[qij] -= eps;
    //   Wp.update_others();
    //   Wm.update_others();

    //   std::cout << ( md.H(p, Wp) - md.H(p, Wm) )/(2.0*eps) << std::endl;
    // }
    for(int i=0; i<Nc; i++){
      for(int j=0; j<Nc; j++){
        const double eps = 1.0e-5;
        Gauge WP = W;
        Gauge WM = W;

        WP.W(i,j) += eps;
        WP.update_others();
        WM.W(i,j) -= eps;
        WM.update_others();

        std::cout << ( md.H(p, WP) - md.H(p, WM) )/(2.0*eps) << std::endl;
        // std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
        std::cout << check(Nc*i + j) << std::endl;
      }
    }
    for(int i=0; i<Nc; i++){
      for(int j=0; j<Nc; j++){
        const double eps = 1.0e-5;
        Gauge WP = W;
        Gauge WM = W;

        WP.W(i,j) += I*eps;
        WP.update_others();
        WM.W(i,j) -= I*eps;
        WM.update_others();

        // std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
        std::cout << ( md.H(p, WP) - md.H(p, WM) )/(2.0*eps) << std::endl;
        std::cout << check(Nc*Nc + Nc*i + j) << std::endl;
      }
    }
    std::cout << std::endl;

  }

  // { // dHdp
  //   Force p = K.gen( W, rng );
  //   const Force dHdp = md.dHdp( p, W );
  //   Force::VG tmp = dHdp.wbasis( W.J() );
  //   const double eps = 1.0e-5;

  //   for(int i=0; i<2*Nc*Nc; i++){
  //     std::cout << tmp[i] << std::endl;

  //     // Force pp = p;
  //     // Force pm = p;
  //     // pp[i] += eps;
  //     // pm[i] -= eps;
  //     Force::VG tmp_p = tmp;
  //     Force::VG tmp_m = tmp;
  //     tmp_p[i] += eps;
  //     tmp_m[i] -= eps;
  //     Force pp( tmp_p, W.J() );
  //     Force pm( tmp_m, W.J() );

  //     std::cout << ( md.H(pp, W)-md.H(pm, W) )/(2.0*eps) << std::endl;
  //   }
  // }



  // {
  //   const double stot = 1.0;

  //   for(int nsteps=10; nsteps<400; nsteps*=2){
  //     Integrator md(S, K, stot, nsteps);
  //     HMC hmc(md, rng, stot, nsteps);

  //     rng.seed( 1 );
  //     W.randomize( rng, 1.0 );

  //     {
  //       Force p = K.gen( W, rng );

  //       const double Hinit = md.H(p,W);
  //       // std::cout << Hinit << std::endl;
  //       for(int i=0; i<md.nsteps; i++) md.onestep( p, W );
  //       const double Hfin = md.H(p,W);
  //       // std::cout << Hfin << std::endl;
  //       const double diff = Hfin-Hinit;
  //       std::cout << hmc.tau << "\t" << diff << std::endl;
  //     }

  //   }
  // }



  // ------------------


  return 0;
}
