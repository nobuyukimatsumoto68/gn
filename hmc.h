#pragma once

template <class Force, class Gauge, class Action>
struct HMC {
  const Action S;
  const double stot;
  const int nsteps;
  const double tau;

  HMC(const Action& S_, const double stot_=1.0, const int nsteps_=10)
    : S(S_)
    , stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
  {}

  double H( const Force& pi, const Gauge& W ) {
    double res = 0.0;
    res += 0.5 * pi.square();
    res += S(W);
    return res;
  }

  void leapfrog_explicit( Force& pi, Gauge& W ) const {
    pi += -0.5*tau * S.d(W);
    W += tau * pi;
    pi += -0.5*tau * S.d(W);
  }

  // void leapfrog_explicit( Force& pi, Gauge& W ) const {
  //   for(int n=0; n<nsteps; n++){
  //     pi += -0.5*tau * S.d(W);
  //     W += tau * pi;
  //     pi += -0.5*tau * S.d(W);
  //   }
  // }

  // void run( Force2& phi0, const Gauge2& U0, double& r, double& dH, bool& is_accept, const bool no_reject = false ) {
  //   Force2 phi = phi0;
  //   Force2 p = gen();
  //   const double h0 = H(phi, U0, p);
  //   leapfrog( phi, U0, p );
  //   const double h1 = H(phi, U0, p);

  //   dH = h1-h0;
  //   r = std::min( 1.0, std::exp(-dH) );
  //   const double a = rnd_master.uniform();
  //   if( a < r || no_reject ){
  //     phi0 = phi;
  //     is_accept=true;
  //   }
  //   else is_accept=false;
  // }
  
};
