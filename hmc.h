#pragma once

#include <random>

template <class Force, class Gauge, class Action>
struct HMC {
  const Action S;
  const double stot;
  const int nsteps;
  const double tau;

  HMC(const Action& S_, const double stot_=1.0, const int nsteps_=10, const int seed_=1)
    : S(S_)
    , stot(stot_)
    , nsteps(nsteps_)
    , tau(stot/nsteps)
  {
  }

  double H( const Force& pi, const Gauge& W ) {
    double res = 0.0;
    res += 0.5 * pi.square();
    res += S(W);
    return res;
  }

  void leapfrog_explicit_singlestep( Force& pi, Gauge& W ) const {
    pi += -0.5*tau * S.d(W);
    W += tau * pi;
    pi += -0.5*tau * S.d(W);
  }

  void leapfrog_explicit( Force& pi, Gauge& W ) const {
    for(int n=0; n<nsteps; n++) leapfrog_explicit_singlestep(pi,W);
  }

  void run( Gauge& W0,
	    double& r,
	    double& dH,
	    bool& is_accept,
	    double (*pi_init)(),
	    double (*uniform)(),
	    const bool no_reject = false ) {
    Force pi(W0.Nc);
    pi.rand( pi_init );
    Gauge W( W0 );
    const double h0 = H(pi, W);
    leapfrog_explicit( pi, W );
    const double h1 = H(pi, W);

    dH = h1-h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = uniform();
    if( a < r || no_reject ){
      W0 = W;
      is_accept=true;
    }
    else is_accept=false;
  }

};
