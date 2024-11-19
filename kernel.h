#pragma once

/*
  Kernel objects should have:

  double operator()( const Force& pi, const Gauge& W ) const; // G(U) * pi
  Force d( const Force& pi, const Gauge& W ) const;
  double det( const Gauge& W ) const; // det G(U)
  Force det_log_d( const Gauge& W ) const; // d [log det G(U)]

  Force gen( Rng& rng ) const ;
*/


struct TrivialKernel {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;
  using Rng = SingleRng;

  const int Nc;

  TrivialKernel( const int Nc ) : Nc(Nc) {};
  double operator()( const Force& f, const Gauge& W ) const { return f.pi.squaredNorm(); }
  Force d( const Force& f, const Gauge& W ) const { return Force(Nc); } // zero
  Force act( const Gauge& W, const Force& f ) const { return f; }
  double det( const Gauge& W ) const { return 1.0; }
  Force det_log_d( const Gauge& W ) const { return Force(Nc); } // zero

  Force gen( Rng& rng ) const {
    Force p(Nc);
    for(int i=0; i<p.size(); i++) p[i] = rng.gaussian();
    return p;
  }
};


struct IdpWHW {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;
  using Rng = SingleRng;

  const int Nc;
  const int N;
  const double alpha;

  IdpWHW( const int Nc, const double alpha )
    : Nc(Nc)
    , N(2*Nc*Nc)
    , alpha(alpha)
  {};

  MR wiwj( const Gauge& W ) const {
    MR res = MR::Zero(N,N);
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) res(i,j) += W[i]*W[j];
    return res;
  }

  MR wiwj_d( const Gauge& W, const int i ) const {
    MR res = MR::Zero(N,N);
    for(int j=0; j<N; j++) res(i,j) += W[j];
    for(int j=0; j<N; j++) res(j,i) += W[j];
    return res;
  }

  MR matrix( const Gauge& W ) const {
    const MR A = MR::Identity(N,N) + alpha*wiwj( W );
    return A.transpose() * A;
  }

  MR matrix_d( const Gauge& W, const int i ) const {
    MR res = MR::Zero(N,N);
    const MR A = MR::Identity(N,N) + alpha*wiwj( W );
    const MR dA = alpha*wiwj_d( W, i );
    return dA.transpose() * A + A.transpose() * dA;
  }

  double operator()( const Force& f, const Gauge& W ) const {
    return f.pi.dot( matrix(W)*f.pi );
  }

  Force d( const Force& f, const Gauge& W ) const {
    Force res(Nc);
    for(int i=0; i<N; i++) res[i] = f.pi.dot( matrix_d(W, i)*f.pi );
    return res;
  }

  Force act( const Gauge& W, const Force& f ) const { return matrix(W)*f; }

  double det( const Gauge& W ) const { return matrix(W).determinant(); }

  Force det_log_d( const Gauge& W ) const {
    Force res(Nc);
    const MR inv =  matrix(W).inverse();
    for(int i=0; i<N; i++) res[i] = ( inv * matrix_d(W, i) ).trace();
    return res;
  }

  Force gen( const Gauge& W, Rng& rng ) const {
    Force z(Nc);
    for(int i=0; i<N; i++) z[i] = rng.gaussian();
    const Force res = matrix( W ).inverse() * z;
    return res;
  }
};
