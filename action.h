#pragma once

/*
  Action objects should have:

  double operator()( const Gauge& W ) const;
  Force d( const Gauge& W ) const;
*/


struct GaussianAction {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;


  const double beta;

  GaussianAction(const double beta_)
    : beta(beta_)
  {}

  double operator()( const Gauge& W ) const {
    double res = 0.0;
    res += W().squaredNorm();
    res *= 0.5*beta/W.Nc;
    return res;
  }

  MR dx( const Gauge& W ) const {
    MR res = W().real();
    res *= beta/W.Nc;
    return res;
  }

  MR dy( const Gauge& W ) const {
    MR res = W().imag();
    res *= beta/W.Nc;
    return res;
  }

  Force d( const Gauge& W ) const {
    const int Nc = W.Nc;

    MR m_dx = dx(W);
    MR m_dy = dy(W);

    VR res = VR::Zero(2*Nc*Nc);
    res.segment(0, Nc*Nc) = Eigen::Map<VR>( m_dx.data(), Nc*Nc );
    res.segment(Nc*Nc, Nc*Nc) = Eigen::Map<VR>( m_dy.data(), Nc*Nc );
    return Force(Nc, res);
  }

};


struct GaussianPhiAction {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;


  const double lambda;

  GaussianPhiAction(const double lambda_)
    : lambda(lambda_)
  {}

  double operator()( const Gauge& W ) const {
    double res = 0.0;
    res += ( W.Phi - W.id() ).squaredNorm();
    res *= 0.5*lambda/W.Nc;
    return res;
  }

  double dphia( const Gauge& W, const int a ) const {
    double res = ( W.Phi * W.t[a] ).trace().real();
    res *= lambda/W.Nc;
    return res;
  }

  double dphi0( const Gauge& W ) const {
    double res = ( W.Phi - W.id() ).trace().real();
    res *= lambda/W.Nc;
    return res;
  }

  Force d( const Gauge& W ) const {
    const int Nc = W.Nc;

    VR dSb = VR::Zero(2*Nc*Nc);
    for(int a=0; a<W.NA; a++) dSb(Nc*Nc+a) = dphia(W,a);
    dSb(2*Nc*Nc-1) = dphi0(W);

    VR res = W.J().inverse() * dSb;
    return Force(Nc, res);
  }

};



struct WilsonGaussianAction {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;


  const double beta;
  const double lambda;

  WilsonGaussianAction(const double beta_,
                       const double lambda_)
    : beta(beta_)
    , lambda(lambda_)
  {}

  double operator()( const Gauge& W ) const {
    double res = 0.0;
    res -= beta/W.Nc * ( W.U ).trace().real();
    res += 0.5*lambda/W.Nc * ( W.Phi - W.id() ).squaredNorm();
    return res;
  }

  double Da( const Gauge& W, const int a ) const {
    double res = beta/W.Nc * ( W.t[a]*W.U ).trace().imag();
    return res;
  }

  double dphia( const Gauge& W, const int a ) const {
    double res = ( W.Phi * W.t[a] ).trace().real();
    res *= lambda/W.Nc;
    return res;
  }

  double dphi0( const Gauge& W ) const {
    double res = ( W.Phi - W.id() ).trace().real();
    res *= lambda/W.Nc;
    return res;
  }

  Force d( const Gauge& W ) const {
    const int Nc = W.Nc;

    VR dSb = VR::Zero(2*Nc*Nc);
    for(int a=0; a<W.NA; a++) dSb(a) = Da(W,a);
    for(int a=0; a<W.NA; a++) dSb(Nc*Nc+a) = dphia(W,a);
    dSb(2*Nc*Nc-1) = dphi0(W);

    VR res = W.J().inverse() * dSb;
    return Force(Nc, res);
  }

};


struct WilsonGaussianAndDet {
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;


  const double beta;
  const double lambda;
  const double kappa;

  WilsonGaussianAndDet(const double beta_,
                       const double lambda_,
                       const double kappa_
                       )
    : beta(beta_)
    , lambda(lambda_)
    , kappa(kappa_)
  {}

  double operator()( const Gauge& W ) const {
    double res = 0.0;
    res -= beta/W.Nc * ( W.U ).trace().real();
    res += 0.5*lambda/W.Nc * ( W.Phi - W.id() ).squaredNorm();
    res -= kappa * std::log( W.Phi.determinant().real() );
    return res;
  }

  double Da( const Gauge& W, const int a ) const {
    double res = beta/W.Nc * ( W.t[a]*W.U ).trace().imag();
    return res;
  }

  double dphia( const Gauge& W, const int a ) const {
    double res = ( W.Phi * W.t[a] ).trace().real();
    res *= lambda/W.Nc;
    res -= kappa * (W.Phi.inverse()*W.t[a]).trace().real();
    return res;
  }

  double dphi0( const Gauge& W ) const {
    double res = ( W.Phi - W.id() ).trace().real();
    res *= lambda/W.Nc;
    res -= kappa * W.Phi.inverse().trace().real();
    return res;
  }

  Force d( const Gauge& W ) const {
    const int Nc = W.Nc;

    VR dSb = VR::Zero(2*Nc*Nc);
    for(int a=0; a<W.NA; a++) dSb(a) = Da(W,a);
    for(int a=0; a<W.NA; a++) dSb(Nc*Nc+a) = dphia(W,a);
    dSb(2*Nc*Nc-1) = dphi0(W);

    VR res = W.J().inverse() * dSb;
    return Force(Nc, res);
  }

};



struct WilsonGaussianAndDet2D {
  using Force = Force2D;
  using Gauge = Dim2Gauge;
  using M = LinkConfig;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;

  const Lattice& lattice;
  const double beta;
  const double lambda;
  const double kappa;

  WilsonGaussianAndDet2D(const Lattice& lattice,
			 const double beta_,
			 const double lambda_,
			 const double kappa_
			 )
    : lattice(lattice)
    , beta(beta_)
    , lambda(lambda_)
    , kappa(kappa_)
  {}

  double operator()( const Gauge& W ) const {
    double res = 0.0;
    for(Lattice::Idx i=0; i<W.lattice.vol; i++){
      const Lattice::Coord x = lattice.get_coord(i);
      const Lattice::Coord xp0 = lattice.cshift(x, 0);
      const Lattice::Coord xp0p1 = lattice.cshift(xp0, 1);
      const Lattice::Coord xp1 = lattice.cshift(x, 1);

      const MC plaq = W(x,0).U * W(xp0,1).U * W(xp1,0).U.adjoint() * W(x,1).U.adjoint();
      res -= beta/W.Nc * plaq.trace().real();

      res += 0.5*lambda/W.Nc * ( W(x,0).Phi - W(x,0).id() ).squaredNorm();
      res -= kappa * std::log( W(x,0).Phi.determinant().real() );

      res += 0.5*lambda/W.Nc * ( W(x,1).Phi - W(x,1).id() ).squaredNorm();
      res -= kappa * std::log( W(x,1).Phi.determinant().real() );
    }
    return res;
  }

  // double Da( const Gauge& W, const int a ) const {
  //   double res = beta/W.Nc * ( W.t[a]*W.U ).trace().imag();
  //   return res;
  // }

  // double dphia( const Gauge& W, const int a ) const {
  //   double res = ( W.Phi * W.t[a] ).trace().real();
  //   res *= lambda/W.Nc;
  //   res -= kappa * (W.Phi.inverse()*W.t[a]).trace().real();
  //   return res;
  // }

  // double dphi0( const Gauge& W ) const {
  //   double res = ( W.Phi - W.id() ).trace().real();
  //   res *= lambda/W.Nc;
  //   res -= kappa * W.Phi.inverse().trace().real();
  //   return res;
  // }

  // Force d( const Gauge& W ) const {
  //   const int Nc = W.Nc;

  //   VR dSb = VR::Zero(2*Nc*Nc);
  //   for(int a=0; a<W.NA; a++) dSb(a) = Da(W,a);
  //   for(int a=0; a<W.NA; a++) dSb(Nc*Nc+a) = dphia(W,a);
  //   dSb(2*Nc*Nc-1) = dphi0(W);

  //   VR res = W.J().inverse() * dSb;
  //   return Force(Nc, res);
  // }

};
