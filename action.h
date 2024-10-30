#pragma once

/*
  Action objects should have:

  double operator()( const Gauge& W ) const;
  Force d( const Gauge& W ) const;
*/


struct GaussianAction { // Force = ForceSingleLink, Gauge = LinkConfig
  const double beta;

  GaussianAction(const double beta_)
    : beta(beta_)
  {}
  
  double operator()( const LinkConfig& ell ) const {
    double res = 0.0;
    res += ell.W.squaredNorm();
    res *= 0.5*beta;
    return res;
  }

  MR dx( const LinkConfig& ell ) const {
    MR res = ell.W.real();
    res *= beta;
    return res;
  }

  MR dy( const LinkConfig& ell ) const {
    MR res = ell.W.imag();
    res *= beta;
    return res;
  }

  ForceSingleLink d( const LinkConfig& W ) const {
    const int Nc = W.Nc;

    MR m_dx = dx(W);
    MR m_dy = dy(W);

    VR res = VR::Zero(2*Nc*Nc);
    res.segment(0, Nc*Nc) = Eigen::Map<VR>( m_dx.data(), Nc*Nc );
    res.segment(Nc*Nc, Nc*Nc) = Eigen::Map<VR>( m_dy.data(), Nc*Nc );
    return ForceSingleLink(res);
  }

};
