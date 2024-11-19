#pragma once

#include <cassert>
#include <cmath>
#include <vector>


double& re(std::complex<double>& c)
{ return reinterpret_cast<double (&)[2]>(c)[0]; }

double& im(std::complex<double>& c)
{ return reinterpret_cast<double (&)[2]>(c)[1]; }

double re(const std::complex<double>& c)
{
  return reinterpret_cast<const double (&)[2]>(c)[0];
}
double im(const std::complex<double>& c)
{
  return reinterpret_cast<const double (&)[2]>(c)[1];
}


/*
  Gauge objects should have:
  Gauge operator+=(const Force& rhs);
  friend Gauge operator-(Gauge v, const Gauge& w);
*/


struct LinkConfig { // Force=ForceSingleLink
  using Gauge=LinkConfig;
  using Force=ForceSingleLink;

  const int Nc;
  const int NA;
  std::vector<MC> t; // generators; tr(TaTb) = delta_{ab}

  MC W;
  double theta;
  MC U;
  MC Phi;

  LinkConfig( const int Nc )
    : Nc( Nc )
    , NA( Nc*Nc-1 )
    , W( id() )
    , theta( 0.0 )
    , U( id() )
    , Phi( id() )
  {
    check_consistency();
    assert(Nc>=2);
    set_generators();
  }

  LinkConfig( const LinkConfig& other )
    : Nc( other.Nc )
    , NA( Nc*Nc-1 )
    , t( other.t )
    , W( other.W )
    , theta( other.theta )
    , U( other.U )
    , Phi( other.Phi )
  {
    check_consistency();
    assert(Nc>=2);
  }

  LinkConfig& operator=(const LinkConfig& other) {
    if (this == &other) return *this;

    assert(Nc==other.Nc);
    W = other.W;
    theta = other.theta;
    U = other.U;
    Phi = other.Phi;
    check_consistency();

    return *this;
  }

  Gauge& operator+=(const Force& f) {
    assert(Nc==f.Nc);
    VR dwr = f.pi.segment(0, Nc*Nc);
    VR dwi = f.pi.segment(Nc*Nc, Nc*Nc);
    W += Eigen::Map<MR>( dwr.data(), Nc, Nc );
    W += I*Eigen::Map<MR>( dwi.data(), Nc, Nc );
    update_others();
    return *this;
  }
  friend Gauge operator+(Gauge v, const Force& w) { v += w; return v; }

  Gauge& operator-=(const Force& f) {
    assert(Nc==f.Nc);
    VR dwr = f.pi.segment(0, Nc*Nc);
    VR dwi = f.pi.segment(Nc*Nc, Nc*Nc);
    W -= Eigen::Map<MR>( dwr.data(), Nc, Nc );
    W -= I*Eigen::Map<MR>( dwi.data(), Nc, Nc );
    update_others();
    return *this;
  }
  friend Gauge operator-(Gauge v, const Force& w) { v -= w; return v; }

  Gauge& operator+=(const Gauge& rhs) {
    assert(Nc==rhs.Nc);
    W += rhs.W;
    update_others();
    return *this;
  }
  friend Gauge operator+(Gauge v, const Gauge& w) { v += w; return v; }

  Gauge& operator-=(const Gauge& rhs) {
    W -= rhs.W;
    update_others();
    return *this;
  }
  friend Gauge operator-(Gauge v, const Gauge& w) { v -= w; return v; }

  inline MC id() const { return MC::Identity(Nc,Nc); }
  inline Complex u1( const double alpha ) const { return std::exp(I*alpha); }
  inline Complex u() const { return u1(theta); }
  inline MC operator()() const { return W; }
  inline Complex operator()(const int i, const int j) const { return W(i,j); }
  inline Complex& operator()(const int i, const int j) { return W(i,j); }

  double norm() { return W.norm(); }

  void get_qij( int& q, int& i, int& j, const int qij ) const {
    if(qij<Nc*Nc){
      q=0;
      j=qij%Nc;
      i=(qij-j+1)/Nc;
    }
    else if(qij<2*Nc*Nc) {
      q=1;
      j=qij%Nc;
      i=(qij-j+1-Nc*Nc)/Nc;
    }
    else assert( false );
  }

  double operator[](const int qij) const {
    int q,i,j;
    get_qij( q,i,j, qij );
    if(q==0) return re(W(i,j));
    else if(q==1) return im(W(i,j));
    else assert( false );
  }

  double& operator[](const int qij) {
    int q,i,j;
    get_qij( q,i,j, qij );
    if(q==0) return re(W(i,j));
    else if(q==1) return im(W(i,j));
    else assert( false );
  }

  double mod2pi( const double alpha ) const {
    double res = alpha + 4.0*M_PI;
    res -= int(std::floor(res/(2.0*M_PI)))*2.0*M_PI;
    return res;
  }

  void randomize( const std::function<double()>& f1,
		  const std::function<double()>& f2){
    for(int i=0; i<Nc; i++){
      for(int j=0; j<Nc; j++){
	W(i, j) = f1() + I*f2();
      }}
    update_others();
  }

  void set_generators(){
    for(int i=0; i<Nc; i++){
      for(int j=i+1; j<Nc; j++){
        {
          MC tmp = MC::Zero(Nc,Nc);
          tmp(i,j) = 1.0;
          tmp(j,i) = 1.0;
          t.push_back(tmp/std::sqrt(2.0));
        }
        {
          MC tmp = MC::Zero(Nc,Nc);
          tmp(i,j) = -I;
          tmp(j,i) =  I;
          t.push_back(tmp/std::sqrt(2.0));
        }
      }}

    for(int m=1; m<Nc; m++){
      MC tmp = MC::Zero(Nc,Nc);
      for(int i=0; i<Nc; i++){
        if(i<m) tmp(i,i) = 1.0;
        else if(i==m) tmp(i,i) = -m;
      }
      t.push_back( tmp/std::sqrt(m*(m+1.0)) );
    }
    // for( auto& elem : t ) elem *= 1.0/std::sqrt(2.0);
  }


  void check_consistency( const double TOL=1.0e-10 ) const {
    const MC check = u()*Phi*U;
    const double norm = (check-W).norm()/(std::sqrt(2.0)*Nc);
    if(norm > TOL) std::clog << "norm = " << norm << std::endl;
    assert( norm<TOL );
  }


  void decomposition(){
    // Eigen::JacobiSVD<MC> svd;
    Eigen::BDCSVD<MC> svd;
    svd.compute(W, Eigen::ComputeFullU | Eigen::ComputeFullV); // U S V^\dagger
    Phi = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixU().adjoint();
    MC Omega = svd.matrixU() * svd.matrixV().adjoint();
    Omega *= u1(-theta);
    const double dtheta = std::arg( Omega.determinant() ) / Nc;
    theta = mod2pi( theta + dtheta );
    U = u1(-dtheta) * Omega;
  }

  void update(){
    W = u()*Phi*U;
    check_consistency();
  }

  void update_from( const MC& Wnew ){
    W = Wnew;
    decomposition();
    check_consistency();
  }

  void update_others(){
    decomposition();
    check_consistency();
  }

  MR J() const {
    MR res = MR::Zero(2*Nc*Nc, 2*Nc*Nc);
    for(int a=0; a<NA; a++){
      const MC mat = ( I*u()*Phi*t[a]*U );
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(a, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(a, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    { // 0th element
      const MC mat = ( I*u()*Phi*U );
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(Nc*Nc-1, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(Nc*Nc-1, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    for(int a=0; a<NA; a++){
      const MC mat = ( u()*t[a]*U );
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(Nc*Nc+a, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(Nc*Nc+a, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    { // 0th element
      const MC mat = ( u()*U );
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(2*Nc*Nc-1, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(2*Nc*Nc-1, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    return res;
  }

  friend std::ostream& operator<<(std::ostream& os, const Gauge& v){ os << v.W; return os; }

};
