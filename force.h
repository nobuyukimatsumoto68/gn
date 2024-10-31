#pragma once

/*
  Force objects should have:

  operator+=(const Force& rhs);
  double square() const;
  void rand();
  Force& operator+=(const Force& rhs);
  friend Force operator*(const double a, Force v);
  // friend Force operator*(Force v, const double a);
*/


struct ForceSingleLink{
  const int Nc;
  const int n;
  VR pi;

  ForceSingleLink(const int Nc_, const int seed_=2)
    : Nc(Nc_)
    , n(2*Nc*Nc)
  {
    pi = VR::Zero(n);
  }

  ForceSingleLink(const VR& pi_)
    : Nc(pi_.size())
    , n(2*Nc*Nc)
  {
    pi = pi_;
  }

  double square() const { return pi.squaredNorm(); }

  void rand( const std::function<double()>& gauss ){
    // pi = VR::Random(2*Nc*Nc);
    for(int i=0; i<pi.size(); i++) pi(i) = gauss();
  } // @@@ make it multivariate Gaussian

  ForceSingleLink& operator+=(const ForceSingleLink& rhs){
    pi += rhs.pi;
    return *this;
  }

  friend ForceSingleLink operator*(ForceSingleLink v, const double a) { v.pi *= a; return v; }
  friend ForceSingleLink operator*(const double a, ForceSingleLink v) { v.pi *= a; return v; }
  friend std::ostream& operator<<(std::ostream& os, const ForceSingleLink& v) { os << v.pi; return os; }
};
