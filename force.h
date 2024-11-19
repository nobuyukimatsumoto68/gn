#pragma OConce

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
  using Force = ForceSingleLink;

  const int Nc;
  const int n;
  VR pi;

  ForceSingleLink(const int Nc_, const int seed_=2)
    : Nc(Nc_)
    , n(2*Nc*Nc)
  {
    pi = VR::Zero(n);
  }

  ForceSingleLink(const int Nc, const VR& pi_)
    : Nc(Nc)
    , n(2*Nc*Nc)
  {
    pi = pi_;
  }

  ForceSingleLink(const ForceSingleLink& other)
    : Nc(other.Nc)
    , n(other.n)
    , pi(other.pi)
  {
    assert(this != &other);
  }

  double square() const { return pi.squaredNorm(); }
  double norm() const { return pi.norm(); }

  // void rand( const Kernel& kernel, const Rng& rng ){
  //   // pi = VR::Random(2*Nc*Nc);
  //   for(int i=0; i<pi.size(); i++) pi(i) = gauss();
  // } // @@@ make it multivariate Gaussian

  Force& operator=(const Force& other) {
    if (this == &other) return *this;

    assert(Nc==other.Nc);
    assert(n==other.n);
    pi = other.pi;

    return *this;
  }

  double operator[](const int i) const { return pi[i]; }
  double& operator[](const int i) { return pi[i]; }
  int size() const { return pi.size(); }

  ForceSingleLink& operator+=(const ForceSingleLink& rhs){
    pi += rhs.pi;
    return *this;
  }
  ForceSingleLink& operator-=(const ForceSingleLink& rhs){
    pi -= rhs.pi;
    return *this;
  }

  // friend ForceSingleLink operator*(ForceSingleLink v, const double a) { v.pi *= a; return v; }
  friend ForceSingleLink operator*(const double a, ForceSingleLink v) { v.pi *= a; return v; }
  friend ForceSingleLink operator*(const MR& mat, ForceSingleLink v) { v.pi = mat*v.pi; return v; }
  friend ForceSingleLink operator-(ForceSingleLink v, const ForceSingleLink& w) {
    assert( v.Nc==w.Nc );
    v.pi -= w.pi;
    return v; }
  friend std::ostream& operator<<(std::ostream& os, const ForceSingleLink& v) { os << v.pi; return os; }
};
