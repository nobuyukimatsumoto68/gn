#pragma OConce

#include <Eigen/Dense>

/*
  Force objects should have:

  operator+=(const Force& rhs);
  double square() const;
  void rand();
  Force& operator+=(const Force& rhs);
  friend Force operator*(const double a, Force v);
  friend Force operator-(Force v, const Force& w);
*/


struct ForceSingleLink{
  using Force = ForceSingleLink;

  using Complex = std::complex<double>;
  using MC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VC = Eigen::VectorXcd;
  using VR = Eigen::VectorXd;

  const int Nc;
  const int n;
  VR pi;

  ForceSingleLink(const int Nc_)
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

  // double square() const { return pi.squaredNorm(); }
  double norm() const { return pi.norm(); }

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

  Force& operator+=(const Force& rhs){
    assert( Nc==rhs.Nc );
    pi += rhs.pi;
    return *this;
  }
  friend Force operator+(Force v, const Force& w) { v += w; return v; }

  Force& operator-=(const Force& rhs){
    assert( Nc==rhs.Nc );
    pi -= rhs.pi;
    return *this;
  }
  friend Force operator-(Force v, const Force& w) { v -= w; return v; }

  friend Force operator*(const double a, Force v) { v.pi *= a; return v; }
  friend Force operator*(const MR& mat, Force v) { v.pi = mat*v.pi; return v; }

  friend std::ostream& operator<<(std::ostream& os, const Force& v) { os << v.pi; return os; }
};




struct Force2D {
  using Force=Force2D;

  const Lattice& lattice;
  std::vector<ForceSingleLink> field;

  Force2D( const Lattice& lattice, const int Nc )
    : lattice( lattice )
    , field( lattice.n_links(), Nc )
  {}

};
