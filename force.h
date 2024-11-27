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
  using MC = Eigen::Matrix<Complex, Nc, Nc, Eigen::RowMajor>;
  using VG = Eigen::Matrix<double, NG, 1>;
  using MG = Eigen::Matrix<double, NG, NG, Eigen::RowMajor>;

  VG pi;

  ForceSingleLink() { pi = VG::Zero(); }

  ForceSingleLink(const VG& pi_) { pi = pi_; }

  ForceSingleLink(const ForceSingleLink& other)
    : pi(other.pi)
  {}

  inline double norm() const { return pi.norm(); }

  Force& operator=(const Force& other) {
    if (this == &other) return *this;
    pi = other.pi;
    return *this;
  }

  // void Lmult(const MR& mat) {
  //   pi = mat*pi;
  // }
  void invert(const MG& mat) {
    const Eigen::PartialPivLU<MG> llt(mat);
    pi = llt.solve(this->pi);
  }

  inline double operator[](const int i) const { return pi[i]; }
  inline double& operator[](const int i) { return pi[i]; }
  inline int size() const { return pi.size(); }

  Force& operator+=(const Force& rhs){
    pi += rhs.pi;
    return *this;
  }
  friend Force operator+(Force v, const Force& w) { v += w; return v; }

  Force& operator-=(const Force& rhs){
    pi -= rhs.pi;
    return *this;
  }
  friend Force operator-(Force v, const Force& w) { v -= w; return v; }

  Force& operator*=(const double a){
    pi *= a;
    return *this;
  }

  friend Force operator*(const double a, Force v) { v.pi *= a; return v; }
  friend Force operator*(const MG& mat, Force v) { v.pi = mat*v.pi; return v; }

  friend std::ostream& operator<<(std::ostream& os, const Force& v) { os << v.pi; return os; }
};




struct Force2D {
  using Force=Force2D;
  using V=ForceSingleLink;

  using Idx = std::size_t;
  using Coord=std::array<int, DIM>;

  const Lattice& lattice;
  std::vector<V> field;

  Force2D( const Lattice& lattice )
    : lattice( lattice )
    , field( lattice.n_links() )
  {}

  inline V operator[](const Idx il) const { return field[il]; }
  inline V& operator[](const Idx il) { return field[il]; }

  V operator()(const Idx ix, const int mu) const {
    assert( mu>=0 );
    return field[ DIM*ix + mu];
  }
  V& operator()(const Idx ix, const int mu) {
    assert( mu>=0 );
    return field[ DIM*ix + mu];
  }

  V operator()(const Coord& x, const int mu) const {
    assert( mu>=0 );
    return field[ DIM*lattice.idx(x) + mu];
  }
  V& operator()(const Coord& x, const int mu) {
    assert( mu>=0 );
    return field[ DIM*lattice.idx(x) + mu];
  }

  Force& operator=(const Force& other) {
    if (this == &other) return *this;
    assert(&lattice==&other.lattice);
    field = other.field;
    return *this;
  }

  Force& operator+=(const Force& rhs){
    assert(&lattice==&rhs.lattice);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(Idx il=0; il<lattice.n_links(); il++) field[il] += rhs.field[il];
    return *this;
  }
  friend Force operator+(Force v, const Force& w) { v += w; return v; }

  Force& operator-=(const Force& rhs){
    assert(&lattice==&rhs.lattice);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(Idx il=0; il<lattice.n_links(); il++) field[il] -= rhs.field[il];
    return *this;
  }
  friend Force operator-(Force v, const Force& w) { v -= w; return v; }

  friend Force operator*(const double a, Force v) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(Idx il=0; il<v.lattice.n_links(); il++) v.field[il] *= a;
    return v;
  }

  double norm() const {
    double res = 0.0;
    // for(Idx il=0; il<lattice.n_links(); il++) res += field[il].norm();

    std::vector<double> tmp( lattice.n_links() );
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(Idx i=0; i<lattice.n_links(); i++) tmp[i] = field[i].norm();

    for(auto elem : tmp ) res += elem;
    res /= std::sqrt( lattice.n_links() );
    return res;
  }




};
