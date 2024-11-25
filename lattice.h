#pragma once

#include <array>

struct Lattice {
  // using Complex = std::complex<double>;
  static constexpr int dim = 2;
  using Idx = std::size_t;
  using Coord=std::array<int, dim>;

  const Coord size;
  const int vol;

  // Lattice( const std::initializer_list<int> size )
  //   : size{size}
  // {};

  Lattice( const Coord& size )
    : size(size)
    , vol(size[0]*size[1])
  {}

  // Lattice( const Idx i )
  //   : x(dim)
  // {
  //   set_x(i);
  // }

  inline Idx n_sites() const { return dim; }
  inline Idx n_links() const { return vol*dim; }

  Coord get_coord( const Idx i ) const {
    Coord x;
    Idx tmp = i;
    for(int mu=dim-1; mu>=0; mu--){
      x[mu] = tmp%size[mu];
      tmp /= size[mu];
    }
    return x;
  }

  int operator[](const int mu) const { return size[mu]; }

  // std::string print() const {
  //   std::stringstream ss;
  //   for(int mu=0; mu<dim; mu++) ss << x[mu] << " ";
  //   return ss.str();
  // }

  Coord cshift( const Coord& x, const int mu ) const {
    Coord xn( x );
    const int sign = 2*(mu>=0)-1;
    const int mu0 = (sign==1) ? mu : -mu-1;
    xn[mu0] = (xn[mu0]+sign+size[mu0]) % size[mu0];
    return xn;
  }

  Idx idx(const Coord& x) const {
    Idx res = 0;
    for(int mu=0; mu<dim; mu++) res = res*size[mu] + x[mu];
    return res;
  }

};
