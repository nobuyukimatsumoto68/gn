#pragma once

#include <array>

struct Lattice {
  using Idx = std::size_t;
  using Coord=std::array<int, DIM>;

  const Coord size;
  const int vol;

  Lattice( const Coord& size )
    : size(size)
    , vol(size[0]*size[1])
  {}

  inline Idx n_sites() const { return vol; }
  inline Idx n_links() const { return vol*DIM; }

  Coord get_coord( const Idx i ) const {
    Coord x;
    Idx tmp = i;
    for(int mu=DIM-1; mu>=0; mu--){
      x[mu] = tmp%size[mu];
      tmp /= size[mu];
    }
    return x;
  }

  inline int operator[](const int mu) const { return size[mu]; }

  Coord cshift( const Coord& x, const int mu ) const {
    Coord xn( x );
    const int sign = 2*(mu>=0)-1;
    const int mu0 = (sign==1) ? mu : -mu-1;
    xn[mu0] = (xn[mu0]+sign+size[mu0]) % size[mu0];
    return xn;
  }

  inline Idx cshift( const Idx& ix, const int mu ) const {
    return idx(cshift( get_coord(ix), mu ));
  }


  Idx idx(const Coord& x) const {
    Idx res = 0;
    for(int mu=0; mu<DIM; mu++) res = res*size[mu] + x[mu];
    return res;
  }

};
