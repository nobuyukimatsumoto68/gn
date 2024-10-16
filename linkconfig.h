#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include "typedefs_and_globalfuncs.h"

struct LinkConfig {
  const int Nc;
  const int NA;
  std::vector<MC> t; // generators; tr(TaTb) = delta_{ab}

  MC W;
  double theta;
  MC U;
  MC Phi;
  int n;

  LinkConfig( const int Nc )
    : Nc(Nc)
    , NA(Nc*Nc-1)
    , W(MC::Identity(Nc,Nc))
    , theta(0.0)
    , U(MC::Identity(Nc,Nc))
    , Phi(MC::Identity(Nc,Nc))
    , n(0)
  {
    assert( check_consistency() );
    assert(Nc>=2);
    set_generators();
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

  inline Complex u() const { return std::exp(I*theta); }

  bool check_consistency( const double TOL=1.0e-13 ) const {
    const MC check = u()*Phi*U;
    const double norm = (check-W).norm()/(std::sqrt(2.0)*Nc);
    return norm<TOL;
  }

  void decomposition(){
    Eigen::JacobiSVD<MC> svd;
    svd.compute(W, Eigen::ComputeFullU | Eigen::ComputeFullV); // U S V^\dagger
    Phi = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixU().adjoint();
    const MC Omega = svd.matrixU() * svd.matrixV().adjoint();
    theta = (1.0/Nc) * std::arg( Omega.determinant() ) + 2.0*M_PI*n/Nc;
    U = 1.0/u() * Omega;
  }

  void update( const MC& Wnew ){
    W = Wnew;
    decomposition();
    assert( check_consistency() );
  }


  MR J() const { // (i,j) -> one dim index in row major format
    MR res = MR::Zero(2*Nc*Nc, 2*Nc*Nc);
    for(int a=0; a<NA; a++){
      const MC mat = ( I*u()*Phi*t[a]*U ).transpose(); // transpose to convert to row major
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(a, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(a, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    { // 0th element
      const MC mat = ( I*u()*Phi*U ).transpose(); // transpose to convert to row major
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(Nc*Nc-1, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(Nc*Nc-1, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    for(int a=0; a<NA; a++){
      const MC mat = ( u()*t[a]*U ).transpose(); // transpose to convert to row major
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(Nc*Nc+a, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(Nc*Nc+a, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    { // 0th element
      const MC mat = ( u()*U ).transpose(); // transpose to convert to row major
      MR Re = mat.real(); // to avoid bugs of Eigen; no const, no .real().data()
      MR Im = mat.imag(); // to avoid bugs of Eigen; no const, no .real().data()
      res.block(2*Nc*Nc-1, 0, 1, Nc*Nc) = Eigen::Map<VR>( Re.data(), Nc*Nc ).transpose();
      res.block(2*Nc*Nc-1, Nc*Nc, 1, Nc*Nc) = Eigen::Map<VR>( Im.data(), Nc*Nc ).transpose();
    }
    return res;
  }

};
