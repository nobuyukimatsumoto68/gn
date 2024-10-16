#include <cassert>
#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

using MC = Eigen::MatrixXcd;
using MR = Eigen::MatrixXd;

using Complex = std::complex<double>;
constexpr Complex I = Complex(0.0, 1.0);


std::ostream& operator<<(std::ostream& os, const MC& W) {
  os << std::scientific << std::setprecision(15);
  for(int i=0; i<W.rows(); i++){
    for(int j=0; j<W.rows(); j++){
      os << std::setw(22) << W(i,j).real() << " "
	 << std::setw(22) << W(i,j).imag() << "\t";
    }
    os << std::endl;
  }
  return os;
}




struct LinkConfig {
  const int Nc;
  const int NA;
  std::vector<MC> t; // generators; tr(TaTb) = delta_{ab}

  MC W;
  double theta;
  MC U;
  MC Phi;

  LinkConfig() = delete;

  LinkConfig( const int Nc )
    : Nc(Nc)
    , NA(Nc*Nc-1)
    , W(MC::Identity(Nc,Nc))
    , theta(0.0)
    , U(MC::Identity(Nc,Nc))
    , Phi(MC::Identity(Nc,Nc))
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
  }

  bool check_consistency( const double TOL=1.0e-13 ) const {
    const MC check = std::exp(I*theta) * Phi*U;
    const double norm = (check-W).norm()/(std::sqrt(2.0)*Nc);
    return norm<TOL;
  }

  void decomposition( const int n=0 ){
    Eigen::JacobiSVD<MC> svd;
    svd.compute(W, Eigen::ComputeFullU | Eigen::ComputeFullV); // U S V^\dagger
    Phi = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixU().adjoint();
    const MC Omega = svd.matrixU() * svd.matrixV().adjoint();
    theta = (1.0/Nc) * std::arg( Omega.determinant() ) + 2.0*M_PI*n/Nc;
    U = std::exp( -I*theta ) * Omega;
  }

  void update( const MC& Wnew, const int n=0 ){
    W = Wnew;
    decomposition( n );
    assert( check_consistency() );
  }


};


int main(){

  const int Nc=3;
  LinkConfig ell(Nc);

  ell.update( MC::Random(Nc, Nc) );
  std::cout << ell.W << std::endl;

  // std::cout << std::exp(I*ell.theta) * ell.Phi*ell.U << std::endl;
  // std::cout << ell.U.determinant() << std::endl;


  return 0;
}
