#pragma once

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <iomanip>


using MC = Eigen::MatrixXcd;
using MR = Eigen::MatrixXd;
using VC = Eigen::VectorXcd;
using VR = Eigen::VectorXd;

using Complex = std::complex<double>;
constexpr Complex I = Complex(0.0, 1.0);

std::ostream& operator<<(std::ostream& os, const MR& W) {
  os << std::scientific << std::setprecision(15);
  for(int i=0; i<W.rows(); i++){
    for(int j=0; j<W.rows(); j++){
      os << std::setw(22) << W(i,j) << "\t";
    }
    os << std::endl;
  }
  return os;
}

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
