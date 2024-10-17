  const int Nc=3;
  LinkVar W(Nc);

  for(auto m : W.t){
    std::cout << m << std::endl;
    std::cout << (m*m).trace() << std::endl;
  }


  ell.update( MC::Random(Nc, Nc) );
  std::cout << ell.W << std::endl;

  // std::cout << std::exp(I*ell.theta) * ell.Phi*ell.U << std::endl;
  // std::cout << ell.U.determinant() << std::endl;


  const MC PhiTaU = (ell.Phi * ell.t[0] * ell.U).transpose(); // for row major
  MR Re = PhiTaU.real();
  // std::cout << Re << std::endl;
  const MR Im = PhiTaU.imag();
  // std::cout << Re.data() << std::endl;
  // MR Re = ell.W.real();
  // MR Im = ell.W.imag();
  Eigen::VectorXd tmp = Eigen::Map<Eigen::VectorXd>( Re.data(), Nc*Nc );
  std::cout << tmp << std::endl;
