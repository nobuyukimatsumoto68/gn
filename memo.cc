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




  // std::cout << std::exp(I*ell.theta) * ell.Phi*ell.U << std::endl;
  // std::cout << ell.U.determinant() << std::endl;

  std::cout << ell.J() << std::endl;
  // std::cout << ell.J().determinant() << std::endl;

  { // a, U
    const double eps = 1.0e-5;

    for(int a=0; a<ell.NA; a++){
      LinkConfig ellP = ell;
      LinkConfig ellM = ell;

      ellP.U = ( I*eps*ell.t[a]).exp() * ell.U;
      ellP.update();

      ellM.U = (-I*eps*ell.t[a]).exp() * ell.U;
      ellM.update();

      MC dW = (ellP.W - ellM.W) / (2.0*eps);
      MR Re = dW.real();
      MR Im = dW.imag();
      VR re = Eigen::Map<VR>( Re.data(), Nc*Nc );
      VR im = Eigen::Map<VR>( Im.data(), Nc*Nc );
      std::cout << re.transpose() << " "
                << im.transpose() << std::endl;
    }
  }
  { // 0, U
    const double eps = 1.0e-5;
    LinkConfig ellP = ell;
    LinkConfig ellM = ell;

    ellP.theta = ell.theta + eps;
    ellP.update();

    ellM.theta = ell.theta - eps;
    ellM.update();

    MC dW = (ellP.W - ellM.W) / (2.0*eps);
    MR Re = dW.real();
    MR Im = dW.imag();
    VR re = Eigen::Map<VR>( Re.data(), Nc*Nc );
    VR im = Eigen::Map<VR>( Im.data(), Nc*Nc );
    std::cout << re.transpose() << " "
              << im.transpose() << std::endl;
  }
  { // a, Phi
    const double eps = 1.0e-5;

    for(int a=0; a<ell.NA; a++){
      LinkConfig ellP = ell;
      LinkConfig ellM = ell;

      ellP.Phi = ell.Phi + eps*ell.t[a];
      ellP.update();

      ellM.Phi = ell.Phi - eps*ell.t[a];
      ellM.update();

      MC dW = (ellP.W - ellM.W) / (2.0*eps);
      MR Re = dW.real();
      MR Im = dW.imag();
      VR re = Eigen::Map<VR>( Re.data(), Nc*Nc );
      VR im = Eigen::Map<VR>( Im.data(), Nc*Nc );
      std::cout << re.transpose() << " "
                << im.transpose() << std::endl;
    }
  }
  { // 0, Phi
    const double eps = 1.0e-5;
    LinkConfig ellP = ell;
    LinkConfig ellM = ell;

    ellP.Phi = ell.Phi + eps * ell.id();
    ellP.update();

    ellM.Phi = ell.Phi - eps * ell.id();
    ellM.update();

    MC dW = (ellP.W - ellM.W) / (2.0*eps) ;
    MR Re = dW.real();
    MR Im = dW.imag();
    VR re = Eigen::Map<VR>( Re.data(), Nc*Nc );
    VR im = Eigen::Map<VR>( Im.data(), Nc*Nc );
    std::cout << re.transpose() << " "
              << im.transpose() << std::endl;
  }




  {
    const double beta = 2.0;
    TrivialAction S(beta);

    std::cout << "S = " << S(ell) << std::endl;

    {
      const int i=0;
      const int j=0;
      std::cout << "dS = " << std::endl
		<< S.d(ell) << std::endl;
    }
    std::cout << "dS (check) = " << std::endl;
    for(int i=0; i<ell.Nc; i++){
      for(int j=0; j<ell.Nc; j++){
	const double eps = 1.0e-5;
	LinkConfig ellP = ell;
	LinkConfig ellM = ell;

	ellP.W(i,j) += eps;
	ellP.update_others();
	ellM.W(i,j) -= eps;
	ellM.update_others();

	std::cout << ( S(ellP)-S(ellM) )/(2.0*eps) << " ";
      }
    }
    for(int i=0; i<ell.Nc; i++){
      for(int j=0; j<ell.Nc; j++){
	const double eps = 1.0e-5;
	LinkConfig ellP = ell;
	LinkConfig ellM = ell;

	ellP.W(i,j) += I*eps;
	ellP.update_others();
	ellM.W(i,j) -= I*eps;
	ellM.update_others();

	std::cout << ( S(ellP)-S(ellM) )/(2.0*eps) << " ";
      }
    }
    std::cout << std::endl;
  }

