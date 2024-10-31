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






  int nsteps = 20;
 if (argc>1){ nsteps = atoi(argv[1]); }

  using Force = ForceSingleLink;
  using Gauge = LinkConfig;
  using Action = GaussianAction;

  const int Nc=2;
  const double beta = 2.0;
  Gauge W(Nc);
  Action S(beta);

  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );
  // std::cout << W << std::endl;

  {
    Force pi(Nc);
    pi.rand();

    const double stot = 1.0;

    // for(int nsteps = 10; nsteps <= 20; nsteps+=2 ){
    HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

    // std::cout << pi << std::endl;
    // std::cout << W << std::endl;

    // std::cout << S.d( W ) << std::endl;
    // std::cout << 0.5 * S.d( W ) << std::endl;
    // pi += 0.5 * S.d( W );

    const double Hinit = hmc.H(pi,W);
    // std::cout << Hinit << std::endl;
    hmc.leapfrog_explicit(pi, W);
    const double Hfin = hmc.H(pi,W);
    // std::cout << Hfin << std::endl;

    const double diff = Hfin-Hinit;
    std::cout << hmc.tau << "\t" << diff << std::endl;
    // }

    // std::cout << pi << std::endl;
    // std::cout << W << std::endl;
  }


  { // a, Phi
    const double eps = 1.0e-5;

    for(int a=0; a<W.NA; a++){
      std::cout << S.dphia(W, a) << std::endl;

      LinkConfig WP = W;
      LinkConfig WM = W;

      WP.Phi = W.Phi + eps*W.t[a];
      WP.update();

      WM.Phi = W.Phi - eps*W.t[a];
      WM.update();

      std::cout << ( S(WP)-S(WM) )/(2.0*eps) << std::endl;
    }
  }
  {
    std::cout << S.dphi0(W) << std::endl;

    const double eps = 1.0e-5;
    LinkConfig WP = W;
    LinkConfig WM = W;

    WP.Phi = W.Phi + eps * W.id();
    WP.update();

    WM.Phi = W.Phi - eps * W.id();
    WM.update();

    std::cout << ( S(WP)-S(WM) )/(2.0*eps) << std::endl;
  }




  std::cout << S.d(W) << std::endl;
  for(int i=0; i<W.Nc; i++){
    for(int j=0; j<W.Nc; j++){
      const double eps = 1.0e-5;
      LinkConfig WP = W;
      LinkConfig WM = W;

      WP(i,j) += eps;
      WP.update_others();
      WM(i,j) -= eps;
      WM.update_others();

      std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
    }
  }
  for(int i=0; i<W.Nc; i++){
    for(int j=0; j<W.Nc; j++){
      const double eps = 1.0e-5;
      LinkConfig WP = W;
      LinkConfig WM = W;

      WP(i,j) += I*eps;
      WP.update_others();
      WM(i,j) -= I*eps;
      WM.update_others();

      std::cout << ( S(WP)-S(WM) )/(2.0*eps) << " ";
    }
  }
  std::cout << std::endl;



  // std::cout << S(W) << std::endl;
  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );
  // std::cout << S(W) << std::endl;

  {
    Force pi(Nc);
    pi.rand();

    const double stot = 1.0;

    HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

    const double Hinit = hmc.H(pi,W);
    // std::cout << Hinit << std::endl;
    hmc.leapfrog_explicit(pi, W);
    const double Hfin = hmc.H(pi,W);
    // std::cout << Hfin << std::endl;

    const double diff = Hfin-Hinit;
    std::cout << hmc.tau << "\t" << diff << std::endl;
  }





  { // a, U
    const double eps = 1.0e-5;

    for(int a=0; a<W.NA; a++){
      std::cout << S.Da(W, a) << std::endl;

      LinkConfig WP = W;
      LinkConfig WM = W;

      WP.U = ( I*eps*W.t[a]).exp() * W.U;
      WP.update();

      WM.U = (-I*eps*W.t[a]).exp() * W.U;
      WM.update();

      std::cout << ( S(WP)-S(WM) )/(2.0*eps) << std::endl;
    }
  }



  std::cout << std::scientific << std::setprecision(15);

  int nsteps = 10;
  if (argc>1){ nsteps = atoi(argv[1]); }

  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  const int Nc=2;
  Gauge W(Nc);

  // using Action = GaussianAction;
  // const double beta = 2.0;
  // Action S(beta);

  // using Action = GaussianPhiAction;
  // const double lambda = 2.0;
  // Action S(lambda);

  using Action = WilsonGaussianAction;
  const double beta = 3.0;
  const double lambda = 2.0;
  Action S(beta, lambda);

  // ------------------


  // std::cout << S(W) << std::endl;
  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );
  // W.theta -= M_PI;
  // W.update_others();
  // std::cout << S(W) << std::endl;

  {
    Force pi(Nc);
    pi.rand();

    const double stot = 1.0;
    HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

    const double Hinit = hmc.H(pi,W);
    hmc.leapfrog_explicit(pi, W);
    const double Hfin = hmc.H(pi,W);

    const double diff = Hfin-Hinit;
    std::cout << hmc.tau << "\t" << diff << std::endl;
  }
