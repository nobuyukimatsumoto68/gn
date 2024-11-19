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


  // srand((unsigned int) time(0));
  W.update_from( MC::Random(Nc, Nc) );

  { // scaling test
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

  // { // branch test
  //   Force pi(Nc);
  //   pi.rand();

  //   const double stot = 400.0;
  //   nsteps = 1000;
  //   HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

  //   for(int n=0; n<nsteps; n++){
  //     hmc.leapfrog_explicit_singlestep(pi, W);
  //     std::cout << W.theta << std::endl;
  //   }
  // }




  // { // scaling test
  //   Force pi(Nc);
  //   pi.rand( gaussian );

  //   const double stot = 1.0;
  //   HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

  //   const double Hinit = hmc.H(pi,W);
  //   hmc.leapfrog_explicit(pi, W);
  //   const double Hfin = hmc.H(pi,W);

  //   const double diff = Hfin-Hinit;
  //   std::cout << hmc.tau << "\t" << diff << std::endl;
  // }

  // { // branch test
  //   Force pi(Nc);
  //   pi.rand( gaussian );

  //   const double stot = 400.0;
  //   nsteps = 1000;
  //   HMC<Force, Gauge, Action> hmc(S, stot, nsteps);

  //   for(int n=0; n<nsteps; n++){
  //     hmc.leapfrog_explicit_singlestep(pi, W);
  //     std::cout << W.theta << std::endl;
  //   }
  // }







  // { // dHdp
  //   Force p = K.gen( rng );
  //   const Force dHdp = hmc.dHdp( p, W );
  //   const double eps = 1.0e-5;

  //   for(int i=0; i<2*Nc*Nc; i++){
  //     std::cout << dHdp[i] << std::endl;

  //     Force pp = p;
  //     Force pm = p;
  //     pp[i] += eps;
  //     pm[i] -= eps;

  //     std::cout << ( hmc.H(pp, W)-hmc.H(pm, W) )/(2.0*eps) << std::endl;
  //   }
  // }

  { // dHdW
    Force p = K.gen( rng );
    const Force dHdW = hmc.dHdW( p, W );
    const double eps = 1.0e-5;

    const Force dS = S.d(W);

    for(int qij=0; qij<2*Nc*Nc; qij++){
      std::cout << dHdW[qij] << std::endl;

      // int q,i,j;
      // W.get_qij( q,i,j, qij );
      // std::cout << "q = " << q << ", i = " << i << ", j = " << j << std::endl;

      Gauge Wp = W;
      Gauge Wm = W;
      Wp[qij] += eps;
      Wm[qij] -= eps;
      Wp.update_others();
      Wm.update_others();

      // std::cout << "Wp-Wm = " << std::endl;
      // std::cout << Wp.W-Wm.W << std::endl;
      // std::cout << "Up-Um = " << std::endl;
      // std::cout << Wp.U-Wm.U << std::endl;
      // std::cout << "Phip-Phim = " << std::endl;
      // std::cout << Wp.Phi-Wm.Phi << std::endl;
      // std::cout << hmc.H(p, Wp) << ", " << hmc.H(p, Wm) << std::endl;
      // std::cout << S(Wp) << ", " << S(Wm) << std::endl;
      // std::cout << "dS =" << dS[qij] << std::endl;
      // std::cout << "dSc=" << (S(Wp) - S(Wm))/(2.0*eps) << std::endl;
      std::cout << ( hmc.H(p, Wp)-hmc.H(p, Wm) )/(2.0*eps) << std::endl;
      // std::cout << ( hmc.H(p, Wp)-hmc.H(p, Wm) )/(2.0*eps) << std::endl;
    }
    // for(int ij=0; ij<Nc*Nc; ij++){ // re
    //   std::cout << dHdW[ij] << std::endl;

    //   int i,j;
    //   W.get_ij( i,j, ij );

    //   Gauge Wp = W;
    //   Gauge Wm = W;
    //   Wp(i,j) += eps;
    //   Wm(i,j) -= eps;

    //   std::cout << "Wp-Wm = " << std::endl;
    //   std::cout << Wp.W-Wm.W << std::endl;

    //   std::cout << ( hmc.H(p, Wp)-hmc.H(p, Wm) )/(2.0*eps) << std::endl;
    // }
    // for(int ij=0; ij<Nc*Nc; ij++){ // im
    //   std::cout << dHdW[Nc*Nc+ij] << std::endl;

    //   int i,j;
    //   W.get_ij( i,j, ij );

    //   Gauge Wp = W;
    //   Gauge Wm = W;
    //   Wp(i,j) += I*eps;
    //   Wm(i,j) -= I*eps;

    //   std::cout << ( hmc.H(p, Wp)-hmc.H(p, Wm) )/(2.0*eps) << std::endl;
    // }
  }



  {
    const double eps = 1.0e-5;

    for(int qij=0; qij<2*Nc*Nc; qij++){
      MR dk = K.wiwj_d(W, qij);

      Gauge Wp = W;
      Gauge Wm = W;
      Wp[qij] += eps;
      Wm[qij] -= eps;
      Wp.update_others();
      Wm.update_others();

      MR diff = K.wiwj(Wp) - K.wiwj(Wm);
      MR check = ( diff )/(2.0*eps);
      std::cout << "norm = " << (dk-check).norm() << std::endl;
    }

    for(int qij=0; qij<2*Nc*Nc; qij++){
      MR dk = K.matrix_d(W, qij);

      Gauge Wp = W;
      Gauge Wm = W;
      Wp[qij] += eps;
      Wm[qij] -= eps;
      Wp.update_others();
      Wm.update_others();

      MR diff = K.matrix(Wp) - K.matrix(Wm);
      MR check = ( diff )/(2.0*eps);
      std::cout << "norm = " << (dk-check).norm() << std::endl;
    }

    {
      const Force f = K.gen(W, rng);
      Force dk = K.d(f, W);

      for(int qij=0; qij<2*Nc*Nc; qij++){
	Gauge Wp = W;
	Gauge Wm = W;
	Wp[qij] += eps;
	Wm[qij] -= eps;
	Wp.update_others();
	Wm.update_others();
	double Kp = K(f, Wp);
	double Km = K(f, Wm);
	std::cout << "norm = " << dk[qij] - (Kp-Km)/(2.0*eps) << std::endl;
      }
    }

    {
      Force dk = K.det_log_d(W);
      double det = K.det(W);

      {
	for(int qij=0; qij<2*Nc*Nc; qij++){
	  Gauge Wp = W;
	  Gauge Wm = W;
	  Wp[qij] += eps;
	  Wm[qij] -= eps;
	  Wp.update_others();
	  Wm.update_others();
	  double detp = K.det(Wp);
	  double detm = K.det(Wm);
	  std::cout << "norm = " << dk[qij] - (detp-detm)/(2.0*eps)/det << std::endl;
	}
      }
    }
  }




  MR corr = MR::Zero(K.N, K.N);
  const int nmax=1000000;
  for(int n=0; n<nmax; n++){
    Force p = K.gen( W, rng );
    for(int i=0; i<K.N; i++) for(int j=0; j<K.N; j++) corr(i,j) += p[i]*p[j];
  }
  corr /= nmax;

  std::cout << "corr = " << std::endl
	    << corr << std::endl << std::endl;

  // std::cout << "K(W) = " << std::endl
  // 	    << K.matrix(W) << std::endl << std::endl;
  std::cout << "K(W)^-1 = " << std::endl
	    << K.matrix(W).inverse() << std::endl << std::endl;

  std::cout << "diff = " << std::endl
	    << K.matrix(W).inverse() - corr << std::endl << std::endl;





  // { // dHdp
  //   Force p = K.gen( W, rng );
  //   const Force dHdp = hmc.dHdp( p, W );
  //   const double eps = 1.0e-5;

  //   for(int i=0; i<2*Nc*Nc; i++){
  //     std::cout << dHdp[i] << std::endl;

  //     Force pp = p;
  //     Force pm = p;
  //     pp[i] += eps;
  //     pm[i] -= eps;

  //     std::cout << ( hmc.H(pp, W)-hmc.H(pm, W) )/(2.0*eps) << std::endl;
  //   }
  // }

  // { // dHdW
  //   Force p = K.gen( W, rng );
  //   const Force dHdW = md.dHdW( p, W );
  //   const double eps = 1.0e-5;

  //   const Force dS = S.d(W);

  //   for(int qij=0; qij<2*Nc*Nc; qij++){
  //     std::cout << dHdW[qij] << std::endl;

  //     Gauge Wp = W;
  //     Gauge Wm = W;
  //     Wp[qij] += eps;
  //     Wm[qij] -= eps;
  //     Wp.update_others();
  //     Wm.update_others();

  //     std::cout << ( md.H(p, Wp) - md.H(p, Wm) )/(2.0*eps) << std::endl;
  //   }

  //   std::cout << "dp = "
  // 	      << (-0.5*md.tau * md.dHdW(p, W)).pi.transpose() << std::endl;
  //   std::cout << "p = "
  // 	      << p.pi.transpose() << std::endl;
  //   p = p - 0.5*md.tau * md.dHdW(p, W);
  //   std::cout << "p = "
  // 	      << p.pi.transpose() << std::endl;
  // }

  // W.check_consistency();


    {
      Force p = K.gen( W, rng );

      const double Hinit = md.H(p,W);
      std::cout << Hinit << std::endl;

      // const double tau = 0.1;

      // std::cout << "step1" << std::endl;
      // std::cout << "p = " << std::endl
      // 	      << p.pi.transpose() << std::endl;
      // p = md.get_phalf( p, W );
      // std::cout << "p = " << std::endl
      // 	      << p.pi.transpose() << std::endl;
      // std::cout << "step2" << std::endl;
      // W += 0.5 * md.tau * md.dHdp( p, W );
      // std::cout << "step3" << std::endl;
      // W = md.get_Wp( p, W );
      // std::cout << "step4" << std::endl;
      // p += -0.5*md.tau * md.dHdW(p, W);
      // std::cout << "step5" << std::endl;

      for(int i=0; i<md.nsteps; i++) md.onestep( p, W );
      // p += -0.5*tau * hmc.dHdW(p, W);
      // W += tau * hmc.dHdp(p, W);
      // p += -0.5*tau * hmc.dHdW(p, W);

      const double Hfin = md.H(p,W);
      std::cout << Hfin << std::endl;

      const double diff = Hfin-Hinit;
      std::cout << hmc.tau << "\t" << diff << std::endl;
    }



  // for(int nsteps=10; nsteps<400; nsteps*=2){
  //   const double stot = 1.0;
  //   //const int nsteps = 40;
  //   Integrator md(S, K, stot, nsteps);
  //   HMC hmc(md, rng, stot, nsteps);

  //   rng.seed( seed );
  //   W.randomize( [&](){ return rng.gaussian(); },
  // 		 [&](){ return rng.gaussian(); }
  // 		 );


  //   {
  //     Force p = K.gen( W, rng );

  //     const double Hinit = md.H(p,W);
  //     // std::cout << Hinit << std::endl;
  //     for(int i=0; i<md.nsteps; i++) md.onestep( p, W );
  //     const double Hfin = md.H(p,W);
  //     // std::cout << Hfin << std::endl;
  //     const double diff = Hfin-Hinit;
  //     std::cout << hmc.tau << "\t" << diff << std::endl;
  //   }
  // }



  // using Integrator = ExplicitLeapfrog<Force,Gauge,Action,Kernel>;
  // if (argc>2){ beta = atof(argv[2]); }
  // using Kernel = TrivialKernel;
  // Kernel K(Nc);
