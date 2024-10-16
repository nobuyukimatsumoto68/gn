  const int Nc=3;
  LinkVar W(Nc);

  for(auto m : W.t){
    std::cout << m << std::endl;
    std::cout << (m*m).trace() << std::endl;
  }
