#pragma once

/*
  Kernel objects should have:

  double operator()( const Force& pi, const Gauge& W ) const; // G(U) * pi
  Force d( const Force& pi, const Gauge& W ) const;
  double det( const Gauge& W ) const; // det G(U)
  Force det_log_d( const Gauge& W ) const; // d [log det G(U)]
*/


struct TrivialKernel { // Force = ForceSingleLink, Gauge = LinkConfig
  using Force = ForceSingleLink;
  using Gauge = LinkConfig;

  TrivialKernel(){};
  double operator()( const Force& f, const Gauge& W ) const { return f.pi.squaredNorm(); }
  Force d( const Force& f, const Gauge& W ) const { return Force(W.Nc); } // zero
  Force act( const Gauge& W, const Force& f ) const { return f; }
  double det( const Gauge& W ) const { return 1.0; }
  Force det_log_d( const Gauge& W ) const { return Force(W.Nc); } // zero
};
