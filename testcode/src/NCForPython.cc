#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.

extern "C" {

  void nctest_getmanyxsvalues( double A1, double b1, double A2, double b2, double Q0, double sigma0, unsigned array_size, const double* ekin_array, double* output_xs_array )
  {
    NCP::PhysicsModel pm( A1,  b1,  A2,  b2,  Q0,  sigma0 );
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }

  void nctest_samplemanyscatmu( double A1, double b1, double A2, double b2, double Q0, double sigma0, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(A1,  b1,  A2,  b2,  Q0,  sigma0);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }

}
