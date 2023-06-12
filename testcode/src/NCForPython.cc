#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.

extern "C" {

  void nctest_getmanyxsvalues(double radius, double sld, double sld_solvent,  unsigned array_size, const double* ekin_array, double* output_xs_array )
  {
    
    NCP::PhysicsModel pm(radius, sld, sld_solvent);
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }

  void nctest_getmanyIQvalues(double radius, double sld, double sld_solvent,  unsigned array_size, const double* q, double* output_IQ_array )
  {
    
    NCP::PhysicsModel pm(radius, sld, sld_solvent);
    for (unsigned i = 0; i < array_size; ++i)
      output_IQ_array[i] = pm.calcIQ(q[i]);
  }

  void nctest_samplemanyscatmu( double radius, double sld, double sld_solvent, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(radius, sld, sld_solvent);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }

}
