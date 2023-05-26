#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.

extern "C" {

  void nctest_getmanyxsvalues( double q, double radius, double thickness, double sld_core, double sld_shell, double sld_solvent, double result,  unsigned array_size, const double* ekin_array, double* output_xs_array )
  {
    
    NCP::PhysicsModel pm(q, radius, thickness, sld_core, sld_shell, sld_solvent, result);
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }

  void nctest_samplemanyscatmu( double q, double radius, double thickness, double sld_core, double sld_shell, double sld_solvent, double result, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(q, radius, thickness, sld_core, sld_shell, sld_solvent, result);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }

}
