#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.

extern "C" {

  void nctest_getmanyxsvalues( const std::vector<double>& I, const std::vector<double>& Q, unsigned array_size, const double* ekin_array, double* output_xs_array )
  {
    NCP::SansIsotropic pm(I,Q);
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }

  void nctest_samplemanyscatmu( const std::vector<double>& I, const std::vector<double>& Q, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::SansIsotropic pm(I,Q);
    auto rng = NC::defaultRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(*rng,ekin).mu;
  }

}
