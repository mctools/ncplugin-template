#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.
extern "C" {

 /*void nctest_getmanyxsvalues( int model,double A, double s, double rg, double m,double p, unsigned array_size, const double* ekin_array, double* output_xs_array )
    {
    NCP::PhysicsModel pm( model,A, s, rg, m, p );
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }*/
  void nctest_getmanyxsvalues( NCrystalPlugin_SANSND::PhysicsModel::Model model,char*filename, unsigned int array_size, const double* ekin_array, double* output_xs_array )
  {
    std::string filename_s(filename);
    NCP::PhysicsModel pm(model, filename_s );
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }/*
  void nctest_getmanyxsvalues( int model,double mono_R, unsigned int array_size, const double* ekin_array, double* output_xs_array )
  {
    NCP::PhysicsModel pm( model,mono_R );
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }*/
  /*void nctest_getmanyxsvalues( int model,char*filename, unsigned int array_size, const double* ekin_array, double* output_xs_array )
  {
    std::string filename_s(filename);
    NCP::PhysicsModel pm(model,filename_s);
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcCrossSection(ekin_array[i]);
  }

  void nctest_samplemanyscatmu( int model,double A, double s, double rg, double m, double p, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(model,A, s, rg, m, p);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }
  void nctest_samplemanyscatmu( char*filename, double ekin, unsigned nvalues, double* output_mu )
  {
    std::string filename_s(filename);
    NCP::PhysicsModel pm(filename_s);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }
  void nctest_samplemanyscatmu( int model, double mono_R,double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(model,mono_R);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }*/
  void nctest_samplemanyscatmu( NCrystalPlugin_SANSND::PhysicsModel::Model model, char*filename,double ekin, unsigned nvalues, double* output_mu )
  {
    std::string filename_s(filename);
    NCP::PhysicsModel pm(model,filename_s);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }


}