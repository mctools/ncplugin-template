#ifndef NCPluginTestCode_ExtraTestUtils_hh
#define NCPluginTestCode_ExtraTestUtils_hh

// In this library should be placed extra functions / classes / etc. which are
// required for testing the plugin, but which are not a part of the plugin
// itself, and which might be used in more than one test application
// (app_*/main.cc) or Python script.
//
// Hooks accessible to Python (via ctypes in this example, for simplicity), are
// added in src/ForPython.cc and wrapped into a python interface in
// python/_wrap.py.
//
// Of course, a lot of the validation from python can simply use the standard
// NCrystal python module, the extra code here is just to show how one can
// directly access the internal helper classes used in the plugin from Python.

#include "NCPhysicsModel.hh"

//For our C++ test applications:

namespace NCPluginTestCode {

  std::vector<double> sampleAngles( const NCP::PhysicsModel&, double ekin, unsigned nvalues );
  std::vector<double> sampleXS( const NCP::PhysicsModel& pm, std::vector<double> ekins, unsigned nvalues );

}


#endif
