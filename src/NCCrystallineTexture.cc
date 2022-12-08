#include "NCCrystallineTexture.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCVector.hh"

//PDDF
double calc_pddf( NCrystal::Vector& preferred_orientation, NCrystal::Vector& tau_hkl, double p1 )
{
  //Calculation of the cylindrically symmetric Pole-Density Distribution Function (PDDF)
  //P_hkl(theta_hkl) which depends on the orientation angle theta_hkl, i.e., angle between
  //the preferred orientation and the plan vectors tau_hkl.
  //The implementation is based on the extinction code provided by D. DiJulio
  //preferred_orientation : preferred orientation of texture
  //tau_hkl : reciprocal lattice vector, or wavevector
  //p1 : coefficient in the March-Dollase texture model
  double Theta_hkl;        // orientation angle
  double pddf_hkl;         // Pole-Density Distribution Function
  Theta_hkl = std::acos(preferred_orientation.dot(tau_hkl) / preferred_orientation.mag() / tau_hkl.mag());
  pddf_hkl = std::pow(p1, 2) * std::pow(std::cos(Theta_hkl), 2) + std::pow(std::sin(Theta_hkl), 2) / p1;
  pddf_hkl = std::pow(pddf_hkl, -1.5);
 
  return pddf_hkl;
}

bool NCP::CrystallineTexture::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::CrystallineTexture NCP::CrystallineTexture::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <preferred orientation (array of 3)> <March-Dollase coefficient> <texture fraction>
  //

  //Verify we have exactly one line and five words:
  //Further extended to several lines with five words
  if ( data.size() != 1 || data.at(0).size()!=5 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be five numbers on a single line");

  //Parse and validate values:
  NCrystal::Vector& preferred_orientation;
  double p1, f1;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), preferred_orientation.at(0) )
       || ! NC::safe_str2dbl( data.at(0).at(1), preferred_orientation.at(1) )
       || ! NC::safe_str2dbl( data.at(0).at(2), preferred_orientation.at(2) )
       || ! NC::safe_str2dbl( data.at(0).at(3), p1 )
       || ! NC::safe_str2dbl( data.at(0).at(4), f1 )
       || !(preferred_orientation.mag()>0) || !(p1>0.0) || !(f1==1.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" PO (should not be [0,0,0]) p1 (should be one positive floating point value) f1 (should be 1, only one PO supported in this moment)" );
    
  //Getting the volume, number of atoms and hkl list
  double volume = info.getStructureInfo().volume;
  int n_atoms = info.getStructureInfo().n_atoms;
  const NCrystal::HKLList& hkl_list = info.hklList();

  //Parsing done! Create and return our model:
  return CrystallineTexture(preferred_orientation,p1,f1,hkl_list);
}

NCP::CrystallineTexture::CrystallineTexture( NCrystal::Vector& preferred_orientation, double p1, double f1,
                                             double volume, int n_atoms, const NCrystal::HKLList& hkl_list )
  : m_preferred_orientation(preferred_orientation),
    m_p1(p1),
    m_f1(f1),
    m_volume(volume),
    m_n_atoms(n_atoms),
    m_hkl_list(hkl_list)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( preferred_orientation.mag() > 0.0 );
  nc_assert( m_p1 > 0.0 );
  nc_assert( m_f1 > 0.0 );
  nc_assert( m_volume > 0.0 );
  nc_assert( m_n_atoms > 0 );
}

double NCP::CrystallineTexture::calcCrossSection( double neutron_ekin ) const
{
  double xs_in_barns;
  const double xsectfact = 0.5 / (m_volume * m_n_atoms);

  for ( auto& hkl : m_hkl_list ) {
    const double wl = 2.0 * hkl.dspacing;
    const double E = NC::wl2ekin(wl);
    if ( E <= neutron_ekin ) {
      //auto rec_lat = NCrystal::getReciprocalLatticeRot( info.getStructureInfo() );
      //NCrystal::Vector& waveVector = rec_lat * hkl;
      double pddf_hkl = calc_pddf( m_preferred_orientation, hkl.waveVector, m_p1); // March-Dollase texture model
      const double fdm = hkl.fsquared * hkl.multiplicity * hkl.dspacing * pddf_hkl;
      xs_in_barns += fdm;
    }
  }
  xs_in_barns *= m_f1 * xsectfact;
    
  return xs_in_barns;
}

NCP::CrystallineTexture::ScatEvent NCP::CrystallineTexture::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  //if ( ! (neutron_ekin > m_cutoffekin) ) {
    //Special case: We are asked to sample a scattering event for a neutron
    //energy where we have zero cross section! Although in a real simulation we
    //would usually not expect this to happen, users with custom code might
    //still generate such calls. The only consistent thing to do when the cross
    //section is zero is to not change the neutron state parameters, which means:
    //result.ekin_final = neutron_ekin;
    //result.mu = 1.0;
    //return result;
  //}

  //Implement our actual model here. Of course it is trivial for the example
  //model. For a more realistic or complicated model, it might be that
  //additional helper classes or functions should be created and used, in order
  //to keep the code here manageable:

  //result.ekin_final = neutron_ekin;//Elastic
  //result.mu = randIsotropicScatterMu(rng).dbl();
    
  //Same as coherent elastic scattering
  result.ekin_final = neutron_ekin;
  result.mu = randIsotropicScatterMu(rng).dbl(); // Take isotropic first for test

  return result;
}


