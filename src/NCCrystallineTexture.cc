#include "NCCrystallineTexture.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCMatrix.hh"

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

  //Verify we have exactly two line and five words:
  //Further extended to several lines with five words
  if ( data.size() != 2 || data.at(0).size()!=5 || data.at(1).size()!=5 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be five numbers on two lines");

  //Parse and validate values:
  NCrystal::Vector preferred_orientation1, preferred_orientation2;
  double p1, f1, p2, f2;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), preferred_orientation1.at(0) )
       || ! NC::safe_str2dbl( data.at(0).at(1), preferred_orientation1.at(1) )
       || ! NC::safe_str2dbl( data.at(0).at(2), preferred_orientation1.at(2) )
       || ! NC::safe_str2dbl( data.at(0).at(3), p1 )
       || ! NC::safe_str2dbl( data.at(0).at(4), f1 )
       || ! NC::safe_str2dbl( data.at(1).at(0), preferred_orientation2.at(0) )
       || ! NC::safe_str2dbl( data.at(1).at(1), preferred_orientation2.at(1) )
       || ! NC::safe_str2dbl( data.at(1).at(2), preferred_orientation2.at(2) )
       || ! NC::safe_str2dbl( data.at(1).at(3), p2 )
       || ! NC::safe_str2dbl( data.at(1).at(4), f2 )
       || !(preferred_orientation1.mag()>0) || !(p1>0.0) || !(f1>0.0)
       || !(preferred_orientation2.mag()>0) || !(p2>0.0) || !(f2>0.0)
       || !(f1+f2==1.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" POs (should not be [0,0,0]) p1,f1,p2,f2 (should be four positive floating point value) f1+f2 (should be 1, only two POs supported in this version)" );
    
  //Getting the strcture info (volume, number of atoms, reciprocal lattice rotation matrix) and hkl list
  const NCrystal::StructureInfo& struct_info = info.getStructureInfo();
  const NCrystal::HKLList& hkl_list = info.hklList();

  //Parsing done! Create and return our model:
  return CrystallineTexture(preferred_orientation1,p1,f1,preferred_orientation2,p2,f2,struct_info,hkl_list);
}

NCP::CrystallineTexture::CrystallineTexture( NCrystal::Vector& preferred_orientation1, double p1, double f1,
                                             NCrystal::Vector& preferred_orientation2, double p2, double f2,
                                             const NCrystal::StructureInfo& struct_info,
                                             const NCrystal::HKLList& hkl_list )
  : m_preferred_orientation1(preferred_orientation1),
    m_p1(p1),
    m_f1(f1),
    m_preferred_orientation2(preferred_orientation2),
    m_p2(p2),
    m_f2(f2),
    m_struct_info(struct_info),
    m_hkl_list(hkl_list)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( preferred_orientation1.mag() > 0.0 );
  nc_assert( m_p1 > 0.0 );
  nc_assert( m_f1 > 0.0 );
  nc_assert( preferred_orientation2.mag() > 0.0 );
  nc_assert( m_p2 > 0.0 );
  nc_assert( m_f2 > 0.0 );
}

double NCP::CrystallineTexture::calcCrossSection( double neutron_ekin ) const
{
  const double pi2hhm = NCrystal::kPiSq * NCrystal::const_hhm;
  double volume = m_struct_info.volume;
  int n_atoms = m_struct_info.n_atoms;
  const NCrystal::RotMatrix& rec_lat = NCrystal::getReciprocalLatticeRot( m_struct_info );
  double xsectfact = pi2hhm / (volume * n_atoms * neutron_ekin);
  double xs_in_barns;
  double xs_in_barns1 = 0.0;
  double xs_in_barns2 = 0.0;

  for ( auto& hkl : m_hkl_list ) {
    double E_hkl = 0.5 * pi2hhm / std::pow(hkl.dspacing, 2);
    if ( E_hkl <= neutron_ekin ) {
      NCrystal::Vector hkl_vector = NCrystal::Vector(hkl.hkl.h, hkl.hkl.k, hkl.hkl.l);
      NCrystal::Vector tau_hkl = rec_lat * hkl_vector;
      double pddf_hkl1 = calc_pddf( m_preferred_orientation1, tau_hkl, m_p1); // March-Dollase texture model
      double pddf_hkl2 = calc_pddf( m_preferred_orientation2, tau_hkl, m_p2);
      double fdm1 = hkl.fsquared * hkl.multiplicity * hkl.dspacing * pddf_hkl1;
      double fdm2 = hkl.fsquared * hkl.multiplicity * hkl.dspacing * pddf_hkl2;
      xs_in_barns1 += fdm1;
      xs_in_barns2 += fdm2;
    }
  }
  xs_in_barns = m_f1 * xs_in_barns1 + m_f2 * xs_in_barns2;
  xs_in_barns *= xsectfact;
    
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


