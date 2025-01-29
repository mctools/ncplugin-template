#include "NCCrystallineTexture.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"

//preferred orientation distribution function (Sato 2011)
double sato_mmd_podf( const NCrystal::Vector& preferred_orientation, NCrystal::Vector vec_hkl,
                      double d_hkl, double R, double wl )
{
  //Calculation of the modified March-Dollase preferred orientation distribution function
  //reported in the paper of Sato et al. 2011
  //P_hkl(lambda, d_hkl, hkl, R)
  //preferred_orientation : preferred orientation of texture
  //vec_hkl : (h,k,l)
  //d_hkl : dspacing for the hkl plan
  //R : coefficient in the modified March-Dollase texture model
  //wl : wavelength, Aa
  //Note: P_hkl is symmetric in (h,k,l), i.e., P_hkl=P_-h-k-l
  double P_hkl = 1.;
  unsigned int num_phis = 1000; //can be changed later

  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta >= -1. && sin_theta <= 1. ) {
    double cos_theta = std::sqrt( 1 - NC::ncsquare(sin_theta) );

    double cos_A = preferred_orientation.dot(vec_hkl) / ( std::sqrt( preferred_orientation.mag2() * vec_hkl.mag2() ) );
    double sin_A = std::sqrt( 1 - NC::ncsquare(cos_A) );

    //trapezoidal integration
    P_hkl = 0.;
    for ( auto phi : NC::linspace( 0, NC::k2Pi * (1-1./num_phis), num_phis ) ) {
      double B = cos_A * sin_theta + sin_A * cos_theta * std::sin(phi); //integrand, to be optimised
      //since P_hkl(0)=P_hkl(2pi)
      P_hkl += std::pow( (NC::ncsquare(R * B) + (1 - NC::ncsquare(B)) / R), -1.5 ) / (num_phis+1);
    }
    //double epsilon = 1.E-9; //precision of the integration bounds
    //double Rm1R = NC::ncsquare(R) - 1. / R;
    //double a = Rm1R * NC::ncsquare(sin_A * cos_theta);
    //double b = -2 * Rm1R * sin_A * cos_A * sin_theta * cos_theta;
    //double c = Rm1R * NC::ncsquare(cos_A * sin_theta) + 1. / R;
    //for ( auto x : NC::linspace( -1.+epsilon, 1.-epsilon, num_phis ) ) {
    //  double x1mx = x * std::sqrt(1 - NC::ncsquare(x));
    //  P_hkl += 1. / std::sqrt( std::pow( 4 * a * NC::ncsquare(x1mx) + 2 * b * x1mx + c, 3 ) * (1 - NC::ncsquare(x)) );
    //}
    //double x1 = -1.+epsilon * std::sqrt(1 - NC::ncsquare(-1.+epsilon));
    //P_hkl -= 0.5 / std::sqrt( std::pow( 4 * a * NC::ncsquare(x1) + 2 * b * x1 + c, 3 ) * (1 - NC::ncsquare(-1.+epsilon)) );
    //P_hkl -= 0.5 / std::sqrt( std::pow( 4 * a * NC::ncsquare(x1) - 2 * b * x1 + c, 3 ) * (1 - NC::ncsquare(-1.+epsilon)) );
    //P_hkl /= NC::kPi * num_phis;
  }

  return P_hkl;
}

bool NCP::CrystallineTexture::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::CrystallineTexture NCP::CrystallineTexture::createFromInfo( const NC::SCOrientation& sco,
                                                                 const NC::Info& info,
                                                                 NC::PlaneProvider * std_pp)
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
  double R1, f1, R2, f2;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), preferred_orientation1.at(0) )
       || ! NC::safe_str2dbl( data.at(0).at(1), preferred_orientation1.at(1) )
       || ! NC::safe_str2dbl( data.at(0).at(2), preferred_orientation1.at(2) )
       || ! NC::safe_str2dbl( data.at(0).at(3), R1 )
       || ! NC::safe_str2dbl( data.at(0).at(4), f1 )
       || ! NC::safe_str2dbl( data.at(1).at(0), preferred_orientation2.at(0) )
       || ! NC::safe_str2dbl( data.at(1).at(1), preferred_orientation2.at(1) )
       || ! NC::safe_str2dbl( data.at(1).at(2), preferred_orientation2.at(2) )
       || ! NC::safe_str2dbl( data.at(1).at(3), R2 )
       || ! NC::safe_str2dbl( data.at(1).at(4), f2 )
       || !(preferred_orientation1.mag()>0) || !(R1>0.0) || !(f1>0.0)
       || !(preferred_orientation2.mag()>0) || !(R2>0.0) || !(f2>0.0)
       || !(f1+f2==1.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" POs (should not be [0,0,0]) R1,f1,R2,f2 (should be four positive floating point value) f1+f2 (should be 1, only two POs supported in this version)" );

  //Getting the strcture info (volume, number of atoms, reciprocal lattice rotation matrix)
  const NCrystal::StructureInfo& struct_info = info.getStructureInfo();

  //Parsing done! Create and return our model:
  return CrystallineTexture(sco,preferred_orientation1,R1,f1,preferred_orientation2,R2,f2,struct_info,std_pp);
}

NCP::CrystallineTexture::CrystallineTexture( const NC::SCOrientation& sco,
                                             const NCrystal::Vector& preferred_orientation1, double R1, double f1,
                                             const NCrystal::Vector& preferred_orientation2, double R2, double f2,
                                             const NCrystal::StructureInfo& struct_info,
                                             NCrystal::PlaneProvider * plane_provider )
  : m_preferred_orientation1(preferred_orientation1),
    m_R1(R1),
    m_f1(f1),
    m_preferred_orientation2(preferred_orientation2),
    m_R2(R2),
    m_f2(f2)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( preferred_orientation1.mag() > 0.0 );
  nc_assert( m_R1 > 0.0 );
  nc_assert( m_f1 > 0.0 );
  nc_assert( preferred_orientation2.mag() > 0.0 );
  nc_assert( m_R2 > 0.0 );
  nc_assert( m_f2 > 0.0 );
  nc_assert( plane_provider->canProvide() );

  m_reclat = getReciprocalLatticeRot( struct_info );
  NCrystal::RotMatrix lattice_rot = NC::getLatticeRot( struct_info.lattice_a, struct_info.lattice_b, struct_info.lattice_c,
                                                       struct_info.alpha*NC::kDeg, struct_info.beta*NC::kDeg, struct_info.gamma*NC::kDeg );
  m_lab2cry = getCrystal2LabRot( sco, m_reclat ).getInv();

  //RotMatrix cry2lab = getCrystal2LabRot( sco, m_reclat );
  double V0numAtom = struct_info.n_atoms * struct_info.volume;
  const double xsectfact = 0.5 / V0numAtom;

  plane_provider->prepareLoop();

  double maxdspacing(0);

  NCrystal::Optional<NC::PlaneProvider::Plane> opt_plane;
  while ( ( opt_plane = plane_provider->getNextPlane() ).has_value() ) {
    auto& pl = opt_plane.value();
    nc_assert( pl.dspacing > 0.0 );
    if ( pl.dspacing > maxdspacing )
      maxdspacing = pl.dspacing;
    m_hklPlanes.push_back( HKLPlane{} );
    auto& e = m_hklPlanes.back();
    e.hkl = lattice_rot * pl.demi_normal;
    e.d_hkl = pl.dspacing;
    e.strength = pl.dspacing * pl.fsq * xsectfact;
  }
}

double NCP::CrystallineTexture::calcCrossSection( NC::NeutronEnergy neutron_ekin, const NC::NeutronDirection& ndirlab ) const
{
  (void)ndirlab;//FIXME: Actually use this!!
  //auto ndir = ( m_lab2cry * ndirlab.as<NC::Vector>() ).unit();
  //ndir[0],ndir[1],ndir[2]
  //auto neutron_HKL = m_reclat * ndir;

  double xs_in_barns = 0.0;

  const double wl = neutron_ekin.wavelength().dbl();
  const double wlsq = NC::ncsquare(wl);
  for ( auto& e: m_hklPlanes ) {
    if ( wl > 2 * e.d_hkl )
      break;
    double P1 = sato_mmd_podf( m_preferred_orientation1, e.hkl, e.d_hkl, m_R1, wl );
    double P2 = sato_mmd_podf( m_preferred_orientation2, e.hkl, e.d_hkl, m_R2, wl );
    xs_in_barns += e.strength * (P1 * m_f1 + P2 * m_f2);
  }
  xs_in_barns *= 2.*wlsq; //consideration of the negative hkl

  return xs_in_barns;
}

NC::ScatterOutcome NCP::CrystallineTexture::sampleScatteringEvent( NC::RNG& rng, NC::NeutronEnergy neutron_ekin, const NC::NeutronDirection& ndirlab ) const
{
  //Don't do anything:
  //return { neutron_ekin, ndirlab };
  //return { neutron_ekin, NC::randIsotropicDirection(rng).as<NeutronDirection>() };

  NC::NeutronDirection outndirlab = ndirlab; //outgoing neutron direction
  const double wl = neutron_ekin.wavelength().dbl();
  const double wlsq = NC::ncsquare(wl);
  const double xs = calcCrossSection( neutron_ekin, ndirlab ) / (2.*wlsq); //calculate xs
  const double rnd = rng.generate(); //random number on [0;1]

  double left_bound = 0.;
  double right_bound = 0.;
  for ( auto& e: m_hklPlanes ) {
    if ( wl > 2 * e.d_hkl )
      break;
    double P1 = sato_mmd_podf( m_preferred_orientation1, e.hkl, e.d_hkl, m_R1, wl );
    double P2 = sato_mmd_podf( m_preferred_orientation2, e.hkl, e.d_hkl, m_R2, wl );
    right_bound += e.strength * (P1 * m_f1 + P2 * m_f2) / xs;
    nc_assert( left_bound < right_bound && right_bound <= 1.0 );

    if ( left_bound <= rnd && right_bound > rnd ) {
      const double E_hkl = 0.5 * NC::kPiSq * NC::const_hhm / NC::ncsquare(e.d_hkl);
      const double mu = 1. - 2 * E_hkl / neutron_ekin.dbl();
      nc_assert( NC::ncabs(mu) <= 1.0 );
      outndirlab.as<NC::Vector>() = NC::randDirectionGivenScatterMu( rng, mu, ndirlab.as<NC::Vector>() );
      break;
    }
    else {
      left_bound = right_bound;
    }
  }

  return { neutron_ekin, outndirlab };

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
  //result.ekin_final = neutron_ekin.dbl();
  //result.mu = randIsotropicScatterMu(rng).dbl(); // Take isotropic first for test

  // return result;
}
