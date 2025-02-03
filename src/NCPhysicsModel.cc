#include "NCCrystallineExtinction.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCRandUtils.hh"
//#include <iostream>

// ************************************************************* //
// ******************  Sabine's model  ************************* //
// ************************************************************* //

//Ref: International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616.

//A, B
struct ABpair {
  
  double A;
  double B;
  
};

ABpair calc_AB( double y ) {
  
  nc_assert( y >= 0. );
  
  ABpair result;
  if ( y <= 1.e-9 ) {
    result.A = 1.;
    result.B = 1.;
  }
  else {
    result.A = NCrystal::exp_negarg_approx(-y) * std::sinh(y) / y;
    result.B = 1. / y - NCrystal::exp_negarg_approx(-y) / std::sinh(y);
  }
  
  return result;
}

//extinction factors
struct ExtinctionFactors {
  
  double E_L; //corresponding to 2theta=0
  double E_B; //corresponding to 2theta=pi
  
};

//primary extinction factors
ExtinctionFactors prim_extn_fact( double x, double y ) {
  
  ExtinctionFactors result;
  double EL = NCrystal::exp_negarg_approx(-y);
  if ( x <= 1. ) {
    EL *= (1. - x / 2. + x * x / 4. - 5. * x * x * x / 48. + 7. * x * x * x * x / 192.);
  }
  else {
    EL *= std::sqrt(2. * NCrystal::kInvPi / x );
    EL *= (1. - 1. / 8. / x - 3. / 128. / x / x - 15. / 1024. / x / x / x);
  }
  
  ABpair AB = calc_AB( y );
  
  result.E_L = EL;
  result.E_B = AB.A / std::sqrt(1. + AB.B * x);
  
  return result;
}

//secondary extinction factors
ExtinctionFactors scnd_extn_fact( double x, double y, int tilt_dist ) {
  
  //tilt_dist : distribution type for the tilts between mosaic blocks
  //0 represents rectangular function, 1 for triangular function
  nc_assert( tilt_dist==0 || tilt_dist==1 );
  
  ExtinctionFactors result;
  ABpair AB = calc_AB( y );
  double Bx = AB.B * x;
  if ( tilt_dist == 0 ) {
    if ( x < 1.e-9 ) {
      result.E_L = NCrystal::exp_negarg_approx(-y);
    }
    else {
      result.E_L = NCrystal::exp_negarg_approx(-y) / 2. / x * (1. - NCrystal::exp_negarg_approx(-2 * x));
    }
    result.E_B = AB.A / (1. + Bx);
  }
  else {
    if ( x < 1.e-9 ) {
      result.E_L = NCrystal::exp_negarg_approx(-y);
      result.E_B = AB.A * AB.B;
    }
    else {
      result.E_L = NCrystal::exp_negarg_approx(-y) / x * (1. - (1. - NCrystal::exp_negarg_approx(-2 * x)) / 2. / x);
      result.E_B = 2 * AB.A / Bx / x * (Bx - std::log1p(Bx));
    }
  }
  
  return result;
}

//uncorrelated block model for calculating the extinction factor
double uncorr_blk_mdl( double Nc, double wl, double F_hkl, double l, 
                       double d_hkl, double mu, double G, double L, int tilt_dist ) {
  
  //Calculation of the extinction factor E_hkl using the uncorrelated block model
  //Nc : number of unit cells per unit volume, Aa^-3
  //wl : wavelength, Aa
  //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
  //l : block size, Aa
  //d_hkl : dspacing for the hkl plan, Aa
  //mu : incoherent and absorption cross section per unit volume, Aa^-1
  //G : integral breadth of the angular distribution of mosaic blocks, dimensionless
  //L : grain size, Aa (A grain is formed by crystallites or blocks.)
  //tilt_dist : distribution type for the tilts between mosaic blocks
  //0 represents rectangular function, 1 for triangular function
  
  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta <= 1. ) {
    double sin_theta_square = NC::ncsquare(sin_theta);
    double cos_theta_square = 1. - sin_theta_square;
    double cos_theta = std::sqrt(cos_theta_square);
    double y = mu * l;
    
    //primary extinction
    double xp = NC::ncsquare(Nc * wl * F_hkl * l);
    ExtinctionFactors EpLB = prim_extn_fact( xp, y );
    double Ep = EpLB.E_L * cos_theta_square + EpLB.E_B * sin_theta_square;
    
    //secondary extinction
    if ( sin_theta != 0. && cos_theta != 0. ) {
      double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / 2. / sin_theta / cos_theta;
      //double xs = Ep * Q_theta * G * l;
      double xs = Ep * Q_theta * G * L; //from equation (6.4.9.1) in "International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616."
      ExtinctionFactors EsLB = scnd_extn_fact( xs, y, tilt_dist );
      double Es = EsLB.E_L * cos_theta_square + EsLB.E_B * sin_theta_square;
      
      return Ep * Es;
    }
    else {
      return 0.;
    }
  }
  else {
    
    return 1.;
  }
}

//correlated block model for calculating the extinction factor
double corr_blk_mdl( double Nc, double wl, double F_hkl, double l,
                     double d_hkl, double mu, double g, double L ) {
  
  //Calculation of the extinction factor E_hkl using the uncorrelated block model
  //Nc : number of unit cells per unit volume, Aa^-3
  //wl : wavelength, Aa
  //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
  //l : block size, Aa
  //d_hkl : dspacing for the hkl plan, Aa
  //mu : incoherent and absorption cross section per unit volume, Aa^-1
  //g : standard deviation of the distribution of tilts * \sqrt(pi)/2, dimensionless
  //L : side of cube of the crystal, Aa (meaning to be clarified)

  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta <= 1. ) {
    double sin_theta_square = NC::ncsquare(sin_theta);
    double cos_theta_square = 1. - sin_theta_square;
    double cos_theta = std::sqrt(cos_theta_square);
    double y = mu * l;
    
    //refine both primary and secondary extinction in this model
    double E;
    if ( l > 0. && g == 0. ) {
      double x = NC::ncsquare(Nc * wl * F_hkl * l); //pure primary
      ExtinctionFactors ELB = prim_extn_fact( x, y );
      E = ELB.E_L * cos_theta_square + ELB.E_B * sin_theta_square;
    }
    else {
      if ( sin_theta != 0. && cos_theta != 0. ) {
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / 2. / sin_theta / cos_theta;
        double x = NC::ncsquare(Nc * wl * F_hkl * l + g * Q_theta * (L - l));
        ExtinctionFactors ELB = prim_extn_fact( x, y );
        E = ELB.E_L * cos_theta_square + ELB.E_B * sin_theta_square;
      }
      else {
        E = 0.;
      }
    }
    
    return E;
  }
  else {
    
    return 1.;
  }
}

// ************************************************************* //
// ***************  Becker & Coppens' model ******************** //
// ************************************************************* //

// Refs: Acta Cryst. (1974). A30, 129
//       Acta Cryst. (1995). A51, 662-667

ABpair calc_AB_theta( double cos_2theta, int opt ) {
  
  //Calculations of A(theta) and B(theta)
  //cos_2theta : cos(2*theta_hkl)
  //opt : 0 for primary extinction, 1, 2, 3 for sencondary extinction following a
  //Gaussian, Lorentzian or Fresnel distribution, respectively
  
  nc_assert( opt <= 3 );
  
  ABpair result;
  
  if ( opt == 0 ) {
    result.A = 0.20 + 0.45 * cos_2theta;
    result.B = 0.22 - 0.12 * NC::ncsquare(0.5 - cos_2theta);
  }
  else if ( opt == 1 ) {
    result.A = 0.58 + 0.48 * cos_2theta + 0.24 * NC::ncsquare(cos_2theta);
    result.B = 0.02 - 0.025 * cos_2theta;
  }
  else if ( opt == 2 ) {
    result.A = 0.025 + 0.285 * cos_2theta;
    if ( cos_2theta >= 0. ) {
      result.B = 0.15 - 0.2 * NC::ncsquare(0.75 - cos_2theta);
    }
    else {
      result.B = -0.45 * cos_2theta;
    }
  }
  else {
    result.A = 0.48 + 0.6 * cos_2theta;
    result.B = 0.20 - 0.06 * NC::ncsquare(0.2 - cos_2theta);
  }
  
  return result;
}

double BC_pure_extn_mdl( double Nc, double wl, double F_hkl, double l,
                         double d_hkl, double g, double L, int tilt_dist ) {

  //Calculation of pure primary or secondary extinction factor y using the model of Becker & Coppens
  // pure primary:  BC_pure  l
  // pure secondary type-I:  BC_pure  g  L
  // pure secondary type-II:  BC_pure  l  L  Gauss/Lorentz/Fresnel
  //Nc : number of unit cells per unit volume, Aa^-3
  //wl : wavelength, Aa
  //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
  //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
  //d_hkl : dspacing for the hkl plan, Aa
  //g : width parameter of the mosaic distribution, dimensionless
  //L : "T-bar", mean path length through a mosaic crystal, Aa
  //tilt_dist : option for the calculation of A(theta) and B(theta),
  //1, 2, 3 for orientation of crystallite following a Gaussian, Lorentzian
  //or Fresnel distribution, respectively

  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta <= 1. ) {
    double cos_theta  = std::sqrt(1. - NC::ncsquare(sin_theta));
    double sin_2theta = 2. * sin_theta * cos_theta;
    double cos_2theta = 1. - 2. * NC::ncsquare(sin_theta);
    //double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta; //same as in Sabine's model
    double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl; //division by sin_2theta to be done later

    double y;
    //pure primary extinction
    if ( l > 0. && g == 0. && L == 0. ) {
      ABpair AB_theta = calc_AB_theta( cos_2theta, 0 );
      //double x = 2. / 3. * Q_theta * l * l * sin_2theta / wl;
      double x = 2. / 3. * Q_theta * l * l / wl;
      y = 1. / std::sqrt(1. + 2. * x + AB_theta.A * NC::ncsquare(x) / (1. + AB_theta.B * x));
    }
    //pure secondary extinction type-I
    else if ( l == 0. && g > 0. && L > 0. ) {
      //double x = 2. * std::sqrt(2.) / 3. * g * Q_theta * L;
      if ( sin_2theta != 0. ) {
        Q_theta /= sin_2theta;
        y = 1. / std::sqrt(1. + 2. * g * Q_theta * L);
      }
      else {
        y = 0.;
      }
    }
    //pure secondary extinction type-II
    else if ( l > 0. && g == 0. && L > 0. ) {
      ABpair AB_theta = calc_AB_theta( cos_2theta, tilt_dist );
      //double x = 2. / 3. * Q_theta * L * l * sin_2theta / wl;
      double x = 2. / 3. * Q_theta * L * l / wl;
      y = 1. / std::sqrt(1. + 2.12 * x + AB_theta.A * NC::ncsquare(x) / (1. + AB_theta.B * x));
    }
    else {
      y = 1.;
    }

    return y;
  }
  else {

    return 1.;
  }

}

double BC_mix_extn_mdl( double Nc, double wl, double F_hkl, double l,
                        double d_hkl, double g, double L, int tilt_dist ) {
  
  //Calculation of mixed primary or secondary extinction factor y using the model of Becker & Coppens
  //Nc : number of unit cells per unit volume, Aa^-3
  //wl : wavelength, Aa
  //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
  //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
  //d_hkl : dspacing for the hkl plan, Aa
  //g : width parameter of the mosaic distribution, dimensionless
  //L : "T-bar", mean path length through a mosaic crystal, Aa
  //tilt_dist : option for the calculation of A(theta) and B(theta),
  //1, 2, 3 for orientation of crystallite following a Gaussian, Lorentzian
  //or Fresnel distribution, respectively
  
  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta >= -1. && sin_theta <= 1. ) {
    double cos_theta  = std::sqrt(1. - NC::ncsquare(sin_theta));
    double sin_2theta = 2. * sin_theta * cos_theta;
    double cos_2theta = 1. - 2. * NC::ncsquare(sin_theta);
    //double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta; //same as in Sabine's model
    double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl; //division by sin_2theta to be done later
    
    //primary
    ABpair AB_theta_p = calc_AB_theta( cos_2theta, 0 );
    //double xp = 2. / 3. * Q_theta * l * l * sin_2theta / wl;
    double xp = 2. / 3. * Q_theta * l * l / wl;
    double yp = 1. / std::sqrt(1. + 2. * xp + AB_theta_p.A * NC::ncsquare(xp) / (1. + AB_theta_p.B * xp));
    
    ABpair AB_theta_s = calc_AB_theta( cos_2theta, tilt_dist );
    double xs, ys;
    if ( l < 1.e-9 ) {
      xs = 0.;
      ys = 1.; //becomes pure primary
    }
    else {
      //xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l / sin_2theta) + 1. / (2. * g * g));
      xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l) + NC::ncsquare(sin_2theta) / (2. * g * g));
      ys = 1. / std::sqrt(1. + 2.12 * xs + AB_theta_s.A * NC::ncsquare(xs) / (1. + AB_theta_s.B * xp));
    }
    
    return yp * ys;
  }
  else {
    
    return 1.;
  }
}

// ************************************************************* //
// ******************  RED model  ****************************** //
// ************************************************************* //

//Ref: Acta Cryst. (1987). A43, 167-173, Acta Cryst. (1988). A44, 283-285
double red_mdl( double Nc, double wl, double F_hkl, double l,
                double d_hkl, double R, double L, double c ) {
  
  //Calculation of the extinction factor E_hkl using the Random Elastic Deformation (RED) model
  //Nc : number of unit cells per unit volume, Aa^-3
  //wl : wavelength, Aa
  //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
  //l : mean free path between elastically deformed regions
  //d_hkl : dspacing for the hkl plan, Aa
  //R : Bending radius parameter
  //L : crystal block size
  //c : deformation gradient parameter

  double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
  if ( sin_theta <= 1. ) {
    double sin_theta_square = NC::ncsquare(sin_theta);
    double cos_theta_square = 1. - sin_theta_square;
    double cos_theta = std::sqrt(cos_theta_square); //always positive
    double sin_2theta = 2. * sin_theta * cos_theta;
    double cos_2theta = cos_theta_square - sin_theta_square;
	    
    //refine both primary and secondary extinction in this model
    double x;
    if ( c >= 0. ) {
      double c_theta = std::sqrt(c + (1. - 2. * c) * cos_theta_square);
      double Delta_g = NC::kPi / Nc / wl / F_hkl; //unit, Aa
      double QRc;
      
      if ( sin_2theta != 0. ) {
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta;
        QRc = NCrystal::exp_negarg_approx(-Q_theta * R / c_theta);
      }
      else {
        QRc = 0;
      }
      x = 2. / 3. * (1. - QRc) * L / std::sqrt(NC::ncsquare(l) + NC::ncsquare(Delta_g));
    }
    else {
      if ( sin_2theta != 0. && cos_theta != 0. ) {
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta;
        x = 2. / 3. * (1. - NCrystal::exp_negarg_approx(-Q_theta * R / cos_theta)) * L / l; //original RED, Acta Cryst. (1987). A43, 167-173
      }
      else {
        x = 2. / 3. * L / l;
      }
    }
      
    //double a = 0.58 + 0.48 * cos_2theta + 0.24 * NC::ncsquare(cos_2theta);
    //double b = 0.02 - 0.025 * cos_2theta;
    //same coefficients as in the BC model for secondary extinction with Gaussian distribution
    ABpair AB_theta_s = calc_AB_theta( cos_2theta, 1 );
      
    double E = 1. / std::sqrt(1. + 2.12 * x + (AB_theta_s.A * NC::ncsquare(x) / (1. + AB_theta_s.B * x)));
      
    return E;
  }
  else {
    
    return 1.;
  }
}


// ************************************************************* //
// ******************  Edge-profile function  ****************** //
// ************************************************************* //

// Refs: Vogel's thesis (2000)
//       J. Imaging 2018, 4, 7
//       GUI-RITS(ver1.4.0) profile shape memo 09/27/2023

double jorgensen( double wl, double d_hkl, double a0, double a1,
                  double b0, double b1, double s02, double s12,
                  double s22, int conv_option ) {
  //profile shape function, converted from a Python function
  //conv_option : convolution option, 0 represents before convolution h(Delta)
  //              1 for after convolution k(Delta)
  
  nc_assert( conv_option == 0 || conv_option == 1 );
  
  double d2 = NC::ncsquare( d_hkl );
  double d4 = NC::ncsquare( d2 );
  double alpha = a0 + a1 / d_hkl;
  double beta = b0 + b1 / d4;
  double sigma2 = s02 + s12 * d2 + s22 * d4;
  double Delta = wl - 2. * d_hkl;
  
  double y = (alpha * sigma2 + Delta) / std::sqrt(2 * sigma2);
  double z = (beta *  sigma2 - Delta) / std::sqrt(2 * sigma2);
  double u = 0.5 * alpha * (alpha * sigma2 + 2 * Delta);
  double v = 0.5 * beta  * (beta  * sigma2 - 2 * Delta);
  
  double uy;
  if ( y > 25 ) {
    uy = std::exp(u - NC::ncsquare( y )) / y / std::sqrt(NC::kPi);
  }
  else {
    uy = std::exp(u) * std::erfc(y);
  }
  
  double vz;
  if ( z > 25 ) {
    vz = std::exp(v - NC::ncsquare( z )) / z / std::sqrt(NC::kPi);
  }
  else {
    vz = std::exp(v) * std::erfc(z);
  }
  
  double hkDelta;
  if ( conv_option == 0 ) {
    hkDelta = 0.5 * alpha * beta / (alpha + beta) * (uy + vz);
  }
  else {
    double w = Delta / std::sqrt(2 * sigma2);
    hkDelta = 0.5 * (std::erfc(w) - (beta * uy - alpha * vz) / (alpha + beta));
  }
  
  return hkDelta;
}


// ************************************************************* //
// *****************  parsing & processing  ******************** //
// ************************************************************* //

bool NCP::CrystallineExtinction::isApplicable( const NC::Info& info ) {
  
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::CrystallineExtinction NCP::CrystallineExtinction::createFromInfo( const NC::Info& info ) {
  
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
  //   Sabine_uncorr  l  G  L  rect/tri
  // or
  //   Sabine_corr  l  g  L
  // or
  //   BC_pure  l  g  L  Gauss/Lorentz/Fresnel
  // or
  //   BC_mix  l  g  L  Gauss/Lorentz/Fresnel
  // or
  //   RED_orig  l  R  T
  // or
  //   RED  l  R  T  c
  //

  //Two cases: one/two lines
  //Verify we have exactly one line and four or five words:
  int res_option;
  if ( data.size() == 1 ) {
    if ( data.at(0).size() != 4 && data.at(0).size() != 5 )
      NCRYSTAL_THROW2(BadInput,"One line in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" section has been detected, which should include four or five numbers, please check.");
    res_option = 0;
  }
  else {
    if ( data.size() != 2 || ( data.at(0).size() != 4 && data.at(0).size() != 5 ) || data.at(1).size() != 7 )
      NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" section considering resolution should be two lines, which include four or five, and seven numbers, respectively.");
    res_option = 1;
  }

  //Parse and validate values:
  int model_option;
  if ( data.at(0).at(0).compare("Sabine_uncorr") == 0 ) {
    model_option = 0;
  }
  else if ( data.at(0).at(0).compare("Sabine_corr") == 0 ) {
    model_option = 1;
  }
  else if ( data.at(0).at(0).compare("BC_pure") == 0 ) {
    model_option = 2;
  }
  else if ( data.at(0).at(0).compare("BC_mix") == 0 ) {
    model_option = 3;
  }
  else if ( data.at(0).at(0).compare("RED_orig") == 0 ) {
    model_option = 4;
  }  
  else if ( data.at(0).at(0).compare("RED") == 0 ) {
    model_option = 5;
  }
  else {
    NCRYSTAL_THROW2(BadInput,"Only the Sabine uncorrelated block model (Sabine_uncorr), the correlated block model (Sabine_corr), the BC models (BC_pure and BC_mix), and the random elastic deformation models (RED_orig and RED), are supported, please check the input file.");
  }
  
  double l, Gg, L, c;
  int tilt_dist_option;
  
  if ( model_option == 0 ) {
    if (   ! NC::safe_str2dbl( data.at(0).at(1), l  )
        || ! NC::safe_str2dbl( data.at(0).at(2), Gg )
        || ! NC::safe_str2dbl( data.at(0).at(3), L  )
        || ! (l  >= 0.0)
        || ! (Gg >= 0.0)
        || ! (L  >= 0.0) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" l, Gg and L should be three positive floating point values." );
    if ( data.at(0).at(4).compare("rect") == 0 ) {
      tilt_dist_option = 0;
    }
    else if ( data.at(0).at(4).compare("tri") == 0 ) {
      tilt_dist_option = 1;
    }
    else {
      NCRYSTAL_THROW2( BadInput,"Invalid strings specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" Dstribution option for tilts should be rect or tri." );
    }
  }
  else if ( model_option == 1 ) {
    if (   ! NC::safe_str2dbl( data.at(0).at(1), l  )
        || ! NC::safe_str2dbl( data.at(0).at(2), Gg )
        || ! NC::safe_str2dbl( data.at(0).at(3), L  )
        || ! (l  >= 0.0)
        || ! (Gg >= 0.0)
        || ! (L  >= 0.0) )
      NCRYSTAL_THROW2( BadInput,"Invalid values/strings specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" l, Gg and L should be three positive floating point values." );
  }
  else if ( model_option == 2 || model_option == 3 ) {
    if (   ! NC::safe_str2dbl( data.at(0).at(1), l  )
        || ! NC::safe_str2dbl( data.at(0).at(2), Gg )
        || ! NC::safe_str2dbl( data.at(0).at(3), L  )
        || ! (l  >= 0.0)
        || ! (Gg >= 0.0)
        || ! (L  >= 0.0) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" l, Gg and L should be three positive floating point values." );
    if ( data.at(0).at(4).compare("Gauss") == 0 ) {
      tilt_dist_option = 1;
    }
    else if ( data.at(0).at(4).compare("Lorentz") == 0 ) {
      tilt_dist_option = 2;
    }
    else if ( data.at(0).at(4).compare("Fresnel") == 0 ) {
      tilt_dist_option = 3;
    }
    else {
      NCRYSTAL_THROW2( BadInput,"Invalid strings specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" Dstribution option for the orientation of crystallite should be Gauss, Lorentz or Fresnel." );
    }
  }
  else if ( model_option == 4 ) {
    if (   ! NC::safe_str2dbl( data.at(0).at(1), l  )
        || ! NC::safe_str2dbl( data.at(0).at(2), Gg )
        || ! NC::safe_str2dbl( data.at(0).at(3), L  )
        || ! (l  >= 0.0)
        || ! (Gg >= 0.0)
        || ! (L  >= 0.0) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" l, Gg and L should be three positive floating point values." );
  }
  else {
    if (   ! NC::safe_str2dbl( data.at(0).at(1), l  )
        || ! NC::safe_str2dbl( data.at(0).at(2), Gg )
        || ! NC::safe_str2dbl( data.at(0).at(3), L  )
        || ! NC::safe_str2dbl( data.at(0).at(4), c  )
        || ! (l  >= 0.0)
        || ! (Gg >= 0.0)
        || ! (L  >= 0.0)
        || ! (c  >= 0.0 && c  <= 1.0 ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values/strings specified in the @CUSTOM_"<<pluginNameUpperCase()
                          <<" l, R, T and c should be four positive floating point values, c should be lie between 0 and 1." );
  }
  
  double a0, a1, b0, b1, s02, s12, s22;
  if ( res_option == 1) {
    if (   ! NC::safe_str2dbl( data.at(1).at(0), a0  )
        || ! NC::safe_str2dbl( data.at(1).at(1), a1  )
        || ! NC::safe_str2dbl( data.at(1).at(2), b0  )
        || ! NC::safe_str2dbl( data.at(1).at(3), b1  )
        || ! NC::safe_str2dbl( data.at(1).at(4), s02 )
        || ! NC::safe_str2dbl( data.at(1).at(5), s12 )
        || ! NC::safe_str2dbl( data.at(1).at(6), s22 )
        || ! (s02 >= 0.0)
        || ! (s12 >= 0.0)
        || ! (s22 >= 0.0) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                        <<" s02, s12 and s22 should be three positive floating point values." );
  }
  else {
    a0  = 0.;
    a1  = 0.;
    b0  = 0.;
    b1  = 0.;
    s02 = 0.;
    s12 = 0.;
    s22 = 0.;
  }
    
  //Getting the strcture info (volume, number of atoms)
  if ( !info.hasStructureInfo() )
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");
  const NCrystal::StructureInfo& struct_info = info.getStructureInfo();
  
  //Getting hklList
  if ( !info.hasHKLInfo() )
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");
  const NCrystal::HKLList& hklList = info.hklList();

  //Parsing done! Create and return our model:
  return CrystallineExtinction(model_option,l,Gg,tilt_dist_option,L,c,res_option,a0,a1,b0,b1,s02,s12,s22,struct_info,hklList);
}

NCP::CrystallineExtinction::CrystallineExtinction( int model_option, double l, double Gg,
                                                   int tilt_dist_option, double L, double c,
                                                   int res_option, double a0, double a1,
                                                   double b0, double b1, double s02,
                                                   double s12, double s22,
                                                   const NCrystal::StructureInfo& struct_info,
                                                   const NCrystal::HKLList& hklList )
  : m_model_option(model_option),
    m_l(l),
    m_Gg(Gg),
    m_tilt_dist_option(tilt_dist_option),
    m_L(L),
    m_c(c),	
    m_res_option(res_option),
    m_a0(a0),
    m_a1(a1),
    m_b0(b0),
    m_b1(b1),
    m_s02(s02),
    m_s12(s12),
    m_s22(s22)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_l   >= 0.0 );
  nc_assert( m_Gg  >= 0.0 );
  nc_assert( m_L   >= 0.0 );
  nc_assert( m_c   >= 0.0 );
  nc_assert( m_s02 >= 0.0 );
  nc_assert( m_s12 >= 0.0 );
  nc_assert( m_s22 >= 0.0 );
  
  m_Nc = 1. / struct_info.volume; //number of unit cells per unit volume, Aa^-3
  m_xsectfact = 0.5 / struct_info.n_atoms / struct_info.volume;
  
  for ( auto& hkl : hklList ) {
    hklDFM.emplace_back();
    hklDFM.back().dspacing     = hkl.dspacing;
    hklDFM.back().fsquared     = hkl.fsquared;
    hklDFM.back().multiplicity = hkl.multiplicity;
  }

}

double NCP::CrystallineExtinction::calcCrossSection( double neutron_ekin ) const {

  double E_hkl, R_hkl;
  double xs_in_barns = 0.0;

  const double wl = NC::ekin2wl( neutron_ekin );
  const double wlsq = NC::ncsquare( wl );
  const double mu = 0.; //contributions of absorption and incoherent scattering are removed, as in Sato et al. 2011
  for ( auto& hkl : hklDFM ) {
    
    double F_hkl = std::sqrt(hkl.fsquared) * 1.e-4; //Aa
    if ( m_model_option == 0 ) {
      E_hkl = uncorr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L, m_tilt_dist_option );
    }
    else if ( m_model_option == 1 ) {
      E_hkl = corr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L );
    }
    else if ( m_model_option == 2 ) {
      E_hkl = BC_pure_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
    }
    else if ( m_model_option == 3 ) {
      E_hkl = BC_mix_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
    }
    else if ( m_model_option == 4 ) {
      E_hkl = red_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, -1 );
    }
	  else{
      E_hkl = red_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_c );
    }	
    
    if ( m_res_option == 1 )
      R_hkl = jorgensen( wl, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 1 );
    else {
      R_hkl = 1.;
      if ( wl > 2 * hkl.dspacing )
        break;
    }

    xs_in_barns += hkl.dspacing * hkl.fsquared * hkl.multiplicity * E_hkl * R_hkl;
  }
  xs_in_barns *= m_xsectfact * wlsq; //lambda^2/(2NV)

  return xs_in_barns;
}

//double NCP::CrystallineExtinction::calcCrossSection( double neutron_ekin ) const {
//
//  //double E_hkl, R_hkl;
//  double ER_hkl; //E_hkl*R_hkl
//  double xs_in_barns = 0.0;
//
//  const double wl = NC::ekin2wl( neutron_ekin );
//  const double wlsq = NC::ncsquare( wl );
//  unsigned int num_wlp = 500; //can be changed later
//  const double mu = 0.; //contributions of absorption and incoherent scattering are removed, as in Sato et al. 2011
//  for ( auto& hkl : hklDFM ) {
//
//    double F_hkl = std::sqrt(hkl.fsquared) * 1.e-4; //Aa
//    if ( m_model_option == 0 ) {
//      if ( m_res_option == 1 ) {
//        double wlp_lower_bound = wl;
//        double wlp_upper_bound = NC::ncmax( wl + 1., 2 * hkl.dspacing + 2.5 );
//
//        ER_hkl = 0.;
//        for ( auto wlp : NC::linspace( wlp_lower_bound, wlp_upper_bound, num_wlp ) ) {
//          double h_wlp = jorgensen( wlp, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//          double E_wlp = uncorr_blk_mdl( m_Nc, wlp, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L, m_tilt_dist_option );
//          ER_hkl += h_wlp * E_wlp;
//        }
//        //double h0 = jorgensen( wlp_lower_bound, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//        //R_hkl -= 0.5 * h0;
//        //double hN = jorgensen( wlp_upper_bound, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//        //R_hkl += 0.5 * hN;
//
//        ER_hkl *= (wlp_upper_bound - wlp_lower_bound) / (num_wlp - 1);
//      }
//      else {
//        if ( wl > 2 * hkl.dspacing )
//          break;
//        ER_hkl = uncorr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L, m_tilt_dist_option );
//      }
//    }
//    else if ( m_model_option == 1 ) {
//      ER_hkl = corr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L );
//    }
//    else if ( m_model_option == 2 ) {
//      ER_hkl = BC_pure_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
//    }
//    else {
//      ER_hkl = BC_mix_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
//    }
//
//    //if ( m_res_option == 1 ) {
//    //  double wlp_lower_bound = wl;
//    //  double wlp_upper_bound = NC::ncmax( wl + 1., 2 * hkl.dspacing + 2.5 );
//    //
//    //  R_hkl = 0.;
//    //  for ( auto wlp : NC::linspace( wlp_lower_bound, wlp_upper_bound, num_wlp ) ) {
//    //    double h_wlp = jorgensen( wlp, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//    //    R_hkl += h_wlp;
//    //  }
//    //  double h0 = jorgensen( wlp_lower_bound, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//    //  R_hkl -= 0.5 * h0;
//    //  double hN = jorgensen( wlp_upper_bound, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 0 );
//    //  R_hkl += 0.5 * hN;
//    //
//    //  R_hkl *= (wlp_upper_bound - wlp_lower_bound) / (num_wlp - 1);
//    //  //R_hkl = jorgensen( wl, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 1 );
//    //}
//    //else {
//    //  R_hkl = 1.;
//    //  if ( wl > 2 * hkl.dspacing )
//    //    break;
//    //}
//
//    xs_in_barns += hkl.dspacing * hkl.fsquared * hkl.multiplicity * ER_hkl;
//  }
//  xs_in_barns *= m_xsectfact * wlsq; //lambda^2/(2NV)
//
//  return xs_in_barns;
//}

NCP::CrystallineExtinction::ScatEvent NCP::CrystallineExtinction::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const {

  ScatEvent result;
  result.ekin_final = neutron_ekin; //elastic scattering
  
  const double wl = NC::ekin2wl(neutron_ekin);
  const double wlsq = NC::ncsquare(wl);
  const double xs = calcCrossSection( neutron_ekin ) / m_xsectfact / wlsq; //calculate xs
  const double mu = 0.; //contributions of absorption and incoherent scattering are removed, as in Sato et al. 2011
  const double rnd = rng.generate(); //random number on [0;1]

  double E_hkl, R_hkl;
  double left_bound = 0.;
  double right_bound = 0.;
  for ( auto& hkl : hklDFM ) {
    
    double F_hkl = std::sqrt(hkl.fsquared) * 1.e-4; //Aa
    if ( m_model_option == 0 ) {
      E_hkl = uncorr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L, m_tilt_dist_option );
    }
    else if ( m_model_option == 1 ) {
      E_hkl = corr_blk_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, mu, m_Gg, m_L );
    }
    else if ( m_model_option == 2 ) {
      E_hkl = BC_pure_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
    }
	  else if ( m_model_option == 3 ) {
      E_hkl = BC_mix_extn_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_tilt_dist_option );
    }
    else if ( m_model_option == 4 ) {
      E_hkl = red_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, -1 );
    }
	  else {
      E_hkl = red_mdl( m_Nc, wl, F_hkl, m_l, hkl.dspacing, m_Gg, m_L, m_c );
    }    
    if ( m_res_option == 1 )
      R_hkl = jorgensen( wl, hkl.dspacing, m_a0, m_a1, m_b0, m_b1, m_s02, m_s12, m_s22, 1 );
    else {
      R_hkl = 1.;
      if ( wl > 2 * hkl.dspacing )
        break;
    }
    right_bound += hkl.dspacing * hkl.fsquared * hkl.multiplicity * E_hkl * R_hkl / xs;
    nc_assert( left_bound < right_bound && right_bound <= 1.0 );
      
    if ( left_bound <= rnd && right_bound > rnd ) {
      const double En_hkl = 0.5 * NC::kPiSq * NC::const_hhm / NC::ncsquare(hkl.dspacing);
      const double mu_n = 1. - 2 * En_hkl / neutron_ekin;
      nc_assert( NC::ncabs(mu_n) <= 1.0 );
      result.mu = mu_n;
      break;
    }
    else {
      left_bound = right_bound;
    }
  }

  return result;

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
