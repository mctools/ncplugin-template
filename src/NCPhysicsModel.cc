#include "NCPhysicsModel.hh"

// Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

#include <algorithm>
#include <sys/stat.h>
#include <fstream>

//Until upstream provides an NCPLUGIN_MSG macro, we define it here (HOWEVER:
//Note that it is better for the scattering process to simply emit some json
//information):
#ifndef NCPLUGIN_MSG
#  define NCPLUGIN_MSG(msg) NCRYSTAL_RAWOUT( "plugin::" NCPLUGIN_NAME_CSTR ": " << msg << '\n' )
#endif

bool NCP::PhysicsModel::isApplicable(const NC::Info &info)
{
  // Accept if input is NCMAT data with @CUSTOM_SANSND section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo(const NC::Info &info)
{
  // Parse the content of our custom section. In case of syntax errors, it
  // raises BadInput exceptions
  // Get the relevant custom section data (and verify that there are not multiple
  // such sections in the input data):
  if (info.countCustomSections(pluginNameUpperCase()) != 1)
    NCRYSTAL_THROW2(BadInput, "Multiple @CUSTOM_" << pluginNameUpperCase() << " sections are not allowed");
  auto data = info.getCustomSection(pluginNameUpperCase());

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_SANSND
  //    v                   # plugin version
  //    model               # which model
  //    parametrs           # paramters for the selected model
  //
  // Verify we have exactly 1-1-x input :
  if (data.at(0).size() != 1 || data.at(1).size() != 1) {
    NCRYSTAL_THROW2(BadInput, "Bad input for the @CUSTOM_" << pluginNameUpperCase()
                    << " section (see the plugin readme for more info).");
  }
  // Parse and validate values:
  double supp_version = 2.0;
  double version;
  NC::VectD param;

  if (!NC::safe_str2dbl(data.at(0).at(0), version)) {
    NCRYSTAL_THROW2(BadInput, "Invalid version input in the @CUSTOM_" << pluginNameUpperCase()
                    << " section (see the plugin readme for more info).");
  }
  // special warning for wrong version
  if (!(version == supp_version)) {
    NCRYSTAL_THROW2(BadInput, "Invalid version specified for the " << pluginNameUpperCase()
                    << " plugin. Only the version " << supp_version << " is supported.");
  }
  if ( data.at(1).at(0) == "FILE" ) {
    NCPLUGIN_MSG("FILE MODE selected");
    if (data.at(2).empty()) {
      NCRYSTAL_THROW2(BadInput, "FILE MODE selected but no filename is specified for the "
                      << pluginNameUpperCase() << " plugin ");
    }

    std::string filename = data.at(2).at(0);
    std::string root_rel = "data/";
    std::string rel_path = root_rel + filename;
    NCPLUGIN_MSG("Input file: " << rel_path);
    // check existence
    struct stat buffer;
    if (!(stat(rel_path.c_str(), &buffer) == 0))
      NCRYSTAL_THROW2(BadInput, "The filename specified for the " << pluginNameUpperCase()
                      << " plugin is invalid or the file could not be found in the data/ directory. ");
    Model model = Model::FILE;
    return PhysicsModel(model, filename);
  } else if (data.at(1).at(0) == "PPF") {
    NCPLUGIN_MSG("Mode PPF selected");
    double A1, b1, A2, b2, Q0, corr;
    if (!NC::safe_str2dbl(data.at(2).at(0), A1)
        || !NC::safe_str2dbl(data.at(2).at(1), b1)
        || !NC::safe_str2dbl(data.at(2).at(2), A2)
        || !NC::safe_str2dbl(data.at(2).at(3), b2)
        || !NC::safe_str2dbl(data.at(2).at(4), Q0)
        || !NC::safe_str2dbl(data.at(2).at(5), corr)) {
      NCRYSTAL_THROW2(BadInput, "Invalid values specified for model PPF in the @CUSTOM_"
                      << pluginNameUpperCase()
                      << " section (see the plugin readme for more info)");
    }
    // CHECK THE INPUT PARAM
    nc_assert_always(A1 > 0);
    nc_assert_always(b1 > 0);
    nc_assert_always(A2 > 0);
    nc_assert_always(b2 > 0);
    nc_assert_always(Q0 > 0);
    nc_assert_always(corr > 0);

    param.insert(param.end(), {A1, b1, A2, b2, Q0, corr});
    Model model = Model::PPF;
    // Parsing done! Create and return our model:
    return PhysicsModel(model, param);
  } else if (data.at(1).at(0) == "GPF") {
    NCPLUGIN_MSG("Mode GPF selected");
    double A, s, rg, m, p, Qmin, Q1;
    if ( !NC::safe_str2dbl(data.at(2).at(0), A)
         || !NC::safe_str2dbl(data.at(2).at(1), s)
         || !NC::safe_str2dbl(data.at(2).at(2), rg)
         || !NC::safe_str2dbl(data.at(2).at(3), m)
         || !NC::safe_str2dbl(data.at(2).at(4), p)
         || !NC::safe_str2dbl(data.at(2).at(5), Qmin) ) {
      NCRYSTAL_THROW2(BadInput, "Invalid values specified for model GPF in the @CUSTOM_"
                      << pluginNameUpperCase()
                      << " section (see the plugin readme for more info)");
    }
    // CHECK THE INPUT PARAM
    nc_assert_always(A > 0);
    nc_assert_always(s > 0);
    nc_assert_always(rg > 0);
    nc_assert_always(m > 0);
    nc_assert_always(p > 0);
    nc_assert_always(Qmin > 0);
    if (data.at(2).size()==6) {
      NCPLUGIN_MSG("Q1 not inserted. Default 0.016 used.");
      Q1 = 0.016;
    } if (data.at(2).size()==7) {
      if(!NC::safe_str2dbl(data.at(2).at(6), Q1))
        NCRYSTAL_THROW2(BadInput, "Invalid Q1 value specified for model GPF in the @CUSTOM_" << pluginNameUpperCase()
                        << " section (see the plugin readme for more info)");
    }
    nc_assert_always(Q1 > 0);
    param.insert(param.end(), {A, s, rg, m, p, Qmin, Q1});
    Model model = Model::GPF;
    // Parsing done! Create and return our model:
    return PhysicsModel(model, param);
  } else if (data.at(1).at(0) == "HSFBA") {
    NCPLUGIN_MSG("Mode HSFBA selected");
    Model model = Model::HSFBA;
    double R;
    if (NC::safe_str2dbl(data.at(2).at(0), R)) {
      nc_assert_always(R > 0);
      param.push_back(R);
      return PhysicsModel(model, param);
    } else {
      std::string filename = data.at(2).at(0);
      std::string root_rel = "data/";
      std::string rel_path = root_rel + filename;
      NCPLUGIN_MSG("Input file: " << rel_path);
      // check existence:
      struct stat buffer;
      if (!(stat(rel_path.c_str(), &buffer) == 0)) {
        NCRYSTAL_THROW2(BadInput, "The filename specified for the " << pluginNameUpperCase()
                        << " plugin is invalid or the file could not be found in the data/ directory. ");
      }
      // CHECK THE INPUT PARAM
      return PhysicsModel(model, filename);
    }
  } else {
    NCRYSTAL_THROW2(BadInput, "Invalid model input in the @CUSTOM_" << pluginNameUpperCase()
                    << " section (see the plugin readme for more info).");
  }
};

NCP::PhysicsModel::PhysicsModel(Model model, std::string filename)
  : m_model(model),
    m_helper(([model, filename]() -> NC::IofQHelper
    {
      NC::VectD q;
      NC::VectD IofQ;
      switch(model) {
      case Model::FILE:
        {
          //Parse the input file and create the vector with the q, Iq info
          std::string root_rel = "data/";
          std::string rel_path = root_rel+filename;
          std::ifstream input_file(rel_path);
          if(input_file) {
            double temp_x, temp_y;
            std::string line;
            while (std::getline(input_file, line)){
              std::istringstream iss(line);
              iss >> temp_x >> temp_y;
              q.push_back(temp_x);
              IofQ.push_back(temp_y);
            }
          } else {
            NCRYSTAL_THROW2( BadInput,rel_path << ": Invalid data file for the "<<pluginNameUpperCase()
                             <<" plugin" );
          }
          break;
        }
      case Model::HSFBA:
        {
          NC::VectD Rs,freq;
          std::string root_rel = "data/";
          std::string rel_path = root_rel + filename;
          std::ifstream input_file(rel_path);
          if (input_file)
            {
              double temp_x, temp_y;
              std::string line;
              while (std::getline(input_file, line))
                {
                  std::istringstream iss(line);
                  iss >> temp_x >> temp_y;
                  Rs.push_back(temp_x);
                  freq.push_back(temp_y);
                }
            }
          else
            {
              NCRYSTAL_THROW2(BadInput, rel_path << ": Invalid data file for the " << pluginNameUpperCase()
                              << " plugin");
            }

          //Generate vector of data q and IofQ
          double q_min = std::log10(1e-6);
          int sampling =  std::abs(1-q_min)*10000;
          q = NC::logspace(q_min,1,sampling);
          IofQ = q;
          double b = 6.646E-05;  // [AA] Carbon coherent scattering length
          double n = 0.1771471666666667; // [at/AA^3] <- Diamond atom density
          double physical_constant = 16*NC::kPi*NC::kPi*std::pow(n*b, 2);  // [1/AA^4]
          std::for_each(IofQ.begin(),IofQ.end(),
                        [Rs,freq,physical_constant](double &x) {
                          double R, f, osc_term, Nc;
                          double I=0;
                          for (size_t i = 0; i < Rs.size(); ++i)
                            {
                              R = Rs.at(i) * 10; // converto to AA
                              nc_assert_always(R>0);
                              f = freq.at(i);
                              osc_term = (sin(x*R) - x*R * cos(x*R)) * (sin(x*R) - x*R * cos(x*R));
                              //Determine the number of atoms in a diamon nanoparticle to normalize per-atoms
                              // Nc = V * n = 4/3*pi*R^3 * n
                              Nc = 0.7420 * R * R * R;
                              I += f * osc_term / Nc;
                            }
                          x = physical_constant * std::pow(x, -6) * I * 1e8 ;//# convert to barn
                        }
                        );
          break;
        }
      }
      NC::IofQHelper helper(q,IofQ);
      return helper;
    })())
{
  // NCPLUGIN_MSG("call to constructor for 0");
  // NCPLUGIN_MSG("helper initialized: " << m_helper.has_value());
};

NCP::PhysicsModel::PhysicsModel(Model model, NC::VectD param)
  : m_model(model),
    m_param(param),
    m_helper(([model, param]() -> NC::IofQHelper
    {
      NC::VectD q;
      NC::VectD IofQ;
      switch(model)
        {
        case Model::GPF:
          { //GPF
            NCPLUGIN_MSG("Intialize GPF");
            nc_assert_always(param.size()==7);

            double A=param.at(0);
            double s=param.at(1);
            double rg=param.at(2);
            double m=param.at(3);
            double p=param.at(4);
            double Qmin=param.at(5);
            double Q1=param.at(6);

            //Generate vector of data q and IofQ

            int sampling =  std::abs(1-std::log10(Qmin))*10000;
            q = NC::logspace(std::log10(Qmin),1,sampling);
            IofQ = q;
            //Q1 is when IofQ stops being evaluated as power law and Guinier starts (maybe parameter?)
            NCrystal::VectD::iterator it_q1;
            if (Q1>Qmin){
              it_q1 = std::lower_bound(IofQ.begin(),IofQ.end(), Q1);
              //Approximation valid as long as we have high sampling. Otherwise interpolation needed
              if (it_q1==IofQ.end())
                NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                                 <<" section (see the plugin readme for more info)" );
              nc_assert_always(s!=3);
              double C = A*std::pow(Q1,p-s)*std::exp((-Q1*Q1*rg*rg)/(3-s));
              std::for_each(IofQ.begin(),it_q1,
                            [C,p](double &x) { x = C*std::pow(x,-p);}
                            );
            } else {
              NCPLUGIN_MSG("Given Q_min greater then Q1. No low-q power law implemented (see the plugin readme for more info)");
              it_q1 = IofQ.begin();
            }
            //evaluate Q2 where IofQ stops being evaluated as Guinier and Porod starts
            double Q2 = 1.0/rg*std::sqrt((m-s)*(3-s)/2);
            //IofQ still filled with q here
            auto it_q2 = std::lower_bound(it_q1,IofQ.end(), Q2);
            //Approximation valid as long as we have high sampling. Otherwise interpolation needed
            if (it_q2==IofQ.end())
              NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q2 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                               <<" section (see the plugin readme for more info)" );
            if (it_q2==it_q1)
              NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1==Q2 in the @CUSTOM_"<<pluginNameUpperCase()
                               <<" section (see the plugin readme for more info)" );
            double B = A*std::pow(Q2,m-s)*std::exp((-Q2*Q2*rg*rg)/(3-s));
            std::for_each(it_q1,it_q2,
                          [A,s,rg](double &x) { x = A*std::pow(x,-s)*std::exp((-x*x*rg*rg)/(3-s));}
                          );
            std::for_each(it_q2,IofQ.end(),
                          [B,m](double &x) { x = B*std::pow(x,-m);}
                          );
            break;
          }
        case Model::PPF:
          {
            // IT IS NOT ACTUALLY USED; WHAT TO DO WITH IT?
            nc_assert_always(param.size()==6);
            double A1=param.at(0);
            double b1=param.at(1);
            double A2=param.at(2);
            double b2=param.at(3);
            double Q0=param.at(4);
            double corr=param.at(5);
            double q_min = std::log10(1e-6);
            int sampling =  std::abs(1-q_min)*10000;
            q = NC::logspace(q_min,1,sampling);
            IofQ = q;
            auto it_q0 = std::lower_bound(IofQ.begin(),IofQ.end(), Q0);
            if (it_q0==IofQ.end())
              NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q0 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                               <<" section (see the plugin readme for more info)" );
            std::for_each(IofQ.begin(),it_q0,
                          [A1,b1](double &x) { x = A1*std::pow(x,-b1);}
                          );
            std::for_each(it_q0,IofQ.end(),
                          [A2,b2](double &x) { x = A2*std::pow(x,-b2);}
                          );
            break;
          }
        case Model::HSFBA:
          {
            nc_assert_always(param.size()==1);
            double mono_R = param.at(0);
            double q_min = std::log10(1e-6);
            int sampling =  std::abs(1-q_min)*10000;
            q = NC::logspace(q_min,1,sampling);
            IofQ = q;
            double b = 6.646E-05;  // [AA] Carbon coherent scattering length
            double n = 0.1771471666666667; // [at/AA^3] <- Diamond atom density
            double physical_constant = 16*NC::kPi*NC::kPi*std::pow(n*b, 2);  // [1/AA^4]
            std::for_each(IofQ.begin(),IofQ.end(),
                          [mono_R,physical_constant](double &x) {
                            double R, osc_term, Nc;
                            double I=0;
                            R = mono_R * 10; // converto to AA
                            nc_assert_always(R>0);
                            osc_term = (sin(x*R) - x*R * cos(x*R)) * (sin(x*R) - x*R * cos(x*R));
                            //Determine the number of atoms in a diamon nanoparticle to normalize per-atoms
                            // Nc = V * n = 4/3*pi*R^3 * n
                            Nc = 0.7420 * R * R * R;
                            I += std::pow(x, -6) * osc_term / Nc * 1e8;//# convert to barn
                            x = physical_constant* I ;
                          }
                          );
            break;
          }
        default:
          NCRYSTAL_THROW2(LogicError,"Bad internal state in costructor for model SANSND plugin");

        }
      //Initialize the helper
      NC::IofQHelper helper(q,IofQ);
      return helper; })())
{
  // NCPLUGIN_MSG("call to constructor for 1 and 2");
  // NCPLUGIN_MSG("helper initialized: " << m_helper.has_value());
};

double NCP::PhysicsModel::calcCrossSection(double neutron_ekin) const
{
  if ( !( neutron_ekin > 0.0 ) )
    return 0.0;

  NC::NeutronEnergy ekin(neutron_ekin);
  double k = NC::k2Pi / NC::ekin2wl(neutron_ekin); // wavevector
  nc_assert_always(k != 0);
  double SANS_xs;

  switch (m_model)
    {
    case Model::PPF:
      {
        if (m_param.has_value())
          {
            double A1 = m_param.value().at(0);
            double b1 = m_param.value().at(1);
            double A2 = m_param.value().at(2);
            double b2 = m_param.value().at(3);
            double Q0 = m_param.value().at(4);
            double corr = m_param.value().at(5);
            SANS_xs = (2 * NC::kPi / (k * k)) * (A1 / (-b1 + 2) * std::pow(Q0, -b1 + 2) + A2 / (-b2 + 2) * std::pow(2 * k, -b2 + 2) - A2 / (-b2 + 2) * std::pow(Q0, -b2 + 2));
            SANS_xs = SANS_xs*corr; // correction added for comparison in thesis
          }
        else
          {
            NCRYSTAL_THROW2(LogicError, "Attempt to use not-initialized parameter in xs sampling in the " << pluginNameUpperCase()
                            << " plugin");
          }
        break;
      }
    case Model::FILE:
    case Model::GPF:
    case Model::HSFBA:
      {
        if (m_helper.has_value())
          {
            nc_assert_always(k != 0);
            SANS_xs = 2 * NC::kPi / (k * k) * m_helper.value().calcQIofQIntegral(ekin);
          }
        else
          {
            NCRYSTAL_THROW2(LogicError, "Attempt to use not-initialized IofQHelper in xs sampling in the " << pluginNameUpperCase()
                            << " plugin");
          }
        break;
      }
    default:
      NCRYSTAL_THROW2(LogicError, "Bad internal state for model in SANSND plugin");
    }

  if ( SANS_xs < 0.0 )
    SANS_xs = 0.0;//Added by TK since the code above seems to give negative results!! Probably this fix hides the real underlying issue.

  return SANS_xs;
}

double NCP::PhysicsModel::sampleScatteringVector(NC::RNG &rng, double neutron_ekin) const
{
  double Q;
  switch (m_model)
    {
    case Model::PPF:
      {
        double rand = rng.generate();
        double k = NC::k2Pi / NC::ekin2wl(neutron_ekin); // wavevector
        // sample a random scattering vector Q from the inverse PPF CDF (see plugin readme)
        if (m_param.has_value())
          {
            double A1 = m_param.value().at(0);
            double b1 = m_param.value().at(1);
            double A2 = m_param.value().at(2);
            double b2 = m_param.value().at(3);
            double Q0 = m_param.value().at(4);
            double corr=m_param.value().at(5);
            double xs = (2 * NC::kPi / (k * k)) * (A1 / (-b1 + 2) * std::pow(Q0, -b1 + 2) + A2 / (-b2 + 2) * std::pow(2 * k, -b2 + 2) - A2 / (-b2 + 2) * std::pow(Q0, -b2 + 2));
            nc_assert_always(k != 0);
            nc_assert_always(xs != 0);
            double ratio_sigma = (2 * NC::kPi / (k * k)) / xs; // cross section over total cross section ratio
            double CDF_Q0 = (A1 * std::pow(Q0, -b1 + 2) / (-b1 + 2)) * ratio_sigma;
            if (rand < CDF_Q0)
              {
                Q = std::pow(((-b1 + 2) * rand / A1) / ratio_sigma, 1 / (-b1 + 2));
              }
            else
              {
                Q = std::pow((rand / ratio_sigma - (A1 / (-b1 + 2)) * std::pow(Q0, -b1 + 2) + (A2 / (-b2 + 2)) * std::pow(Q0, -b2 + 2)) * (-b2 + 2) / A2, 1 / (-b2 + 2));
              }
          }
        else
          {
            NCRYSTAL_THROW2(LogicError, "Attempt to use not-initialized parameter in xs sampling in the " << pluginNameUpperCase()
                            << " plugin");
          }
        break;
      }
    case Model::FILE:
    case Model::HSFBA:
    case Model::GPF:
      {
        if (m_helper.has_value())
          {
            NC::NeutronEnergy ekin(neutron_ekin);
            Q = m_helper.value().sampleQValue(rng, ekin);
          }
        else
          {
            NCRYSTAL_THROW2(LogicError, "Attempt to use not-initialized IofQHelper in Q sampling in the " << pluginNameUpperCase()
                            << " plugin");
          }
        break;
      }
    default:
      NCRYSTAL_THROW2(LogicError, "Bad internal state for model in SANSND plugin");
    }
  return Q;
}
NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent(NC::RNG &rng, double neutron_ekin) const
{
  ScatEvent result;

  // section for energy limits
  /*if ( ! (neutron_ekin > 1) ) {
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
    }*/

  // Implement our actual model here:
  result.ekin_final = neutron_ekin; // Elastic
  double ksquared = NC::k4PiSq * NC::ekin2wlsqinv(neutron_ekin);
  if ( ! (ksquared>0.0) ) {
    result.mu = 1;
    return result;
  }
  double Q = sampleScatteringVector(rng, neutron_ekin);
  result.mu = NC::ncclamp( 1 - 0.5 * (Q * Q / ksquared), -1.0, 1.0 );
  nc_assert_always(!std::isnan(result.mu));

  return result;
}
