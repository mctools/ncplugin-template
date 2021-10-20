#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <iostream>

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_SANSND section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, it 
  //raises BadInput exceptions
  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

 // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_SANSND
  //    v                   # plugin version
  //    filename            # read file mode
  //   
  //    or
  //
  //    v                   # plugin version
  //    model               # which model 
  //    parametrs           # paramters for the selected model
  //              

  bool file_mode=false;
  //Verify we have exactly 1-2 (file_mode on) or 1-1 input (fit mode) :
  if ( data.size() == 2 ) {
    if ( data.at(0).size()!=1 || data.at(1).size()!=1 ) {
    NCRYSTAL_THROW2(BadInput,"Bad input for file mode in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");
    } else { file_mode=true; }   
  } else if ( data.size() == 3 ) {
    if ( data.at(0).size()!=1 || data.at(1).size()!=1 )
    NCRYSTAL_THROW2(BadInput,"Bad input for fit mode in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).")          
  } else {
    NCRYSTAL_THROW2(BadInput,"Bad input in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");    
  }
  //Parse and validate values:
  double supp_version = 2.0;
  double version;
  //NC::VectD param;
  double p0,p1,p2,p3,p4;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), version ))
    NCRYSTAL_THROW2( BadInput,"Invalid version input in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");  
  //special warning for wrong version
  if ( ! (version==supp_version) ) 
    NCRYSTAL_THROW2( BadInput,"Invalid version specified for the "<<pluginNameUpperCase()
                      <<" plugin. Only the version "<<supp_version<<" is supported." );
  
  if (file_mode){
    std::string filename= data.at(1).at(0);
    std::string root_rel = "data/";
    std::string rel_path = root_rel+filename;
    std::cout << rel_path << std::endl;
    //check existence
    struct stat buffer;   
    if (!(stat (rel_path.c_str(), &buffer) == 0))
      NCRYSTAL_THROW2( BadInput,"The filename specified for the "<<pluginNameUpperCase()
                      <<" plugin is invalid or the file could not be found in the data/ directory. " );
    //Checks done! Create and return our model:
    return PhysicsModel(filename);
  } else {
      std::string model = data.at(1).at(0);
      if (model == "GP") {
        if ( ! NC::safe_str2dbl( data.at(2).at(0), p0 )
        || ! NC::safe_str2dbl( data.at(2).at(1), p1 )
        || ! NC::safe_str2dbl( data.at(2).at(2), p2 )
        || ! NC::safe_str2dbl( data.at(2).at(3), p3 )
         ) {
        NCRYSTAL_THROW2( BadInput,"Invalid values specified for " << model << " model in the @CUSTOM_"<<pluginNameUpperCase()
                        <<" section (see the plugin readme for more info)" ); 
                    } else {
              //CHECK THE INPUT PARAM
               
              //param.insert(param.end(), {A,s,rg,m});
            }
      } else if (model == "PPF") {
        if ( ! NC::safe_str2dbl( data.at(2).at(0), p0 )
            || ! NC::safe_str2dbl( data.at(2).at(1), p1 )
            || ! NC::safe_str2dbl( data.at(2).at(2), p2 )
            || ! NC::safe_str2dbl( data.at(2).at(3), p3 )
            || ! NC::safe_str2dbl( data.at(2).at(4), p4 )
            //|| ! NC::safe_str2dbl( data.at(2).at(5), sigma0 )
         ) {
        NCRYSTAL_THROW2( BadInput,"Invalid values specified for " << model << " model in the @CUSTOM_"<<pluginNameUpperCase()
                        <<" section (see the plugin readme for more info)" );
            } else {
              //CHECK THE INPUT PARAM
              //nc_assert(Q0>0);
              //nc_assert(sigma0>0);
              //param.insert(param.end(), {A1,b1,A2,b2,Q0,sigma0});
            }
      } else {
        NCRYSTAL_THROW2( BadInput,"Invalid model specified in the @CUSTOM_"<<pluginNameUpperCase()
                        <<" section (see the plugin readme for more info)" );
      }
      //Parsing done! Create and return our model:
      return PhysicsModel(model,p0,p1,p2,p3,p4);
  }
  
}

NCP::PhysicsModel::PhysicsModel( std::string model, double p0, double p1, double p2, double p3, double p4 )
  : m_model(model),
    m_param({p0,p1,p2,p3,p4}),
    m_helper(([this]() -> NC::IofQHelper { 
      
      //Generate vector of data q and IofQ
      NC::VectD q = NC::logspace(-4,1,10000);
      NC::VectD IofQ = q;
      if (m_model == "GP") {
        double A=m_param.at(0);
        double s=m_param.at(1);
        double rg=m_param.at(2);
        double m=m_param.at(3);
        //evaluate Q1 where IofQ stops being evaluated as Guinier and Porod starts
        double Q1 = 1.0/rg*std::sqrt((m-s)*(3-s)/2);
        //IofQ still filled with q here
        auto it_q1 = std::lower_bound(IofQ.begin(),IofQ.end(), Q1);
        //Approximation valid as long as we have high sampling. Otherwise interpolation needed
        if (it_q1==IofQ.end()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        double B = std::pow(Q1,m-s)*std::exp((-Q1*Q1*rg*rg)/(3-s));
        std::for_each(IofQ.begin(),it_q1,
                      [A,s,rg](double &x) { x = A*std::pow(x,-s)*std::exp((-x*x*rg*rg)/(3-s));}
                      ); 
        std::advance(it_q1,1);
        std::for_each(it_q1,IofQ.end(),
                      [A,B,m](double &x) { x = A*B*std::pow(x,-m);}
                      );
      } else if (m_model == "PPF") {
        double A1=m_param.at(0);
        double b1=m_param.at(1);
        double A2=m_param.at(2);
        double b2=m_param.at(3);
        double Q0=m_param.at(4);
        //double sigma0=m_param.at(5);

        auto it_q0 = std::lower_bound(IofQ.begin(),IofQ.end(), Q0);
        if (it_q0==IofQ.end()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q0 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        std::for_each(IofQ.begin(),it_q0,
                      [A1,b1](double &x) { x = A1*std::pow(x,-b1);}
                      ); 
        std::advance(it_q0,1);
        std::for_each(it_q0,IofQ.end(),
                      [A2,b2](double &x) { x = A2*std::pow(x,-b2);}
                      );    
      }     
      //Initialize the helper           
      NC::IofQHelper helper(q,IofQ);
      return helper;
    })())
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  
}

NCP::PhysicsModel::PhysicsModel(std::string filename)
  : m_model(),
    m_param(),
    m_helper(([filename]() -> NC::IofQHelper {
      
      //Parse the input file and create the vector with the q, Iq info
      NC::VectD q;
      NC::VectD IofQ;
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
      /*for(auto i :q){
        std::cout<<i<<std::endl;;
      }*/ 
      NC::IofQHelper helper(q,IofQ);
      return helper;
        })())
{   
  
};

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{
  NC::NeutronEnergy ekin(neutron_ekin);
  double k =  NC::k2Pi/ NC::ekin2wl(neutron_ekin); //wavevector
  double SANS_xs = 5.551/(2*k*k)*m_helper.calcQIofQIntegral(ekin);
  return SANS_xs;
}

double NCP::PhysicsModel::sampleScatteringVector( NC::RNG& rng, double neutron_ekin ) const 
{
  double rand = rng.generate();
  //sample a random scattering vector Q from the inverse CDF (see plugin readme)
  return rand;
}
NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  //section for energy limits
  /*if ( ! (neutron_ekin > 1) ) {
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
  }*/

  //Implement our actual model here:
  result.ekin_final = neutron_ekin;//Elastic
  double Q = sampleScatteringVector(rng, neutron_ekin);
  double ksquared = NC::k4PiSq*NC::ekin2wlsqinv(neutron_ekin);
  result.mu = 1-0.5*(Q*Q/ksquared);

  return result;
}


