#ifndef NCrystal_LookUpTable_hh
#define NCrystal_LookUpTable_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <functional>
#include <vector>

#include "NCrystal/internal/NCMath.hh"


namespace NCPluginNamespace {
  //A table for linear interplation of a curve.
  //Integration of curve can also be applied.
  class LookUpTable {
  public:
    //extrapolation methods below lower and beyond upper bounds
    enum Extrapolate {
      kConst_Zero,
      kZero_Zero,
      kZero_Const,
      kOverX_Zero,
      kOverSqrtX_Zero,
      kOverSqrtX_OverSqrtX,
      kConst_OverSqrtX,
    };
  public:
    LookUpTable(const std::vector<double>& x, const std::vector<double>& f, Extrapolate extrap=Extrapolate::kConst_Zero);
    LookUpTable();
    virtual ~LookUpTable();
    bool empty() const;
    void sanityCheck() const;
    void set(const std::vector<double>& x, const std::vector<double>& f, Extrapolate extrap=Extrapolate::kConst_Zero);
    double get(double x) const;
    void print() const;
    virtual double integrate(double lower_x, double upper_x);

  private:
    void init();
    std::vector<double> m_x, m_f;
    std::function<double(double)>  m_func_extrapLower, m_func_extrapUpper;

    double extrapZero(double );
    double extrapConstUpper(double );
    double extrapConstLower(double );
    double extrapOverSqrtXLower(double x);
    double extrapOverSqrtXUpper(double x);
    double extrapOverXLower(double x);
    double extrapOverXUpper(double x);
  };
  #include "NCLookUpTable.icc"
}

#endif
