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

#include "NCrystal/NCPluginBoilerplate.hh"

#include "NCLookUpTable.hh"
#include  <algorithm>
namespace NC = NCrystal;

NCP::LookUpTable::LookUpTable() = default;
NCP::LookUpTable::~LookUpTable() = default;

NCP::LookUpTable::LookUpTable(const std::vector<double>& x, const std::vector<double>& f, Extrapolate extrap)
{
  set(x,f,extrap);
}

void NCP::LookUpTable::set(const std::vector<double>& x, const std::vector<double>& f, Extrapolate extrap)
{
  switch(extrap) {
    case kConst_Zero:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapConstLower, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      break;
    case kZero_Zero:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      break;
    case kZero_Const:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapConstUpper, this, std::placeholders::_1);
      break;
    case kOverX_Zero:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapOverXLower, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      break;
    case kOverSqrtX_Zero:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapOverSqrtXLower, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapZero, this, std::placeholders::_1);
      break;
    case kOverSqrtX_OverSqrtX:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapOverSqrtXLower, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapOverSqrtXUpper, this, std::placeholders::_1);
      break;
    case kConst_OverSqrtX:
      m_func_extrapLower = std::bind( &NCP::LookUpTable::extrapConstLower, this, std::placeholders::_1);
      m_func_extrapUpper = std::bind( &NCP::LookUpTable::extrapOverSqrtXUpper, this, std::placeholders::_1);
      break;
    default :
      NCRYSTAL_THROW(CalcError, "extrapolation functions are not defined");
  }

  m_x=x;
  m_f=f;
  sanityCheck();
}

void NCP::LookUpTable::sanityCheck() const
{
  if(m_x.size()!=m_f.size())
    NCRYSTAL_THROW(CalcError, "x and f have different size");

  if(m_x.empty())
    NCRYSTAL_THROW(CalcError, "empty input array");

  if(!std::is_sorted(std::begin(m_x), std::end(m_x)))
    NCRYSTAL_THROW(CalcError, "x is not sorted");
}

void NCP::LookUpTable::print() const
{
  sanityCheck();
  printf("Look-up table content:\n");
  for(unsigned i=0;i<m_x.size();i++)
  {
    printf("%e %e\n", m_x[i], m_f[i]);
  }
}

bool NCP::LookUpTable::empty() const
{
  return m_x.empty() && m_f.empty();
}

double NCP::LookUpTable::integrate(double , double )
{
  NCRYSTAL_THROW(CalcError, "NCP::LookUpTable::integral to be implemented");
  return -1;
}

double NCP::LookUpTable::get(double x) const
{
  nc_assert(!NC::ncisnan(x));

  auto it = std::lower_bound(m_x.begin(),m_x.end(),x);
  if (it == m_x.begin())
    return m_func_extrapLower(x);
  else if( it == m_x.end())
    return m_func_extrapUpper(x);

  size_t idx = it-m_x.begin()-1;
  double result = m_f[idx] + (x-m_x[idx])*(m_f[idx+1]-m_f[idx])/(m_x[idx+1]-m_x[idx]);
  return result;
}
