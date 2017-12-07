/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "tools/Grid.h"
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "GridIntegrationWeights.h"

#include "tools/Keywords.h"

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST MATHEVAL_DIST
/*
Target distribution given by a matheval parsed function (static or dynamic).

Use as a target distribution the distribution defined by
\f[
p(\mathbf{s}) =
\frac{f(\mathbf{s})}{\int d\mathbf{s} \, f(\mathbf{s})}
\f]
where \f$f(\mathbf{s})\f$ is some arbitrary mathematical function that
is parsed by the matheval library.

The function \f$f(\mathbf{s})\f$ is given by the FUNCTION keywords by
using _s1_,_s2_,..., as variables for the arguments
\f$\mathbf{s}=(s_1,s_2,\ldots,s_d)\f$.
If one variable is not given the target distribution will be
taken as uniform in that argument.

It is also possible to include the free energy surface \f$F(\mathbf{s})\f$
in the target distribution by using the _FE_ variable. In this case the
target distribution is dynamic and needs to be updated with current
best estimate of \f$F(\mathbf{s})\f$, similarly as for the
\ref WELL_TEMPERED "well-tempered target distribution".
Furthermore, the inverse temperature \f$\beta = (k_{\mathrm{B}}T)^{-1}\f$ and
the thermal energy \f$k_{\mathrm{B}}T\f$ can be included
by using the _beta_ and _kBT_ variables.

The target distribution will be automatically normalized over the region on
which it is defined on. Therefore, the function given in
FUNCTION needs to be non-negative and normalizable. The
code will perform checks to make sure that this is indeed the case.


\attention
The MATHEVAL_DIST only works if libmatheval is installed on the system and
PLUMED has been linked to it.

\par Examples

Here we use as shifted [Maxwell-Boltzmann distribution](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution)
as a target distribution in one-dimension.
Note that it is not need to include the normalization factor as the distribution will be
automatically normalized.
\verbatim
TARGET_DISTRIBUTION={MATHEVAL_DIST
                     FUNCTION=(s1+20)^2*exp(-(s1+20)^2/(2*10.0^2))}
\endverbatim

Here we have a two dimensional target distribution where we
use a [generalized normal distribution](https://en.wikipedia.org/wiki/Generalized_normal_distribution)
for argument \f$s_2\f$ while the distribution for \f$s_1\f$ is taken as
uniform as the variable _s1_ is not included in the function.
\verbatim
TARGET_DISTRIBUTION={MATHEVAL_DIST
                     FUNCTION=exp(-(abs(s2-20.0)/5.0)^4.0)}
\endverbatim

By using the _FE_ variable the target distribution can depend on
the free energy surface \f$F(\mathbf{s})\f$. For example,
the following input is identical to using \ref WELL_TEMPERED with
BIASFACTOR=10.
\verbatim
TARGET_DISTRIBUTION={MATHEVAL_DIST
                    FUNCTION=exp(-(beta/10.0)*FE)}
\endverbatim
Here the inverse temperature is automatically obtained by using the _beta_
variable. It is also possible to use the _kBT_ variable. The following
syntax will give the exact same results as the syntax above
\verbatim
TARGET_DISTRIBUTION={MATHEVAL_DIST
                    FUNCTION=exp(-(1.0/(kBT*10.0))*FE)}
\endverbatim










*/
//+ENDPLUMEDOC

class TD_Matheval : public TargetDistribution {
private:
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  //
  void* evaluator_pntr_;
  //
  std::vector<unsigned int> cv_var_idx_;
  std::vector<std::string> cv_var_str_;
  //
  std::string cv_var_prefix_str_;
  std::string fes_var_str_;
  std::string kbt_var_str_;
  std::string beta_var_str_;
  //
  bool use_fes_;
  bool use_kbt_;
  bool use_beta_;
public:
  static void registerKeywords( Keywords&);
  explicit TD_Matheval( const TargetDistributionOptions& to );
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_Matheval();
};

#ifdef __PLUMED_HAS_MATHEVAL
VES_REGISTER_TARGET_DISTRIBUTION(TD_Matheval,"MATHEVAL_DIST")


void TD_Matheval::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","FUNCTION","The function you wish to use for the target distribution where you should use the variables _s1_,_s2_,... for the arguments. You can also use the current estimate of the FES by using the variable _FE_ and the temperature by using the _kBT_ and _beta_ variables.");
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
}


TD_Matheval::~TD_Matheval(){
  evaluator_destroy(evaluator_pntr_);
}



TD_Matheval::TD_Matheval(const TargetDistributionOptions& to):
TargetDistribution(to),
evaluator_pntr_(NULL),
//
cv_var_idx_(0),
cv_var_str_(0),
//
cv_var_prefix_str_("s"),
fes_var_str_("FE"),
kbt_var_str_("kBT"),
beta_var_str_("beta"),
//
use_fes_(false),
use_kbt_(false),
use_beta_(false)
{
  std::string func_str;
  parse("FUNCTION",func_str);
  checkRead();
  //
  evaluator_pntr_=evaluator_create(const_cast<char*>(func_str.c_str()));
  if(evaluator_pntr_==NULL) plumed_merror(getName()+": there was some problem in parsing matheval formula "+func_str);
  //
  char** var_names;
  int var_count;
  evaluator_get_variables(evaluator_pntr_,&var_names,&var_count);
  //
  for(int i=0; i<var_count; i++){
    std::string curr_var = var_names[i];
    unsigned int cv_idx;
    if(curr_var.substr(0,cv_var_prefix_str_.size())==cv_var_prefix_str_ && Tools::convert(curr_var.substr(cv_var_prefix_str_.size()),cv_idx) && cv_idx>0){
      cv_var_idx_.push_back(cv_idx-1);
    }
    else if(curr_var==fes_var_str_){
      use_fes_=true;
      setDynamic();
      setFesGridNeeded();
    }
    else if(curr_var==kbt_var_str_){
      use_kbt_=true;
    }
    else if(curr_var==beta_var_str_){
      use_beta_=true;
    }
    else {
      plumed_merror(getName()+": problem with parsing matheval formula, cannot recognise the variable "+curr_var);
    }
  }
  //
  std::sort(cv_var_idx_.begin(),cv_var_idx_.end());
  cv_var_str_.resize(cv_var_idx_.size());
  for(unsigned int j=0; j<cv_var_idx_.size(); j++){
    std::string str1; Tools::convert(cv_var_idx_[j]+1,str1);
    cv_var_str_[j] = cv_var_prefix_str_+str1;
  }
}


void TD_Matheval::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins){
  if(cv_var_idx_.size()>0 && cv_var_idx_[cv_var_idx_.size()-1]>getDimension()){
    plumed_merror(getName()+": mismatch between CVs given in FUNC and the dimension of the target distribution");
  }
}


double TD_Matheval::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_Matheval");
  return 0.0;
}


void TD_Matheval::updateGrid(){
  std::vector<char*> var_char(cv_var_str_.size());
  std::vector<double> var_values(cv_var_str_.size());
  for(unsigned int j=0; j<cv_var_str_.size(); j++){
    var_char[j] = const_cast<char*>(cv_var_str_[j].c_str());
  }
  if(use_fes_){
    plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to the free energy in the target distribution");
    var_char.push_back(const_cast<char*>(fes_var_str_.c_str()));
    var_values.push_back(0.0);
  }
  if(use_kbt_){
    var_char.push_back(const_cast<char*>(kbt_var_str_.c_str()));
    var_values.push_back(1.0/getBeta());
  }
  if(use_beta_){
    var_char.push_back(const_cast<char*>(beta_var_str_.c_str()));
    var_values.push_back(getBeta());
  }
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  //
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    std::vector<double> point = targetDistGrid().getPoint(l);
    for(unsigned int k=0; k<cv_var_idx_.size() ; k++){
      var_values[k] = point[cv_var_idx_[k]];
    }
    if(use_fes_){
      var_values[cv_var_idx_.size()] = getFesGridPntr()->getValue(l);
    }
    double value = evaluator_evaluate(evaluator_pntr_,var_char.size(),&var_char[0],&var_values[0]);

    if(value<0.0 && !isTargetDistGridShiftedToZero()){plumed_merror(getName()+": The target distribution function gives negative values. You should change the definition of the function used for the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");}
    targetDistGrid().setValue(l,value);
    norm += integration_weights[l]*value;
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  if(norm>0.0){
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  }
  else if(!isTargetDistGridShiftedToZero()){
    plumed_merror(getName()+": The target distribution function cannot be normalized proberly. You should change the definition of the function used for the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");
  }
  logTargetDistGrid().setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
    std::vector<char*> rw_var_char(cv_var_str_.size());
    std::vector<double> rw_var_values(cv_var_str_.size());
    for(unsigned int j=0; j<cv_var_str_.size(); j++){
      rw_var_char[j] = const_cast<char*>(cv_var_str_[j].c_str());
    }
    if(use_fes_){
      plumed_massert(getFesRWGridPntr()!=NULL,"the FES reweight grid has to be linked to the free energy in the target distribution");
      rw_var_char.push_back(const_cast<char*>(fes_var_str_.c_str()));
      rw_var_values.push_back(0.0);
    }
    if(use_kbt_){
      rw_var_char.push_back(const_cast<char*>(kbt_var_str_.c_str()));
      rw_var_values.push_back(1.0/getBeta());
    }
    if(use_beta_){
      rw_var_char.push_back(const_cast<char*>(beta_var_str_.c_str()));
      rw_var_values.push_back(getBeta());
    }
    std::vector<double> rw_integration_weights = GridIntegrationWeights::getIntegrationWeights(getReweightGridPntr());
    norm = 0.0;
    //
    for(Grid::index_t l=0; l<reweightGrid().getSize(); l++){
      std::vector<double> point = reweightGrid().getPoint(l);
      for(unsigned int k=0; k<cv_var_idx_.size() ; k++){
        rw_var_values[k] = point[cv_var_idx_[k]];
      }
      if(use_fes_){
        rw_var_values[cv_var_idx_.size()] = getFesRWGridPntr()->getValue(l);
      }
      double value = evaluator_evaluate(evaluator_pntr_,rw_var_char.size(),&rw_var_char[0],&rw_var_values[0]);

      if(value<0.0 && !isTargetDistGridShiftedToZero()){plumed_merror(getName()+": The reweight grid function gives negative values. You should change the definition of the function used for the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");}
      reweightGrid().setValue(l,value);
      norm += rw_integration_weights[l]*value;
      logReweightGrid().setValue(l,-std::log(value));
    }
    if(norm>0.0){
      reweightGrid().scaleAllValuesAndDerivatives(1.0/norm);
    }
    else if(!isTargetDistGridShiftedToZero()){
      plumed_merror(getName()+": The reweight grid function cannot be normalized proberly. You should change the definition of the function used for the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");
    }
    logReweightGrid().setMinToZero();
  }
}


#endif


}
}
