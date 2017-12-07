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

#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "GridIntegrationWeights.h"

#include "tools/Keywords.h"
#include "tools/Grid.h"


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST WELL_TEMPERED
/*
Well-tempered target distribution (dynamic).

Use as a target distribution the well-tempered distribution \cite Barducci:2008
given by
\f[
p(\mathbf{s}) =
\frac{e^{-(\beta/\gamma) F(\mathbf{s})}}
{\int d\mathbf{s}\, e^{-(\beta/\gamma) F(\mathbf{s})}} =
\frac{[P_{0}(\mathbf{s})]^{1/\gamma}}
{\int d\mathbf{s}\, [P_{0}(\mathbf{s})]^{1/\gamma}}
\f]
where \f$\gamma\f$ is a so-called bias factor and \f$P_{0}(\mathbf{s})\f$ is the
unbiased canonical distribution of the CVs. This target distribution thus
correponds to a biased ensemble where, as compared to the unbiased one,
the probability peaks have been broaden and the fluctations of the CVs are
enhanced.
The value of the bias factor \f$\gamma\f$ determines by how much the fluctations
are enhanced.

The well-tempered distribution can be view as sampling on
an effective free energy surface \f$\tilde{F}(\mathbf{s}) = (1/\gamma) F(\mathbf{s})\f$
which has largely the same metastable states as the original \f$F(\mathbf{s})\f$
but with barriers that have been reduced by a factor of \f$\gamma\f$.
Generally one should use a value of \f$\gamma\f$ that results in
effective barriers on the order of few \f$k_{\mathrm{B}}T\f$
such that thermal fluctuations can easily induce transitions
between different metastable states.

At convergence the relationship between the bias potential and the free
energy surface is given by
\f[
F(\mathbf{s}) = - \left(\frac{1}{1-\gamma^{-1}} \right) V(\mathbf{s})
\f]

This target distribution depends directly on the free energy surface
\f$F(\mathbf{s})\f$ which is quantity that we do not know a-priori and
want to obtain. Therefore, this target distribution
is iteratively updated \cite Valsson-JCTC-2015 according to
\f[
p^{(m+1)}(\mathbf{s}) =
\frac{e^{-(\beta/\gamma) F^{(m+1)}(\mathbf{s})}}
{\int d\mathbf{s}\, e^{-(\beta/\gamma) F^{(m+1)}(\mathbf{s})}}
\f]
where \f$F^{(m+1)}(\mathbf{s})\f$ is the current best estimate of the
free energy surface obtained according to
\f[
F^{(m+1)}(\mathbf{s}) =
- V^{(m+1)}(\mathbf{s}) - \frac{1}{\beta} \log p^{(m)}(\mathbf{s}) =
- V^{(m+1)}(\mathbf{s}) + \frac{1}{\gamma} F^{(m)}(\mathbf{s})
\f]
The frequency of performing this update needs to be set in the
optimizer used in the calculation. Normally it is sufficient
to do it every 100-1000 bias update iterations.

\par Examples
Employ a well-tempered target distribution with a bias factor of 10
\verbatim
TARGET_DISTRIBUTION={WELL_TEMPERED BIASFACTOR=10}
\endverbatim

*/
//+ENDPLUMEDOC

class TD_WellTempered: public TargetDistribution {
private:
  double bias_factor_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_WellTempered(const TargetDistributionOptions& to);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_WellTempered(){}
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_WellTempered,"WELL_TEMPERED")


void TD_WellTempered::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","BIASFACTOR","The bias factor used for the well-tempered distribution.");
  keys.use("BIAS_CUTOFF");
}


TD_WellTempered::TD_WellTempered( const TargetDistributionOptions& to ):
TargetDistribution(to),
bias_factor_(0.0)
{
  parse("BIASFACTOR",bias_factor_);
  if(bias_factor_<=1.0){
    plumed_merror("WELL_TEMPERED target distribution: the value of the bias factor doesn't make sense, it should be larger than 1.0");
  }
  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_WellTempered::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_WellTempered");
  return 0.0;
}


void TD_WellTempered::updateGrid(){
  double beta_prime = getBeta()/bias_factor_;
  plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_WellTempered!");
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    double value = beta_prime * getFesGridPntr()->getValue(l);
    logTargetDistGrid().setValue(l,value);
    value = exp(-value);
    norm += integration_weights[l]*value;
    targetDistGrid().setValue(l,value);
  }
  targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  logTargetDistGrid().setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
    std::vector<double> rw_integration_weights = GridIntegrationWeights::getIntegrationWeights(getReweightGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<reweightGrid().getSize(); l++){
      double value = beta_prime * getFesRWGridPntr()->getValue(l);
      logReweightGrid().setValue(l,value);
      value = exp(-value);
      norm += rw_integration_weights[l]*value;
      reweightGrid().setValue(l,value);
    }
    reweightGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logReweightGrid().setMinToZero();
  }
}


}
}
