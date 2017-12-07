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

#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "GridIntegrationWeights.h"


namespace PLMD{

// class Grid;
class Action;

namespace ves{

//+PLUMEDOC VES_TARGETDIST PRODUCT_COMBINATION
/*
Target distribution given by product combination of distributions (static or dynamic).

Employ a target distribution that is a product combination of the other
distributions, defined as
\f[
p(\mathbf{s}) =
\frac{\prod_{i} p_{i}(\mathbf{s})}
{\int d \mathbf{s} \prod_{i} p_{i}(\mathbf{s})}
\f]
where the distributions \f$p_{i}(\mathbf{s})\f$ are in full dimensional space
of the arguments used.

Note the difference between this target distribution and the one defined in
\ref PRODUCT_DISTRIBUTION. Here we have a non-separable distribution given
as a product of distribution \f$p_{i}(\mathbf{s})\f$ which are in full dimensional
space of the arguments used.

The distributions \f$p_{i}(\mathbf{s})\f$ are given by using a separate numbered
DISTRIBUTION keyword for each distribution. The keywords for each distribution
should be enclosed within curly brackets.

The target distribution resulting from the product combination will be
automatically normalized. Therefore, the product combination needs to
be a proper distribution that is non-negative and normalizable. The
code will perform checks to make sure that this is indeed the case.

The product combination will be a dynamic target distribution if one or more
of the distributions used is a dynamic distribution. Otherwise it will be a
static distribution.

\par Examples

In the following example the overall interval on which the
target distribution is defined is from -4 to 4.
We employ a product of a Gaussian distribution with two centers
and distribution that is uniform on the interval -3 to 3 and
then smoothly decays to zero outside that interval.
The overall effect will then be to cut off the tails of the
Gaussian distribution
\verbatim
TARGET_DISTRIBUTION={PRODUCT_COMBINATION
                     DISTRIBUTION1={GAUSSIAN
                                   CENTER1=-2.9 SIGMA1=1.0
                                   CENTER2=+2.9 SIGMA2=0.4}
                     DISTRIBUTION2={UNIFORM
                                    MINIMA=-3.0 SIGMA_MINIMA=0.20
                                    MAXIMA=+3.0 SIGMA_MAXIMA=0.15}}
\endverbatim

*/
//+ENDPLUMEDOC

class VesBias;

class TD_ProductCombination: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  // Added by Y. Isaac Yang to calculate the reweighting factor
  std::vector<Grid*> rw_grid_pntrs_;
  //
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  // Added by Y. Isaac Yang to calculate the reweighting factor
  void setupAdditionalReweightGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
public:
  static void registerKeywords(Keywords&);
  explicit TD_ProductCombination(const TargetDistributionOptions& to);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_ProductCombination();
  //
  void linkVesBias(VesBias*);
  void linkAction(Action*);
  //
  void linkBiasGrid(Grid*);
  void linkBiasWithoutCutoffGrid(Grid*);
  void linkFesGrid(Grid*);
  // Added by Y. Isaac Yang to calculate the reweighting factor
  void linkBiasRWGrid(Grid*);
  void linkBiasWithoutCutoffRWGrid(Grid*);
  void linkFesRWGrid(Grid*);
  //
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_ProductCombination,"PRODUCT_COMBINATION")


void TD_ProductCombination::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","The target distributions to be used in the product combination, each given within a separate numbered DISTRIBUTION keyword and enclosed in curly brackets {}.");
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
}


TD_ProductCombination::TD_ProductCombination( const TargetDistributionOptions& to ):
TargetDistribution(to),
distribution_pntrs_(0),
grid_pntrs_(0),
rw_grid_pntrs_(0),
ndist_(0)
{
  for(unsigned int i=1;; i++) {
    std::string keywords;
    if(!parseNumbered("DISTRIBUTION",i,keywords) ){break;}
    std::vector<std::string> words = Tools::getWords(keywords);
    TargetDistribution* dist_pntr_tmp = targetDistributionRegister().create( (words) );
    if(dist_pntr_tmp->isDynamic()){setDynamic();}
    if(dist_pntr_tmp->fesGridNeeded()){setFesGridNeeded();}
    if(dist_pntr_tmp->biasGridNeeded()){setBiasGridNeeded();}
    distribution_pntrs_.push_back(dist_pntr_tmp);
  }
  ndist_ = distribution_pntrs_.size();
  grid_pntrs_.assign(ndist_,NULL);
  // Added by Y. Isaac Yang to calculate the reweighting factor
  rw_grid_pntrs_.assign(ndist_,NULL);
  //
  checkRead();
}


TD_ProductCombination::~TD_ProductCombination(){
  for(unsigned int i=0; i<ndist_; i++){
    delete distribution_pntrs_[i];
  }
}


double TD_ProductCombination::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_ProductCombination");
  return 0.0;
}


void TD_ProductCombination::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
  }
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TD_ProductCombination::setupAdditionalReweightGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupReweightGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    rw_grid_pntrs_[i]=distribution_pntrs_[i]->getReweightGridPntr();
  }
}

void TD_ProductCombination::updateGrid(){
  for(unsigned int i=0; i<ndist_; i++){
	// Added by Y. Isaac Yang to calculate the reweighting factor
	if(isReweightGridActive())
	  distribution_pntrs_[i]->setReweightGridActive();
	//
    distribution_pntrs_[i]->update();
  }
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    double value = 1.0;
    for(unsigned int i=0; i<ndist_; i++){
      value *= grid_pntrs_[i]->getValue(l);
    }
    if(value<0.0 && !isTargetDistGridShiftedToZero()){plumed_merror(getName()+": The target distribution function gives negative values. You should change the definition of the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");}
    norm += integration_weights[l]*value;
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }

  if(norm>0.0){
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  }
  else if(!isTargetDistGridShiftedToZero()){
    plumed_merror(getName()+": The target distribution function cannot be normalized proberly. You should change the definition of the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");
  }
  logTargetDistGrid().setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
	norm = 0.0;
    for(Grid::index_t l=0; l<reweightGrid().getSize(); l++){
      double value = 1.0;
      for(unsigned int i=0; i<ndist_; i++){
        value *= rw_grid_pntrs_[i]->getValue(l);
    }
    if(value<0.0 && !isTargetDistGridShiftedToZero()){plumed_merror(getName()+": The reweight grid function gives negative values. You should change the definition of the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");}
      norm += integration_weights[l]*value;
      reweightGrid().setValue(l,value);
      logReweightGrid().setValue(l,-std::log(value));
	}

	if(norm>0.0){
      reweightGrid().scaleAllValuesAndDerivatives(1.0/norm);
	}
	else if(!isTargetDistGridShiftedToZero()){
      plumed_merror(getName()+": The reweight grid function cannot be normalized proberly. You should change the definition of the target distribution to avoid this. You can also use the SHIFT_TO_ZERO keyword to avoid this problem.");
	}
  }
}


void TD_ProductCombination::linkVesBias(VesBias* vesbias_pntr_in){
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void TD_ProductCombination::linkAction(Action* action_pntr_in){
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void TD_ProductCombination::linkBiasGrid(Grid* bias_grid_pntr_in){
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void TD_ProductCombination::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in){
  TargetDistribution::linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  }
}


void TD_ProductCombination::linkFesGrid(Grid* fes_grid_pntr_in){
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TD_ProductCombination::linkBiasRWGrid(Grid* bias_rwgrid_pntr_in){
  TargetDistribution::linkBiasRWGrid(bias_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasRWGrid(bias_rwgrid_pntr_in);
  }
}
void TD_ProductCombination::linkBiasWithoutCutoffRWGrid(Grid* bias_withoutcutoff_rwgrid_pntr_in){
  TargetDistribution::linkBiasWithoutCutoffRWGrid(bias_withoutcutoff_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasWithoutCutoffRWGrid(bias_withoutcutoff_rwgrid_pntr_in);
  }
}
void TD_ProductCombination::linkFesRWGrid(Grid* fes_rwgrid_pntr_in){
  TargetDistribution::linkFesRWGrid(fes_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesRWGrid(fes_rwgrid_pntr_in);
  }
}


}
}
