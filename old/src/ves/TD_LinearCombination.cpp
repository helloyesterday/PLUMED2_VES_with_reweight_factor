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


namespace PLMD{

// class Grid;
class Action;

namespace ves{

//+PLUMEDOC VES_TARGETDIST LINEAR_COMBINATION
/*
Target distribution given by linear combination of distributions (static or dynamic).

Employ a target distribution that is a linear combination of the other
distributions, defined as
\f[
p(\mathbf{s}) = \sum_{i} w_{i} \, p_{i}(\mathbf{s})
\f]
where the weights \f$w_{i}\f$ are normalized to 1, \f$\sum_{i}w_{i}=1\f$.

The distributions \f$p_{i}(\mathbf{s})\f$ are given by using a separate numbered
DISTRIBUTION keyword for each distribution. The keywords for each distribution
should be enclosed within curly brackets.

The weights \f$w_{i}\f$ can be given using
the WEIGHTS keyword. The distributions are weighted equally if no weights are given.

It is assumed that all the distributions \f$p_{i}(\mathbf{s})\f$ given with the numbered
DISTRIBUTION keywords are normalized. If that is not the case you should
normalize each distribution separately by using the NORMALIZE
keyword within the curly brackets of each separate DISTRIBUTION keyword.
Note that normalizing the overall
linear combination will generally lead to different results than normalizing
each distribution separately.

The linear combination will be a dynamic target distribution if one or more
of the distributions used is a dynamic distribution. Otherwise it will be a
static distribution.

\par Examples

Here we employ a linear combination of a uniform and a Gaussian distribution.
No weights are given so the two distributions will be weighted equally.
\verbatim
TARGET_DISTRIBUTION={LINEAR_COMBINATION
                     DISTRIBUTION1={UNIFORM}
                     DISTRIBUTION2={GAUSSIAN
                                    CENTER=-2.0
                                    SIGMA=0.5}}
\endverbatim

Here we employ a linear combination of a uniform and two Gaussian distribution.
The weights are automatically normalized to 1 such that giving
WEIGHTS=1.0,1.0,2.0 as we do here is equal to giving WEIGHTS=0.25,0.25,0.50.
\verbatim
TARGET_DISTRIBUTION={LINEAR_COMBINATION
                     DISTRIBUTION1={UNIFORM}
                     DISTRIBUTION2={GAUSSIAN
                                    CENTER=-2.0,-2.0
                                    SIGMA=0.5,0.3}
                     DISTRIBUTION3={GAUSSIAN
                                    CENTER=+2.0,+2.0
                                    SIGMA=0.3,0.5}
                     WEIGHTS=1.0,1.0,2.0}
\endverbatim

In the above example the two Gaussians are given using two separate
DISTRIBUTION keywords. As the \ref GAUSSIAN target distribution allows multiple
centers is it also possible to use just one DISTRIBUTION keyword for the two
Gaussians. This is shown in the following example which will give the
exact same result as the one above as the weights have been appropriately
adjusted
\verbatim
TARGET_DISTRIBUTION={LINEAR_COMBINATION
                     DISTRIBUTION1={UNIFORM}
                     DISTRIBUTION2={GAUSSIAN
                                    CENTER1=-2.0,-2.0
                                    SIGMA1=0.5,0.3
                                    CENTER2=+2.0,+2.0
                                    SIGMA2=0.3,0.5
                                    WEIGHTS=1.0,2.0}
                     WEIGHTS=0.25,0.75}
\endverbatim

*/
//+ENDPLUMEDOC

class VesBias;

class TD_LinearCombination: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  // Added by Y. Isaac Yang to calculate the reweighting factor
  std::vector<Grid*> rw_grid_pntrs_;
  //
  std::vector<double> weights_;
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  // Added by Y. Isaac Yang to calculate the reweighting factor
  void setupAdditionalReweightGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
public:
  static void registerKeywords(Keywords&);
  explicit TD_LinearCombination(const TargetDistributionOptions& to);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_LinearCombination();
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


VES_REGISTER_TARGET_DISTRIBUTION(TD_LinearCombination,"LINEAR_COMBINATION")


void TD_LinearCombination::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","The target distributions to be used in the linear combination, each given within a separate numbered DISTRIBUTION keyword and enclosed in curly brackets {}.");
  keys.add("optional","WEIGHTS","The weights of target distributions. Have to be as many as the number of target distributions given with the numbered DISTRIBUTION keywords. If no weights are given the distributions are weighted equally. The weights are automatically normalized to 1.");
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  //keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_LinearCombination::TD_LinearCombination( const TargetDistributionOptions& to ):
TargetDistribution(to),
distribution_pntrs_(0),
grid_pntrs_(0),
rw_grid_pntrs_(0),
weights_(0),
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
  if(ndist_==0){plumed_merror(getName()+ ": no distributions are given.");}
  if(ndist_==1){plumed_merror(getName()+ ": giving only one distribution does not make sense.");}
  //
  if(!parseVector("WEIGHTS",weights_,true)){weights_.assign(distribution_pntrs_.size(),1.0);}
  if(distribution_pntrs_.size()!=weights_.size()){
    plumed_merror(getName()+ ": there has to be as many weights given in WEIGHTS as numbered DISTRIBUTION keywords");
  }
  //
  double sum_weights=0.0;
  for(unsigned int i=0;i<weights_.size();i++){sum_weights+=weights_[i];}
  for(unsigned int i=0;i<weights_.size();i++){weights_[i]/=sum_weights;}
  checkRead();
}


TD_LinearCombination::~TD_LinearCombination(){
  for(unsigned int i=0; i<ndist_; i++){
    delete distribution_pntrs_[i];
  }
}


double TD_LinearCombination::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_LinearCombination");
  return 0.0;
}


void TD_LinearCombination::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
  }
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TD_LinearCombination::setupAdditionalReweightGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupReweightGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    rw_grid_pntrs_[i]=distribution_pntrs_[i]->getReweightGridPntr();
  }
}
//


void TD_LinearCombination::updateGrid(){
  for(unsigned int i=0; i<ndist_; i++){
	// Added by Y. Isaac Yang to calculate the reweighting factor
	if(isReweightGridActive())
	  distribution_pntrs_[i]->setReweightGridActive();
	//
    distribution_pntrs_[i]->update();
  }
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    double value = 0.0;
    for(unsigned int i=0; i<ndist_; i++){
      value += weights_[i]*grid_pntrs_[i]->getValue(l);
    }
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  logTargetDistGrid().setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
	for(Grid::index_t l=0; l<reweightGrid().getSize(); l++){
      double value = 0.0;
      for(unsigned int i=0; i<ndist_; i++){
        value += weights_[i]*rw_grid_pntrs_[i]->getValue(l);
      }
      reweightGrid().setValue(l,value);
      logReweightGrid().setValue(l,-std::log(value));
    }
    logReweightGrid().setMinToZero();
  }
  //
}


void TD_LinearCombination::linkVesBias(VesBias* vesbias_pntr_in){
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void TD_LinearCombination::linkAction(Action* action_pntr_in){
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void TD_LinearCombination::linkBiasGrid(Grid* bias_grid_pntr_in){
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void TD_LinearCombination::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in){
  TargetDistribution::linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  }
}


void TD_LinearCombination::linkFesGrid(Grid* fes_grid_pntr_in){
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TD_LinearCombination::linkBiasRWGrid(Grid* bias_rwgrid_pntr_in){
  TargetDistribution::linkBiasRWGrid(bias_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasRWGrid(bias_rwgrid_pntr_in);
  }
}
void TD_LinearCombination::linkBiasWithoutCutoffRWGrid(Grid* bias_withoutcutoff_rwgrid_pntr_in){
  TargetDistribution::linkBiasWithoutCutoffRWGrid(bias_withoutcutoff_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasWithoutCutoffRWGrid(bias_withoutcutoff_rwgrid_pntr_in);
  }
}
void TD_LinearCombination::linkFesRWGrid(Grid* fes_rwgrid_pntr_in){
  TargetDistribution::linkFesRWGrid(fes_rwgrid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesRWGrid(fes_rwgrid_pntr_in);
  }
}


}
}
