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
#include "TargetDistModifer.h"
#include "TargetDistributionRegister.h"
#include "VesBias.h"
#include "GridIntegrationWeights.h"
#include "VesTools.h"

#include "core/Value.h"
#include "tools/Grid.h"
#include "tools/File.h"
#include "tools/Keywords.h"

#include "GridProjWeights.h"

#include <iostream>

namespace PLMD{
namespace ves{

TargetDistributionOptions::TargetDistributionOptions( const std::vector<std::string>& input):
words(input)
{}


void TargetDistribution::registerKeywords( Keywords& keys ){
  keys.reserve("hidden","BIAS_CUTOFF","Add a bias cutoff to the target distribution.");
  keys.reserve("optional","WELLTEMPERED_FACTOR","Broaden the target distribution such that it is taken as [p(s)]^(1/g) where g is the well tempered factor given here. If this option is active the distribution will be automatically normalized.");
  keys.reserveFlag("SHIFT_TO_ZERO",false,"Shift the minimum value of the target distribution to zero. This can for example be used to avoid negative values in the target distribution. If this option is active the distribution will be automatically normalized.");
  keys.reserveFlag("NORMALIZE",false,"Renormalized the target distribution over the intervals on which it is defined to make sure that it is properly normalized to 1. In most cases this should not be needed as the target distributions should be normalized. The code will issue a warning (but still run) if this is needed for some reason.");
}


TargetDistribution::TargetDistribution( const TargetDistributionOptions& to):
name_(to.words[0]),
input(to.words),
type_(static_targetdist),
force_normalization_(false),
check_normalization_(true),
check_nonnegative_(true),
shift_targetdist_to_zero_(false),
dimension_(0),
grid_args_(0),
targetdist_grid_pntr_(NULL),
log_targetdist_grid_pntr_(NULL),
reweight_grid_active_(false),
reweight_grid_pntr_(NULL),
log_reweight_grid_pntr_(NULL),
targetdist_modifer_pntrs_(0),
action_pntr_(NULL),
vesbias_pntr_(NULL),
needs_bias_grid_(false),
needs_bias_withoutcutoff_grid_(false),
needs_fes_grid_(false),
bias_grid_pntr_(NULL),
bias_withoutcutoff_grid_pntr_(NULL),
fes_grid_pntr_(NULL),
bias_rwgrid_pntr_(NULL),
bias_withoutcutoff_rwgrid_pntr_(NULL),
fes_rwgrid_pntr_(NULL),
static_grid_calculated(false),
bias_cutoff_active_(false),
bias_cutoff_value_(0.0),
keywords(targetDistributionRegister().getKeywords(name_))
{
  input.erase( input.begin() );
  //
  parse("BIAS_CUTOFF",bias_cutoff_value_,true);
  if(bias_cutoff_value_<0.0){
    plumed_merror(getName()+": negative value in BIAS_CUTOFF does not make sense");
  }
  if(bias_cutoff_value_>0.0){
    setupBiasCutoff();
  }
  //
  if(keywords.exists("WELLTEMPERED_FACTOR")){
    double welltempered_factor=0.0;
    parse("WELLTEMPERED_FACTOR",welltempered_factor,true);
    //
    if(welltempered_factor>0.0){
      if(bias_cutoff_active_){plumed_merror(getName()+": using WELLTEMPERED_FACTOR with bias cutoff is not allowed.");}
      TargetDistModifer* pntr = new WellTemperedModifer(welltempered_factor);
      targetdist_modifer_pntrs_.push_back(pntr);
    }
    else if(welltempered_factor<0.0){
      plumed_merror(getName()+": negative value in WELLTEMPERED_FACTOR does not make sense");
    }
  }
  //
  if(keywords.exists("SHIFT_TO_ZERO")){
    parseFlag("SHIFT_TO_ZERO",shift_targetdist_to_zero_);
    if(shift_targetdist_to_zero_){
      if(bias_cutoff_active_){plumed_merror(getName()+": using SHIFT_TO_ZERO with bias cutoff is not allowed.");}
      check_nonnegative_=false;
    }
  }
  //
  if(keywords.exists("NORMALIZE")){
    bool force_normalization=false;
    parseFlag("NORMALIZE",force_normalization);
    if(force_normalization){
      if(shift_targetdist_to_zero_){plumed_merror(getName()+": using NORMALIZE with SHIFT_TO_ZERO is not needed, the target distribution will be automatically normalized.");}
      if(bias_cutoff_active_){plumed_merror(getName()+": using NORMALIZE with bias cutoff is not allowed, the target distribution will be automatically normalized.");}
      setForcedNormalization();
    }
  }

}


TargetDistribution::~TargetDistribution() {
  if(targetdist_grid_pntr_!=NULL){
    delete targetdist_grid_pntr_;
  }
  if(log_targetdist_grid_pntr_!=NULL){
    delete log_targetdist_grid_pntr_;
  }
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(reweight_grid_pntr_!=NULL){
    delete reweight_grid_pntr_;
  }
  if(log_reweight_grid_pntr_!=NULL){
    delete log_reweight_grid_pntr_;
  }
  //
  for(unsigned int i=0; i<targetdist_modifer_pntrs_.size(); i++){
    delete targetdist_modifer_pntrs_[i];
  }
}


double TargetDistribution::getBeta() const {
  plumed_massert(vesbias_pntr_!=NULL,"The VesBias has to be linked to use TargetDistribution::getBeta()");
  return vesbias_pntr_->getBeta();
}


void TargetDistribution::setDimension(const unsigned int dimension){
  plumed_massert(dimension_==0,"setDimension: the dimension of the target distribution has already been set");
  dimension_=dimension;
}


void TargetDistribution::linkVesBias(VesBias* vesbias_pntr_in){
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void TargetDistribution::linkAction(Action* action_pntr_in){
  action_pntr_ = action_pntr_in;
}


void TargetDistribution::linkBiasGrid(Grid* bias_grid_pntr_in){
  bias_grid_pntr_ = bias_grid_pntr_in;
}


void TargetDistribution::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in){
  bias_withoutcutoff_grid_pntr_ = bias_withoutcutoff_grid_pntr_in;
}


void TargetDistribution::linkFesGrid(Grid* fes_grid_pntr_in){
  fes_grid_pntr_ = fes_grid_pntr_in;
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TargetDistribution::linkBiasRWGrid(Grid* bias_rwgrid_pntr_in){
  bias_rwgrid_pntr_ = bias_rwgrid_pntr_in;
}
void TargetDistribution::linkBiasWithoutCutoffRWGrid(Grid* bias_withoutcutoff_rwgrid_pntr_in){
  bias_withoutcutoff_rwgrid_pntr_ = bias_withoutcutoff_rwgrid_pntr_in;
}
void TargetDistribution::linkFesRWGrid(Grid* fes_rwgrid_pntr_in){
  fes_rwgrid_pntr_ = fes_rwgrid_pntr_in;
}
//

void TargetDistribution::setupBiasCutoff() {
  if(!keywords.exists("BIAS_CUTOFF")){
    plumed_merror(getName()+": this target distribution does not support a bias cutoff");
  }
  bias_cutoff_active_=true;
  setBiasWithoutCutoffGridNeeded();
  setDynamic();
  // as the p(s) includes the derivative factor so normalization
  // check can be misleading
  check_normalization_=false;
  force_normalization_=false;
}


void TargetDistribution::parseFlag(const std::string& key, bool& t) {
  Tools::parseFlag(input,key,t);
}


void TargetDistribution::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from the target distribution input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg);
  }
}


std::string TargetDistribution::description() {
  std::string str="Type: " + name_;
  return str;
}


void TargetDistribution::setupGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  if(getDimension()==0){
    setDimension(arguments.size());
  }
  unsigned int dimension = getDimension();
  plumed_massert(arguments.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(min.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(max.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(nbins.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  grid_args_=arguments;
  targetdist_grid_pntr_ =     new Grid("targetdist",arguments,min,max,nbins,false,false);
  log_targetdist_grid_pntr_ = new Grid("log_targetdist",arguments,min,max,nbins,false,false);
  setupAdditionalGrids(arguments,min,max,nbins);
}

// Added by Y. Isaac Yang to calculate the reweighting factor
void TargetDistribution::setupReweightGrids(const std::vector<Value*>& arguments,const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  unsigned int dimension = getDimension();
  plumed_massert(arguments.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(min.size()==dimension,"TargetDistribution::setupReweightGrids: mismatch between number of values given for grid parameters");
  plumed_massert(max.size()==dimension,"TargetDistribution::setupReweightGrids: mismatch between number of values given for grid parameters");
  plumed_massert(nbins.size()==dimension,"TargetDistribution::setupReweightGrids: mismatch between number of values given for grid parameters");
  reweight_grid_pntr_ =     new Grid("reweight",arguments,min,max,nbins,false,false);
  log_reweight_grid_pntr_ = new Grid("log_reweight",arguments,min,max,nbins,false,false);
  setReweightGridActive();
  setupAdditionalReweightGrids(arguments,min,max,nbins);
}
//

void TargetDistribution::calculateStaticDistributionGrid(){
  if(static_grid_calculated && !bias_cutoff_active_){return;}
  // plumed_massert(isStatic(),"this should only be used for static distributions");
  plumed_massert(targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  plumed_massert(log_targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++)
  {
   std::vector<double> argument = targetdist_grid_pntr_->getPoint(l);
   double value = getValue(argument);
   targetdist_grid_pntr_->setValue(l,value);
   log_targetdist_grid_pntr_->setValue(l,-std::log(value));
  }
  log_targetdist_grid_pntr_->setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
    plumed_massert(reweight_grid_pntr_!=NULL,"the grids have not been setup using setupReweightGrids");
    plumed_massert(log_reweight_grid_pntr_!=NULL,"the grids have not been setup using setupReweightGrids");
    for(Grid::index_t l=0; l<reweight_grid_pntr_->getSize(); l++)
    {
      std::vector<double> rw_argument = reweight_grid_pntr_->getPoint(l);
      double value = getValue(rw_argument);
      reweight_grid_pntr_->setValue(l,value);
      log_reweight_grid_pntr_->setValue(l,-std::log(value));
    }
    log_reweight_grid_pntr_->setMinToZero();
  }
  //
  static_grid_calculated = true;
}


double TargetDistribution::integrateGrid(const Grid* grid_pntr){
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(grid_pntr);
  double sum = 0.0;
  for(Grid::index_t l=0; l<grid_pntr->getSize(); l++){
    sum += integration_weights[l]*grid_pntr->getValue(l);
  }
  return sum;
}


double TargetDistribution::normalizeGrid(Grid* grid_pntr){
  double normalization = TargetDistribution::integrateGrid(grid_pntr);
  grid_pntr->scaleAllValuesAndDerivatives(1.0/normalization);
  return normalization;
}


Grid TargetDistribution::getMarginalDistributionGrid(Grid* grid_pntr, const std::vector<std::string>& args) {
  plumed_massert(grid_pntr->getDimension()>1,"doesn't make sense calculating the marginal distribution for a one-dimensional distribution");
  plumed_massert(args.size()<grid_pntr->getDimension(),"the number of arguments for the marginal distribution should be less than the dimension of the full distribution");
  //
  std::vector<std::string> argnames = grid_pntr->getArgNames();
  std::vector<unsigned int> args_index(0);
  for(unsigned int i=0; i<argnames.size(); i++){
    for(unsigned int l=0; l<args.size(); l++){
      if(argnames[i]==args[l]){args_index.push_back(i);}
    }
  }
  plumed_massert(args.size()==args_index.size(),"getMarginalDistributionGrid: problem with the arguments of the marginal");
  //
  MarginalWeight* Pw = new MarginalWeight();
  Grid proj_grid = grid_pntr->project(args,Pw);
  delete Pw;
  //
  // scale with the bin volume used for the integral such that the
  // marginals are proberly normalized to 1.0
  double intVol = grid_pntr->getBinVolume();
  for(unsigned int l=0; l<args_index.size(); l++){
    intVol/=grid_pntr->getDx()[args_index[l]];
  }
  proj_grid.scaleAllValuesAndDerivatives(intVol);
  //
  return proj_grid;
}


Grid TargetDistribution::getMarginal(const std::vector<std::string>& args){
  return TargetDistribution::getMarginalDistributionGrid(targetdist_grid_pntr_,args);
}


void TargetDistribution::update() {
  //
  updateGrid();
  //
  for(unsigned int i=0; i<targetdist_modifer_pntrs_.size(); i++){
    applyTargetDistModiferToGrid(targetdist_modifer_pntrs_[i]);
  }
  //
  if(bias_cutoff_active_){updateBiasCutoffForTargetDistGrid();}
  //
  if(shift_targetdist_to_zero_ && !(bias_cutoff_active_)){setMinimumOfTargetDistGridToZero();}
  if(force_normalization_ && !(bias_cutoff_active_) ){normalizeTargetDistGrid();}
  //
  // if(check_normalization_ && !force_normalization_ && !shift_targetdist_to_zero_){
  if(check_normalization_ && !(bias_cutoff_active_)){
    double normalization = integrateGrid(targetdist_grid_pntr_);
    const double normalization_thrshold = 0.1;
    if(normalization < 1.0-normalization_thrshold || normalization > 1.0+normalization_thrshold){
      std::cerr << "PLUMED WARNING - the target distribution grid in " + getName() + " is not proberly normalized, integrating over the grid gives: " << normalization << " - You can avoid this problem by using the NORMALIZE keyword\n";
    }
  }
  
  //
  if(check_nonnegative_){
    const double nonnegative_thrshold = -0.02;
    double grid_min_value = targetdist_grid_pntr_->getMinValue();
    if(grid_min_value<nonnegative_thrshold){
      std::cerr << "PLUMED WARNING - the target distribution grid in " + getName() + " has negative values, the lowest value is: " << grid_min_value << " - You can avoid this problem by using the SHIFT_TO_ZERO keyword\n";
    }
  }
  //
}


void TargetDistribution::updateBiasCutoffForTargetDistGrid() {
  plumed_massert(vesbias_pntr_!=NULL,"The VesBias has to be linked to use updateBiasCutoffForTargetDistGrid()");
  plumed_massert(vesbias_pntr_->biasCutoffActive(),"updateBiasCutoffForTargetDistGrid() should only be used if the bias cutoff is active");
  // plumed_massert(targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  // plumed_massert(log_targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  plumed_massert(getBiasWithoutCutoffGridPntr()!=NULL,"the bias without cutoff grid has to be linked");
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(targetdist_grid_pntr_);
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++)
  {
   double value = targetdist_grid_pntr_->getValue(l);
   double bias = getBiasWithoutCutoffGridPntr()->getValue(l);
   double deriv_factor_swf = 0.0;
   double swf = vesbias_pntr_->getBiasCutoffSwitchingFunction(bias,deriv_factor_swf);
   // this comes from the p(s)
   value *= swf;
   norm += integration_weights[l]*value;
   // this comes from the derivative of V(s)
   value *= deriv_factor_swf;
   targetdist_grid_pntr_->setValue(l,value);
   // double log_value = log_targetdist_grid_pntr_->getValue(l) - std::log(swf);
   // log_targetdist_grid_pntr_->setValue(l,log_value);
  }
  targetdist_grid_pntr_->scaleAllValuesAndDerivatives(1.0/norm);
  // log_targetdist_grid_pntr_->setMinToZero();
  
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
    std::vector<double> rw_integration_weights = GridIntegrationWeights::getIntegrationWeights(reweight_grid_pntr_);
	double norm = 0.0;
	for(Grid::index_t l=0; l<reweight_grid_pntr_->getSize(); l++)
	{
	  double value = reweight_grid_pntr_->getValue(l);
	  double bias = getBiasWithoutCutoffRWGridPntr()->getValue(l);
	  double deriv_factor_swf = 0.0;
	  double swf = vesbias_pntr_->getBiasCutoffSwitchingFunction(bias,deriv_factor_swf);
	  // this comes from the p(s)
	  value *= swf;
	  norm += rw_integration_weights[l]*value;
	  // this comes from the derivative of V(s)
	  value *= deriv_factor_swf;
	  reweight_grid_pntr_->setValue(l,value);
	}
	reweight_grid_pntr_->scaleAllValuesAndDerivatives(1.0/norm);
  }
}

void TargetDistribution::applyTargetDistModiferToGrid(TargetDistModifer* modifer_pntr) {
  // plumed_massert(targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  // plumed_massert(log_targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids");
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(targetdist_grid_pntr_);
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++)
  {
   double value = targetdist_grid_pntr_->getValue(l);
   std::vector<double> cv_values = targetdist_grid_pntr_->getPoint(l);
   value = modifer_pntr->getModifedTargetDistValue(value,cv_values);
   norm += integration_weights[l]*value;
   targetdist_grid_pntr_->setValue(l,value);
   log_targetdist_grid_pntr_->setValue(l,-std::log(value));
  }
  targetdist_grid_pntr_->scaleAllValuesAndDerivatives(1.0/norm);
  log_targetdist_grid_pntr_->setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
    std::vector<double> rw_integration_weights = GridIntegrationWeights::getIntegrationWeights(reweight_grid_pntr_);
    norm = 0.0;
    for(Grid::index_t l=0; l<reweight_grid_pntr_->getSize(); l++)
    {
      double value = reweight_grid_pntr_->getValue(l);
      std::vector<double> rw_cv_values = reweight_grid_pntr_->getPoint(l);
      value = modifer_pntr->getModifedTargetDistValue(value,rw_cv_values);
      norm += rw_integration_weights[l]*value;
      reweight_grid_pntr_->setValue(l,value);
      log_reweight_grid_pntr_->setValue(l,-std::log(value));
    }
    reweight_grid_pntr_->scaleAllValuesAndDerivatives(1.0/norm);
    log_reweight_grid_pntr_->setMinToZero();
  }
  //
}


void TargetDistribution::updateLogTargetDistGrid() {
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++)
  {
    log_targetdist_grid_pntr_->setValue(l,-std::log(targetdist_grid_pntr_->getValue(l)));
  }
  log_targetdist_grid_pntr_->setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
  {
	for(Grid::index_t l=0; l<reweight_grid_pntr_->getSize(); l++)
    {
      log_reweight_grid_pntr_->setValue(l,-std::log(reweight_grid_pntr_->getValue(l)));
    }
    log_reweight_grid_pntr_->setMinToZero();
  }
}


void TargetDistribution::setMinimumOfTargetDistGridToZero(){
  targetDistGrid().setMinToZero();
  // Added by Y. Isaac Yang to calculate the reweighting factor
  if(isReweightGridActive())
    reweightGrid().setMinToZero();
  //
  normalizeTargetDistGrid();
  updateLogTargetDistGrid();
}


void TargetDistribution::readInRestartTargetDistGrid(const std::string& grid_fname){
  plumed_massert(isDynamic(),"this should only be used for dynamically updated target distributions!");
  IFile gridfile;
  if(!gridfile.FileExist(grid_fname)){
    plumed_merror(getName()+": problem with reading previous target distribution when restarting, cannot find file " + grid_fname);
  }
  gridfile.open(grid_fname);
  Grid* restart_grid = Grid::create("targetdist",grid_args_,gridfile,false,false,false);
  if(restart_grid->getSize()!=targetdist_grid_pntr_->getSize()){
    plumed_merror(getName()+": problem with reading previous target distribution when restarting, the grid is not of the correct size!");
  }
  VesTools::copyGridValues(restart_grid,targetdist_grid_pntr_);
  updateLogTargetDistGrid();
  delete restart_grid;
}

void TargetDistribution::clearLogTargetDistGrid(){
  log_targetdist_grid_pntr_->clear();
}

void TargetDistribution::clearLogReweightGrid(){
  log_reweight_grid_pntr_->clear();
}

}
}
