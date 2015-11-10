#ifndef NEUTRINO_SOLUTIONS_H
#define NEUTRINO_SOLUTIONS_H

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/Polynomial.h"

// Helper Class

class neutrinoSolutions
{
private:
  
  double b1E, b1px, b1py, b1pz, b2E, b2px, b2py, b2pz, l1E, l1px, l1py, l1pz, l2E, l2px, l2py, l2pz, METx, METy, mW0, mt0;
  double mW1,mW2,mt1,mt2,mtBar,mtDelta;
  double mW1Temp,mW2Temp,mt1Temp,mt2Temp;
  
  std::vector<std::complex<double> > nu1E, nu1px, nu1py, nu1pz, nu2E, nu2px, nu2py, nu2pz;

public:
  neutrinoSolutions(const double& b1E_, const double& b1px_, const double& b1py_, const double& b1pz_, 
		    const double& l1E_, const double& l1px_, const double& l1py_, const double& l1pz_, 
		    const double& b2E_, const double& b2px_, const double& b2py_, const double& b2pz_, 
		    const double& l2E_, const double& l2px_, const double& l2py_, const double& l2pz_, 
		    const double& METx_, const double& METy_, const double& mW0_, const double& mt0_);
  neutrinoSolutions();
  ~neutrinoSolutions();
  void setupMeasurements(const double b1E_, const double B1px_, const double b1py_, const double b1pz_, 
			 const double l1E_, const double l1px_, const double l1py_, const double l1pz_, 
			 const double b2E_, const double b2px_, const double b2py_, const double b2pz_, 
			 const double l2E_, const double l2px_, const double l2py_, const double l2pz_, 
			 const double METx_,const double METy_, 
			 const double mW1_, const double mW2_,
			 const double mt1_, const double mt2_);
  //void getMasses(double& mt1, double& mt2, double& mW1, double& mW2);
  void getRealNeutrinoVectors(std::vector<double>& nu1E_real, std::vector<double>& nu1px_real, std::vector<double>& nu1py_real, std::vector<double>& nu1pz_real, 
			      std::vector<double>& nu2E_real, std::vector<double>& nu2px_real, std::vector<double>& nu2py_real, std::vector<double>& nu2pz_real);
  void getNeutrinoVectors();

  //void setTopMassRange(const double& low, const double& high){topMass.setRange(low,high);}
  //RooPlot* getTopMassFrame(){return topMass.frame();}
  //RooAbsReal* getImSolution(){return imSolution;}
  double imagness(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_);
  double imagness2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_);
  double nu1E2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_);
  double nu2E2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_);
  double topPairMass(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_, const int& i);
  double topPairMassHigh(double* mt, double *p);
  double topPairMassLow(double* mt, double *p);
  //double topPairMassDifference(double* mt, double *p);
  double imagnessTopMass(double* mt, double *p);
  double imagnessTopSplitting(double* mtDelta, double *p);
  double imagnessTopSplitting2D(double* mt, double *p);
  double imagnessTopMass1(double* mt, double *p);
  double imagnessTopMass2(double* mt, double *p);
  double imagnessTopMass2D(double* mt, double *p);
  double nu1E2TopMass(double* mt, double *p);
  double nu1E2TopSplitting(double* mt, double *p);
  double topPairPzFunction(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_, const int& i);
  double topPairPzMass(double* mt, double *p);
  double topPairPzMass1(double* mt, double *p);
  double topPairPzMass2(double* mt, double *p);
  double topPairPzSplitting(double* mt, double *p);
  double splittingRightEdge(double* mt, double *p);
  double splittingLeftEdge(double* mt, double *p);
  double splittingEdgeDifference(double* mt, double *p);
  void mW1Set(const double& mW1_){mW1 = mW1_;};
  void mW2Set(const double& mW2_){mW2 = mW2_;};
  void mt1Set(const double& mt1_){mt1 = mt1_; mtDelta = 0.5*(mt1 - mt2); mtBar = 0.5*(mt1 + mt2);};
  void mt2Set(const double& mt2_){mt2 = mt2_; mtDelta = 0.5*(mt1 - mt2); mtBar = 0.5*(mt1 + mt2);};
  void setMasses(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_){mt1 = mt1_;mt2 = mt2_;mW1 = mW1_;mW2 = mW2_;mtDelta = 0.5*(mt1 - mt2);mtBar = 0.5*(mt1 + mt2);};
};

#endif
