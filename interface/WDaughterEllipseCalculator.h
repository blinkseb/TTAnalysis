// C++ implementation of arXiv:1305.1872
// Analytic solutions for constrained momentum of second W daughter in decay of top quarks
// Gala Nicolas Kaufman (gnn4@cornell.edu)

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <list>
#include <utility>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TMath.h>
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

#include "Math/RootFinderAlgorithms.h"
#include "Math/Polynomial.h"

class WDaughterEllipseCalculator{

 private:

  const double& bJetPx_;
  const double& bJetPy_;
  const double& bJetPz_;
  const double& bJetE_ ;

  const double& WDaughter1Px_;
  const double& WDaughter1Py_;
  const double& WDaughter1Pz_;
  const double& WDaughter1E_ ;

  double bJetP2_  , bJetP_    , bJetMass2_;
  double bJetBeta_, bJetBeta2_, bJetGamma_, bJetGamma2_    ;
  double WDaughter1Beta_, WDaughter1Beta2_, WDaughter1Gamma_, WDaughter1Gamma2_;
  double WDaughter1P2_  , WDaughter1P_    , WDaughter1Mass2_;
  double WDaughter1Phi_ , WDaughter1Theta_;

  //particle masses
  double mW_ , mTop_ , WDaughter2Mass_ ;
  double mW2_, mTop2_, WDaughter2Mass2_;
  double mTopAssumed_;

  //parameters
  double x0_, x0p_;
  double Sx_, Sy_;
  double epsilon2_;

  double c_, s_,c2_,s2_; //cosine and sine of theta_{b,l}

  double omega_;
  double Omega_;
  double x1_, y1_;
  double Z2_;

  //matrices
  TMatrixD Ab_;
  TMatrixD AWDaughter1_;

  TMatrixD Htilde_;
  TMatrixD H_;
  TMatrixD Hperp_;
  TMatrixD HperpInv_;

  TMatrixD Nperp_;

  TVectorD WDaughterPerp_;
  TVectorD pWDaughter_;

  double mTopEdgeLow_, mTopEdgeHigh_;

 public:

  void setBJetFactors();
  void setMeasuredWDaughterFactors();

  double getZ2(double, double, double);
  double getZ2();

  void setAngles();

  void initializeMatrices();

  TMatrixD rotationMatrix(int, double);

  void Wsurface();
  void bJetEllipsoid();
  void measuredWDaughterEllipsoid();
  void calcZ2();
  void WDaughterSolution();
  void labSystemTransform();

  bool errorFlag_;

  WDaughterEllipseCalculator(const double& , const double& , const double& , const double& , const double& , const double& , const double& , const double& );
  WDaughterEllipseCalculator(const double& , const double& , const double& , const double& , const double& , const double& , const double& , const double& , double , double , double );

  ~WDaughterEllipseCalculator();

  bool badPoint(){return errorFlag_;};

  void setupEllipse(double , double , double );

  //void setBJet(const double , const double, const double , const double );
  //void setLepton(const double , const double, const double , const double );

  void setTopMass(double& mTop){ mTop_=mTop; mTop2_=mTop_*mTop_; };
  void setWBosonMass(double& mW){ mW_=mW; mW2_=mW_*mW_; };
  void setWDaughterMass(double& WDaughter2Mass){ WDaughter2Mass_=WDaughter2Mass; WDaughter2Mass2_=WDaughter2Mass_*WDaughter2Mass_; };
  void setMasses(double& , double& , double& );

  double getTopMass() { return mTop_; } ;
  double getTopMassEdgeLow(){ return mTopEdgeLow_; };
  double getTopMassEdgeHigh(){ return mTopEdgeHigh_; };
  
  TMatrixD* getHomogeneousWDaughterEllipse();
  TMatrixD* getExtendedWDaughterEllipse();

  TVectorD* getWDaughterMomentum(double theta);

  void calcWDaughterEllipse();
  void calcExtendedWDaughterEllipse();
  
  //void calcTopMassCorrection();

  //double getBJetEnergy() { return bJetE_; } ;
  void printFactors();

};

