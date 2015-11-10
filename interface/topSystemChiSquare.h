#ifndef TOPSYSTEMCHISQUARE
#define TOPSYSTEMCHISQUARE

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

#include <cp3_llbb/TTAnalysis/interface/Types.h>

#include "WDaughterEllipseCalculator.h"

class topSystemChiSquare
{

  // private:
   protected:
  // public:

  double bJetPx_;
  double bJetPy_;
  double bJetPz_;
  double bJetE_ ;

  double WDaughter1Px_;
  double WDaughter1Py_;
  double WDaughter1Pz_;
  double WDaughter1E_ ;

  double mTop_;
  double mW_;
  //double mWDaughter1_;
  double mWDaughter2_;

  //widths
  double bJetPtWidth_ ;
  double bJetEtaWidth_;
  double bJetPhiWidth_;

  double WDaughter1PtWidth_ ;
  double WDaughter1EtaWidth_;
  double WDaughter1PhiWidth_;

  double sigmaMTop_;
  double sigmaMW_;

  //delta parameters for measured vs true values
  double bJetPtDelta_ ;
  double bJetEtaDelta_; 
  double bJetPhiDelta_;

  double WDaughter1PtDelta_ ;
  double WDaughter1EtaDelta_;
  double WDaughter1PhiDelta_;

  double deltaMTop_;
  double deltaMW_;

  bool rangeFlag_;
  bool hasHighEdge_;
  double mTopEdgeLow_, mTopEdgeHigh_;
  double deltaMTopRangeLow_, deltaMTopRangeHigh_;

  //reconstructed values
  PxPyPzELorentzVector reconstructed_bJetLorentzVector_;

  double reconstructed_bJetPt_ ;
  double reconstructed_bJetPhi_;
  double reconstructed_bJetEta_;
  double reconstructed_bJetMass2_ ;

  PxPyPzELorentzVector reconstructed_WDaughter1LorentzVector_;

  double reconstructed_WDaughter1Pt_ ;
  double reconstructed_WDaughter1Phi_;
  double reconstructed_WDaughter1Eta_;
  double reconstructed_WDaughter1Mass2_ ;


  
  WDaughterEllipseCalculator WDaughter2Calc_;

  double WDaughter2Px_;
  double WDaughter2Py_;
  double WDaughter2Pz_;
  double WDaughter2E_ ;
		  
  double WDaughter2Pt_ ;
  double WDaughter2Phi_;
  double WDaughter2Eta_;

  double topPx_;
  double topPy_;
  double topPz_;
  double topE_ ;


  void setBJet(double, double, double, double);
  void setObservedBJet(double, double, double, double);
  void setWDaughter1(double, double, double, double);
  void setObservedWDaughter1(double, double, double, double);

  void setTopMass(double m){ mTop_=m; };
  void setWMass(double m){ mW_=m; };
  void setWDaughter2Mass(double m){ mWDaughter2_=m; };

  void resetBJet();
  void resetWDaughter1();
  void resetWDaughter2(double);
  //void resetAll();

  void setBJetDeltas(const double, const double, const double );
  void setWDaughter1Deltas(const double, const double, const double );

  void setBJetWidths(double , double, double );
  void setWDaughter1Widths(double , double , double );
  void setTopMassWidth(double sigmaMTop){ sigmaMTop_=sigmaMTop; };
  void setWMassWidth(double sigmaMW){ sigmaMW_=sigmaMW; };

  void setWMassDelta(double delta) { deltaMW_=delta; } ;


  double breitWignerError(double&,double&,const double&);

 public:

  topSystemChiSquare(const PhysicsObject& bjet, const PhysicsObject& lepton,
		     const double& , const double& , 
		     const double& , const double& ,
		     const double& );

  topSystemChiSquare(const topSystemChiSquare& other);

  virtual ~topSystemChiSquare();

  virtual void printTopConstituents() =0;

  //void printTopInfo();

  void calcTopMassRange();
  void setTopMassDelta(double delta) { deltaMTop_=delta; } ;
  void setDeltas(double, double, double, 
		 double, double, double,
		 //double, 
		 double);
  TMatrixD* getHomogeneousWDaughterEllipse();
  
  double theta_;

  void getBJet(double& , double& , double& , double& );
  void getWDaughter1(double& , double& , double& , double& );
  void getWDaughter2(double& , double& , double& , double& );

  void setupWDaughter2Ellipse();
  void calcWDaughter2Ellipse();
  void getWDaughter2Momentum(double& , double& , double& );
  virtual void setEllipseAngle(double) =0;
  double getEllipseAngle() { return theta_; } ;
  virtual void setWDaughter2(double, double, double) =0;
  virtual void calcWDaughter2Deltas() =0;
  virtual void getWDaughter2Deltas(double& , double& , double& ) =0;
  virtual void printWDaughter2() =0;
  double getZ2(double, double, double);

  bool hasTopMassRange() { return rangeFlag_; } ;
  void getTopMassRange(bool& , double& , double& );
  void getTopMassDeltaRange(bool& , double& , double& );
  double getTopMass() { return mTop_+sigmaMTop_*deltaMTop_; };

  double getWMass() { return mW_+sigmaMW_*deltaMW_; };

  virtual double getChiSquare() =0 ;
  virtual double getHadronicChiSquare() =0;
  
  double getTopMassChiSquare();

  void calcTopMomentum();
  void getTopMomentum(double& , double& , double& , double& );

};

//void topSystemChiSquare::printTopInfo()
//{
//  cout << "In the parent class:" << endl;
//  cout << "Second W daughter momentum:"
//       << "\npx = " << WDaughter2Px_
//       << "\npy = " << WDaughter2Py_
//       << "\npz = " << WDaughter2Pz_
//       << endl;
//}

#endif

