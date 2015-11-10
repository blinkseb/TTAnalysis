#ifndef TOPEVENTMINIMIZER
#define TOPEVENTMINIMIZER

#include "neutrinoSolutions.h"
#include "topSystemChiSquare.h"
#include "leptonicTopSystemChiSquare.h"
#include "lightJetChiSquareMinimumSolver.h"

class topEventMinimizer
{
 private:
  int nTops_;

  std::vector<int> bJets_;
  std::vector<int> firstWDaughters_;
  std::vector<int> secondWDaughters_;

  std::vector<PhysicsObject> nonTopObjects_;

  double nonTopPx_;
  double nonTopPy_;
  double nonTopPz_;

  double mTop_;
  double sigmaMTop_;

  double mW_;
  double sigmaMW_;


  std::vector<double> bJets_PtDeltas_;
  std::vector<double> bJets_PhiDeltas_;
  std::vector<double> bJets_EtaDeltas_;
  std::vector<double> firstWDaughters_PtDeltas_;
  std::vector<double> firstWDaughters_PhiDeltas_;
  std::vector<double> firstWDaughters_EtaDeltas_;
  std::vector<double> topMassDeltas_;
  std::vector<double> WMassDeltas_;

  std::vector<double> bJets_PtDeltasBest_;
  std::vector<double> bJets_PhiDeltasBest_;
  std::vector<double> bJets_EtaDeltasBest_;
  std::vector<double> firstWDaughters_PtDeltasBest_;
  std::vector<double> firstWDaughters_PhiDeltasBest_;
  std::vector<double> firstWDaughters_EtaDeltasBest_;
  std::vector<double> topMassDeltasBest_;
  std::vector<double> WMassDeltasBest_;

  std::vector<double> secondWDaughters_PtDeltasBest_;
  std::vector<double> secondWDaughters_PhiDeltasBest_;
  std::vector<double> secondWDaughters_EtaDeltasBest_;

  std::vector<double> topMassDeltasCurrent_;
  std::vector<double> topMassDeltasInnerBest_;


  std::vector<double> nonTopObjects_PxDeltasBest_;
  std::vector<double> nonTopObjects_PyDeltasBest_;
  //vector<double> nonTopObjects_PzDeltasBest_;


  //vector<pair<double, double> > WDaughterMasses_;

  double chi2_;

  std::vector<std::pair<topSystemChiSquare*, bool> > topSystemChiSquares_;

  double nonTopChi2_;
  double dx_, dy_, dz_;

  lightJetChiSquareMinimumSolver nonTopChiSquare_;

  std::vector<double> ellipseAngles_;
  std::vector<double> ellipseAnglesBest_;
  std::vector<double> ellipseAnglesCurrent_;
  std::vector<double> ellipseAnglesInnerBest_;

  double hadChi2_;
  double topChi2_;
  double topMassChi2_;

  double chi2Best_;
  double innerChi2Best_;
  double topChi2Best_;
  double hadChi2Best_;
  double topMassChi2Best_;
  double nonTopChi2Best_;

  double thisInnerChi2Best_;
  double thisTopMassChi2Best_;
  double thisNonTopChi2Best_;
  double thisHadChi2Best_;

  double maxConsideredChiSquareRoot_;

  int innerMinStatus_;
  int outerMinStatus_;

  //ROOT::Math::Minimizer* ellipseAngleMin_;
  //ROOT::Math::Minimizer* topMassMin_;
  ROOT::Math::Minimizer* innerMin_;
  ROOT::Math::Minimizer* outerMin_;

  bool checkInputSizes();

  //void setBJets();
  //void setWDaughters();

  void initializeDeltas();
  void initializeChiSquares();

  void setRecoil(double , double, double ); 

  void calcWDaughterEllipses();

  void getDxDyFromEllipses();
  void calcNonTopMomentum();

  void buildBestNonTopObjects();

  void setupNonTopChiSquare();

  void calcTopMassChiSquare();
  void calcTopChiSquare();

  double getTopChiSquare();
  double getTopMassChiSquare();
  double getNonTopChiSquare();

  void setBestValues();

 public:

  topEventMinimizer(const std::vector<PhysicsObject>&,
		    double , double ,
		    double , double );

  ~topEventMinimizer();


  topSystemChiSquare* makeLeptonicTop(const PhysicsObject& bjet, const PhysicsObject& lepton);

  void addLeptonicTop(const PhysicsObject& bjet, const PhysicsObject& lepton);

  void printTopConstituents();
  void printNonTopObjects();

  void calcTopMassRanges();

  //double ellipseAngleMinimizationOperator(const double* );
  double outerMinimizationOperator(const double* );
  double innerMinimizationOperator(const double* );

  void findStartingValues(int);
  void minimizeNonTopChiSquare();
  void minimizeTotalChiSquare();

  void calcTotalChiSquare();
  double getChiSquare();

  int getInnerMinimizerStatus() { return innerMinStatus_; } ;
  int getOuterMinimizerStatus() { return outerMinStatus_; } ;
  double getOuterMinimizerEdm() { if(outerMin_) return outerMin_->Edm(); return -1.; };

  double getBestTotalChiSquare() { return chi2Best_; } ;
  double getBestTopSystemChiSquare() { return topChi2Best_; } ;
  double getBestTopMassChiSquare() { return topMassChi2Best_; } ;
  double getBestNonTopChiSquare() { return nonTopChi2Best_; } ;

  void getBestDeltas(std::vector<double>& , std::vector<double>& , std::vector<double>& ,
		     std::vector<double>& , std::vector<double>& , std::vector<double>& ,
		     std::vector<double>& , std::vector<double>& , std::vector<double>& ,
		     std::vector<double>& , std::vector<double>& ,
		     std::vector<double>& , std::vector<double>& );

  void getBJet(int, double& , double& , double& , double& );
  void getWDaughter1(int, double& , double& , double& , double& );
  void getWDaughter2(int, double& , double& , double& , double& );
  void getTop(int, double& , double& , double& , double& );
  void getW(int, double& , double& , double& , double& );
  void getNonTopObject(int, double& , double& );

  void plotEllipses(TString);
  
};

#endif

