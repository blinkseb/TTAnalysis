#include <cp3_llbb/TTAnalysis/interface/leptonicTopSystemChiSquare.h>

leptonicTopSystemChiSquare::leptonicTopSystemChiSquare(const PhysicsObject& bjet, const PhysicsObject& lepton, const double& mTop, const double& sigmaMTop, const double& mW, const double& sigmaMW ) :
  topSystemChiSquare   (bjet, lepton,
                        mTop, sigmaMTop,
                        mW, sigmaMW,
                        0.
                        ),
  chi2_                 (0.)
{

}

leptonicTopSystemChiSquare::~leptonicTopSystemChiSquare()
{
  //cout << "destructor" << endl;
}

void leptonicTopSystemChiSquare::printTopConstituents()
{
    std::cout << "Leptonic top decay products:" << std::endl;
    std::cout << "b-jet: "
       << "\npx = " << bJetPx_
       << "\npy = " << bJetPy_
       << "\npz = " << bJetPz_
       << "\ne  = " << bJetE_
//       << "\nm  = " << sqrt(max(0.,bJetE_*bJetE_-bJetPx_*bJetPx_-bJetPy_*bJetPy_-bJetPz_*bJetPz_)) 
       << std::endl;

    std::cout << "lepton: "
       << "\npx = " << WDaughter1Px_
       << "\npy = " << WDaughter1Py_
       << "\npz = " << WDaughter1Pz_
       << "\ne  = " << WDaughter1E_
//       << "\nm  = " << sqrt(max(0.,WDaughter1E_*WDaughter1E_-WDaughter1Px_*WDaughter1Px_-WDaughter1Py_*WDaughter1Py_-WDaughter1Pz_*WDaughter1Pz_)) 
       << std::endl;

//  cout << "top mass: " << getTopMass() << endl;
//
//  cout << "W mass: " << getWMass() << endl;
//
////
////  cout << "neutrino: "
////       << "\npx = " << WDaughter2Px_
////       << "\npy = " << WDaughter2Py_
////       << "\npz = " << WDaughter2Pz_
////       << "\ne  = " << WDaughter2E_
////       << "\nm  = " << sqrt(max(0.,WDaughter2E_*WDaughter2E_-WDaughter2Px_*WDaughter2Px_-WDaughter2Py_*WDaughter2Py_-WDaughter2Pz_*WDaughter2Pz_))
////       << "\nreco m = " << sqrt(reconstructed_WDaughter2Mass2_) << endl;
//  double low, high;
//  getTopMassRange(low,high);
}

void leptonicTopSystemChiSquare::setEllipseAngle(double theta)
{
  theta_=theta;
  resetWDaughter2(theta_);
  calcTopMomentum();
}

void leptonicTopSystemChiSquare::getWDaughter2Deltas(double& ptDelta, double& phiDelta, double& etaDelta)
{
  ptDelta =0.;
  phiDelta=0.;
  etaDelta=0.;
}

void leptonicTopSystemChiSquare::calcChiSquare()
{
  chi2_=bJetPtDelta_*bJetPtDelta_
    +bJetPhiDelta_*bJetPhiDelta_
    +bJetEtaDelta_*bJetEtaDelta_
    +WDaughter1PtDelta_*WDaughter1PtDelta_
    +WDaughter1PhiDelta_*WDaughter1PhiDelta_
    +WDaughter1EtaDelta_*WDaughter1EtaDelta_
//    +breitWignerError(mTop_,sigmaMTop_,deltaMTop_)   //moved to inner minimization
    +breitWignerError(mW_,sigmaMW_,deltaMW_);
//    +deltaMW_*deltaMW_;
}

double leptonicTopSystemChiSquare::getChiSquare()
{
  calcChiSquare();
  return chi2_;
}
