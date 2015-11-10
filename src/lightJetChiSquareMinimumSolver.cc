#include <cp3_llbb/TTAnalysis/interface/lightJetChiSquareMinimumSolver.h>
#include "TDecompSVD.h"
#include "TDecompLU.h"


lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(const std::vector<PhysicsObject>& jets, double& dx, double& dy, double& dz, bool do3D ) :
  do3D_                (do3D),
  jetPxWidths2_        (std::vector<double>(jets.size(),0.)),
  jetPyWidths2_        (std::vector<double>(jets.size(),0.)),
  jetPzWidths2_        (std::vector<double>(jets.size(),0.)),
  jetPxPyWidths_       (std::vector<double>(jets.size(),0.)),
  jetPxPzWidths_       (std::vector<double>(jets.size(),0.)),
  jetPyPzWidths_       (std::vector<double>(jets.size(),0.)),
  dx_                  (dx),
  dy_                  (dy),
  dz_                  (dz),
  dxCheck_             (0.),
  dyCheck_             (0.),
  dzCheck_             (0.),
  jetSigmas2D_         (std::vector<TMatrixD>(jets.size(),TMatrixD(2,2))),
  jetSigmas3D_         (std::vector<TMatrixD>(jets.size(),TMatrixD(3,3))),
  inverter2D_          (new TDecompLU(2)),
  inverter3D_          (new TDecompLU(3)),
  inverseSumSigmas2D_  (2,2),
  inverseSumSigmas3D_  (3,3),
  minDeltasX_          (jets.size(),0.),
  minDeltasY_          (jets.size(),0.),
  minDeltasZ_          (jets.size(),0.),
  chi2_                (0.),
  nJets_               (jets.size())
{
  setCartesianWidths(jets);
  calcSigmas();
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(int nObjects,
							       double& dx, double& dy, double& dz, bool do3D ) :
  do3D_                (do3D),
  jetPxWidths2_        (std::vector<double>(nObjects,0.)),
  jetPyWidths2_        (std::vector<double>(nObjects,0.)),
  jetPzWidths2_        (std::vector<double>(nObjects,0.)),
  jetPxPyWidths_       (std::vector<double>(nObjects,0.)),
  jetPxPzWidths_       (std::vector<double>(nObjects,0.)),
  jetPyPzWidths_       (std::vector<double>(nObjects,0.)),
  dx_                  (dx),
  dy_                  (dy),
  dz_                  (dz),
  dxCheck_             (0.),
  dyCheck_             (0.),
  dzCheck_             (0.),
  jetSigmas2D_         (std::vector<TMatrixD>(nObjects,TMatrixD(2,2))),
  jetSigmas3D_         (std::vector<TMatrixD>(nObjects,TMatrixD(3,3))),
  inverter2D_          (new TDecompLU(2)),
  inverter3D_          (new TDecompLU(3)),
  inverseSumSigmas2D_  (2,2),
  inverseSumSigmas3D_  (3,3),
  minDeltasX_          (nObjects,0.),
  minDeltasY_          (nObjects,0.),
  minDeltasZ_          (nObjects,0.),
  chi2_                (0.),
  nJets_               (nObjects)
{
  //std::cout << "Light jet chi square constructor" << std::endl;
  //std::cout << "dx is " << dx_ << std::endl;
  //std::cout << "dy is " << dy_ << std::endl;
  //std::cout << "dz is " << dz_ << std::endl;
  //std::cout << "dxCheck is " << dxCheck_ << std::endl;
  //std::cout << "dyCheck is " << dyCheck_ << std::endl;
  //std::cout << "dzCheck is " << dzCheck_ << std::endl;
  //std::cout << "End light jet chi square constructor" << std::endl;
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver& other) :
  do3D_                (other.do3D_),
  jetPxWidths2_        (other.jetPxWidths2_),
  jetPyWidths2_        (other.jetPyWidths2_),
  jetPzWidths2_        (other.jetPzWidths2_),
  jetPxPyWidths_       (other.jetPxPyWidths_),
  jetPxPzWidths_       (other.jetPxPzWidths_),
  jetPyPzWidths_       (other.jetPxPzWidths_),
  dx_                  (other.dx_),
  dy_                  (other.dy_),
  dz_                  (other.dz_),
  dxCheck_             (other.dxCheck_),
  dyCheck_             (other.dyCheck_),
  dzCheck_             (other.dzCheck_),
  jetSigmas2D_         (other.jetSigmas2D_),
  jetSigmas3D_         (other.jetSigmas3D_),
  inverter2D_          (other.inverter2D_),
  inverter3D_          (other.inverter3D_),
  inverseSumSigmas2D_  (other.inverseSumSigmas2D_),
  inverseSumSigmas3D_  (other.inverseSumSigmas3D_),
  minDeltasX_          (other.minDeltasX_),
  minDeltasY_          (other.minDeltasY_),
  minDeltasZ_          (other.minDeltasZ_),
  chi2_                (other.chi2_),
  nJets_               (other.nJets_)
{
}

lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
  delete inverter2D_;
  delete inverter3D_;
}

//void lightJetChiSquareMinimumSolver::setRecoil(double x, double y, double z)
//{
//  std::cout << "dx is " << dx_ << std::endl;
//  std::cout << "dy is " << dy_ << std::endl;
//  std::cout << "dz is " << dz_ << std::endl;
//  std::cout << "dxCheck is " << dxCheck_ << std::endl;
//  std::cout << "dyCheck is " << dyCheck_ << std::endl;
//  std::cout << "dzCheck is " << dzCheck_ << std::endl;
//
//  //std::cout << "Setting recoil" << std::endl;
////  dxCheck_=x;
////  dyCheck_=y;
////  if(do3D_) dzCheck_=z;
//
////  std::cout << "dx is " << dx_ << std::endl;
////  std::cout << "dy is " << dy_ << std::endl;
////  std::cout << "dz is " << dz_ << std::endl;
////  std::cout << "dxCheck is " << dxCheck_ << std::endl;
////  std::cout << "dyCheck is " << dyCheck_ << std::endl;
////  std::cout << "dzCheck is " << dzCheck_ << std::endl;
//}

void lightJetChiSquareMinimumSolver::setupEquations(const std::vector<PhysicsObject>& jets)
{
  if(jets.size() != jetPxWidths2_.size())
    {
      std::cout << "Unequal number of cartesian and radial jets!" << std::endl;
      return;
    }
  setCartesianWidths(jets);
  calcSigmas();
}

void lightJetChiSquareMinimumSolver::setCartesianWidths(const std::vector<PhysicsObject>& jets) 
{
  for( unsigned int i = 0 ; i < nJets_;  i++)
    {
        std::cout << jets[i].pt_width << std::endl;

      double halfPt2 = 0.5*jets[i].p4.Pt()*jets[i].p4.Pt();
      double sigmaPt2 = log(1+jets[i].pt_width);
      sigmaPt2*=sigmaPt2;
      double sigmaPhi2 = jets[i].phi_width*jets[i].phi_width;
      double sigmaEta = jets[i].eta_width;
      double sigmaEta2 = sigmaEta*sigmaEta;

      double expSigmaPt2 = exp(sigmaPt2);
      double expTwoSigmaPt2 = expSigmaPt2*expSigmaPt2;
      double expMinusSigmaPhi2 = exp(-sigmaPhi2);
      double expMinusTwoSigmaPhi2 = expMinusSigmaPhi2*expMinusSigmaPhi2;
      double expMinusSigmaPhi2Over2 = exp(-0.5*sigmaPhi2);
      double expSigmaEta2 = exp(sigmaEta2);
      double expTwoSigmaEta2 = expSigmaEta2*expSigmaEta2;
      double expSigmaEta2OverTwo = exp(-0.5*sigmaEta2);
      double cosPhi = cos(jets[i].p4.Phi());
      double sinPhi = sin(jets[i].p4.Phi());
      double cosTwoPhi = cos(2.*jets[i].p4.Phi());
      double sinTwoPhi = sin(2.*jets[i].p4.Phi());
      double coshTwoEta = cosh(2.*jets[i].p4.Eta());
      double sinhEta  = sinh(jets[i].p4.Eta());
      double sinhEta2 = sinhEta*sinhEta;
      jetPxWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 + cosTwoPhi*expMinusTwoSigmaPhi2)
				      - expSigmaPt2 * (1 + cosTwoPhi) * expMinusSigmaPhi2) ;
      jetPyWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 - cosTwoPhi*expMinusTwoSigmaPhi2)
				      - expSigmaPt2 * (1 - cosTwoPhi) * expMinusSigmaPhi2) ;
      jetPzWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (expTwoSigmaEta2*coshTwoEta - 1) //FIXMEis there a typo here? expMinusTwoSigmaEta instead?
					  - 2.*expSigmaPt2*expSigmaEta2*sinhEta2 ); 
      jetPxPyWidths_[i] = halfPt2 * sinTwoPhi *( expTwoSigmaPt2 *expMinusTwoSigmaPhi2
						- expSigmaPt2 * expMinusSigmaPhi2);
      jetPxPzWidths_[i] = 2.*halfPt2*cosPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2 * ( expTwoSigmaPt2 - expSigmaPt2 );
      jetPyPzWidths_[i] = 2.*halfPt2*sinPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2 * ( expTwoSigmaPt2 - expSigmaPt2 );

//      std::cout << "calculating widths:\n"
//           << "pt  is " << jets.at(i).Pt()  << " with width of " << log(1+jetPtWidths .at(i)) << "\n"
//           << "phi is " << jets.at(i).Phi() << " with width of " << log(1+jetPhiWidths .at(i)) << "\n"
//           << "px  is " << jets.at(i).Pt()*cos(jets.at(i).Phi()) << " with width of " << sqrt(jetPxWidths2_.at(i)) << "\n"
//           << "py  is " << jets.at(i).Pt()*sin(jets.at(i).Phi()) << " with width of " << sqrt(jetPyWidths2_.at(i)) << "\n"
//           << "correlation coefficient is " << jetPxPyWidths_.at(i)/(sqrt(jetPxWidths2_.at(i))*sqrt(jetPyWidths2_.at(i))) << std::endl;
    }
}

void lightJetChiSquareMinimumSolver::calcSigmas()
{
  inverseSumSigmas2D_.Zero();
  inverseSumSigmas3D_.Zero();

  for( unsigned int i = 0 ; i < nJets_ ; i++ )
    {
      double px2(0.), py2(0.), pz2(0.), pxpy(0.), pxpz(0.), pypz(0.);
      px2 = jetPxWidths2_.at(i);
      py2 = jetPyWidths2_.at(i);
      pxpy = jetPxPyWidths_.at(i);
      if( do3D_ ) 
	{
	  pz2 = jetPzWidths2_.at(i);
	  pxpz = jetPxPzWidths_.at(i);
	  pypz = jetPyPzWidths_.at(i);
	}
      if( do3D_ )
	{
	  double array3D[9]= {0.};
	  array3D[0] = px2;
	  array3D[1] = pxpy;
	  array3D[2] = pxpz;
	  array3D[3] = pxpy;
	  array3D[4] = py2;
	  array3D[5] = pypz;
	  array3D[6] = pxpz;
	  array3D[7] = pypz;
	  array3D[8] = pz2;
	  jetSigmas3D_.at(i) = TMatrixD(3,3,array3D);
	  inverseSumSigmas3D_ += jetSigmas3D_.at(i);
	}
      else
	{
	  double array2D[4] = {0.};
	  array2D[0] = px2;
	  array2D[1] = pxpy;
	  array2D[2] = pxpy;
	  array2D[3] = py2;
	  jetSigmas2D_.at(i) = TMatrixD(2,2,array2D);
	  inverseSumSigmas2D_ += jetSigmas2D_.at(i);
	}
    }

  if( do3D_ )
    {
      inverseSumSigmas3D_.Print();
      dynamic_cast<TDecompLU*>(inverter3D_)->SetMatrix(TMatrixD(inverseSumSigmas3D_));
      inverter3D_->Decompose();
      dynamic_cast<TDecompLU*>(inverter3D_)->Invert(inverseSumSigmas3D_);
    }
  else
    {
      inverseSumSigmas2D_.Print();
      dynamic_cast<TDecompLU*>(inverter2D_)->SetMatrix(TMatrixD(inverseSumSigmas2D_));
      inverter2D_->Decompose();
      dynamic_cast<TDecompLU*>(inverter2D_)->Invert(inverseSumSigmas2D_);
      //inverseSumSigmas2D_.Print();
    }
}

void lightJetChiSquareMinimumSolver::calcMin()
{
  //std::cout << "dx is " << dx_ << std::endl;
  //std::cout << "dy is " << dy_ << std::endl;
  //std::cout << "dz is " << dz_ << std::endl;
  //std::cout << "dxCheck is " << dxCheck_ << std::endl;
  //std::cout << "dyCheck is " << dyCheck_ << std::endl;
  //std::cout << "dzCheck is " << dzCheck_ << std::endl;

  //if(dxCheck_ == dx_ && dyCheck_ == dy_ && dzCheck_ == dz_) return;
  if(do3D_)
    {
      if(dxCheck_ == dx_ && dyCheck_ == dy_ && dzCheck_ == dz_) return;
    }
  else
    {
      if(dxCheck_ == dx_ && dyCheck_ == dy_) return;
    }
  
  //std::cout << "Calculating minimum chi^2" << std::endl;

  dxCheck_ = dx_;
  dyCheck_ = dy_; 
  if( do3D_ ) dzCheck_ = dz_;

  //chi2_ = dx_*dx_ + dy_*dy_;
  //if(do3D_) chi2_ += dz_*dz_;
  //return;
  chi2_ = 0;

  double dArray3D[3] = { dx_ , dy_ , dz_ };
  TVectorD dVec3D(3,dArray3D);

  double dArray2D[2] = { dx_ , dy_ };
  TVectorD dVec2D(2,dArray2D);

  for( unsigned int i = 0 ; i < nJets_ ; i++ )
    {
      if( do3D_ )
	{
	  TMatrixD thisJetB(jetSigmas3D_.at(i));
	  thisJetB *= inverseSumSigmas3D_;
	  TVectorD thisJetDelta = thisJetB*dVec3D;
	  minDeltasX_.at(i) = thisJetDelta[0];
	  minDeltasY_.at(i) = thisJetDelta[1];
	  minDeltasZ_.at(i) = thisJetDelta[2];
	}
      else
	{
	  TMatrixD thisJetB(jetSigmas2D_.at(i));
          thisJetB *= inverseSumSigmas2D_;
          TVectorD thisJetDelta = thisJetB*dVec2D;
          minDeltasX_.at(i) = thisJetDelta[0];
          minDeltasY_.at(i) = thisJetDelta[1];
	  //minDeltasZ_.at(i) = 0.;
	}
    }  

  if( do3D_ )
    {
      chi2_ = dVec3D * ( inverseSumSigmas3D_ * dVec3D ) ;
    }
  else
    {
      chi2_ = dVec2D * ( inverseSumSigmas2D_ * dVec2D ) ;
    }

}

void lightJetChiSquareMinimumSolver::printResults()
{
  for( unsigned int i = 0 ; i < nJets_ ; i++ )
    {
      std::cout << "delta px " << i+1 << " = " << minDeltasX_.at(i) << std::endl;
      std::cout << "delta py " << i+1 << " = " << minDeltasY_.at(i) << std::endl;
      if( do3D_ )
	{
	  std::cout << "delta pz " << i+1 << " = " << minDeltasZ_.at(i) << std::endl;
	}
    }

}

double lightJetChiSquareMinimumSolver::getChiSquare()
{
  calcMin();
//  std::vector<double>::iterator thisDeltaX = minDeltasX_.begin();
//  double deltaXCheck(0.);
//  double deltaYCheck(0.);
//  for(std::vector<double>::iterator thisDeltaY = minDeltasY_.begin(); thisDeltaY != minDeltasY_.end(); thisDeltaX++, thisDeltaY++)
//    {
//      deltaXCheck+=*thisDeltaX;
//      deltaYCheck+=*thisDeltaY;
//    }
//  std::cout << "delta x = " << dx_ << " and delta x check = " << deltaXCheck << std::endl;
//  std::cout << "delta y = " << dy_ << " and delta y check = " << deltaYCheck << std::endl;
  //printResults();
  return chi2_;
}
