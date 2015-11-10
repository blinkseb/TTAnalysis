#include <cp3_llbb/TTAnalysis/interface/WDaughterEllipseCalculator.h>

#include <complex>

WDaughterEllipseCalculator::WDaughterEllipseCalculator(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
						       const double& WDaughter1Px, const double& WDaughter1Py,const double& WDaughter1Pz, const double& WDaughter1E
						       ):
  bJetPx_(bJetPx),
  bJetPy_(bJetPy),
  bJetPz_(bJetPz),
  bJetE_ (bJetE),
  WDaughter1Px_(WDaughter1Px),
  WDaughter1Py_(WDaughter1Py),
  WDaughter1Pz_(WDaughter1Pz),
  WDaughter1E_ (WDaughter1E),
  Ab_(4,4),
  AWDaughter1_(4,4),
  Htilde_(3,3),
  H_(3,3),
  Hperp_(3,3),
  HperpInv_(3,3),
  Nperp_(3,3),
  WDaughterPerp_(3),
  pWDaughter_(3),
  mTopEdgeLow_(-1.),
  mTopEdgeHigh_(-1.),
  errorFlag_(false)
{
  initializeMatrices();
}

WDaughterEllipseCalculator::WDaughterEllipseCalculator(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
						       const double& WDaughter1Px, const double& WDaughter1Py, const double& WDaughter1Pz, const double& WDaughter1E,
						       double mTop, double mW, double WDaughter2Mass)  :
  bJetPx_(bJetPx),
  bJetPy_(bJetPy),
  bJetPz_(bJetPz),
  bJetE_ (bJetE),
  WDaughter1Px_(WDaughter1Px),
  WDaughter1Py_(WDaughter1Py),
  WDaughter1Pz_(WDaughter1Pz),
  WDaughter1E_ (WDaughter1E),
  mTopAssumed_(mTop),
  Ab_(4,4),
  AWDaughter1_(4,4),
  Htilde_(3,3),
  H_(3,3),
  Hperp_(3,3),
  HperpInv_(3,3),
  Nperp_(3,3),
  WDaughterPerp_(3),
  pWDaughter_(3),
  mTopEdgeLow_(-1.),
  mTopEdgeHigh_(-1.),
  errorFlag_(false)
{
  //cout << "constructor" << endl;

  setBJetFactors();
  //printFactors();
  setMeasuredWDaughterFactors();
  //printFactors();

  setMasses(mTop,mW,WDaughter2Mass);
  //printFactors();
  setAngles();
  //printFactors();

  initializeMatrices();

  //printFactors();

  //cout << "end constructor" << endl;
}

WDaughterEllipseCalculator::~WDaughterEllipseCalculator()
{
  //cout << "destructor" << endl;
}


void WDaughterEllipseCalculator::setupEllipse(double mTop, double mW, double WDaughter2Mass)
{
  //cout << "setting up the measured W daughter ellipse" << endl;
  setBJetFactors();
  setMeasuredWDaughterFactors();
  setMasses(mTop,mW,WDaughter2Mass);
  setAngles();
  initializeMatrices();
  Wsurface();
  measuredWDaughterEllipsoid();
  bJetEllipsoid();
  calcZ2();
  //calcTopMassCorrection();
}

double WDaughterEllipseCalculator::getZ2(double mTop, double mW, double WDaughter2Mass)
{
  //double tempmTop=mTop_, tempmW=mW_, tempWDaughter2Mass=WDaughter2Mass_;
//  setBJetFactors();
//  setMeasuredWDaughterFactors();
//  setMasses(mTop,mW,WDaughter2Mass);
//  setAngles();
//  initializeMatrices();
//  Wsurface();
//  measuredWDaughterEllipsoid();
//  bJetEllipsoid();
//  calcZ2();
  setupEllipse(mTop,mW,WDaughter2Mass);
  //cout << "Top mass is " <<mTop_ << endl;
  //cout << "W mass is " << mW_ << endl;
  //cout << "second W daughter mass is " << WDaughter2Mass_ << endl;
  return Z2_;
}

double WDaughterEllipseCalculator::getZ2()
{
  setBJetFactors();
  setMeasuredWDaughterFactors();
  setMasses(mTop_,mW_,WDaughter2Mass_);
  //cout << "Top mass is " <<mTop_ << endl;
  //cout << "W mass is " << mW_ << endl;
  setAngles();
  initializeMatrices();
  Wsurface();
  measuredWDaughterEllipsoid();
  bJetEllipsoid();
  calcZ2();
  return Z2_;
}

void WDaughterEllipseCalculator::setMasses(double& mTop, double& mW, double& WDaughter2Mass)
{
  //cout << "setting masses" << endl;
  mTop_=mTop;
  mW_=mW;
  WDaughter2Mass_=WDaughter2Mass; 
  mW2_=mW_*mW_;
  mTop2_=mTop_*mTop_;
  WDaughter2Mass2_=WDaughter2Mass_*WDaughter2Mass_;
}

void WDaughterEllipseCalculator::setBJetFactors()
{
  //cout << "setting b-jet factors" << endl;
  bJetP2_ = (bJetPx_*bJetPx_ + bJetPy_*bJetPy_ + bJetPz_*bJetPz_);
  bJetP_ = sqrt(bJetP2_);
  double bJetE2 =  bJetE_*bJetE_;
  bJetMass2_ = bJetE2 - bJetP2_;
  bJetBeta2_=bJetP2_/bJetE2;
  bJetBeta_=sqrt(bJetBeta2_);
  bJetGamma2_=1.0/(1.0 - bJetBeta2_);
  bJetGamma_=sqrt(bJetGamma2_);
}

void WDaughterEllipseCalculator::setMeasuredWDaughterFactors()
{
  //cout << "setting measured W daughter factors" << endl;
  WDaughter1P2_ = (WDaughter1Px_*WDaughter1Px_ + WDaughter1Py_*WDaughter1Py_ + WDaughter1Pz_*WDaughter1Pz_);
  WDaughter1P_ = sqrt(WDaughter1P2_);
  double WDaughter1E2 =  WDaughter1E_*WDaughter1E_;
  WDaughter1Mass2_ = std::max(0.,WDaughter1E2 - WDaughter1P2_);
  WDaughter1Beta2_ = WDaughter1P2_/WDaughter1E2;
  WDaughter1Beta_ = sqrt(WDaughter1Beta2_);
  WDaughter1Gamma2_ = 1.0/(1.0 - WDaughter1Beta2_);
  WDaughter1Gamma_ = sqrt(WDaughter1Gamma2_);
  WDaughter1Phi_ = (WDaughter1Px_ == 0.0 && WDaughter1Py_ == 0.0) ? 0.0 : atan2(WDaughter1Py_,WDaughter1Px_);
  WDaughter1Theta_ = (WDaughter1Px_ == 0.0 && WDaughter1Py_ == 0.0 && WDaughter1Pz_ == 0.0) ? 0.0 : atan2(sqrt(WDaughter1Px_*WDaughter1Px_ + WDaughter1Py_*WDaughter1Py_),WDaughter1Pz_);
}

void WDaughterEllipseCalculator::setAngles()
{
  //cout << "setting angles" << endl;
  double WDaughter1DotbJet = (WDaughter1Px_*bJetPx_ + WDaughter1Py_*bJetPy_ + WDaughter1Pz_*bJetPz_);
  c2_ = WDaughter1DotbJet*WDaughter1DotbJet/(bJetP2_*WDaughter1P2_);
  c_=WDaughter1DotbJet/sqrt(bJetP2_*WDaughter1P2_);
  s2_ = 1.-c2_;
  s_=sqrt(s2_);
  //cout << "Cosine: " << c_ << endl;
  //cout << "Sine  : " << s_ << endl;
}

void WDaughterEllipseCalculator::initializeMatrices()
{
  //cout << "initializing matrices" << endl;

  Ab_.Zero();
  AWDaughter1_.Zero();
  Htilde_.Zero();
  H_.Zero();
  Hperp_.Zero();
  HperpInv_.Zero();
  Nperp_.Zero();

  WDaughterPerp_.Zero();
  pWDaughter_.Zero();

  //cout << "the b-jet matrix has " << Ab_.GetNrows() << " rows and " << Ab_.GetNcols() << " columns"  << endl;
  //cout << "the measuredWDaughter matrix has " << AWDaughter1_.GetNrows() << " rows and " << AWDaughter1_.GetNcols() << " columns"  << endl;
}

void WDaughterEllipseCalculator::Wsurface()
{
  //cout << "calculating the W surface" << endl;
  x0p_=-(0.5/  bJetE_)*(mTop2_- mW2_-bJetMass2_);
  x0_ =-(0.5/WDaughter1E_)*(mW2_-WDaughter2Mass2_-WDaughter1Mass2_);
  Sx_=x0_/WDaughter1Beta_-WDaughter1P_/WDaughter1Beta2_+WDaughter1P_;
  epsilon2_=(mW2_-WDaughter2Mass2_)-WDaughter1Beta2_*(mW2_-WDaughter2Mass2_);
  //cout << "b-jet energy is " << bJetE_ << endl;
  //cout << "b-jet mass^2 is " << bJetMass2_ << endl;
  //cout << "mW^2 is " << mW2_ << endl;
  //cout << "mTop^2 is " << mTop2_ << endl;
  //cout << "x0p is " << x0p_ << endl;
  //cout << "1st W daughter energy is " << WDaughter1E_ << endl;
  //cout << "1st W daughter mass^2 is " << WDaughter1Mass2_ << endl;
  //cout << "2nd W daughter mass^2 is " << WDaughter2Mass2_ << endl;
  //cout << "x0 is " << x0_ << endl;
  //cout << "c is " << c_ << endl;
  //cout << "s is " << s_ << endl;
  //cout << "Sx is " << Sx_ << endl;
  //cout << "epsilon^2 is " << epsilon2_ << endl;
}

void WDaughterEllipseCalculator::bJetEllipsoid()
{
  //cout << "calculating the b-jet ellipsoid" << endl;

  Ab_[0][0]=1-c2_*bJetBeta2_;
  Ab_[1][0]=-c_*s_*bJetBeta2_;
  Ab_[2][0]=0;
  Ab_[3][0]=c_*x0p_*bJetBeta_;
    
  Ab_[0][1]=-c_*s_*bJetBeta2_;
  Ab_[1][1]=1-s2_*bJetBeta2_;
  Ab_[2][1]=0;
  Ab_[3][1]=s_*x0p_*bJetBeta_;
    
  Ab_[0][2]=0;
  Ab_[1][2]=0;
  Ab_[2][2]=1;
  Ab_[3][2]=0;
    
  Ab_[0][3]=c_*x0p_*bJetBeta_;
  Ab_[1][3]=s_*x0p_*bJetBeta_;
  Ab_[2][3]=0;
  Ab_[3][3]=mW2_-x0p_*x0p_;

  //cout << "Measured b-jet ellipsoid:" << endl;
  //Ab_.Print();
}

void WDaughterEllipseCalculator::measuredWDaughterEllipsoid()
{
  //cout << "calculating the measured W daughter ellipsoid" << endl;

  AWDaughter1_[0][0]=1.-WDaughter1Beta2_;
  AWDaughter1_[1][0]=0;
  AWDaughter1_[2][0]=0;
  AWDaughter1_[3][0]=Sx_*WDaughter1Beta2_;
    
  AWDaughter1_[0][1]=0;
  AWDaughter1_[1][1]=1;
  AWDaughter1_[2][1]=0;
  AWDaughter1_[3][1]=0;
    
  AWDaughter1_[0][2]=0;
  AWDaughter1_[1][2]=0;
  AWDaughter1_[2][2]=1;
  AWDaughter1_[3][2]=0;
    
  AWDaughter1_[0][3]=Sx_*WDaughter1Beta2_;
  AWDaughter1_[1][3]=0;
  AWDaughter1_[2][3]=0;
  AWDaughter1_[3][3]=mW2_-x0_*x0_-epsilon2_;

  //cout << "Measured W daughter ellipsoid:" << endl;
  //AWDaughter1_.Print();
}

void WDaughterEllipseCalculator::calcZ2()
{
  //cout << "Calculating Z^2" << endl;
  Sy_=(1./s_)*(x0p_/bJetBeta_-c_*Sx_);
  //cout << "Sy is " << Sy_ << endl;
  omega_=(1./s_)*(WDaughter1Beta_/bJetBeta_-c_); //only the positive slope
  //cout << "omega is " << omega_ << endl;
  double Omega2=omega_*omega_+1.-WDaughter1Beta2_;
  //cout << "Omega^2 is " << Omega2 << endl;
  Omega_=sqrt(Omega2);
  //cout << "Omega is " << Omega_ << endl;
  x1_=Sx_-(1./Omega2)*(Sx_+omega_*Sy_);
  y1_=Sy_-(1./Omega2)*omega_*(Sx_+omega_*Sy_);
  Z2_=x1_*x1_*Omega2-(Sy_-omega_*Sx_)*(Sy_-omega_*Sx_)-(mW2_-x0_*x0_-epsilon2_);
  //if(Z2_ > 0) cout << "Z^2 is positive: " << Z2_ << endl;
  //else if(Z2_ < 0) cout << "Z^2 is negative: " << Z2_ << endl;
  //else if(Z2_ == 0) cout << "Z^2 is exactly zero" << endl;
  //else cout << "Z^2 is nan" << endl;
  //if(Z2_<0.) Htilde_.Print();
}

void WDaughterEllipseCalculator::WDaughterSolution()
{
  //cout << "calculating second W daughter ellipse" << endl;
  //cout << "Current top mass is " << mTop_ << endl;

  calcZ2();
  
  //cout << "checking sign of Z^2: " << Z2_ << endl;

  if( Z2_<0 ) 
    {
      //cout << "Z^2 is negative !!" << endl;
      //cout << "Z^2 = " << Z2_ << endl;
      
      // Vary the top mass to force Z^2>0
      // x0p_ depends on mTop_, so Wsurface and bJetEllipsoid need to be called to propagate the changes to the top mass

      errorFlag_ = true;
      Z2_ = 0;

//      calcTopMassCorrection(); //calls to Wsurface and bJetEllipsoid
//
//      if( Z2_<0 ) Z2_=0; //for now reject events 
      
    }
  else 
    {
      //cout << "Z^2 is positive, no need to vary the top mass" << endl;
      errorFlag_ = false;
    }
  double Z=sqrt(Z2_);

  Htilde_[0][0]=Z/Omega_;
  Htilde_[0][1]=0;
  Htilde_[0][2]=x1_-WDaughter1P_;
	
  Htilde_[1][0]=omega_*Z/Omega_;
  Htilde_[1][1]=0;
  Htilde_[1][2]=y1_;
  	
  Htilde_[2][0]=0;
  Htilde_[2][1]=Z;
  Htilde_[2][2]=0;

  //Htilde_.Print();

}

TMatrixD WDaughterEllipseCalculator::rotationMatrix(int axis, double angle)
{
  TMatrixD r(3,3);
  r.Zero();
  if (axis!=0 && axis!=1 && axis!=2) return r;
  
  for( int i=0; i<3; i++ ) 
    {
      r[i][i]=cos(angle);
    }

  for( int i=-1; i<=1; i++ )  
    {
      double row=(axis-i)%3; if(row<0) row+=3;
      double col=(axis+i)%3; if(col<0) col+=3;
      r[row][col]=i*sin(angle)+(1-i*i);
    }

  return r;

}

void WDaughterEllipseCalculator::labSystemTransform()
{
  //cout << "Boosting back into the lab frame" << endl;

  //rotate Htilde to H
  TMatrixD Rz=rotationMatrix(2,-WDaughter1Phi_);
  TMatrixD Ry=rotationMatrix(1,0.5*M_PI-WDaughter1Theta_);
  double bJetP[3]={bJetPx_,bJetPy_, bJetPz_};
  TMatrixD bJet_xyz(3,1,bJetP);
  TMatrixD rM(Ry,TMatrixD::kMult,TMatrixD(Rz,TMatrixD::kMult,bJet_xyz));
  double* rA=rM.GetMatrixArray();
  double phi=-TMath::ATan2(rA[2],rA[1]);
  TMatrixD Rx=rotationMatrix(0,phi);

  H_ = TMatrixD(TMatrixD::kTransposed,Rz);
  TMatrixD RyT(TMatrixD::kTransposed,Ry);
  TMatrixD RxT(TMatrixD::kTransposed,Rx);

  H_*=RyT;
  H_*=RxT;
  H_*=Htilde_;

  //calculate Hperp
  double Hvalues[9]={H_[0][0],H_[0][1],H_[0][2],H_[1][0],H_[1][1],H_[1][2],0,0,1};
  TArrayD Harray(9,Hvalues);
  Hperp_.SetMatrixArray(Harray.GetArray());

  //Hperp_.Print();
}

void WDaughterEllipseCalculator::calcWDaughterEllipse()
{
  //cout << "Calculating W daughter ellipse" << endl;
  //setBJetFactors();
  //setMeasuredWDaughterFactors(); 
  //setAngles();
  //printFactors();
  Wsurface();
  //printFactors();
  measuredWDaughterEllipsoid();
  bJetEllipsoid();
  //printFactors();
  WDaughterSolution();
  labSystemTransform();
  //Ab_.Print();
  //AWDaughter1_.Print();
  //printFactors();
}

void WDaughterEllipseCalculator::calcExtendedWDaughterEllipse()
{
  //avoid non-invertible matrices
  if( Z2_<=0 ) 
    {
      //cout << "Z^2 is negative: " << Z2_ << ", cannot calculate extended representation" << endl;
      Nperp_.Zero();
      return; 
    }
  //cout << "Z^2 is positive: " << Z2_ << ", now calculating extended representation" << endl;
  //cout << "Z^2 = " << Z2_ << endl;
  //Hperp_.Print();
  HperpInv_=Hperp_;
  HperpInv_.Invert();
  TMatrixD U(3,3);
  U.Zero();
  U[0][0]=1;
  U[1][1]=1;
  U[2][2]=-1;
  Nperp_=TMatrixD(HperpInv_,TMatrixD::kTransposeMult,TMatrixD(U,TMatrixD::kMult,HperpInv_));
}

TMatrixD* WDaughterEllipseCalculator::getExtendedWDaughterEllipse()
{
  return &Nperp_;
}

TMatrixD* WDaughterEllipseCalculator::getHomogeneousWDaughterEllipse()
{
  return &Hperp_;
}

TVectorD* WDaughterEllipseCalculator::getWDaughterMomentum(double theta)
{
  //calcExtendedWDaughterEllipse();
  //Hperp_.Print();
  double tArray[3]={cos(theta),sin(theta),1.};
  TVectorD t(3,tArray);
  WDaughterPerp_=t;
  WDaughterPerp_*=Hperp_;
  pWDaughter_=WDaughterPerp_;
  pWDaughter_*=HperpInv_;
  pWDaughter_*=H_;
  return &pWDaughter_;
}

//void WDaughterEllipseCalculator::calcTopMassCorrection()
//{ //Find the ranges of squared top mass values in which Z^2>=0:
//  //Z^2=0 is equivalent to a polynomial of degree 2 in mTop^2
//
//  //adapted from the old topSystemChiSquare::getStartingTopMassRange()
//
//  //cout << "Current Z^2 is " << Z2_ << endl;
//
//  if(Z2_ > 0) return;
//
//  if(Z2_ != Z2_) 
//    {
//      cout << "Z2 is nan!" << endl;
//      return;
//    }
//
//  //cout << "Z^2 is negative!" << endl;
//
//  cout << "Calculating allowed values for the top mass" << endl;
//  //cout << "Current top mass is " << mTop_ << endl;
//
//  double WDaughter1P4 = WDaughter1P2_*WDaughter1P2_;
//
//  double bJetE2 = bJetE_*bJetE_;
//
//  double WDaughter1E2 = WDaughter1E_*WDaughter1E_;
//  double WDaughter1E3 = WDaughter1E2*WDaughter1E_;
//  double WDaughter1E4 = WDaughter1E2*WDaughter1E2;
//
//  double c = (WDaughter1Px_*bJetPx_+WDaughter1Py_*bJetPy_+WDaughter1Pz_*bJetPz_);
//  double s2 = 1.-c*c/(bJetP2_*WDaughter1P2_);
//  c/=(bJetP_*WDaughter1P_);
//
//  double mTopEdge1 = (bJetE_*WDaughter1E3 - bJetE_*WDaughter1E_*(WDaughter2Mass2_ - mW2_ + WDaughter1P2_) + WDaughter1E2*(bJetE2 + mW2_ - bJetP_*(bJetP_ + c*WDaughter1P_)) +
//		      WDaughter1P_*((-bJetE2 - mW2_ + bJetP2_)*WDaughter1P_ + c*bJetP_*(WDaughter2Mass2_ - mW2_ + WDaughter1P2_)) -
//		      sqrt((WDaughter1E4 + pow(WDaughter2Mass2_ - mW2_,2) + 2*(WDaughter2Mass2_ + mW2_)*WDaughter1P2_ + WDaughter1P4 - 2*WDaughter1E2*(WDaughter2Mass2_ + mW2_ + WDaughter1P2_))*
//			   (pow(c*WDaughter1E_*bJetP_ - bJetE_*WDaughter1P_,2) + bJetP2_*(WDaughter1E_ - WDaughter1P_)*(WDaughter1E_ + WDaughter1P_)*s2)))/((WDaughter1E_ - WDaughter1P_)*(WDaughter1E_ + WDaughter1P_));
//  double mTopEdge2 = (bJetE_*WDaughter1E3 - bJetE_*WDaughter1E_*(WDaughter2Mass2_ - mW2_ + WDaughter1P2_) + WDaughter1E2*(bJetE2 + mW2_ - bJetP_*(bJetP_ + c*WDaughter1P_)) +
//		      WDaughter1P_*((-bJetE2 - mW2_ + bJetP2_)*WDaughter1P_ + c*bJetP_*(WDaughter2Mass2_ - mW2_ + WDaughter1P2_)) +
//		      sqrt((WDaughter1E4 + pow(WDaughter2Mass2_ - mW2_,2) + 2*(WDaughter2Mass2_ + mW2_)*WDaughter1P2_ + WDaughter1P4 - 2*WDaughter1E2*(WDaughter2Mass2_ + mW2_ + WDaughter1P2_))*
//			   (pow(c*WDaughter1E_*bJetP_ - bJetE_*WDaughter1P_,2) + bJetP2_*(WDaughter1E_ - WDaughter1P_)*(WDaughter1E_ + WDaughter1P_)*s2)))/((WDaughter1E_ - WDaughter1P_)*(WDaughter1E_ + WDaughter1P_));
//
//
//  //cout << "Top mass squared edges: (" << mTopEdge2 << ", " << mTopEdge1 << ")" << endl;
//
//  mTopEdge1 = max(0.,mTopEdge1);
//  mTopEdge2 = max(0.,mTopEdge2);
//
//  mTopEdgeLow_ = min(sqrt(mTopEdge1),sqrt(mTopEdge2));
//  mTopEdgeHigh_ = max(sqrt(mTopEdge1),sqrt(mTopEdge2));
//
//  cout << "Top mass low and high edges: (" << mTopEdgeLow_ << ", " << mTopEdgeHigh_ << ")" << endl;
//
//
//  double mTop=0;
//
//  if(mTop_ < mTopEdgeLow_)
//    {
//      cout << "Case 2" << endl;
//      cout << "Starting at the low edge and incrementing upwards" << endl;
//      mTop=mTopEdgeLow_;
//      setTopMass(mTop);
//      Wsurface();
//      bJetEllipsoid();
//      calcZ2();
//      //cout << "Here Z^2 = " << Z2_ << endl;
//      while(Z2_ < 0.)
//	{
//	  //mTop += 1.e-3*(mTopEdgeHigh_-mTopEdgeLow_);
//	  mTop += 1.e-2*mTopAssumed_;
//	  cout << "Testing top mass: " << mTop << endl;
//	  if(mTop>=mTopEdgeHigh_) 
//	    {
//	      cout << "Reached the high edge !!!" << endl;
//	      break;
//	    }
//	  setTopMass(mTop);
//	  Wsurface();
//	  bJetEllipsoid();
//	  calcZ2();
//	}
//    }
//  else if(mTop_ > mTopEdgeLow_ && mTop_ < mTopEdgeHigh_)
//    {
//      cout << "Case 1" << endl;
//      double inc = (mTop_-mTopEdgeLow_ < mTopEdgeHigh_-mTop_) ? -1. : 1.;
//      if(mTopEdgeLow_ == 0.) inc=1.; //avoid negative mass values
//      if(inc == 1.) 
//	{
//	  cout << "Starting at the high edge and incrementing upwards" << endl;
//	  mTop=mTopEdgeHigh_;
//	}
//      else 
//	{
//	  cout << "Starting at the low edge and incrementing downwards" << endl;
//	  mTop=mTopEdgeLow_;
//	}
//      setTopMass(mTop);
//      Wsurface();
//      bJetEllipsoid();
//      calcZ2();
//      //cout << "Here Z^2 = " << Z2_ << endl;
//      while(Z2_ < 0.)
//	{
//	  //mTop += inc*1.e-3*(mTopEdgeHigh_-mTopEdgeLow_); //no bound (upper or lower) on mTop in this case
//	  mTop += inc*1.e-2*mTopAssumed_;
//	  if(mTop < 0.)
//	    {
//	      cout << "Switching to starting at the high edge and incrementing upwards" << endl;
//	      mTop=mTopEdgeHigh_;
//	      inc=1.;
//	    }
//	  cout << "Testing top mass: " << mTop << endl;
//	  setTopMass(mTop);
//	  Wsurface();
//	  bJetEllipsoid();
//	  calcZ2();
//	}
//    }
//  else
//    {
//      cout << "Case 3" << endl;
//      cout << "Starting at the high edge and incrementing downwards" << endl;
//      mTop=mTopEdgeHigh_;
//      setTopMass(mTop);
//      Wsurface();
//      bJetEllipsoid();
//      calcZ2();
//      //cout << "Here Z^2 = " << Z2_ << endl;
//      while(Z2_ < 0.)
//        {
//          //mTop -= 1.e-3*(mTopEdgeHigh_-mTopEdgeLow_);
//	  mTop -= 1.e-2*mTopAssumed_;
//	  cout << "Testing top mass: " << mTop << endl;
//          if(mTop<=mTopEdgeLow_)
//            {
//              cout << "Reached the low edge !!!" << endl;
//              break;
//            }
//          setTopMass(mTop);
//          Wsurface();
//          bJetEllipsoid();
//          calcZ2();
//        }
//    }
//
//  cout << "Corrected top mass is: " << mTop_ << endl;
//  cout << "New Z^2 is: " << Z2_ << endl;
//
//  return;
//}

void WDaughterEllipseCalculator::printFactors()
{
  //cout << "x0 = " << x0_ << endl;
  //cout << "x0p = " << x0p_ << endl;
  //cout << "Sx = " << Sx_ << endl;
  //cout << "Sy = " << Sy_ << endl;
  //cout << "eps^2 = " << epsilon2_ << endl;
  //cout << "top mass = " << mTop_ << endl;
  //cout << "top mass squared = " << mTop2_ << endl;
  //cout << "W mass squared = " << mW2_ << endl;
  //cout << "b-jet mass squared = " << bJetMass2_ << endl;
  //cout << "measured W daughter E  = " << WDaughter1E_ << endl;
  //cout << "measured W daughter m2 = " << WDaughter1E_*WDaughter1E_-WDaughter1Px_*WDaughter1Px_-WDaughter1Py_*WDaughter1Py_-WDaughter1Pz_*WDaughter1Pz_ << endl;
  //cout << "using m2 = " << WDaughter1Mass2_ << endl;
  //cout << "using second W daughter m2 = "<< WDaughter2Mass2_ << endl;
    std::cout << "b-jet momentum:"
       << "\npx = " << bJetPx_ 
       << "\npy = " << bJetPy_ 
       << "\npz = " << bJetPz_ 
       << "\nE  = " << bJetE_  << std::endl;
    std::cout << "first W daughter momentum:"
       << "\npx = " << WDaughter1Px_ 
       << "\npy = " << WDaughter1Py_ 
       << "\npz = " << WDaughter1Pz_ 
       << "\nE  = " << WDaughter1E_  << std::endl;
}
