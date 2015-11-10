#include <cp3_llbb/TTAnalysis/interface/neutrinoSolutions.h>

using namespace ROOT::Math;
using namespace std;
double neutrinoSolutions::splittingRightEdge(double* mt, double *)
{
  if(*mt != mt1 && *mt != mt2) setMasses(*mt,*mt,mW1,mW2);
  TF1 neutrinoEImPartSplitting("neutrinoEImPartSplitting",this,&neutrinoSolutions::imagnessTopSplitting,0,1,0,"neutrinoEImPartSplitting","neutrinoEImPartSplitting");
  if(neutrinoEImPartSplitting.Eval(0)<0) return -1;
  double massEdgeLow = 0;
  double massEdgeHigh = 0;
  for(int i = 1; i<101; i++)
    {
      massEdgeHigh = *mt*(double)i*0.01;
      if(neutrinoEImPartSplitting.Eval(massEdgeHigh) < 0) break;
      massEdgeLow = massEdgeHigh;
    }
  //cout << "the edge is between " << massEdgeLow << " with reality " << neutrinoEImPartSplitting.Eval(massEdgeLow) << endl;
  //cout << " and " << massEdgeHigh  << " with reality " << neutrinoEImPartSplitting.Eval(massEdgeHigh) << endl;
  if(neutrinoEImPartSplitting.Eval(massEdgeHigh) > 0) return -1.;
  //cout << "found imaginary part" << endl;
  //if(*mt != mt1 && *mt != mt2) setMasses(*mt,*mt,mW1,mW2);
  //cout << "checking 0: " << neutrinoEImPartSplitting.Eval(0) << endl;
  //cout << "checking " << *mt << ": " << neutrinoEImPartSplitting.Eval(*mt) << endl;
  //if((neutrinoEImPartSplitting.Eval(0)<0 && neutrinoEImPartSplitting.Eval(*mt)<0) || (neutrinoEImPartSplitting.Eval(0)>0 && neutrinoEImPartSplitting.Eval(*mt)>0))
  //  {
  //    return -1.;
  //  }
  WrappedTF1 WrappedNeutrinoImPartSplitting(neutrinoEImPartSplitting);
  Roots::Bisection neutrinoEZeros;
  neutrinoEZeros.SetFunction( WrappedNeutrinoImPartSplitting, massEdgeLow , massEdgeHigh );
  bool hasRightRoot = neutrinoEZeros.Solve();
  double splittingRightEdgeValue = !hasRightRoot ? -1. : neutrinoEZeros.Root();
  //cout << "actual splitting value is " << splittingRightEdgeValue << endl;
  return splittingRightEdgeValue;
}

double neutrinoSolutions::splittingLeftEdge(double* mt, double *)
{
  if(*mt != mt1 && *mt != mt2) setMasses(*mt,*mt,mW1,mW2);
  TF1 neutrinoEImPartSplitting("neutrinoEImPartSplitting",this,&neutrinoSolutions::imagnessTopSplitting,0,1,0,"neutrinoEImPartSplitting","neutrinoEImPartSplitting");
  if(neutrinoEImPartSplitting.Eval(0)<0) return 1;
  double massEdgeLow = 0;
  double massEdgeHigh = 0;
  for(int i = 1; i<101; i++)
    {
      massEdgeLow = -*mt*(double)i*0.01;
      if(neutrinoEImPartSplitting.Eval(massEdgeLow) < 0) break;
      massEdgeHigh = massEdgeLow;
    }
  //cout << "the edge is between " << massEdgeLow << " with reality " << neutrinoEImPartSplitting.Eval(massEdgeLow) << endl;
  //cout << " and " << massEdgeHigh  << " with reality " << neutrinoEImPartSplitting.Eval(massEdgeHigh) << endl;
  if(neutrinoEImPartSplitting.Eval(massEdgeLow) > 0) return 1.;
  //cout << "found imaginary part" << endl;
  //if(*mt != mt1 && *mt != mt2) setMasses(*mt,*mt,mW1,mW2);
  //cout << "checking 0: " << neutrinoEImPartSplitting.Eval(0) << endl;
  //cout << "checking " << -*mt << ": " << neutrinoEImPartSplitting.Eval(-*mt) << endl;
  //if((neutrinoEImPartSplitting.Eval(0)<0 && neutrinoEImPartSplitting.Eval(-*mt)<0) || (neutrinoEImPartSplitting.Eval(0)>0 && neutrinoEImPartSplitting.Eval(-*mt)>0))
  //  {
  //    return 1.;
  //  }
  WrappedTF1 WrappedNeutrinoImPartSplitting(neutrinoEImPartSplitting);
  Roots::Bisection neutrinoEZeros;
  neutrinoEZeros.SetFunction( WrappedNeutrinoImPartSplitting, massEdgeLow , massEdgeHigh );

  bool hasLeftRoot = neutrinoEZeros.Solve(1000);
  double splittingLeftEdgeValue = !hasLeftRoot ? 1. : neutrinoEZeros.Root();
  //cout << "actual splitting value is " << splittingLeftEdgeValue << endl;

  return splittingLeftEdgeValue;
}


double neutrinoSolutions::splittingEdgeDifference(double* mt, double *)
{
  if(*mt != mt1 && *mt != mt2) setMasses(*mt,*mt,mW1,mW2);
  //cout << "setting masses to " << *mt << endl;
  double reality = imagness(*mt,*mt,mW1,mW2);
  if(reality < 1)
    {
      cout << "No real solution for top mass " << *mt << endl; 
      double realityLow = imagness(*mt-*mt/10,*mt-*mt/10,mW1,mW2);
      cout << "realness of 10% low: " << realityLow << endl;
      double realityHigh = imagness(*mt+*mt/10,*mt+*mt/10,mW1,mW2);
      cout << "realness of 10% high: " << realityLow << endl;
      setMasses(*mt,*mt,mW1,mW2);
      if((realityLow<0 && realityHigh>0) || (realityLow>0 && realityHigh<0))
	{
	  TF1 neutrinoEImPartMass("neutrinoEImPartMass",this,&neutrinoSolutions::imagnessTopMass,0,1,0,"neutrinoEImPartMass","neutrinoEImPartMass");
	  WrappedTF1 WrappedNeutrinoImPart(neutrinoEImPartMass);
	  Roots::Bisection neutrinoEZeros;
	  neutrinoEZeros.SetFunction( WrappedNeutrinoImPart, *mt-*mt/10 , *mt+*mt/10 );
	  bool hasEdge = neutrinoEZeros.Solve(1000);
	  double edge = !hasEdge ? 0 : neutrinoEZeros.Root();
	  if(!hasEdge) return numeric_limits<double>::quiet_NaN();
	  double newMass = edge-(edge - *mt)/2;
	  return splittingEdgeDifference(&newMass, NULL);
	}
      else return 0;//numeric_limits<double>::quiet_NaN();
    }
  double splittingRightEdgeValue = splittingRightEdge(mt,NULL);
  double splittingLeftEdgeValue  = splittingLeftEdge(mt,NULL);

  //cout << "right edge at " << splittingRightEdgeValue << endl;
  //cout << "left edge at " << splittingLeftEdgeValue << endl;

  //cout << "high edge is " << splittingRightEdgeValue << endl;
  //cout << "low edge is " << splittingLeftEdgeValue << endl;

  if(splittingLeftEdgeValue > 0 ) splittingLeftEdgeValue = 0;
  if(splittingRightEdgeValue < 0 ) splittingRightEdgeValue = 0;

  return splittingRightEdgeValue+splittingLeftEdgeValue;
}

neutrinoSolutions::neutrinoSolutions(const double& b1E_, const double& b1px_, const double& b1py_, const double& b1pz_, 
				     const double& l1E_, const double& l1px_, const double& l1py_, const double& l1pz_, 
				     const double& b2E_, const double& b2px_, const double& b2py_, const double& b2pz_, 
				     const double& l2E_, const double& l2px_, const double& l2py_, const double& l2pz_, 
				     const double& METx_, const double& METy_, const double& mW0_, const double& mt0_):
                                     b1E (b1E_ ),
  				     b1px(b1px_),
  				     b1py(b1py_),
  				     b1pz(b1pz_),
  				     b2E (b2E_ ),
  				     b2px(b2px_),
  				     b2py(b2py_),
  				     b2pz(b2pz_),
  				     l1E (l1E_ ),
  				     l1px(l1px_),
  				     l1py(l1py_),
  				     l1pz(l1pz_),
  				     l2E (l2E_ ),
  				     l2px(l2px_),
  				     l2py(l2py_),
  				     l2pz(l2pz_),
  				     METx(METx_),
  				     METy(METy_),
  				     mW0 (mW0_ ),
  				     mt0 (mt0_ ),
				     mW1(mW0_),
				     mW2(mW0_),
				     mt1(mt0_),
				     mt2(mt0_),
				     mtBar(mt0_),
				     mtDelta(0),
				     mW1Temp(mW0_+1),
				     mW2Temp(mW0_+1),
				     mt1Temp(mt0_+1),
				     mt2Temp(mt0_+1)
{
  getNeutrinoVectors();
}

neutrinoSolutions::neutrinoSolutions():
                                     b1E (0),
  				     b1px(0),
  				     b1py(0),
  				     b1pz(0),
  				     b2E (0),
  				     b2px(0),
  				     b2py(0),
  				     b2pz(0),
  				     l1E (0),
  				     l1px(0),
  				     l1py(0),
  				     l1pz(0),
  				     l2E (0),
  				     l2px(0),
  				     l2py(0),
  				     l2pz(0),
  				     METx(0),
  				     METy(0),
  				     mW0 (0),
  				     mt0 (0),
				     mW1(0),
				     mW2(0),
				     mt1(0),
				     mt2(0),
				     mtBar(0),
				     mtDelta(0),
				     mW1Temp(1.),
				     mW2Temp(1.),
				     mt1Temp(1.),
				     mt2Temp(1.)
{
  getNeutrinoVectors();
}

void neutrinoSolutions::setupMeasurements(const double b1px_, const double b1py_, const double b1pz_, const double b1E_, 
					  const double l1px_, const double l1py_, const double l1pz_, const double l1E_, 
					  const double b2px_, const double b2py_, const double b2pz_, const double b2E_, 
					  const double l2px_, const double l2py_, const double l2pz_, const double l2E_, 
					  const double METx_, const double METy_, 
					  const double mW1_,  const double mW2_ ,
					  const double mt1_,  const double mt2_ )
{
  b1E  = b1E_ ;
  b1px = b1px_;
  b1py = b1py_;
  b1pz = b1pz_;
  b2E  = b2E_ ;
  b2px = b2px_;
  b2py = b2py_;
  b2pz = b2pz_;
  l1E  = l1E_ ;
  l1px = l1px_;
  l1py = l1py_;
  l1pz = l1pz_;
  l2E  = l2E_ ;
  l2px = l2px_;
  l2py = l2py_;
  l2pz = l2pz_;
  METx = METx_;
  METy = METy_;
  mW1 = mW1_;
  mW2 = mW2_;
  mt1 = mt1_;
  mt2 = mt2_;
  mW1Temp = mW1_+1;
  mW2Temp = mW2_+1;
  mt1Temp = mt1_+1;
  mt2Temp = mt2_+1;
  getNeutrinoVectors();
}

neutrinoSolutions::~neutrinoSolutions()
{
  //delete imSolution;
}

double neutrinoSolutions::topPairMass(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_,const int& i)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();

  vector<double> pairMasses;

  for(unsigned int j = 1; j < nu1E.size(); j++) if(abs( nu1E[i].real()) > 1e9 * abs( nu1E[i].imag()) /*nu1E[j].imag() == 0*/) {
      double topPairE  = ( nu1E[j].real() +  nu2E[j].real() + b1E  + l1E  + b2E  + l2E );
      double topPairPx = (nu1px[j].real() + nu2px[j].real() + b1px + l1px + b2px + l2px);
      double topPairPy = (nu1py[j].real() + nu2py[j].real() + b1py + l1py + b2py + l2py);
      double topPairPz = (nu1pz[j].real() + nu2pz[j].real() + b1pz + l1pz + b2pz + l2pz);

      pairMasses.push_back(topPairE*topPairE - topPairPx*topPairPx - topPairPy*topPairPy - topPairPz*topPairPz);

      //cout <<"pair mass squared is " << topPairE*topPairE - topPairPx*topPairPx - topPairPy*topPairPy - topPairPz*topPairPz << endl;
    }

  //return pairMasses.size();

  sort(pairMasses.begin(),pairMasses.end());
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  if(pairMasses.size() == 0) return 0;

  if(i < 0 ) return (sqrt(pairMasses[0]) - 2 * mtBar);

  unsigned int index = abs(i);

  if(!(pairMasses.size() > index)) return 0;

  return (sqrt(pairMasses[pairMasses.size() - 1 - index]) - 2 * mtBar);
}

double neutrinoSolutions::topPairMassHigh(double* mt, double* )
{
  return topPairMass(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2,0);
}

double neutrinoSolutions::topPairMassLow(double* mt, double* )
{
  return topPairMass(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2,-1);
}

double neutrinoSolutions::topPairPzFunction(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_,const int& i)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();

  vector<double> pairPzs;

  for(unsigned int j = 1; j < nu1E.size(); j++) if(/*abs( nu1E[i].real()) > 1e9 * abs( nu1E[i].imag())*/ nu1E[j].imag() == 0) {
      double topPairPz = (nu1pz[j].real() + nu2pz[j].real() + b1pz + l1pz + b2pz + l2pz);

      pairPzs.push_back(topPairPz);

      //cout <<"pair mass squared is " << topPairE*topPairE - topPairPx*topPairPx - topPairPy*topPairPy - topPairPz*topPairPz << endl;
    }
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  //return pairPzs.size();

  sort(pairPzs.begin(),pairPzs.end());

  if(pairPzs.size() == 0) return 0;

  if(i < 0 ) return (pairPzs[0]);

  unsigned int index = abs(i);

  if(!(pairPzs.size() > index)) return 0;

  return (pairPzs[pairPzs.size() - 1 - index]);
}

double neutrinoSolutions::topPairPzMass(double* mt, double* )
{
  return topPairPzFunction(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2,0);
}

double neutrinoSolutions::topPairPzMass1(double* mt, double* )
{
  return topPairPzFunction(*mt,mt2,mW1,mW2,0);
}

double neutrinoSolutions::topPairPzMass2(double* mt, double* )
{
  return topPairPzFunction(mt1,*mt,mW1,mW2,0);
}

double neutrinoSolutions::topPairPzSplitting(double* mt, double* )
{
  return topPairPzFunction(mtBar + 0.5* *mt, mtBar - 0.5 * *mt ,mW1,mW2,0);
}

//double neutrinoSolutions::imagness(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_)
//{
//  setMasses(mt1_,mt2_,mW1_,mW2_);
//  getNeutrinoVectors();
//  double imagpart = nu1E[0].imag();
//  for(unsigned int i = 1; i < nu1E.size(); i++)
//    {
//      if(abs(nu1E[i].imag())<=abs(imagpart)) imagpart = nu1E[i].imag();
//    }
//  return imagpart;
//}

double neutrinoSolutions::nu1E2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();
  double imagpart2 = nu1E[0].imag()*nu1E[0].imag();
  double realpart2 = nu1E[0].real()*nu1E[0].real();
  for(unsigned int i = 1; i < nu1E.size(); i++)
    {
      double imag2=nu1E[i].imag()*nu1E[i].imag();
      if(imag2<imagpart2) {imagpart2 = imag2;realpart2 = nu1E[0].real()*nu1E[0].real();}
    }
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  return realpart2;
}

double neutrinoSolutions::nu2E2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();
  double imagpart2 = nu2E[0].imag()*nu2E[0].imag();
  double realpart2 = nu2E[0].real()*nu2E[0].real();
  for(unsigned int i = 1; i < nu1E.size(); i++)
    {
      double imag2=nu2E[i].imag()*nu2E[i].imag();
      if(imag2<imagpart2) {imagpart2 = imag2;realpart2 = nu2E[0].real()*nu2E[0].real();}
    }
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  return realpart2 + imagpart2;
}

double neutrinoSolutions::imagness2(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();
  double imagpart2 = nu1E[0].imag()*nu1E[0].imag();
  for(unsigned int i = 1; i < nu1E.size(); i++)
    {
      double imag2=nu1E[i].imag()*nu1E[i].imag();
      if(imag2<imagpart2) imagpart2 = imag2;
    }
  //return (imagpart2>0) ? imagpart2 : -1;
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  return (imagpart2>0) ? 1 : -1;
}

double neutrinoSolutions::imagness(const double& mt1_, const double& mt2_, const double& mW1_, const double& mW2_)
{
  double temp_mt1(mt1),temp_mt2(mt2),temp_mW1(mW1),temp_mW2(mW2);
  setMasses(mt1_,mt2_,mW1_,mW2_);
  getNeutrinoVectors();
  double realness = -1;
  for(unsigned int i = 0; i < nu1E.size(); i++)
    {
      //if(nu1E[i].imag() ==0){realness = 1;}
      if( abs( nu1E[i].real()) > 1e9 * abs( nu1E[i].imag()) && abs( nu2E[i].real()) > 1e9 * abs( nu2E[i].imag()) &&
	  abs(nu1px[i].real()) > 1e9 * abs(nu1px[i].imag()) && abs(nu2px[i].real()) > 1e9 * abs(nu2px[i].imag()) &&
	  abs(nu1py[i].real()) > 1e9 * abs(nu1py[i].imag()) && abs(nu2py[i].real()) > 1e9 * abs(nu2py[i].imag()) &&
	  abs(nu1pz[i].real()) > 1e9 * abs(nu1pz[i].imag()) && abs(nu2pz[i].real()) > 1e9 * abs(nu2pz[i].imag())  )realness = 1;
      //cout << nu1E[i].imag() << endl;
    }
  setMasses(temp_mt1,temp_mt2,temp_mW1,temp_mW2);
  return realness;
}

double neutrinoSolutions::imagnessTopMass(double* mt, double* )
{
  return imagness(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2);
}

double neutrinoSolutions::nu1E2TopMass(double* mt, double* )
{
  return nu1E2(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2);
}

double neutrinoSolutions::nu1E2TopSplitting(double* mt, double* )
{
  return nu1E2(mtBar+0.5* *mt,mtBar-0.5* *mt,mW1,mW2);
}

//double neutrinoSolutions::nu2E2TopMass(double* mt, double* )
//{
//  return nu2E2(*mt+0.5*mtDelta,*mt-0.5*mtDelta,mW1,mW2);
//}

double neutrinoSolutions::imagnessTopSplitting(double* mt, double* )
{
  return imagness(mtBar+0.5* *mt,mtBar-0.5* *mt,mW1,mW2);
}

double neutrinoSolutions::imagnessTopSplitting2D(double* mt, double* )
{
  return imagness(mt[0]+mt[1],mt[0]-mt[1],mW1,mW2);
}

double neutrinoSolutions::imagnessTopMass1(double* mt, double* )
{
  return imagness(*mt,mt2,mW1,mW2);
}

double neutrinoSolutions::imagnessTopMass2(double* mt, double* )
{
  return imagness(mt1,*mt,mW1,mW2);
}

double neutrinoSolutions::imagnessTopMass2D(double* mt, double* )
{
  return imagness(mt[0],mt[1],mW1,mW2);
}

void neutrinoSolutions::getRealNeutrinoVectors(vector<double>& nu1E_real, vector<double>& nu1px_real, vector<double>& nu1py_real, vector<double>& nu1pz_real, 
					       vector<double>& nu2E_real, vector<double>& nu2px_real, vector<double>& nu2py_real, vector<double>& nu2pz_real)
{
  getNeutrinoVectors();
  for(unsigned int i = 0; i < nu1E.size(); i++)
    {
      //cout << "Solution is " << nu1E[i] << endl;
      if( abs( nu1E[i].real()) > 1e9 * abs( nu1E[i].imag()) && abs( nu2E[i].real()) > 1e9 * abs( nu2E[i].imag()) &&
	  abs(nu1px[i].real()) > 1e9 * abs(nu1px[i].imag()) && abs(nu2px[i].real()) > 1e9 * abs(nu2px[i].imag()) &&
	  abs(nu1py[i].real()) > 1e9 * abs(nu1py[i].imag()) && abs(nu2py[i].real()) > 1e9 * abs(nu2py[i].imag()) &&
	  abs(nu1pz[i].real()) > 1e9 * abs(nu1pz[i].imag()) && abs(nu2pz[i].real()) > 1e9 * abs(nu2pz[i].imag()) &&
	  //nu1E [i].imag() == 0. && 
	  //nu1px[i].imag() == 0. &&
	  //nu1py[i].imag() == 0. &&
	  //nu1pz[i].imag() == 0. &&
	  //nu2E [i].imag() == 0. &&
	  //nu2px[i].imag() == 0. &&
	  //nu2py[i].imag() == 0. &&
	  //nu2pz[i].imag() == 0. &&
	  nu1E[i].real() > 0.0 && nu2E[i].real() > 0.0 )
	{
	  nu1E_real.push_back(  nu1E[i].real() );
	 nu1px_real.push_back( nu1px[i].real() );
	 nu1py_real.push_back( nu1py[i].real() );
	 nu1pz_real.push_back( nu1pz[i].real() );
          nu2E_real.push_back(  nu2E[i].real() );
	 nu2px_real.push_back( nu2px[i].real() );
	 nu2py_real.push_back( nu2py[i].real() );
	 nu2pz_real.push_back( nu2pz[i].real() );
	}
    }
}

void neutrinoSolutions::getNeutrinoVectors()
{  
  if( mW1Temp == mW1 &&
      mW2Temp == mW2 &&
      mt1Temp == mt1 &&
      mt2Temp == mt2  ) return;

  //cout << "State of the mass system: " << endl;
  //cout << "mt1: " << mt1 << endl;
  //cout << "mt2: " << mt2 << endl;
  //cout << "mW1: " << mW1 << endl;
  //cout << "mW2: " << mW2 << endl;

  mW1Temp = mW1;
  mW2Temp = mW2;
  mt1Temp = mt1;
  mt2Temp = mt2;

  nu1E.clear();
  nu2E.clear();
  nu1px.clear();
  nu1py.clear();
  nu1pz.clear();
  nu2px.clear();
  nu2py.clear();
  nu2pz.clear();  

  //First solve equations
  //    (p_nu1 + p_l1) = (m_W)^2 - (m_l1)^2                   |  (p_nu2 + p_l2) = (m_W)^2 - (m_l2)^2                  
  //    (p_nu1 + p_l1 + p_b1) = (m_t)^2 - (m_W)^2 - (m_l1)^2  |  (p_nu2 + p_l2 + p_b1) = (m_t)^2 - (m_W)^2 - (m_l2)^2 
  //    p_nu1_x + p_nu2_x = METx                              |  p_nu1_y + p_nu2_y = METy                               
  
  //double mW = 80.4;
  //double mt = 173.;

  //Here we could put in mass for a neutrino!
  //Solving involve having ml_i^2 -> ml_i^2 + mnu_i^2

  double mW12 = mW1*mW1;
  double mt12 = mt1*mt1;
  double mW22 = mW2*mW2;
  double mt22 = mt2*mt2;
  double ml12 = l1E*l1E - l1px*l1px - l1py*l1py - l1pz*l1pz;
  double ml22 = l2E*l2E - l2px*l2px - l2py*l2py - l2pz*l2pz;
  double mb12 = b1E*b1E - b1px*b1px - b1py*b1py - b1pz*b1pz;
  double mb22 = b2E*b2E - b2px*b2px - b2py*b2py - b2pz*b2pz;
  //cout << "mb12: " << mb12 << endl;
  //cout << "mb22: " << mb22 << endl;
  //cout << "ml12: " << ml12 << endl;
  //cout << "ml22: " << ml22 << endl;
  
  double nuDenom = 2*(b1pz*b2pz*l1py*l2px - b1py*b2pz*l1pz*l2px - b1pz*b2pz*l1px*l2py + b1px*b2pz*l1pz*l2py + b1pz*b2py*l1px*l2pz - b1pz*b2px*l1py*l2pz + b1py*b2px*l1pz*l2pz - b1px*b2py*l1pz*l2pz);
  
  //cout << "nuDenom is " << nuDenom << endl;
  
  //nuPi = nuPi0 + nuPi1*nuE1 + nuPi2
  
  double nu1Px0 = 2*b1E*b2pz*l1E*l1pz*l2py - 2*b1px*b2pz*l1px*l1pz*l2py - 2*b1py*b2pz*l1py*l1pz*l2py - 2*b1pz*b2pz*pow(l1pz,2)*l2py - 2*b1E*b2py*l1E*l1pz*l2pz + 
    2*b1px*b2py*l1px*l1pz*l2pz + 2*b1py*b2py*l1py*l1pz*l2pz + 2*b1pz*b2py*pow(l1pz,2)*l2pz + 2*b1pz*b2E*l1py*l2E*l2pz - 
    2*b1py*b2E*l1pz*l2E*l2pz - 2*b1pz*b2px*l1py*l2px*l2pz + 2*b1py*b2px*l1pz*l2px*l2pz - 2*b1pz*b2py*l1py*l2py*l2pz + 2*b1py*b2py*l1pz*l2py*l2pz - 
    2*b1pz*b2pz*l1py*pow(l2pz,2) + 2*b1py*b2pz*l1pz*pow(l2pz,2) + b2pz*l1pz*l2py*mb12 - b2py*l1pz*l2pz*mb12 + b1pz*l1py*l2pz*mb22 - 
    b1py*l1pz*l2pz*mb22 + 2*b1pz*b2pz*l1py*l2px*METx - 2*b1py*b2pz*l1pz*l2px*METx - 2*b1pz*b2px*l1py*l2pz*METx + 2*b1py*b2px*l1pz*l2pz*METx + 
    2*b1pz*b2pz*l1py*l2py*METy - 2*b1py*b2pz*l1pz*l2py*METy - 2*b1pz*b2py*l1py*l2pz*METy + 2*b1py*b2py*l1pz*l2pz*METy - b1pz*b2pz*l2py*ml12 + 
    b1pz*b2py*l2pz*ml12 - b1pz*b2pz*l1py*ml22 + b1py*b2pz*l1pz*ml22 - b2pz*l1pz*l2py*mt12 + b2py*l1pz*l2pz*mt12 - b1pz*l1py*l2pz*mt22 + 
    b1py*l1pz*l2pz*mt22 + b1pz*b2pz*l2py*mW12 + b2pz*l1pz*l2py*mW12 - b1pz*b2py*l2pz*mW12 - b2py*l1pz*l2pz*mW12 + b1pz*b2pz*l1py*mW22 - 
    b1py*b2pz*l1pz*mW22 + b1pz*l1py*l2pz*mW22 - b1py*l1pz*l2pz*mW22;
  double nu1Px1 = -2*b1pz*b2pz*l1E*l2py + 2*b1E*b2pz*l1pz*l2py + 2*b1pz*b2py*l1E*l2pz - 2*b1E*b2py*l1pz*l2pz;
  double nu1Px2 = -2*b1pz*b2pz*l1py*l2E + 2*b1py*b2pz*l1pz*l2E + 2*b1pz*b2E*l1py*l2pz - 2*b1py*b2E*l1pz*l2pz;
  
  double nu1Py0 = -2*b1E*b2pz*l1E*l1pz*l2px + 2*b1px*b2pz*l1px*l1pz*l2px + 2*b1py*b2pz*l1py*l1pz*l2px + 2*b1pz*b2pz*pow(l1pz,2)*l2px + 2*b1E*b2px*l1E*l1pz*l2pz - 
   2*b1px*b2px*l1px*l1pz*l2pz - 2*b1py*b2px*l1py*l1pz*l2pz - 2*b1pz*b2px*pow(l1pz,2)*l2pz - 2*b1pz*b2E*l1px*l2E*l2pz + 
   2*b1px*b2E*l1pz*l2E*l2pz + 2*b1pz*b2px*l1px*l2px*l2pz - 2*b1px*b2px*l1pz*l2px*l2pz + 2*b1pz*b2py*l1px*l2py*l2pz - 2*b1px*b2py*l1pz*l2py*l2pz + 
   2*b1pz*b2pz*l1px*pow(l2pz,2) - 2*b1px*b2pz*l1pz*pow(l2pz,2) - b2pz*l1pz*l2px*mb12 + b2px*l1pz*l2pz*mb12 - b1pz*l1px*l2pz*mb22 + 
   b1px*l1pz*l2pz*mb22 - 2*b1pz*b2pz*l1px*l2px*METx + 2*b1px*b2pz*l1pz*l2px*METx + 2*b1pz*b2px*l1px*l2pz*METx - 2*b1px*b2px*l1pz*l2pz*METx - 
   2*b1pz*b2pz*l1px*l2py*METy + 2*b1px*b2pz*l1pz*l2py*METy + 2*b1pz*b2py*l1px*l2pz*METy - 2*b1px*b2py*l1pz*l2pz*METy + b1pz*b2pz*l2px*ml12 - 
   b1pz*b2px*l2pz*ml12 + b1pz*b2pz*l1px*ml22 - b1px*b2pz*l1pz*ml22 + b2pz*l1pz*l2px*mt12 - b2px*l1pz*l2pz*mt12 + b1pz*l1px*l2pz*mt22 - 
   b1px*l1pz*l2pz*mt22 - b1pz*b2pz*l2px*mW12 - b2pz*l1pz*l2px*mW12 + b1pz*b2px*l2pz*mW12 + b2px*l1pz*l2pz*mW12 - b1pz*b2pz*l1px*mW22 + 
   b1px*b2pz*l1pz*mW22 - b1pz*l1px*l2pz*mW22 + b1px*l1pz*l2pz*mW22;
  double nu1Py1 = 2*b1pz*b2pz*l1E*l2px - 2*b1E*b2pz*l1pz*l2px - 2*b1pz*b2px*l1E*l2pz + 2*b1E*b2px*l1pz*l2pz;
  double nu1Py2 = 2*b1pz*b2pz*l1px*l2E - 2*b1px*b2pz*l1pz*l2E - 2*b1pz*b2E*l1px*l2pz + 2*b1px*b2E*l1pz*l2pz;
  
  double nu1Pz0 = 2*b1E*b2pz*l1E*l1py*l2px - 2*b1px*b2pz*l1px*l1py*l2px - 2*b1py*b2pz*pow(l1py,2)*l2px - 2*b1pz*b2pz*l1py*l1pz*l2px - 2*b1E*b2pz*l1E*l1px*l2py + 
   2*b1px*b2pz*pow(l1px,2)*l2py + 2*b1py*b2pz*l1px*l1py*l2py + 2*b1pz*b2pz*l1px*l1pz*l2py + 2*b1E*b2py*l1E*l1px*l2pz - 
   2*b1px*b2py*pow(l1px,2)*l2pz - 2*b1E*b2px*l1E*l1py*l2pz + 2*b1px*b2px*l1px*l1py*l2pz - 2*b1py*b2py*l1px*l1py*l2pz + 
   2*b1py*b2px*pow(l1py,2)*l2pz - 2*b1pz*b2py*l1px*l1pz*l2pz + 2*b1pz*b2px*l1py*l1pz*l2pz + 2*b1py*b2E*l1px*l2E*l2pz - 
   2*b1px*b2E*l1py*l2E*l2pz - 2*b1py*b2px*l1px*l2px*l2pz + 2*b1px*b2px*l1py*l2px*l2pz - 2*b1py*b2py*l1px*l2py*l2pz + 2*b1px*b2py*l1py*l2py*l2pz - 
   2*b1py*b2pz*l1px*pow(l2pz,2) + 2*b1px*b2pz*l1py*pow(l2pz,2) + b2pz*l1py*l2px*mb12 - b2pz*l1px*l2py*mb12 + b2py*l1px*l2pz*mb12 - 
   b2px*l1py*l2pz*mb12 + b1py*l1px*l2pz*mb22 - b1px*l1py*l2pz*mb22 + 2*b1py*b2pz*l1px*l2px*METx - 2*b1px*b2pz*l1py*l2px*METx - 
   2*b1py*b2px*l1px*l2pz*METx + 2*b1px*b2px*l1py*l2pz*METx + 2*b1py*b2pz*l1px*l2py*METy - 2*b1px*b2pz*l1py*l2py*METy - 
   2*b1py*b2py*l1px*l2pz*METy + 2*b1px*b2py*l1py*l2pz*METy - b1py*b2pz*l2px*ml12 + b1px*b2pz*l2py*ml12 + b1py*b2px*l2pz*ml12 - 
   b1px*b2py*l2pz*ml12 - b1py*b2pz*l1px*ml22 + b1px*b2pz*l1py*ml22 - b2pz*l1py*l2px*mt12 + b2pz*l1px*l2py*mt12 - b2py*l1px*l2pz*mt12 + 
   b2px*l1py*l2pz*mt12 - b1py*l1px*l2pz*mt22 + b1px*l1py*l2pz*mt22 + b1py*b2pz*l2px*mW12 + b2pz*l1py*l2px*mW12 - b1px*b2pz*l2py*mW12 - 
   b2pz*l1px*l2py*mW12 - b1py*b2px*l2pz*mW12 + b1px*b2py*l2pz*mW12 + b2py*l1px*l2pz*mW12 - b2px*l1py*l2pz*mW12 + b1py*b2pz*l1px*mW22 - 
   b1px*b2pz*l1py*mW22 + b1py*l1px*l2pz*mW22 - b1px*l1py*l2pz*mW22;
  double nu1Pz1 = -2*b1py*b2pz*l1E*l2px + 2*b1E*b2pz*l1py*l2px + 2*b1px*b2pz*l1E*l2py - 2*b1E*b2pz*l1px*l2py + 2*b1py*b2px*l1E*l2pz - 2*b1px*b2py*l1E*l2pz + 2*b1E*b2py*l1px*l2pz - 2*b1E*b2px*l1py*l2pz;
  double nu1Pz2 = -2*b1py*b2pz*l1px*l2E + 2*b1px*b2pz*l1py*l2E + 2*b1py*b2E*l1px*l2pz - 2*b1px*b2E*l1py*l2pz;
  
  double nu2Px0 = -2*b1E*b2pz*l1E*l1pz*l2py + 2*b1px*b2pz*l1px*l1pz*l2py + 2*b1py*b2pz*l1py*l1pz*l2py + 2*b1pz*b2pz*pow(l1pz,2)*l2py + 2*b1E*b2py*l1E*l1pz*l2pz - 
   2*b1px*b2py*l1px*l1pz*l2pz - 2*b1py*b2py*l1py*l1pz*l2pz - 2*b1pz*b2py*pow(l1pz,2)*l2pz - 2*b1pz*b2E*l1py*l2E*l2pz + 
   2*b1py*b2E*l1pz*l2E*l2pz + 2*b1pz*b2px*l1py*l2px*l2pz - 2*b1py*b2px*l1pz*l2px*l2pz + 2*b1pz*b2py*l1py*l2py*l2pz - 2*b1py*b2py*l1pz*l2py*l2pz + 
   2*b1pz*b2pz*l1py*pow(l2pz,2) - 2*b1py*b2pz*l1pz*pow(l2pz,2) - b2pz*l1pz*l2py*mb12 + b2py*l1pz*l2pz*mb12 - b1pz*l1py*l2pz*mb22 + 
   b1py*l1pz*l2pz*mb22 - 2*b1pz*b2pz*l1px*l2py*METx + 2*b1px*b2pz*l1pz*l2py*METx + 2*b1pz*b2py*l1px*l2pz*METx - 2*b1px*b2py*l1pz*l2pz*METx - 
   2*b1pz*b2pz*l1py*l2py*METy + 2*b1py*b2pz*l1pz*l2py*METy + 2*b1pz*b2py*l1py*l2pz*METy - 2*b1py*b2py*l1pz*l2pz*METy + b1pz*b2pz*l2py*ml12 - 
   b1pz*b2py*l2pz*ml12 + b1pz*b2pz*l1py*ml22 - b1py*b2pz*l1pz*ml22 + b2pz*l1pz*l2py*mt12 - b2py*l1pz*l2pz*mt12 + b1pz*l1py*l2pz*mt22 - 
   b1py*l1pz*l2pz*mt22 - b1pz*b2pz*l2py*mW12 - b2pz*l1pz*l2py*mW12 + b1pz*b2py*l2pz*mW12 + b2py*l1pz*l2pz*mW12 - b1pz*b2pz*l1py*mW22 + 
   b1py*b2pz*l1pz*mW22 - b1pz*l1py*l2pz*mW22 + b1py*l1pz*l2pz*mW22;
  double nu2Px1 = 2*b1pz*b2pz*l1E*l2py - 2*b1E*b2pz*l1pz*l2py - 2*b1pz*b2py*l1E*l2pz + 2*b1E*b2py*l1pz*l2pz;
  double nu2Px2 = 2*b1pz*b2pz*l1py*l2E - 2*b1py*b2pz*l1pz*l2E - 2*b1pz*b2E*l1py*l2pz + 2*b1py*b2E*l1pz*l2pz;
  
  double nu2Py0 = 2*b1E*b2pz*l1E*l1pz*l2px - 2*b1px*b2pz*l1px*l1pz*l2px - 2*b1py*b2pz*l1py*l1pz*l2px - 2*b1pz*b2pz*pow(l1pz,2)*l2px - 2*b1E*b2px*l1E*l1pz*l2pz + 
   2*b1px*b2px*l1px*l1pz*l2pz + 2*b1py*b2px*l1py*l1pz*l2pz + 2*b1pz*b2px*pow(l1pz,2)*l2pz + 2*b1pz*b2E*l1px*l2E*l2pz - 
   2*b1px*b2E*l1pz*l2E*l2pz - 2*b1pz*b2px*l1px*l2px*l2pz + 2*b1px*b2px*l1pz*l2px*l2pz - 2*b1pz*b2py*l1px*l2py*l2pz + 2*b1px*b2py*l1pz*l2py*l2pz - 
   2*b1pz*b2pz*l1px*pow(l2pz,2) + 2*b1px*b2pz*l1pz*pow(l2pz,2) + b2pz*l1pz*l2px*mb12 - b2px*l1pz*l2pz*mb12 + b1pz*l1px*l2pz*mb22 - 
   b1px*l1pz*l2pz*mb22 + 2*b1pz*b2pz*l1px*l2px*METx - 2*b1px*b2pz*l1pz*l2px*METx - 2*b1pz*b2px*l1px*l2pz*METx + 2*b1px*b2px*l1pz*l2pz*METx + 
   2*b1pz*b2pz*l1py*l2px*METy - 2*b1py*b2pz*l1pz*l2px*METy - 2*b1pz*b2px*l1py*l2pz*METy + 2*b1py*b2px*l1pz*l2pz*METy - b1pz*b2pz*l2px*ml12 + 
   b1pz*b2px*l2pz*ml12 - b1pz*b2pz*l1px*ml22 + b1px*b2pz*l1pz*ml22 - b2pz*l1pz*l2px*mt12 + b2px*l1pz*l2pz*mt12 - b1pz*l1px*l2pz*mt22 + 
   b1px*l1pz*l2pz*mt22 + b1pz*b2pz*l2px*mW12 + b2pz*l1pz*l2px*mW12 - b1pz*b2px*l2pz*mW12 - b2px*l1pz*l2pz*mW12 + b1pz*b2pz*l1px*mW22 - 
   b1px*b2pz*l1pz*mW22 + b1pz*l1px*l2pz*mW22 - b1px*l1pz*l2pz*mW22;
  double nu2Py1 = -2*b1pz*b2pz*l1E*l2px + 2*b1E*b2pz*l1pz*l2px + 2*b1pz*b2px*l1E*l2pz - 2*b1E*b2px*l1pz*l2pz;
  double nu2Py2 = -2*b1pz*b2pz*l1px*l2E + 2*b1px*b2pz*l1pz*l2E + 2*b1pz*b2E*l1px*l2pz - 2*b1px*b2E*l1pz*l2pz;
  
  double nu2Pz0 = -2*b1E*b2py*l1E*l1pz*l2px + 2*b1px*b2py*l1px*l1pz*l2px + 2*b1py*b2py*l1py*l1pz*l2px + 2*b1pz*b2py*pow(l1pz,2)*l2px + 2*b1pz*b2E*l1py*l2E*l2px - 
   2*b1py*b2E*l1pz*l2E*l2px - 2*b1pz*b2px*l1py*pow(l2px,2) + 2*b1py*b2px*l1pz*pow(l2px,2) + 2*b1E*b2px*l1E*l1pz*l2py - 
   2*b1px*b2px*l1px*l1pz*l2py - 2*b1py*b2px*l1py*l1pz*l2py - 2*b1pz*b2px*pow(l1pz,2)*l2py - 2*b1pz*b2E*l1px*l2E*l2py + 
   2*b1px*b2E*l1pz*l2E*l2py + 2*b1pz*b2px*l1px*l2px*l2py - 2*b1pz*b2py*l1py*l2px*l2py - 2*b1px*b2px*l1pz*l2px*l2py + 2*b1py*b2py*l1pz*l2px*l2py + 
   2*b1pz*b2py*l1px*pow(l2py,2) - 2*b1px*b2py*l1pz*pow(l2py,2) - 2*b1pz*b2pz*l1py*l2px*l2pz + 2*b1py*b2pz*l1pz*l2px*l2pz + 
   2*b1pz*b2pz*l1px*l2py*l2pz - 2*b1px*b2pz*l1pz*l2py*l2pz - b2py*l1pz*l2px*mb12 + b2px*l1pz*l2py*mb12 + b1pz*l1py*l2px*mb22 - 
   b1py*l1pz*l2px*mb22 - b1pz*l1px*l2py*mb22 + b1px*l1pz*l2py*mb22 - 2*b1pz*b2py*l1px*l2px*METx + 2*b1px*b2py*l1pz*l2px*METx + 
   2*b1pz*b2px*l1px*l2py*METx - 2*b1px*b2px*l1pz*l2py*METx - 2*b1pz*b2py*l1py*l2px*METy + 2*b1py*b2py*l1pz*l2px*METy + 
   2*b1pz*b2px*l1py*l2py*METy - 2*b1py*b2px*l1pz*l2py*METy + b1pz*b2py*l2px*ml12 - b1pz*b2px*l2py*ml12 + b1pz*b2py*l1px*ml22 - 
   b1pz*b2px*l1py*ml22 + b1py*b2px*l1pz*ml22 - b1px*b2py*l1pz*ml22 + b2py*l1pz*l2px*mt12 - b2px*l1pz*l2py*mt12 - b1pz*l1py*l2px*mt22 + 
   b1py*l1pz*l2px*mt22 + b1pz*l1px*l2py*mt22 - b1px*l1pz*l2py*mt22 - b1pz*b2py*l2px*mW12 - b2py*l1pz*l2px*mW12 + b1pz*b2px*l2py*mW12 + 
   b2px*l1pz*l2py*mW12 - b1pz*b2py*l1px*mW22 + b1pz*b2px*l1py*mW22 - b1py*b2px*l1pz*mW22 + b1px*b2py*l1pz*mW22 + b1pz*l1py*l2px*mW22 - 
   b1py*l1pz*l2px*mW22 - b1pz*l1px*l2py*mW22 + b1px*l1pz*l2py*mW22;
  double nu2Pz1 = 2*b1pz*b2py*l1E*l2px - 2*b1E*b2py*l1pz*l2px - 2*b1pz*b2px*l1E*l2py + 2*b1E*b2px*l1pz*l2py;
  double nu2Pz2 = 2*b1pz*b2py*l1px*l2E - 2*b1pz*b2px*l1py*l2E + 2*b1py*b2px*l1pz*l2E - 2*b1px*b2py*l1pz*l2E + 2*b1pz*b2E*l1py*l2px - 
    2*b1py*b2E*l1pz*l2px - 2*b1pz*b2E*l1px*l2py + 2*b1px*b2E*l1pz*l2py;
  
  //cout << "nu1Px = " << nu1Px0/nuDenom << " + " << nu1Px1/nuDenom << " * Nu1E + " << nu1Px2/nuDenom << " * Nu2E" << endl; 
  //cout << "nu1Py = " << nu1Py0/nuDenom << " + " << nu1Py1/nuDenom << " * Nu1E + " << nu1Py2/nuDenom << " * Nu2E" << endl; 
  //cout << "nu1Pz = " << nu1Pz0/nuDenom << " + " << nu1Pz1/nuDenom << " * Nu1E + " << nu1Pz2/nuDenom << " * Nu2E" << endl; 
  //cout << "nu2Px = " << nu2Px0/nuDenom << " + " << nu2Px1/nuDenom << " * Nu1E + " << nu2Px2/nuDenom << " * Nu2E" << endl; 
  //cout << "nu2Py = " << nu2Py0/nuDenom << " + " << nu2Py1/nuDenom << " * Nu1E + " << nu2Py2/nuDenom << " * Nu2E" << endl; 
  //cout << "nu2Pz = " << nu2Pz0/nuDenom << " + " << nu2Pz1/nuDenom << " * Nu1E + " << nu2Pz2/nuDenom << " * Nu2E" << endl; 
  
  //Now we want to solve the equations (p_nu1)^2 = 0 and (p_nu2)^2 = 0
  // use the form p = p00 + p10 * nu1E + p01 * nu2E + p11 * nu1E*nu2E + p20 * nu1E^2 + p02 * nu2E^2 
  //We could also solve for a neutrino with mass.  This would mean that p00->p00 - nuDenom*nuDenom * mnu_1^2
  //We could also solve for a neutrino with mass.  This would mean that q00->q00 - nuDenom*nuDenom * mnu_2^2
  
  double p00 = -pow(nu1Px0,2) - pow(nu1Py0,2) - pow(nu1Pz0,2);
  double p10 = -2*nu1Px0*nu1Px1 - 2*nu1Py0*nu1Py1 - 2*nu1Pz0*nu1Pz1;
  double p01 = -2*nu1Px0*nu1Px2 - 2*nu1Py0*nu1Py2 - 2*nu1Pz0*nu1Pz2;
  double p11 = -2*nu1Px1*nu1Px2 - 2*nu1Py1*nu1Py2 - 2*nu1Pz1*nu1Pz2;
  double p20 = -pow(nu1Px1,2) - pow(nu1Py1,2) - pow(nu1Pz1,2) + pow(nuDenom,2);
  double p02 = -pow(nu1Px2,2) - pow(nu1Py2,2) - pow(nu1Pz2,2);
    
  //cout << p00/(nuDenom*nuDenom) << " + " << p10/(nuDenom*nuDenom) << " * nu1E + " << p01/(nuDenom*nuDenom) << " * nu2E + " << p11/(nuDenom*nuDenom) << " * nu1E*nu2E + " << p20/(nuDenom*nuDenom) << " * nu1E^2 + " << p02/(nuDenom*nuDenom) << " * nu2E^2 " << endl; 
  
  double q00 = -pow(nu2Px0,2) - pow(nu2Py0,2) - pow(nu2Pz0,2);
  double q10 = -2*nu2Px0*nu2Px1 - 2*nu2Py0*nu2Py1 - 2*nu2Pz0*nu2Pz1;
  double q01 = -2*nu2Px0*nu2Px2 - 2*nu2Py0*nu2Py2 - 2*nu2Pz0*nu2Pz2;
  double q11 = -2*nu2Px1*nu2Px2 - 2*nu2Py1*nu2Py2 - 2*nu2Pz1*nu2Pz2;
  double q20 = -pow(nu2Px1,2) - pow(nu2Py1,2) - pow(nu2Pz1,2);
  double q02 = -pow(nu2Px2,2) - pow(nu2Py2,2) - pow(nu2Pz2,2) + pow(nuDenom,2);
  
  //cout << q00/(nuDenom*nuDenom) << " + " << q10/(nuDenom*nuDenom) << " * nu1E + " << q01/(nuDenom*nuDenom) << " * nu2E + " << q11/(nuDenom*nuDenom) << " * nu1E*nu2E + " << q20/(nuDenom*nuDenom) << " * nu1E^2 + " << q02/(nuDenom*nuDenom) << " * nu2E^2 " << endl; 
  
  //Reduce these to one quartic equation
  // p * q02 - q * p02 gives linear in nu2E which then yields a quartic in nu1E
  
  double nu2EDenom0 = p01*q02 - p02*q01;
  double nu2EDenom1 = p11*q02 - p02*q11;
  
  double nu2ENum0 = p02*q00 - p00*q02;
  double nu2ENum1 = p02*q10 - p10*q02;
  double nu2ENum2 = p02*q20 - p20*q02;

  //cout << "p * q02 - q * p02" << endl;
  //cout << -p02*q00+p00*q02 << " + " << p10*q02-p02*q10 << " * nu1E + " << p20*q02-p02*q20 << " * nu1E^2 + " << -p02*q01+p01*q02 << " * nu2E + " << p11*q02-p02*q11 << " * nu1E*nu2E" << endl;
  //
  //cout << "numerator:" << endl;
  //cout << nu2ENum0 << " + " << nu2ENum1 << " * nu1E + " << nu2ENum2 << " * nu1E^2" << endl;
  //cout << "denominator:" << endl;
  //cout << nu2EDenom0 << " + " << nu2EDenom1 << " * nu1E" << endl;
  
  double r0 = pow(p02,3)*pow(q00,2) - p01*pow(p02,2)*q00*q01 + p00*pow(p02,2)*pow(q01,2) + pow(p01,2)*p02*q00*q02 - 2*p00*pow(p02,2)*q00*q02 - 
   p00*p01*p02*q01*q02 + pow(p00,2)*p02*pow(q02,2);
  double r1 = -(pow(p02,2)*p11*q00*q01) + pow(p02,2)*p10*pow(q01,2) - 2*pow(p02,2)*p10*q00*q02 + 2*p01*p02*p11*q00*q02 - p01*p02*p10*q01*q02 - 
   p00*p02*p11*q01*q02 + 2*p00*p02*p10*pow(q02,2) + 2*pow(p02,3)*q00*q10 - p01*pow(p02,2)*q01*q10 + pow(p01,2)*p02*q02*q10 - 
   2*p00*pow(p02,2)*q02*q10 - p01*pow(p02,2)*q00*q11 + 2*p00*pow(p02,2)*q01*q11 - p00*p01*p02*q02*q11;
  double r2 = pow(p02,2)*p20*pow(q01,2) + p02*pow(p11,2)*q00*q02 - 2*pow(p02,2)*p20*q00*q02 - p02*p10*p11*q01*q02 - p01*p02*p20*q01*q02 + 
   p02*pow(p10,2)*pow(q02,2) + 2*p00*p02*p20*pow(q02,2) - pow(p02,2)*p11*q01*q10 - 2*pow(p02,2)*p10*q02*q10 + 2*p01*p02*p11*q02*q10 + 
   pow(p02,3)*pow(q10,2) - pow(p02,2)*p11*q00*q11 + 2*pow(p02,2)*p10*q01*q11 - p01*p02*p10*q02*q11 - p00*p02*p11*q02*q11 - 
   p01*pow(p02,2)*q10*q11 + p00*pow(p02,2)*pow(q11,2) + 2*pow(p02,3)*q00*q20 - p01*pow(p02,2)*q01*q20 + pow(p01,2)*p02*q02*q20 - 
   2*p00*pow(p02,2)*q02*q20;
  double r3 = -(p02*p11*p20*q01*q02) + 2*p02*p10*p20*pow(q02,2) + p02*pow(p11,2)*q02*q10 - 2*pow(p02,2)*p20*q02*q10 + 2*pow(p02,2)*p20*q01*q11 - 
   p02*p10*p11*q02*q11 - p01*p02*p20*q02*q11 - pow(p02,2)*p11*q10*q11 + pow(p02,2)*p10*pow(q11,2) - pow(p02,2)*p11*q01*q20 - 
   2*pow(p02,2)*p10*q02*q20 + 2*p01*p02*p11*q02*q20 + 2*pow(p02,3)*q10*q20 - p01*pow(p02,2)*q11*q20;
  double r4 = p02*pow(p20,2)*pow(q02,2) - p02*p11*p20*q02*q11 + pow(p02,2)*p20*pow(q11,2) + p02*pow(p11,2)*q02*q20 - 2*pow(p02,2)*p20*q02*q20 - 
   pow(p02,2)*p11*q11*q20 + pow(p02,3)*pow(q20,2);
  
  //cout << "the quartic inputs: \n" 
  //     << "r4 = " << r4 << "\n" 
  //     << "r3 = " << r3 << "\n" 
  //     << "r2 = " << r2 << "\n" 
  //     << "r1 = " << r1 << "\n" 
  //     << "r0 = " << r0 << "\n" ;
  
  Polynomial* quartic = NULL;
  vector<complex<double> > nu1EValues;
  if(r4!=0)         quartic = new Polynomial(r4,r3,r2,r1,r0);
  else if (r3 != 0) quartic = new Polynomial(r3,r2,r1,r0);
  else if (r2 != 0) quartic = new Polynomial(r2,r1,r0);
  else if (r1 != 0) quartic = new Polynomial(r1,r0);

  //if(!quartic) cout << "degenerate neutrino solutions" << endl;
  if(quartic) 
    {
      nu1EValues = quartic->FindRoots();
      delete quartic;
      quartic = NULL;
    }
  for(vector<complex<double> >::iterator iNu1E = nu1EValues.begin(); iNu1E!=nu1EValues.end(); iNu1E++)
    {
      complex<double> thisNu1E  = *iNu1E;
      complex<double> thisNu2E  = (nu2ENum0 + nu2ENum1 * thisNu1E + nu2ENum2 * thisNu1E*thisNu1E) / (nu2EDenom0 + nu2EDenom1 * thisNu1E);
      complex<double> thisNu1Px = (nu1Px0 + nu1Px1 * thisNu1E + nu1Px2 * thisNu2E)/nuDenom;
      complex<double> thisNu1Py = (nu1Py0 + nu1Py1 * thisNu1E + nu1Py2 * thisNu2E)/nuDenom;
      complex<double> thisNu1Pz = (nu1Pz0 + nu1Pz1 * thisNu1E + nu1Pz2 * thisNu2E)/nuDenom;
      complex<double> thisNu2Px = (nu2Px0 + nu2Px1 * thisNu1E + nu2Px2 * thisNu2E)/nuDenom;
      complex<double> thisNu2Py = (nu2Py0 + nu2Py1 * thisNu1E + nu2Py2 * thisNu2E)/nuDenom;
      complex<double> thisNu2Pz = (nu2Pz0 + nu2Pz1 * thisNu1E + nu2Pz2 * thisNu2E)/nuDenom;

      //cout << "nu1E  = " << thisNu1E << endl; 
      //cout << "nu1Px = " << thisNu1Px << endl; 
      //cout << "nu1Py = " << thisNu1Py << endl; 
      //cout << "nu1Pz = " << thisNu1Pz << endl; 
      //cout << "nu2E  = " << thisNu2E << endl; 
      //cout << "nu2Px = " << thisNu2Px << endl; 
      //cout << "nu2Py = " << thisNu2Py << endl; 
      //cout << "nu2Pz = " << thisNu2Pz << endl; 
      
      nu1E.push_back(  thisNu1E  );
      nu2E.push_back(  thisNu2E  );
      nu1px.push_back( thisNu1Px );
      nu1py.push_back( thisNu1Py );
      nu1pz.push_back( thisNu1Pz );
      nu2px.push_back( thisNu2Px );
      nu2py.push_back( thisNu2Py );
      nu2pz.push_back( thisNu2Pz );	     
    }

} 

//vector<pair<double, double> > getIntersectionPoints(TGraph& one, TGraph& two)
//{
//  // Return a TGraph with the points of intersection
//  vector<pair<double, double> > intersectionPoints;
//  //TGraph *interPoint = new TGraph();
//  Int_t i = 0;
//  
//  // Loop over all points in this TGraph
//  for(size_t a_i = 0; a_i < one.GetN()-1; ++a_i)
//    {
//      // Loop over all points in the other TGraph
//      for(size_t b_i = 0; b_i < two.GetN()-1; ++b_i)
//      {
//         
//	// Get the current point, and the next point for each of the objects
//	Double_t x1, y1, x2, y2 = 0;
//	Double_t ax1, ay1, ax2, ay2 = 0;
//	one.GetPoint(a_i, x1, y1);
//	one.GetPoint(a_i+1, x2, y2);
//	two.GetPoint(b_i, ax1, ay1);
//	two.GetPoint(b_i+1, ax2, ay2);
//
//	// Calculate the intersection between two straight lines, x axis
//	Double_t x = (ax1 *(ay2 *(x1-x2)+x2 * y1 - x1 * y2 )+ ax2 * (ay1 * (-x1+x2)- x2 * y1+x1 * y2)) 
//	  / (-(ay1-ay2) * (x1-x2)+(ax1-ax2)* (y1-y2));
//	
//	// Calculate the intersection between two straight lines, y axis
//	Double_t y = (ax1 * ay2 * (y1-y2)+ax2 * ay1 * (-y1+y2)+(ay1-ay2) * (x2 * y1-x1 * y2))/(-(ay1-ay2) * (x1-x2)+(ax1-ax2) * (y1-y2));
//	
//	// Find the tightest interval along the x-axis defined by the four points
//	Double_t xrange_min = max(min(x1, x2), min(ax1, ax2));
//	Double_t xrange_max = min(max(x1, x2), max(ax1, ax2));
//	
//	// If points from the two lines overlap, they are trivially intersecting
//	if ((x1 == ax1 && y1 == ay1) or (x2 == ax2 && y2 == ay2)){               
//	  intersectionPoints.push_back( make_pair( (x1 == ax1 && y1 == ay1) ? x1 : x2, (x1 == ax1 && y1 == ay1) ? y1 : y2) ) ; //SetPoint(i, (x1 == ax1 and y1 == ay1) ? x1 : x2, (x1 == ax1 and y1 == ay1) ? y1 : y2);
//	  i++;            
//	} 
//	// If the intersection between the two lines is within the tight range, add it to the list of intersections.
//	else if(x > xrange_min && x < xrange_max)
//	  {
//            intersectionPoints.push_back(make_pair(x, y)) ; //interPoint->SetPoint(i,x, y);
//            i++;
//	  }
//      }
//    }
//  return intersectionPoints;
//}

