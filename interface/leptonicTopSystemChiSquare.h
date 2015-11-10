#ifndef LEPTONICTOPSYSTEMCHISQUARE
#define LEPTONICTOPSYSTEMCHISQUARE 

#include "topSystemChiSquare.h"

class leptonicTopSystemChiSquare : public topSystemChiSquare
{
 private:
  double chi2_;

 public:
  leptonicTopSystemChiSquare(const PhysicsObject& bjet, const PhysicsObject& lepton,
                             const double& , const double& ,
                             const double& , const double& );

  ~leptonicTopSystemChiSquare();

  virtual void printTopConstituents();

  virtual void setEllipseAngle(double);

  virtual void setWDaughter2(double, double, double) { return;};
  virtual void calcWDaughter2Deltas() { return;} ;
  virtual void getWDaughter2Deltas(double& , double& , double& );
  virtual void printWDaughter2() {return;};

  void calcChiSquare();
  double getChiSquare();

  virtual double getHadronicChiSquare() {return 0;};

};

#endif
