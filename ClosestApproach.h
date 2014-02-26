// Calculate the distance of closest approach between
// two ellipses representing the neutrino momenta in
// dileptonic ttbar decays.
// The inputs are the Hperp matrices defined in Eq.7
// of arXiv:1305.1872.


#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

#include <TMatrix.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TArrayD.h>
#include <TMath.h>
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/GSLMinimizer1D.h"

#include "Math/IFunction.h"

//#include "Ellipse.cxx"

using namespace std;
using namespace ROOT::Math;
 
class distanceDerivative: public ROOT::Math::IGradientFunctionOneDim
{

 private:

  double b111_, b112_, b122_;
  double b211_, b212_, b222_;
  double xc_, yc_; 

 public:

  void setB1(double b11, double b12, double b22)
  {
    b111_ = b11;
    b112_ = b12;
    b122_ = b22;
  }

  void setB2(double b11, double b12, double b22)
  {
    b211_ = b11;
    b212_ = b12;
    b222_ = b22;
  }

  void setCenter(double x,double y)
  {
    xc_ = x;
    yc_ = y;
  }

   distanceDerivative():
     b111_(0),
     b112_(0),
     b122_(0),
     b211_(0),
     b212_(0),
     b222_(0),
     xc_(0),
     yc_(0)
       {
       }

   distanceDerivative(const distanceDerivative& other):
     b111_(other.b111_),
     b112_(other.b112_),
     b122_(other.b122_),
     b211_(other.b211_),
     b212_(other.b212_),
     b222_(other.b222_),
     xc_(other.xc_),
     yc_(other.yc_)
       {
       }

   double DoEval(double theta) const
   {
     double cosTheta=cos(theta);
     double sinTheta=sin(theta);
     double cosTwoTheta=cos(2.*theta);
     double sinTwoTheta=sin(2.*theta);

     double b111=b111_, b112=b112_, b122=b122_;
     double b211=b211_, b212=b212_, b222=b222_;

     double d1p=

       (((b111*xc_ + b112*yc_)*cosTheta + (b112*xc_ + b122*yc_)*sinTheta)*(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5) + 
       (2*b112*(b111 + b122)*cosTwoTheta - 2*(pow(b111,2) + pow(b112,2))*cosTheta*sinTheta + 2*(pow(b112,2) + pow(b122,2))*cosTheta*sinTheta)/
       (b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) - 
       (2*((b112*xc_ + b122*yc_)*cosTheta - (b111*xc_ + b112*yc_)*sinTheta))/sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) - 
       ((2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta)*
	((pow(b111,2) + pow(b112,2))*pow(cosTheta,2) + (pow(b112,2) + pow(b122,2))*pow(sinTheta,2) + b112*(b111 + b122)*sinTwoTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,2) + 
       (((b211*xc_ + b212*yc_)*cosTheta + (b212*xc_ + b222*yc_)*sinTheta)*(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5) + 
       (2*b212*(b211 + b222)*cosTwoTheta - 2*(pow(b211,2) + pow(b212,2))*cosTheta*sinTheta + 2*(pow(b212,2) + pow(b222,2))*cosTheta*sinTheta)/
       (b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) - 
       (2*((b212*xc_ + b222*yc_)*cosTheta - (b211*xc_ + b212*yc_)*sinTheta))/sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) - 
       ((2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta)*
	((pow(b211,2) + pow(b212,2))*pow(cosTheta,2) + (pow(b212,2) + pow(b222,2))*pow(sinTheta,2) + b212*(b211 + b222)*sinTwoTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,2) + 
       (2*((b111 + b122)*b212 + b112*(b211 + b222))*cosTwoTheta - 2*(b111*b211 - b122*b222)*sinTwoTheta)/
       (sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta))
       - ((2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta)*
	  (b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (2.*sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*
	pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5)) - 
       ((2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta)*
	(b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (2.*pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5)*
	sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta))   ;

     return d1p;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new distanceDerivative(*this);
   }
 
   double DoDerivative(double theta) const
   {
     double cosTheta=cos(theta);
     double sinTheta=sin(theta);
     double cosTwoTheta=cos(2.*theta);
     double sinTwoTheta=sin(2.*theta);

     double b111=b111_, b112=b112_, b122=b122_;
     double b211=b211_, b212=b212_, b222=b222_;

     double d2p=

       (-3*((b111*xc_ + b112*yc_)*cosTheta + (b112*xc_ + b122*yc_)*sinTheta)*pow(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta,2))/
       (2.*pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,2.5)) - 
       (2*(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta)*
	(2*b112*(b111 + b122)*cosTwoTheta - 2*(pow(b111,2) + pow(b112,2))*cosTheta*sinTheta + 2*(pow(b112,2) + pow(b122,2))*cosTheta*sinTheta)
	)/pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,2) + 
       (2*((b112*xc_ + b122*yc_)*cosTheta - (b111*xc_ + b112*yc_)*sinTheta)*(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5) + 
       (((b111*xc_ + b112*yc_)*cosTheta + (b112*xc_ + b122*yc_)*sinTheta)*
	(-2*b111*pow(cosTheta,2) + 2*b122*pow(cosTheta,2) + 2*b111*pow(sinTheta,2) - 2*b122*pow(sinTheta,2) - 4*b112*sinTwoTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5) - 
       (2*(-((b111*xc_ + b112*yc_)*cosTheta) - (b112*xc_ + b122*yc_)*sinTheta))/sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) + 
       (-2*(pow(b111,2) + pow(b112,2))*pow(cosTheta,2) + 2*(pow(b112,2) + pow(b122,2))*pow(cosTheta,2) + 
	2*(pow(b111,2) + pow(b112,2))*pow(sinTheta,2) - 2*(pow(b112,2) + pow(b122,2))*pow(sinTheta,2) - 
	4*b112*(b111 + b122)*sinTwoTheta)/(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) + 
       (2*pow(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta,2)*
	((pow(b111,2) + pow(b112,2))*pow(cosTheta,2) + (pow(b112,2) + pow(b122,2))*pow(sinTheta,2) + b112*(b111 + b122)*sinTwoTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,3) - 
       ((-2*b111*pow(cosTheta,2) + 2*b122*pow(cosTheta,2) + 2*b111*pow(sinTheta,2) - 2*b122*pow(sinTheta,2) - 4*b112*sinTwoTheta)*
	((pow(b111,2) + pow(b112,2))*pow(cosTheta,2) + (pow(b112,2) + pow(b122,2))*pow(sinTheta,2) + b112*(b111 + b122)*sinTwoTheta))/
       pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,2) - 
       (3*((b211*xc_ + b212*yc_)*cosTheta + (b212*xc_ + b222*yc_)*sinTheta)*pow(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta,2))/
       (2.*pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,2.5)) - 
       (2*(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta)*
	(2*b212*(b211 + b222)*cosTwoTheta - 2*(pow(b211,2) + pow(b212,2))*cosTheta*sinTheta + 2*(pow(b212,2) + pow(b222,2))*cosTheta*sinTheta)
	)/pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,2) + 
       (2*((b212*xc_ + b222*yc_)*cosTheta - (b211*xc_ + b212*yc_)*sinTheta)*(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5) + 
       (((b211*xc_ + b212*yc_)*cosTheta + (b212*xc_ + b222*yc_)*sinTheta)*
	(-2*b211*pow(cosTheta,2) + 2*b222*pow(cosTheta,2) + 2*b211*pow(sinTheta,2) - 2*b222*pow(sinTheta,2) - 4*b212*sinTwoTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5) - 
       (2*(-((b211*xc_ + b212*yc_)*cosTheta) - (b212*xc_ + b222*yc_)*sinTheta))/sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) + 
       (-2*(pow(b211,2) + pow(b212,2))*pow(cosTheta,2) + 2*(pow(b212,2) + pow(b222,2))*pow(cosTheta,2) + 
	2*(pow(b211,2) + pow(b212,2))*pow(sinTheta,2) - 2*(pow(b212,2) + pow(b222,2))*pow(sinTheta,2) - 
	4*b212*(b211 + b222)*sinTwoTheta)/(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) + 
       (2*pow(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta,2)*
	((pow(b211,2) + pow(b212,2))*pow(cosTheta,2) + (pow(b212,2) + pow(b222,2))*pow(sinTheta,2) + b212*(b211 + b222)*sinTwoTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,3) - 
       ((-2*b211*pow(cosTheta,2) + 2*b222*pow(cosTheta,2) + 2*b211*pow(sinTheta,2) - 2*b222*pow(sinTheta,2) - 4*b212*sinTwoTheta)*
	((pow(b211,2) + pow(b212,2))*pow(cosTheta,2) + (pow(b212,2) + pow(b222,2))*pow(sinTheta,2) + b212*(b211 + b222)*sinTwoTheta))/
       pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,2) - 
       ((2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta)*
	(2*((b111 + b122)*b212 + b112*(b211 + b222))*cosTwoTheta - 2*(b111*b211 - b122*b222)*sinTwoTheta))/
       (sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*
	pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5)) - 
       ((2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta)*
	(2*((b111 + b122)*b212 + b112*(b211 + b222))*cosTwoTheta - 2*(b111*b211 - b122*b222)*sinTwoTheta))/
       (pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5)*
	sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta)) + 
       (-4*(b111*b211 - b122*b222)*cosTwoTheta - 4*((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta)/
       (sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta))
       + (3*pow(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta,2)*
	  (b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (4.*sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*
	pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,2.5)) + 
       ((2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta)*(2*b212*cosTwoTheta - 2*b211*cosTheta*sinTheta + 2*b222*cosTheta*sinTheta)*
	(b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (2.*pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5)*
	pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5)) - 
       ((-2*b211*pow(cosTheta,2) + 2*b222*pow(cosTheta,2) + 2*b211*pow(sinTheta,2) - 2*b222*pow(sinTheta,2) - 4*b212*sinTwoTheta)*
	(b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (2.*sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*
	pow(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta,1.5)) + 
       (3*pow(2*b112*cosTwoTheta - 2*b111*cosTheta*sinTheta + 2*b122*cosTheta*sinTheta,2)*
	(b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (4.*pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,2.5)*
	sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta)) - 
       ((-2*b111*pow(cosTheta,2) + 2*b122*pow(cosTheta,2) + 2*b111*pow(sinTheta,2) - 2*b122*pow(sinTheta,2) - 4*b112*sinTwoTheta)*
	(b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta))/
       (2.*pow(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta,1.5)*
	sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta))   ;

     return d2p;
   }
}; 


class ClosestApproach
{

 private:

  TMatrixD B1_, B2_;

  double b111_, b112_, b122_;
  double b211_, b212_, b222_;
  double x1c_, x2c_, xc_;
  double y1c_, y2c_, yc_;
  double theta0_;
  double theta_, d_;

  double dist(double);

  TMatrixD translateEllipse(const TMatrixD);

  distanceDerivative derivative_;


 public:
  ClosestApproach(TMatrixD, TMatrixD);
  ~ClosestApproach();

  void setNeutrinoOneEllipse(TMatrixD );
  void setNeutrinoTwoEllipse(TMatrixD );

  void calculate();

  void getClosestApproach(double& ,double& );

};
  

// Go from HEP paper parameterization (Hperp: 3x3, non-symmetric) 
// to liquid crystal paper parameterization (A: 2x2, symmetric)
TMatrixD ClosestApproach::translateEllipse(const TMatrixD Hperp)
{
  TMatrixD H2(2,2);
  H2[0][0]=Hperp[0][0]; H2[0][1]=Hperp[0][1];
  H2[1][0]=Hperp[1][0]; H2[1][1]=Hperp[1][1];
  //TMatrixD B=TMatrixD(H2,TMatrixD::kMultTranspose,H2);
  TMatrixD B(H2); B.Invert();
  return B;
}





double ClosestApproach::dist(double theta)
{
  double cosTheta=cos(theta);
  double sinTheta=sin(theta);
  double cosTwoTheta=cos(2.*theta);
  double sinTwoTheta=sin(2.*theta);

  double b111=b111_, b112=b112_, b122=b122_;
  double b211=b211_, b212=b212_, b222=b222_;

  double d2= 

    pow(xc_,2) + pow(yc_,2) - (2*((b111*xc_ + b112*yc_)*cosTheta + (b112*xc_ + b122*yc_)*sinTheta))/
    sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) + 
    ((pow(b111,2) + pow(b112,2))*pow(cosTheta,2) + (pow(b112,2) + pow(b122,2))*pow(sinTheta,2) + b112*(b111 + b122)*sinTwoTheta)/
    (b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta) - 
    (2*((b211*xc_ + b212*yc_)*cosTheta + (b212*xc_ + b222*yc_)*sinTheta))/sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) + 
    ((pow(b211,2) + pow(b212,2))*pow(cosTheta,2) + (pow(b212,2) + pow(b222,2))*pow(sinTheta,2) + b212*(b211 + b222)*sinTwoTheta)/
    (b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta) + 
    (b111*b211 + 2*b112*b212 + b122*b222 + (b111*b211 - b122*b222)*cosTwoTheta + ((b111 + b122)*b212 + b112*(b211 + b222))*sinTwoTheta)/
    (sqrt(b111*pow(cosTheta,2) + b122*pow(sinTheta,2) + b112*sinTwoTheta)*sqrt(b211*pow(cosTheta,2) + b222*pow(sinTheta,2) + b212*sinTwoTheta))  ;

  return sqrt(d2);
}
