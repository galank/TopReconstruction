#include "ClosestApproach.h"


ClosestApproach::ClosestApproach(TMatrixD nu1, TMatrixD nu2) :
  B1_(TMatrixD(2,2)),
  B2_(TMatrixD(2,2)),
  theta_(-99.),
  d_(-1.)
{

  //2x2 representations of the ellipses
  //cout << "calculating the symmetric 2x2 representation of the neutrino ellipses" << endl;

  B1_=translateEllipse(nu1);
  b111_=B1_[0][0]; b112_=B1_[0][1]; b122_=B1_[1][1];
  cout << b111_ << " " << b112_ << " " << b122_ << endl;

  B2_=translateEllipse(nu2);
  b211_=B2_[0][0]; b212_=B2_[0][1]; b222_=B2_[1][1];
  cout << b211_ << " " << b212_ << " " << b222_ << endl;

  //Centers of the ellipses
  x1c_=nu1[0][2]; y1c_=nu1[1][2];
  x2c_=nu2[0][2]; y2c_=nu2[1][2];
  xc_=x2c_-x1c_;
  yc_=y2c_-y1c_;
  cout << xc_ << " " << yc_ << endl;

  derivative_.setB1(b111_,b112_,b122_);
  derivative_.setB2(b211_,b212_,b222_);
  derivative_.setCenter(xc_,yc_);


  //Initial guess
  //theta0_=atan2(yc_,xc_);
  theta0_=atan2(yc_,xc_)+TMath::Pi()/3.;
  
}

ClosestApproach::~ClosestApproach()
{
}

void ClosestApproach::setNeutrinoOneEllipse(TMatrixD nu)
{
  B1_=translateEllipse(nu);
  b111_=B1_[0][0]; b112_=B1_[0][1]; b122_=B1_[1][1];

  x1c_=nu[0][2]; y1c_=nu[1][2];
  xc_=x2c_-x1c_;
  yc_=y2c_-y1c_;

  derivative_.setB1(b111_,b112_,b122_);
}

void ClosestApproach::setNeutrinoTwoEllipse(TMatrixD nu)
{
  B2_=translateEllipse(nu);
  b211_=B2_[0][0]; b212_=B2_[0][1]; b222_=B1_[1][1];

  x2c_=nu[0][2]; y2c_=nu[1][2];
  xc_=x2c_-x1c_;
  yc_=y2c_-y1c_;

  derivative_.setB2(b211_,b212_,b222_);
}

void ClosestApproach::calculate()
{
  //Find the minimum of the distance of closest approach between the neutrinos
  //Three possible cases depending on the type of neutrino solution: point or ellipse


  //Case 1: Two points
  //The closest approach is the distance between the points
  //No minimization necessary
  if( b111_==0 && b112_==0 && b122_==0 && b211_==0 && b212_==0 && b222_==0)
    {
      d_=sqrt(xc_*xc_+yc_*yc_);
      theta_=0;
      return;
    }


  //Case 2: One point and one ellipse
  //Case 3: Two ellipses
  //Minimization is necessary

  else
    {
      
      //Method 1:
      //using the function itself and a minimum finder
      //ROOT::Math::Functor1D func(&dist);
      //ROOT::Math::GSLMinimizer1D* min = new ROOT::Math::GSLMinimizer1D;
      //min->SetFunction(func,TMath::Pi(),0.,2.*TMath::Pi());
      //int maxIter=10000;
      //double absTol=0.001, relTol=0.001;
      //bool hasSol=min->Minimize(maxIter,absTol,relTol);
      //cout << "Number of iterations: " << min->Iterations() << endl;
      //cout << "Status: " << min->Status() << endl;
      //if(hasSol)
      //  {
      //    double thetaMin=min->XMinimum();
      //    cout << "The minimum is: " << thetaMin << endl;
      //    return dist(thetaMin);
      //  }
      //return -1;
      
      
      //Method 2:
      //using the derivative and a root finder
      ROOT::Math::Roots::Steffenson rf;
      rf.SetFunction(derivative_,theta0_);
      int maxIter=10000; 
      //double absTol=1E-8, relTol=1E-10;
      bool hasSol=rf.Solve(maxIter);
      cout << "Number of iterations: " << rf.Iterations() << endl;
      //cout << "Status: " << rf.Status() << endl;
      if(hasSol) 
	{
	  //cout << "The root is: " << rf.Root() << endl;
	  //cout << "The distance of closest approach is " << dist(rf.Root()) << endl;
	  //return dist(rf.Root());
	  theta_=rf.Root();
	  d_=dist(theta_);
	  return;
	}
      return;

    }
}

void ClosestApproach::getClosestApproach(double& d, double& theta)
{
  calculate();
  d=d_;
  theta=theta_;
  return;
}
