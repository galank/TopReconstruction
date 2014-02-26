#include "NeutrinoEllipseCalculator.h"

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator()
{
  //cout << "constructor" << endl;
}

NeutrinoEllipseCalculator::~NeutrinoEllipseCalculator()
{
  //cout << "destructor" << endl;
  delete bJet_;
  delete lepton_;
  delete Ab_;
  delete Al_;
  delete Htilde_;
  delete H_;
  delete Hperp_;
  delete Nperp_;
}

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator(//const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& bJetLorentzVector,
						     //const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& leptonLorentzVector,
						     //TLorentzVector& bJetLorentzVector, TLorentzVector& leptonLorentzVector,
						     //math::XYZTLorentzVector& bJetLorentzVector, 
						     double bJetPx, double bJetPy, double bJetPz, double bJetE,
						     //math::XYZTLorentzVector& leptonLorentzVector,
						     double leptonPx, double leptonPy, double leptonPz, double leptonE,
						     double& mTop, double& mW, double& mNu) :
  bJet_  (new math::XYZTLorentzVector(0,0,0,0)),
  lepton_(new math::XYZTLorentzVector(0,0,0,0)),
  Ab_(new TMatrixD(4,4)),
  Al_(new TMatrixD(4,4)),
  Htilde_(new TMatrixD(3,3)),
  H_(new TMatrixD(3,3)),
  Hperp_(new TMatrixD(3,3)),
  Nperp_(new TMatrixD(3,3))
{
  //setbJet(bJetLorentzVector);
  //setLepton(leptonLorentzVector);
  setBJet(bJetPx,bJetPy,bJetPz,bJetE);
  setLepton(leptonPx,leptonPy,leptonPz,leptonE);

  setBJetRelativisticFactors();
  setLeptonRelativisticFactors();

  setMasses(mTop,mW,mNu);

  setAngles();

  initializeMatrices();

  bJetSF_=1.;
}

void NeutrinoEllipseCalculator::setMasses(double& mTop, double& mW, double& mNu)
{
  setTopMass(mTop);
  setWBosonMass(mW);
  setNeutrinoMass(mNu);

  mW2_=mW_*mW_;
  mt2_=mt_*mt_;
  mnu2_=mnu_*mnu_;
}

//void NeutrinoEllipseCalculator::setbJet(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& bJetLorentzVector)
//void NeutrinoEllipseCalculator::setbJet(const TLorentzVector& bJetLorentzVector)
//void NeutrinoEllipseCalculator::setbJet(const math::XYZTLorentzVector& bJetLorentzVector)
void NeutrinoEllipseCalculator::setBJet(const double px, const double py, const double pz, const double E)
{
  //bJet_=bJetLorentzVector;
  bJet_->SetPxPyPzE(px,py,pz,E);
  setBJetRelativisticFactors();
  setAngles();
}

//void NeutrinoEllipseCalculator::setLepton(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& leptonLorentzVector)
//void NeutrinoEllipseCalculator::setLepton(const TLorentzVector& leptonLorentzVector)
//void NeutrinoEllipseCalculator::setLepton(const math::XYZTLorentzVector& leptonLorentzVector)
void NeutrinoEllipseCalculator::setLepton(const double px, const double py, const double pz, const double E)
{
  //lepton_=leptonLorentzVector;
  lepton_->SetPxPyPzE(px,py,pz,E);
  setLeptonRelativisticFactors();
  setAngles();
}

void NeutrinoEllipseCalculator::setBJetRelativisticFactors()
{
  bJetBeta_=bJet_->Beta();
  bJetBeta2_=pow(bJetBeta_,2);
  bJetGamma_=bJet_->Gamma();
  bJetGamma2_=pow(bJetGamma_,2);
}

void NeutrinoEllipseCalculator::setLeptonRelativisticFactors()
{
  leptonBeta_=lepton_->Beta();
  leptonBeta2_=pow(leptonBeta_,2);
  leptonGamma_=lepton_->Gamma();
  leptonGamma2_=pow(leptonGamma_,2);
}

void NeutrinoEllipseCalculator::setAngles()
{
  c_=ROOT::Math::VectorUtil::CosTheta(*lepton_,*bJet_);
  s_=sqrt(1.-c_*c_);
}

void NeutrinoEllipseCalculator::initializeMatrices()
{
  Ab_->Zero();
  Al_->Zero();
  Htilde_->Zero();
  H_->Zero();
  Hperp_->Zero();
  Nperp_->Zero();
}

void NeutrinoEllipseCalculator::Wsurface()
{
  x0p_=-(0.5/  bJet_->E())*(mt2_- mW2_-bJet_  ->M2());
  x0_ =-(0.5/lepton_->E())*(mW2_-mnu2_-lepton_->M2());
  Sx_=(1./leptonBeta2_)*(x0_*leptonBeta_-lepton_->P()*(1.-leptonBeta2_));
  epsilon2_=(1.-leptonBeta2_)*(mW2_-mnu2_);
}

void NeutrinoEllipseCalculator::bJetEllipsoid()
{
  Ab_[0][0]=1-pow(c_,2)*bJetBeta2_;
  Ab_[1][0]=-c_*s_*bJetBeta2_;
  Ab_[2][0]=0;
  Ab_[3][0]=c_*x0p_*bJetBeta_;
    
  Ab_[0][1]=-c_*s_*bJetBeta2_;
  Ab_[1][1]=1-pow(s_,2)*bJetBeta2_;
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
}

void NeutrinoEllipseCalculator::leptonEllipsoid()
{
  Al_[0][0]=1.-leptonBeta2_;
  Al_[1][0]=0;
  Al_[2][0]=0;
  Al_[3][0]=Sx_*leptonBeta2_;
    
  Al_[0][1]=0;
  Al_[1][1]=1;
  Al_[2][1]=0;
  Al_[3][1]=0;
    
  Al_[0][2]=0;
  Al_[1][2]=0;
  Al_[2][2]=1;
  Al_[3][2]=0;
    
  Al_[0][3]=Sx_*leptonBeta2_;
  Al_[1][3]=0;
  Al_[2][3]=0;
  Al_[3][3]=mW2_-x0_*x0_-epsilon2_;
}

void NeutrinoEllipseCalculator::neutrinoSolution()
{
  Sy_=(1./s_)*(x0p_/bJetBeta_-c_*Sx_);
  omega_=(1./s_)*(leptonBeta_/bJetBeta_-c_); //only the positive slope
  Omega_=sqrt(max(0.,omega_*omega_+1.-leptonBeta2_));
  double Omega2=Omega_*Omega_;
  x1_=Sx_-(1./Omega2)*(Sx_+omega_*Sy_);
  y1_=Sy_-(1./Omega2)*omega_*(Sx_+omega_*Sy_);
  Z2_=x1_*x1_*Omega2-(Sy_-omega_*Sx_)*(Sy_-omega_*Sx_)-(mW2_-x0_*x0_-epsilon2_);
  double Z=sqrt(max(0.,Z2_));

  Htilde_[0][0]=Z/Omega_;
  Htilde_[0][1]=0;
  Htilde_[0][2]=x1_-lepton_->P();
	
  Htilde_[1][0]=omega_*Z/Omega_;
  Htilde_[1][1]=0;
  Htilde_[1][2]=y1_;
  	
  Htilde_[2][0]=0;
  Htilde_[2][1]=Z;
  Htilde_[2][2]=0;

  if( Z2_<0 ) 
    {
      //cout << "Z^2 is <0" << endl;
      //calcbJetCorrection();
    }

}

TMatrixD NeutrinoEllipseCalculator::rotationMatrix(int axis, double angle)
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

void NeutrinoEllipseCalculator::labSystemTransform()
{
  //rotate Htilde to H
  TMatrixD R(3,3);
  R.Zero();
  TMatrixD Rz=rotationMatrix(2,-lepton_->Phi());
  TMatrixD Ry=rotationMatrix(1,0.5*M_PI-lepton_->Theta());
  double bJetP[3]={bJet_->Px(),bJet_->Py(), bJet_->Pz()};
  TMatrixD bJet_xyz(3,1,bJetP);
  TMatrixD rM(Ry,TMatrixD::kMult,TMatrixD(Rz,TMatrixD::kMult,bJet_xyz));
  double* rA=rM.GetMatrixArray();
  double phi=-TMath::ATan2(rA[2],rA[1]);
  TMatrixD Rx=rotationMatrix(0,phi);
  R=TMatrixD(Rz,TMatrixD::kTransposeMult,TMatrixD(Ry,TMatrixD::kTransposeMult,Rx.T()));
  //double* Ha=(TMatrixD(R,TMatrixD::kMult,*Htilde_)).GetMatrixArray();
  //H_->SetMatrixArray(Ha);
  (*H_)=TMatrixD(R,TMatrixD::kMult,*Htilde_);


  //calculate Hperp
  double Hvalues[9]={(*H_)[0][0],(*H_)[0][1],(*H_)[0][2],(*H_)[1][0],(*H_)[1][1],(*H_)[1][2],0,0,1};
  TArrayD Harray(9,Hvalues);
  Hperp_->SetMatrixArray(Harray.GetArray());


  //calculate Nperp
  //TMatrixD HperpInv(Hperp_);
  //HperpInv.Invert();
  //TMatrixD U(3,3);
  //U.Zero();
  //U[0][0]=1;
  //U[1][1]=1;
  //U[2][2]=-1;
  //Nperp_=TMatrixD(HperpInv,TMatrixD::kTransposeMult,TMatrixD(U,TMatrixD::kMult,HperpInv));
}

void NeutrinoEllipseCalculator::calcNeutrinoEllipse()
{
  Wsurface();
  leptonEllipsoid();
  bJetEllipsoid();
  neutrinoSolution();
  labSystemTransform();
}

TMatrixD NeutrinoEllipseCalculator::getNeutrinoEllipse()
{
  calcNeutrinoEllipse();
  TMatrixD neutrinoEllipse=*Nperp_;
  return neutrinoEllipse;
}

TMatrixD NeutrinoEllipseCalculator::getHomogeneousNeutrinoEllipse()
{
  calcNeutrinoEllipse();
  TMatrixD neutrinoEllipse=*Hperp_;
  return neutrinoEllipse;
}


void NeutrinoEllipseCalculator::calcBJetCorrection()
{
  if( Z2_>0. ) return;

  //If Z^2<0, solve for the b-jet energy necessary to make Z^2 exactly zero
  //In this configuration H11=H12=H21=H22=0 and the neutrino solution is a
  //point (H13,H23)

  //Z^2=0 is equivalent to a polynomial of degree 4 in Eb
  //bJetBeta_=bJetP/bJetE does not change when the energy correction is applied 
  //(Cancels in numerator and denominator)
  //All other parameters are independent of the b-jet energy and momentum

  double r0=-pow(mt2_,2) + 2.*mt2_*mW2_ - pow(mW2_,2) - pow(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_,2)/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))  ;
  double r1=-4*mt2_*Sx_*leptonBeta_ + 4*mW2_*Sx_*leptonBeta_ - (2*(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_)*
								(2*c_*Sx_*bJetBeta_*leptonBeta_ - 2*Sx_*pow(leptonBeta_,2) + 2*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2)))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1 + pow(leptonBeta_,2)))   ;
  double r2=2.*mt2_ - 2.*mW2_ - 2.*mt2_*pow(bJetBeta_,2) + 2.*mW2_*pow(bJetBeta_,2) - 4.*pow(Sx_,2)*pow(leptonBeta_,2) + 
    (-2.*(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_)*(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_) - 
     pow(2.*c_*Sx_*bJetBeta_*leptonBeta_ - 2.*Sx_*pow(leptonBeta_,2) + 2.*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2),2))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2))) - 4.*pow(s_,2)*pow(bJetBeta_,2)*(mW2_ - pow(x0_,2) - epsilon2_)   ;
  double r3=4.*Sx_*leptonBeta_ - 4.*Sx_*pow(bJetBeta_,2)*leptonBeta_ - (2.*(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_)*
									(2.*c_*Sx_*bJetBeta_*leptonBeta_ - 2.*Sx_*pow(leptonBeta_,2) + 2.*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2)))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))   ;
  double r4=-1. + 2.*pow(bJetBeta_,2) - pow(bJetBeta_,4) - pow(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_,2)/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))   ;

  
  //Solve for the roots of the polynomial
  //Consider only the real positive roots

  Polynomial* quartic = new Polynomial(r4,r3,r2,r1,r0);
  vector<double > bJetEValues = quartic->FindRealRoots();

  for( vector<double>::iterator ibJetE = bJetEValues.begin(); ibJetE != bJetEValues.end(); ibJetE++ )
    {
      double bJetE=*ibJetE;
      if( bJetE<0. ) 
	{
	  cout << "This root has a negative value: " << bJetE << "/nRemoving it from the vector of roots" << endl;
	  bJetEValues.erase(ibJetE);
	}
    }


  //Find the root closest to the measured b-jet energy
  
  double minDelta=1.e9;
  double closestEnergy=-1.;

  for( vector<double>::iterator ibJetE = bJetEValues.begin(); ibJetE != bJetEValues.end(); ibJetE++ )
    {
      double bJetE=*ibJetE;
      double thisDelta=abs(bJetE-bJet_->E());
      if( thisDelta<minDelta )
	{
	  minDelta=thisDelta;
	  closestEnergy=bJetE;
	}
    }


  //Energy scale factor
  if( closestEnergy>0 ) bJetSF_=closestEnergy/bJet_->E();

}


void NeutrinoEllipseCalculator::print3By3Matrix(const TMatrixD& m)
{
  cout << m[0][0] << " " << m[0][1] << " " << m[0][2] << "\n"
       << m[1][0] << " " << m[1][1] << " " << m[1][2] << "\n"
       << m[2][0] << " " << m[2][1] << " " << m[2][2] << endl;
}

void NeutrinoEllipseCalculator::print4By4Matrix(const TMatrixD& m)
{
  cout << m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << m[0][3] << "\n"
       << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << m[1][3] << "\n"
       << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << m[2][3] << "\n"
       << m[3][0] << " " << m[3][1] << " " << m[3][2] << " " << m[3][3] << endl;
}
