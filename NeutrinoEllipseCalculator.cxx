#include "NeutrinoEllipseCalculator.h"

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator()
{
  //cout << "constructor" << endl;
}

NeutrinoEllipseCalculator::~NeutrinoEllipseCalculator()
{
  //cout << "destructor" << endl;
}

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator(double bJetPx, double bJetPy, double bJetPz, double bJetE,
						     double leptonPx, double leptonPy, double leptonPz, double leptonE,
						     double mTop, double mW, double mNu) 
{
  setBJet(bJetPx,bJetPy,bJetPz,bJetE);
  setLepton(leptonPx,leptonPy,leptonPz,leptonE);

  setBJetRelativisticFactors();
  setLeptonRelativisticFactors();

  setMasses(mTop,mW,mNu);

  setAngles();

  initializeMatrices();
}

void NeutrinoEllipseCalculator::setupEllipse(double bJetPx, double bJetPy, double bJetPz, double bJetE,
					     double leptonPx, double leptonPy, double leptonPz, double leptonE,
					     double mTop, double mW, double mNu) 
{
  setBJet(bJetPx,bJetPy,bJetPz,bJetE);
  setLepton(leptonPx,leptonPy,leptonPz,leptonE);

  setBJetRelativisticFactors();
  setLeptonRelativisticFactors();

  setMasses(mTop,mW,mNu);

  setAngles();

  initializeMatrices();
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

void NeutrinoEllipseCalculator::setBJet(const double px, const double py, const double pz, const double E)
{
  bJet_.SetPxPyPzE(px,py,pz,E);
  setBJetRelativisticFactors();
  setAngles();
}

void NeutrinoEllipseCalculator::setLepton(const double px, const double py, const double pz, const double E)
{
  lepton_.SetPxPyPzE(px,py,pz,E);
  setLeptonRelativisticFactors();
  setAngles();
}

void NeutrinoEllipseCalculator::setBJetRelativisticFactors()
{
  bJetBeta_=bJet_.Beta();
  bJetBeta2_=pow(bJetBeta_,2);
  bJetGamma_=bJet_.Gamma();
  bJetGamma2_=pow(bJetGamma_,2);
}

void NeutrinoEllipseCalculator::setLeptonRelativisticFactors()
{
  leptonBeta_=lepton_.Beta();
  leptonBeta2_=pow(leptonBeta_,2);
  leptonGamma_=lepton_.Gamma();
  leptonGamma2_=pow(leptonGamma_,2);
}

void NeutrinoEllipseCalculator::setAngles()
{
  c_=ROOT::Math::VectorUtil::CosTheta(lepton_,bJet_);
  s_=sqrt(1.-c_*c_);
}

void NeutrinoEllipseCalculator::initializeMatrices()
{
  Ab_=TMatrixD(4,4);
  Al_=TMatrixD(4,4);
  Htilde_=TMatrixD(3,3);
  H_=TMatrixD(3,3);
  Hperp_=TMatrixD(3,3);
  Nperp_=TMatrixD(3,3);

  Ab_.Zero();
  Al_.Zero();
  Htilde_.Zero();
  H_.Zero();
  Hperp_.Zero();
  Nperp_.Zero();
}

void NeutrinoEllipseCalculator::Wsurface()
{
  x0p_=-(0.5/  bJet_.E())*(mt2_- mW2_-bJet_  .M2());
  x0_ =-(0.5/lepton_.E())*(mW2_-mnu2_-lepton_.M2());
  Sx_=(1./leptonBeta2_)*(x0_*leptonBeta_-lepton_.P()*(1.-leptonBeta2_));
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
  if( Z2_<0 ) calcBJetCorrection();
  double Z=sqrt(max(0.,Z2_));

  Htilde_[0][0]=Z/Omega_;
  Htilde_[0][1]=0;
  Htilde_[0][2]=x1_-lepton_.P();
	
  Htilde_[1][0]=omega_*Z/Omega_;
  Htilde_[1][1]=0;
  Htilde_[1][2]=y1_;
  	
  Htilde_[2][0]=0;
  Htilde_[2][1]=Z;
  Htilde_[2][2]=0;

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
  TMatrixD Rz=rotationMatrix(2,-lepton_.Phi());
  TMatrixD Ry=rotationMatrix(1,0.5*M_PI-lepton_.Theta());
  double bJetP[3]={bJet_.Px(),bJet_.Py(), bJet_.Pz()};
  TMatrixD bJet_xyz(3,1,bJetP);
  TMatrixD rM(Ry,TMatrixD::kMult,TMatrixD(Rz,TMatrixD::kMult,bJet_xyz));
  double* rA=rM.GetMatrixArray();
  double phi=-TMath::ATan2(rA[2],rA[1]);
  TMatrixD Rx=rotationMatrix(0,phi);
  R=TMatrixD(Rz,TMatrixD::kTransposeMult,TMatrixD(Ry,TMatrixD::kTransposeMult,Rx.T()));
  H_=TMatrixD(R,TMatrixD::kMult,Htilde_);


  //calculate Hperp
  double Hvalues[9]={H_[0][0],H_[0][1],H_[0][2],H_[1][0],H_[1][1],H_[1][2],0,0,1};
  TArrayD Harray(9,Hvalues);
  Hperp_.SetMatrixArray(Harray.GetArray());


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

TMatrixD* NeutrinoEllipseCalculator::getNeutrinoEllipse()
{
  //calcNeutrinoEllipse();
  return &Nperp_;
}

TMatrixD* NeutrinoEllipseCalculator::getHomogeneousNeutrinoEllipse()
{
  //calcNeutrinoEllipse();
  return &Hperp_;
}

void NeutrinoEllipseCalculator::calcBJetCorrection()
{
  //Find the ranges of b-jet energies in which Z^2>=0:
  //Z^2=0 is equivalent to a polynomial of degree 4 in Eb
  //bJetBeta_=bJetP/bJetE does not change when the energy correction is applied
  //(Cancels in numerator and denominator)
  //All other parameters are independent of the b-jet energy and momentum
  //if Z^2>0 already this will simply return the interval in which the current b-jet energy is

  double a0=-pow(mt2_,2) + 2.*mt2_*mW2_ - pow(mW2_,2) - pow(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_,2)/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))  ;

  double a1=-4*mt2_*Sx_*leptonBeta_ + 4*mW2_*Sx_*leptonBeta_ - (2*(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_)*
								(2*c_*Sx_*bJetBeta_*leptonBeta_ - 2*Sx_*pow(leptonBeta_,2) + 2*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2)))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1 + pow(leptonBeta_,2)))   ;

  double a2=2.*mt2_ - 2.*mW2_ - 2.*mt2_*pow(bJetBeta_,2) + 2.*mW2_*pow(bJetBeta_,2) - 4.*pow(Sx_,2)*pow(leptonBeta_,2) + 
    (-2.*(c_*mt2_*bJetBeta_ - c_*mW2_*bJetBeta_ - mt2_*leptonBeta_ + mW2_*leptonBeta_)*(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_) - 
     pow(2.*c_*Sx_*bJetBeta_*leptonBeta_ - 2.*Sx_*pow(leptonBeta_,2) + 2.*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2),2))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2))) - 4.*pow(s_,2)*pow(bJetBeta_,2)*(mW2_ - pow(x0_,2) - epsilon2_)   ;

  double a3=4.*Sx_*leptonBeta_ - 4.*Sx_*pow(bJetBeta_,2)*leptonBeta_ - (2.*(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_)*
									(2.*c_*Sx_*bJetBeta_*leptonBeta_ - 2.*Sx_*pow(leptonBeta_,2) + 2.*pow(s_,2)*Sx_*pow(bJetBeta_,2)*pow(leptonBeta_,2)))/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))   ;

  double a4=-1. + 2.*pow(bJetBeta_,2) - pow(bJetBeta_,4) - pow(-(c_*bJetBeta_) + c_*pow(bJetBeta_,3) + leptonBeta_ - pow(bJetBeta_,2)*leptonBeta_,2)/
    (-pow(c_*bJetBeta_ - leptonBeta_,2) + pow(s_,2)*pow(bJetBeta_,2)*(-1. + pow(leptonBeta_,2)))   ;

  
  //Solve for the roots of the polynomial
  //Consider only the real positive roots

  Polynomial quartic(a4,a3,a2,a1,a0);
  vector<double > bJetEs = quartic.FindRealRoots();

  for( vector<double>::iterator ibJetE = bJetEs.begin(); ibJetE != bJetEs.end(); ibJetE++ )
    {
      if( *ibJetE<0 ) 
	{
	  //cout << "This root has a negative value: " << *ibJetE << "/nRemoving it from the vector of roots" << endl;
	  bJetEs.erase(ibJetE);
	}
    }

  //Store the interval closest to the measured jet energy where Z^2 is positive as a vector of pairs of doubles and bools: 

  //bJetLogSFRange_.clear();
  double  lowEdge=0; bool  hasLowEdge=false;
  double highEdge=0; bool hasHighEdge=false;
  nRanges_=0;
  if( bJetEs.size() > 0 ) 
    {
      vector<bool> isPositiveSection;
      vector<double> edgesLogScale;
      double e(0.);//, deltaEMin=abs(bJetEs.at(0)-bJet_.E());
      e = 0.5*bJetEs.at(0);
      (a0+a1*e+a2*pow(e,2)+a3*pow(e,3)+a4*pow(e,4) > 0) ? isPositiveSection.push_back(true) : isPositiveSection.push_back(false);
      unsigned int middleBin = 0;
      for( unsigned int i=0; i<bJetEs.size(); i++ )
	{
	  ( i == bJetEs.size()-1 ) ? e = 2.*bJetEs.at(bJetEs.size()-1) : e = 0.5*(bJetEs.at(i)+bJetEs.at(i+1));
	  if(a0+a1*e+a2*pow(e,2)+a3*pow(e,3)+a4*pow(e,4) > 0)
	    {
	      if(*(isPositiveSection.rbegin())) continue;
	      isPositiveSection.push_back(true);
	    }
	  else
	    {
	      if(!*(isPositiveSection.rbegin())) continue;
	      isPositiveSection.push_back(false);
	    }
	  double edge = log(bJetEs.at(i)/bJet_.E());
	  edgesLogScale.push_back(edge);
	  if(edge < 0 ) middleBin++;
	}
      if(isPositiveSection[middleBin])
	{
	  if(middleBin != isPositiveSection.size()-1)
	    {
	      hasHighEdge=true;
	      highEdge=edgesLogScale[middleBin];
	    }
	  if(middleBin != 0)
	    {
	      hasLowEdge=true;
	      lowEdge=edgesLogScale[middleBin-1];
	    }
	}
      else if(isPositiveSection.size() > 1) 
	{
	  if(middleBin == 0)
	    {
	      if(2 != isPositiveSection.size())
		{
		  hasHighEdge=true;
		  highEdge=edgesLogScale[1];
		}
	      hasLowEdge=true;
	      lowEdge=edgesLogScale[0];
	    }
	  else if(middleBin == isPositiveSection.size() - 1)
	    {
	      hasHighEdge=true;
	      highEdge=edgesLogScale[middleBin-1];
	      if(middleBin-1 != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-2];
		}
	    }
	  else if(abs(edgesLogScale[middleBin-1])>edgesLogScale[middleBin])
	    {
	      middleBin-=1;
	      if(middleBin != isPositiveSection.size()-1)
		{
		  hasHighEdge=true;
		  highEdge=edgesLogScale[middleBin];
		}
	      if(middleBin != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-1];
		}
	    }
	  else
	    {
	      middleBin+=1;
	      if(middleBin != isPositiveSection.size()-1)
		{
		  hasHighEdge=true;
		  highEdge=edgesLogScale[middleBin];
		}
	      if(middleBin != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-1];
		}
	    }
	}
      for(vector<bool>::iterator thisRange = isPositiveSection.begin(); thisRange != isPositiveSection.end(); thisRange++)
	{
	  if(*thisRange) nRanges_++;
	}
      
    }

  bJetLogSFRange_.first =  make_pair( lowEdge, hasLowEdge) ;
  bJetLogSFRange_.second = make_pair(highEdge,hasHighEdge) ;

  return;
}

pair<pair<double,bool>, pair<double,bool> > NeutrinoEllipseCalculator::getBJetLogSFRange(int& nRanges)
{
  nRanges = nRanges_;
  return bJetLogSFRange_;
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
