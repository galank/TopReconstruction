#include "NeutrinoEllipseCalculator.h"
#include <complex>

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
						     const double& bJetPtWidth, const double& bJetPhiWidth, const double& bJetEtaWidth,
						     const double& leptonPx,const  double& leptonPy,const double& leptonPz, const double& leptonE,
						     const double& unscaled_bJetPx, const double& unscaled_bJetPy, const double& unscaled_bJetPz, 
						     const double& unscaled_bJetP, const double& unscaled_bJetE)  :
  bJetPx_(bJetPx),
  bJetPy_(bJetPy),
  bJetPz_(bJetPz),
  bJetE_ (bJetE),
  bJetPtWidth_(bJetPtWidth),
  bJetPhiWidth_(bJetPhiWidth),
  bJetEtaWidth_(bJetEtaWidth),
  leptonPx_(leptonPx),
  leptonPy_(leptonPy),
  leptonPz_(leptonPz),
  leptonE_ (leptonE),
  unscaled_bJetPx_(unscaled_bJetPx),
  unscaled_bJetPy_(unscaled_bJetPy),
  unscaled_bJetPz_(unscaled_bJetPz),
  unscaled_bJetE_ (unscaled_bJetE),
  unscaled_bJetP_ (unscaled_bJetP),
  Ab_(4,4),
  Al_(4,4),
  Htilde_(3,3),
  H_(3,3),
  Hperp_(3,3),
  HperpInv_(3,3),
  Nperp_(3,3),
  nuPerp_(3),
  pNu_(3),
  errorFlag_(false)
{
  initializeMatrices();
}

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
						     const double& bJetPtWidth, const double& bJetPhiWidth, const double& bJetEtaWidth,
						     const double& leptonPx, const double& leptonPy, const double& leptonPz, const double& leptonE,
						     const double& unscaled_bJetPx, const double& unscaled_bJetPy, const double& unscaled_bJetPz, 
						     const double& unscaled_bJetP, const double& unscaled_bJetE,
						     double mTop, double mW, double mNu)  :
  bJetPx_(bJetPx),
  bJetPy_(bJetPy),
  bJetPz_(bJetPz),
  bJetE_ (bJetE),
  bJetPtWidth_(bJetPtWidth),
  bJetPhiWidth_(bJetPhiWidth),
  bJetEtaWidth_(bJetEtaWidth),
  leptonPx_(leptonPx),
  leptonPy_(leptonPy),
  leptonPz_(leptonPz),
  leptonE_ (leptonE),
  unscaled_bJetPx_(unscaled_bJetPx),
  unscaled_bJetPy_(unscaled_bJetPy),
  unscaled_bJetPz_(unscaled_bJetPz),
  unscaled_bJetE_ (unscaled_bJetE),
  unscaled_bJetP_ (unscaled_bJetP),
  Ab_(4,4),
  Al_(4,4),
  Htilde_(3,3),
  H_(3,3),
  Hperp_(3,3),
  HperpInv_(3,3),
  Nperp_(3,3),
  nuPerp_(3),
  pNu_(3),
  errorFlag_(false)
{
  setBJetFactors();
  setLeptonFactors();

  setMasses(mTop,mW,mNu);

  setAngles();

  initializeMatrices();
}

NeutrinoEllipseCalculator::~NeutrinoEllipseCalculator()
{
  //cout << "destructor" << endl;
}


void NeutrinoEllipseCalculator::setupEllipse(double mTop, double mW, double mNu)
{
  //cout << "set up the neutrino ellipse" << endl;
  setBJetFactors();
  setLeptonFactors();
  setMasses(mTop,mW,mNu);
  setAngles();
  initializeMatrices();
  Wsurface();
  leptonEllipsoid();
  bJetEllipsoid();
  calcBJetCorrection();
  //neutrinoSolution();
  //labSystemTransform();
}

double NeutrinoEllipseCalculator::getZ2(double mTop, double mW, double mNu)
{
  setBJetFactors();
  setLeptonFactors();
  setMasses(mTop,mW,mNu);
  setAngles();
  initializeMatrices();
  Wsurface();
  leptonEllipsoid();
  bJetEllipsoid();
  calcZ2();
  return Z2_;
}

void NeutrinoEllipseCalculator::setMasses(double& mTop, double& mW, double& mNu)
{
  mt_=mTop;
  mW_=mW;
  mnu_=mNu; 
  mW2_=mW_*mW_;
  mt2_=mt_*mt_;
  mnu2_=mnu_*mnu_;
}

//void NeutrinoEllipseCalculator::setBJet(const double px, const double py, const double pz, const double E)
//{
//  //bJet_.SetPxPyPzE(px,py,pz,E);
//  setBJetRelativisticFactors();
//  //setAngles();
//}
//
//void NeutrinoEllipseCalculator::setLepton(const double px, const double py, const double pz, const double E)
//{
//  //lepton_.SetPxPyPzE(px,py,pz,E);
//  setLeptonRelativisticFactors();
//  //setAngles();
//}

void NeutrinoEllipseCalculator::setBJetFactors()
{
  bJetP2_ = (bJetPx_*bJetPx_ + bJetPy_*bJetPy_ + bJetPz_*bJetPz_);
  bJetP_ = sqrt(bJetP2_);
  double bJetE2 =  bJetE_*bJetE_;
  bJetE_nonConst_ = bJetE_;
  bJetMass2_ = bJetE2 - bJetP2_;
  bJetBeta2_=bJetP2_/bJetE2;
  bJetBeta_=sqrt(bJetBeta2_);
  bJetGamma2_=1.0/(1.0 - bJetBeta2_);
  bJetGamma_=sqrt(bJetGamma2_);
}

void NeutrinoEllipseCalculator::setTempBJetFactors(double bJetPx_temp, double bJetPy_temp, double bJetPz_temp, double bJetE_temp)
{
  bJetP2_ = (bJetPx_temp*bJetPx_temp + bJetPy_temp*bJetPy_temp + bJetPz_temp*bJetPz_temp);
  bJetP_ = sqrt(bJetP2_);
  double bJetE2 =  bJetE_temp*bJetE_temp;
  bJetE_nonConst_ = bJetE_temp;
  bJetMass2_ = bJetE2 - bJetP2_;
  bJetBeta2_=bJetP2_/bJetE2;
  bJetBeta_=sqrt(bJetBeta2_);
  bJetGamma2_=1.0/(1.0 - bJetBeta2_);
  bJetGamma_=sqrt(bJetGamma2_);
}

void NeutrinoEllipseCalculator::setLeptonFactors()
{
  leptonP2_ = (leptonPx_*leptonPx_ + leptonPy_*leptonPy_ + leptonPz_*leptonPz_);
  leptonP_ = sqrt(leptonP2_);
  double leptonE2 =  leptonE_*leptonE_;
  leptonMass2_ = leptonE2 - leptonP2_;
  leptonBeta2_=leptonP2_/leptonE2;
  leptonBeta_=sqrt(leptonBeta2_);
  leptonGamma2_=1.0/(1.0 - leptonBeta2_);
  leptonGamma_=sqrt(leptonGamma2_);
  leptonPhi_ = (leptonPx_ == 0.0 && leptonPy_ == 0.0) ? 0.0 : atan2(leptonPy_,leptonPx_);
  leptonTheta_ = (leptonPx_ == 0.0 && leptonPy_ == 0.0 && leptonPz_ == 0.0) ? 0.0 : atan2(sqrt(leptonPx_*leptonPx_ + leptonPy_*leptonPy_),leptonPz_);
}

void NeutrinoEllipseCalculator::setAngles()
{
  double leptonDotbJet = (leptonPx_*bJetPx_ + leptonPy_*bJetPy_ + leptonPz_*bJetPz_);
  c2_ = leptonDotbJet*leptonDotbJet/(bJetP2_*leptonP2_);
  c_=leptonDotbJet/sqrt(bJetP2_*leptonP2_);
  s2_ = 1.-c2_;
  s_=sqrt(s2_);
}

void NeutrinoEllipseCalculator::setTempAngles(double bJetPx_temp, double bJetPy_temp, double bJetPz_temp)
{
  double leptonDotbJet = (leptonPx_*bJetPx_temp + leptonPy_*bJetPy_temp + leptonPz_*bJetPz_temp);
  c_=leptonDotbJet/sqrt(bJetP2_*leptonP2_);
  s_=sqrt(1.-leptonDotbJet*leptonDotbJet/(bJetP2_*leptonP2_));
}

void NeutrinoEllipseCalculator::initializeMatrices()
{
  //cout << "initializing matrices" << endl;

  Ab_.Zero();
  Al_.Zero();
  Htilde_.Zero();
  H_.Zero();
  Hperp_.Zero();
  HperpInv_.Zero();
  Nperp_.Zero();

  nuPerp_.Zero();
  pNu_.Zero();

  //cout << "the b-jet matrix has " << Ab_.GetNrows() << " rows and " << Ab_.GetNcols() << " columns"  << endl;
  //cout << "the lepton matrix has " << Al_.GetNrows() << " rows and " << Al_.GetNcols() << " columns"  << endl;
}

void NeutrinoEllipseCalculator::Wsurface()
{
  //cout << "calculating the W surface" << endl;
  x0p_=-(0.5/  bJetE_nonConst_)*(mt2_- mW2_-bJetMass2_);
  x0_ =-(0.5/leptonE_)*(mW2_-mnu2_-leptonMass2_);
  Sx_=x0_/leptonBeta_-leptonP_/leptonBeta2_+leptonP_;
  epsilon2_=(mW2_-mnu2_)-leptonBeta2_*(mW2_-mnu2_);
}

void NeutrinoEllipseCalculator::bJetEllipsoid()
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
}

void NeutrinoEllipseCalculator::leptonEllipsoid()
{
  //cout << "calculating the lepton ellipsoid" << endl;

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

void NeutrinoEllipseCalculator::calcZ2()
{
  Sy_=(1./s_)*(x0p_/bJetBeta_-c_*Sx_);
  omega_=(1./s_)*(leptonBeta_/bJetBeta_-c_); //only the positive slope
  double Omega2=omega_*omega_+1.-leptonBeta2_;
  Omega_=sqrt(Omega2);
  x1_=Sx_-(1./Omega2)*(Sx_+omega_*Sy_);
  y1_=Sy_-(1./Omega2)*omega_*(Sx_+omega_*Sy_);
  Z2_=x1_*x1_*Omega2-(Sy_-omega_*Sx_)*(Sy_-omega_*Sx_)-(mW2_-x0_*x0_-epsilon2_);
}

void NeutrinoEllipseCalculator::neutrinoSolution()
{
  //cout << "calculating neutrino ellipse" << endl;

  calcZ2();
  if( Z2_<0 ) 
    {
      
      cout << "Z2 is " << Z2_ << " and bjet energy is " << bJetE_nonConst_ << "\n"
	   << "b jet energy or momentum out of range for real solutions!\n"
	   << "energy scaling factor is " << bJetE_nonConst_/unscaled_bJetE_ << "\n"
	   << "momentum scaling factor is " << bJetP_/unscaled_bJetP_ << endl;
      cout << "Range from inside: (";
      if(!bJetLogSFRange_.first.second)  cout << "0,";
      else cout << exp(bJetLogSFRange_.first.first) << ",";
      if(!bJetLogSFRange_.second.second)  cout << "inf)" << endl;
      else cout << exp(bJetLogSFRange_.second.first) << ")" << endl;
      cout << "Log range from inside: (";
      if(!bJetLogSFRange_.first.second)  cout << "-inf,";
      else cout << bJetLogSFRange_.first.first << ",";
      if(!bJetLogSFRange_.second.second)  cout << "inf)" << endl;
      else cout << bJetLogSFRange_.second.first << ")" << endl;

      cout << "bJet Scaled Information:\n"
	   << "Px: " << bJetPx_
	   << "\nPy: " << bJetPy_
	   << "\nPz: " << bJetPz_
	   << "\nP : " << bJetP_
	   << "\nE : " << bJetE_ << endl;

      cout << "bJet Unscaled Information:\n"
	   << "Px: " << unscaled_bJetPx_
	   << "\nPy: " << unscaled_bJetPy_
	   << "\nPz: " << unscaled_bJetPz_
	   << "\nP : " << unscaled_bJetP_
	   << "\nE : " << unscaled_bJetE_ << endl;

      cout << "bJet Ratio Information:\n"
	   << "Px: " <<   bJetPx_/unscaled_bJetPx_
	   << "\nPy: " << bJetPy_/unscaled_bJetPy_
	   << "\nPz: " << bJetPz_/unscaled_bJetPz_
	   << "\nP : " << bJetP_ /unscaled_bJetP_
	   << "\nE : " << bJetE_ /unscaled_bJetE_ << endl;

      errorFlag_ = true;
      Z2_ = 0;

      double sf = bJetP_/unscaled_bJetP_;
      setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
      setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
      Wsurface();
      leptonEllipsoid();
      bJetEllipsoid();
      calcZ2();

      cout << "Z2 is " << Z2_ << " using internal functions and scaling factor " << sf << endl;

      cout << "New bJet Unscaled Information:\n"
	   << "Px: " <<   sf* unscaled_bJetPx_
	   << "\nPy: " << sf* unscaled_bJetPy_
	   << "\nPz: " << sf* unscaled_bJetPz_
	   << "\nP : " << sf* unscaled_bJetP_ << " or " << bJetP_ << " or  " << sqrt(bJetPx_*bJetPx_ + bJetPy_*bJetPy_ + bJetPz_*bJetPz_)
	   << "\nE : " << sf* unscaled_bJetE_ << " or " << bJetE_ << endl;

      cout << "bJet Difference Information:\n"
	   << "Px: " <<   sf*unscaled_bJetPx_ - bJetPx_
	   << "\nPy: " << sf*unscaled_bJetPy_ - bJetPy_
	   << "\nPz: " << sf*unscaled_bJetPz_ - bJetPz_
	   << "\nP : " << bJetP_ - sqrt(bJetPx_*bJetPx_ + bJetPy_*bJetPy_ + bJetPz_*bJetPz_)
	   << "\nE : " << bJetE_nonConst_ - bJetE_ << endl;
  
      double tempZ2 = Z2_;
      setBJetFactors();
      setAngles();
      Wsurface();
      calcZ2();
      Z2_ = tempZ2;
    }
  else errorFlag_ = false;
  double Z=sqrt(Z2_);

  Htilde_[0][0]=Z/Omega_;
  Htilde_[0][1]=0;
  Htilde_[0][2]=x1_-leptonP_;
	
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
  ////rotate Htilde to H
  //TMatrixD R(3,3);
  //R.Zero();
  //TMatrixD Rz=rotationMatrix(2,-leptonPhi_);
  //TMatrixD Ry=rotationMatrix(1,0.5*M_PI-leptonTheta_);
  //double bJetP[3]={bJetPx_,bJetPy_, bJetPz_};
  //TMatrixD bJet_xyz(3,1,bJetP);
  //TMatrixD rM(Ry,TMatrixD::kMult,TMatrixD(Rz,TMatrixD::kMult,bJet_xyz));
  //double* rA=rM.GetMatrixArray();
  //double phi=-TMath::ATan2(rA[2],rA[1]);
  //TMatrixD Rx=rotationMatrix(0,phi);
  //R=TMatrixD(Rz,TMatrixD::kTransposeMult,TMatrixD(Ry,TMatrixD::kTransposeMult,Rx.T()));
  //H_=TMatrixD(R,TMatrixD::kMult,Htilde_);

  TMatrixD Rz=rotationMatrix(2,-leptonPhi_);
  TMatrixD Ry=rotationMatrix(1,0.5*M_PI-leptonTheta_);
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
  setBJetFactors();
  //setLeptonFactors();  //only needed if lepton uncertainty included
  setAngles();
  Wsurface();
  leptonEllipsoid();
  bJetEllipsoid();
  neutrinoSolution();
  labSystemTransform();
}

void NeutrinoEllipseCalculator::calcExtendedNeutrinoEllipse()
{
  HperpInv_=Hperp_;
  HperpInv_.Invert();
  TMatrixD U(3,3);
  U.Zero();
  U[0][0]=1;
  U[1][1]=1;
  U[2][2]=-1;
  Nperp_=TMatrixD(HperpInv_,TMatrixD::kTransposeMult,TMatrixD(U,TMatrixD::kMult,HperpInv_));
}

TMatrixD* NeutrinoEllipseCalculator::getExtendedNeutrinoEllipse()
{
  //calcNeutrinoEllipse();
  return &Nperp_;
}

TMatrixD* NeutrinoEllipseCalculator::getHomogeneousNeutrinoEllipse()
{
  //calcNeutrinoEllipse();
  return &Hperp_;
}

TVectorD* NeutrinoEllipseCalculator::getNeutrinoMomentum(double theta)
{
  calcExtendedNeutrinoEllipse();
  double tArray[3]={cos(theta),sin(theta),1.};
  TVectorD t(3,tArray);
  nuPerp_=t;
  nuPerp_*=Hperp_;
  pNu_=nuPerp_;
  pNu_*=HperpInv_;
  pNu_*=H_;
  return &pNu_;
}

void NeutrinoEllipseCalculator::calcBJetCorrection()
{

  //Find the ranges of b-jet energies in which Z^2>=0:
  //Z^2=0 is equivalent to a polynomial of degree 4 in Eb
  //bJetBeta_=bJetP/bJetE does not change when the energy correction is applied
  //(Cancels in numerator and denominator)
  //All other parameters are independent of the b-jet energy and momentum

  //cout << "Finding b-jet energy ranges where Z^2>0" << endl;
 
  double a0 = -(pow(mt2_ - mW2_,2))*((leptonE_ - leptonP_)*(leptonE_ + leptonP_)/(4.*(pow(c_*leptonE_*bJetP_ - bJetE_*leptonP_,2) + pow(bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_)*s2_)));
  
  double a1 = ((mt2_ - mW2_))*(((bJetE_*leptonE_ - c_*bJetP_*leptonP_)*(pow(leptonE_,2) - mnu2_ + mW2_ - pow(leptonP_,2)))/
			      (2.*(pow(c_*leptonE_*bJetP_ - bJetE_*leptonP_,2) + pow(bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_)*s2_)));

  double a2 = (8*c_*bJetE_*leptonE_*mW2_*bJetP_*leptonP_ - pow(bJetE_,2)*(pow(leptonE_,4) + pow(mnu2_ - mW2_,2) + 2*(mnu2_ + mt2_)*pow(leptonP_,2) + pow(leptonP_,4) - 
						    2*pow(leptonE_,2)*(mnu2_ + mt2_ - 2*mW2_ + pow(leptonP_,2))) + 
	       pow(bJetP_,2)*(-2*(mt2_ - mW2_)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_) + 
			    c2_*(pow(leptonE_,4) + pow(mnu2_ - mW2_ + pow(leptonP_,2),2) - 2*pow(leptonE_,2)*(mnu2_ + mW2_ + pow(leptonP_,2))) + 
			    (pow(leptonE_,4) + pow(mnu2_ - mW2_,2) + 2*(mnu2_ + mW2_)*pow(leptonP_,2) + pow(leptonP_,4) - 2*pow(leptonE_,2)*(mnu2_ + mW2_ + pow(leptonP_,2)))*
			    s2_))/(4.*(pow(c_*leptonE_*bJetP_ - bJetE_*leptonP_,2) + pow(bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_)*s2_));

  double a3 = -((bJetE_ - bJetP_)*(bJetE_ + bJetP_)*(bJetE_*leptonE_ - c_*bJetP_*leptonP_)*(pow(leptonE_,2) - mnu2_ + mW2_ - pow(leptonP_,2)))/
    (2.*(pow(c_*leptonE_*bJetP_ - bJetE_*leptonP_,2) + pow(bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_)*s2_));

  double a4 = -(pow(bJetE_ - bJetP_,2)*pow(bJetE_ + bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_))/
   (4.*(pow(c_*leptonE_*bJetP_ - bJetE_*leptonP_,2) + pow(bJetP_,2)*(leptonE_ - leptonP_)*(leptonE_ + leptonP_)*s2_));
  
  //Solve for the roots of the polynomial
  //Consider only the real positive roots

  //cout << "top mass: " << mt_ << " and top mass squared: " << mt2_ << endl;
  //cout << "W   mass: " << mW_ << " and W   mass squared: " << mW2_ << endl;
  //
  //cout << "a0 : " << a0 << endl;
  //cout << "a1 : " << a1 << endl;
  //cout << "a2 : " << a2 << endl;
  //cout << "a3 : " << a3 << endl;

  Polynomial quartic(a4,a3,a2,a1,a0);
  vector<complex<double> > realBJetSFs = quartic.FindNumRoots();
  vector<double> bJetSFs;

  //cout << "iterating over possible energy edges" << endl;

  for( vector<complex<double> >::iterator ibJetSF = realBJetSFs.begin(); ibJetSF != realBJetSFs.end(); ibJetSF++ )
    {
      //cout << "This root is: " << ibJetSF->real() << " + " << ibJetSF->imag() << "i" << endl;
      //cout << "This root is: " << *ibJetSF << endl;
      if( ibJetSF->real() <=0 || abs(ibJetSF->imag()) > 0.0001*abs(ibJetSF->real())) 
	{
	  //cout << "we claim it is negative or too complex" << endl;
	  //cout << "This root has a negative value: " << *ibJetE/bJetE_ << ". Discarding it." << endl;
	  continue;
	}
      else
	{
	  //cout << "This root has a positive value!" << endl;
	  if(bJetSFs.size() > 0 && abs(ibJetSF->real() - *(bJetSFs.rbegin())) < 1e-9 ) continue;
	  //cout << "we claim it's unique" << endl;
	  bJetSFs.push_back(ibJetSF->real());
	  double sf(ibJetSF->real());
	  double sf2 = sf*sf;
	  double Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
	  double Z2up = (a0+a1*sf+a2*pow(1.02*sf,2)+a3*pow(1.02*sf,3)+a4*pow(1.02*sf,4))/(1.02*sf*1.02*sf);
	  double Z2down = (a0+a1*sf+a2*pow(0.98*sf,2)+a3*pow(0.98*sf,3)+a4*pow(0.98*sf,4))/(0.98*sf*0.98*sf);
	  //cout << "Testing scale factor " << sf << endl;
	  //cout << "Z2(sf) is " << Z2 << endl;
	  //cout << "two percent up is: " << Z2up << endl;
	  //cout << "two percent down is: " << Z2down << endl;
	}
    }
  sort(bJetSFs.begin(),bJetSFs.end() );
  //cout << "There are " << bJetSFs.size() << " real positive roots" << endl;

  //Store the interval closest to the measured jet energy where Z^2 is positive as a pair of pairs of doubles and bools: 

  double  lowEdge=0; bool  hasLowEdge=false;
  double highEdge=0; bool hasHighEdge=false;
  nRanges_=0;
  if( bJetSFs.size() > 0 ) 
    {
      //cout << "more than one range" << endl;
      vector<bool> isPositiveSection;
      vector<double> edgesLogScale;
      double sf(1.);
      if(bJetSFs.at(0) < 1) sf = 0.5*bJetSFs.at(0);
      double sf2 = sf*sf;
      double Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
      //cout << "Testing scale factor " << sf << endl;
      //cout << "between 0 and " << bJetSFs.at(0) << endl;
      //cout << "Z2 is " << Z2 << endl;
      (Z2 > 0) ? isPositiveSection.push_back(true) : isPositiveSection.push_back(false);
      unsigned int middleBin = 0;
      for( unsigned int i=0; i<bJetSFs.size(); i++ )
	{
	  ( i == bJetSFs.size()-1 ) ? (bJetSFs.at(bJetSFs.size()-1)<1 ? sf = 1 : sf = 1.1*bJetSFs.at(bJetSFs.size()-1)) : (bJetSFs.at(i)<1&&bJetSFs.at(i+1)>1 ? sf = 1 :sf = 0.5*(bJetSFs.at(i)+bJetSFs.at(i+1)));
	  sf2 = sf*sf;
	  Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
	  //cout << "Testing scale factor " << sf << endl;
	  //( i == bJetSFs.size()-1 ) ? cout << "between " << bJetSFs.at(i) << " and inf" << endl : cout << "between " << bJetSFs.at(i) << " and " << bJetSFs.at(i+1) << endl;
	  //cout << "Z2 is " << Z2 << endl;
	  
	  if(Z2 > 0)
	    {
	      //cout << "Z^2(sf) is positive" << endl;
	      if(*(isPositiveSection.rbegin())) continue;
	      isPositiveSection.push_back(true);
	    }
	  else
	    {
	      //cout << "Z^2(sf) is negative" << endl;
	      if(!*(isPositiveSection.rbegin())) continue;
	      isPositiveSection.push_back(false);
	    }
	  double edge = log(bJetSFs.at(i));
	  //cout << "Edge value: " << edge << endl;
	  edgesLogScale.push_back(edge);
	  if(edge < 0 ) middleBin++;
	}
      //cout << "range including zero in bin: " << middleBin << endl;
      if(isPositiveSection[middleBin])
	{
	  if(middleBin != isPositiveSection.size()-1)
	    {
	      hasHighEdge=true;
	      highEdge=edgesLogScale[middleBin];
	      //cout << "found a high edge: " << highEdge << endl;
	    }
	  if(middleBin != 0)
	    {
	      hasLowEdge=true;
	      lowEdge=edgesLogScale[middleBin-1];
	      //cout << "found a low edge: " << lowEdge << endl;
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
		  //cout << "found a high edge: " << highEdge << endl;
		}
	      hasLowEdge=true;
	      lowEdge=edgesLogScale[0];
	      //cout << "found a low edge: " << lowEdge << endl;
	    }
	  else if(middleBin == isPositiveSection.size() - 1)
	    {
	      hasHighEdge=true;
	      highEdge=edgesLogScale[middleBin-1];
	      //cout << "found a high edge: " << highEdge << endl;
	      if(middleBin-1 != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-2];
		  //cout << "found a low edge: " << lowEdge << endl;
		}
	    }
	  else if(abs(edgesLogScale[middleBin-1])>edgesLogScale[middleBin])
	    {
	      middleBin-=1;
	      if(middleBin != isPositiveSection.size()-1)
		{
		  hasHighEdge=true;
		  highEdge=edgesLogScale[middleBin];
		  //cout << "found a high edge: " << highEdge << endl;
		}
	      if(middleBin != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-1];
		  //cout << "found a low edge: " << lowEdge << endl;
		}
	    }
	  else
	    {
	      middleBin+=1;
	      if(middleBin != isPositiveSection.size()-1)
		{
		  hasHighEdge=true;
		  highEdge=edgesLogScale[middleBin];
		  //cout << "found a high edge: " << highEdge << endl;
		}
	      if(middleBin != 0)
		{
		  hasLowEdge=true;
		  lowEdge=edgesLogScale[middleBin-1];
		  //cout << "found a low edge: " << lowEdge << endl;
		}
	    }
	}
      for(vector<bool>::iterator thisRange = isPositiveSection.begin(); thisRange != isPositiveSection.end(); thisRange++)
	{
	  if(*thisRange) nRanges_++;
	} 
    }
  else
    {
      //cout << "only one range" << endl;
      double Z2 = (a0+a1+a2+a3+a4);
      if(Z2 > 0)
	{
	  nRanges_ = 1;
	  hasLowEdge = false;
	  hasHighEdge = false;
	}
      else 
	{
	  nRanges_ = 0;
	  hasLowEdge = true;
	  hasHighEdge = true;
	}
    }

  //double Z2atOne = (a0+a1+a2+a3+a4);
  //cout << "Z2 in range calculation is " << Z2atOne << endl;

  //make sure edges are definitely positive
  lowEdge /=bJetPtWidth_;
  highEdge /=bJetPtWidth_;
  
  if(hasLowEdge)
    {
      double sf(exp(lowEdge*bJetPtWidth_));
      setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
      setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
      Wsurface();
      leptonEllipsoid();
      bJetEllipsoid();
      calcZ2();
      //double sf2 = sf*sf;
      //double Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
      while(Z2_ < 0)
	{
	  if(hasHighEdge) lowEdge+=0.001*(highEdge - lowEdge);
	  else lowEdge+=0.001;
	  sf = exp(lowEdge*bJetPtWidth_);
	  setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
	  setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
	  Wsurface();
	  leptonEllipsoid();
	  bJetEllipsoid();
	  calcZ2();
	  //sf2 = sf*sf;
	  //Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
	  if(hasHighEdge && lowEdge >= highEdge) break;
	}
      //if(hasHighEdge) lowEdge+=0.05*(highEdge - lowEdge);
      //else lowEdge+=0.05;
      //cout << "lower energy edge is " << bJetP_ << " with scaling factor " << sf << " = exp(" << lowEdge << ")" << endl;
      //cout << "Z2 is " << Z2_ << endl;

    }

  if(hasHighEdge)
    {
      double sf(exp(highEdge*bJetPtWidth_));
      setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
      setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
      Wsurface();
      leptonEllipsoid();
      bJetEllipsoid();
      calcZ2();
      //double sf2 = sf*sf;
      //double Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
      while(Z2_<0)
	{
	  if(hasLowEdge) highEdge-=0.001*(highEdge - lowEdge);
	  else highEdge-=0.001;
	  sf = exp(highEdge*bJetPtWidth_);
	  setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
	  setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
	  Wsurface();
	  leptonEllipsoid();
	  bJetEllipsoid();
	  calcZ2();
	  //sf2 = sf*sf;
	  //Z2 = a0/sf2+a1/sf+a2+a3*sf+a4*sf2;
	  if(hasLowEdge && lowEdge >= highEdge) break;
	}
      //if(hasLowEdge) highEdge-=0.05*(highEdge - lowEdge);
      //else highEdge-=0.05;
      //cout << "upper energy edge is " << bJetP_ << " with scaling factor " << sf  << " = exp(" << highEdge << ")" << endl;
      //cout << "Z2 is " << Z2_ << endl;
    }

  if(hasHighEdge || hasLowEdge)
    {
      setBJetFactors();
      setAngles();
      Wsurface();
      leptonEllipsoid();
      bJetEllipsoid();
      calcZ2();
    }

  bJetLogSFRange_.first =  make_pair( lowEdge, hasLowEdge) ;
  bJetLogSFRange_.second = make_pair(highEdge,hasHighEdge) ;

  if(hasHighEdge && hasLowEdge && lowEdge >= highEdge ) nRanges_ = 0;

  //cout << "bJet Range from inside: (";
  //if(!bJetLogSFRange_.first.second)  cout << "-inf,";
  //else cout << bJetLogSFRange_.first.first << ",";
  //if(!bJetLogSFRange_.second.second)  cout << "inf)" << endl;
  //else cout << bJetLogSFRange_.second.first << ")" << endl;
  //
  //cout << "Z2 at edges: (";
  //if(!bJetLogSFRange_.first.second)  cout << "inf,";
  //else{
  //  double sf(exp(bJetLogSFRange_.first.first*bJetPtWidth_));
  //  setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
  //  setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
  //  Wsurface();
  //  leptonEllipsoid();
  //  bJetEllipsoid();
  //  calcZ2();
  //  cout << Z2_ << ",";
  //}
  //if(!bJetLogSFRange_.second.second)  cout << "inf)" << endl;
  //else{
  //  double sf(exp(bJetLogSFRange_.second.first*bJetPtWidth_));
  //  setTempBJetFactors(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_, sf*unscaled_bJetE_);
  //  setTempAngles(sf*unscaled_bJetPx_, sf*unscaled_bJetPy_, sf*unscaled_bJetPz_);
  //  Wsurface();
  //  leptonEllipsoid();
  //  bJetEllipsoid();
  //  calcZ2();
  //  cout << Z2_ << ")" << endl;
  //}

  //cout << "Low edge? " << boolalpha << hasLowEdge << " " << lowEdge << endl;
  //cout << "High edge? " << boolalpha << hasHighEdge << " " << highEdge << endl;

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
