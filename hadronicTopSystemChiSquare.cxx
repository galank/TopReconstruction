#include "hadronicTopSystemChiSquare.h"

hadronicTopSystemChiSquare::hadronicTopSystemChiSquare(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
						       const double& bJetPtWidth, const double& bJetEtaWidth, const double& bJetPhiWidth,
						       const double& measJet1Px, const double& measJet1Py, const double& measJet1Pz, const double& measJet1E,
						       const double& measJet1PtWidth, const double& measJet1EtaWidth, const double& measJet1PhiWidth,
						       const double& measJet2Px, const double& measJet2Py, const double& measJet2Pz, const double& measJet2E,
						       const double& measJet2PtWidth, const double& measJet2EtaWidth, const double& measJet2PhiWidth,
						       const double& mTop, const double& sigmaMTop,
						       const double& mW, const double& sigmaMW ) :
  topSystemChiSquare         (bJetPx, bJetPy, bJetPz, bJetE,
			      bJetPtWidth, bJetEtaWidth, bJetPhiWidth,
			      measJet1Px, measJet1Py, measJet1Pz, measJet1E,
			      measJet1PtWidth, measJet1EtaWidth, measJet1PhiWidth,
			      mTop, sigmaMTop,
			      mW, sigmaMW,
			      sqrt(max(0.,measJet2E*measJet2E-measJet2Px*measJet2Px-measJet2Py*measJet2Py-measJet2Pz*measJet2Pz))
			      ),
  WDaughter2PtWidth_         (measJet2PtWidth),
  WDaughter2PhiWidth_        (measJet2PhiWidth),
  WDaughter2EtaWidth_        (measJet2EtaWidth),
  WDaughter2PtDelta_         (0.),
  WDaughter2PhiDelta_        (0.),
  WDaughter2EtaDelta_        (0.),
  chi2_                      (0.)
{
  //cout << "constructor" << endl;

  WDaughter2Px_=measJet2Px;
  WDaughter2Py_=measJet2Py;
  WDaughter2Pz_=measJet2Pz;
  WDaughter2E_ =measJet2E ;
  WDaughter2Pt_=sqrt(WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_);
  WDaughter2Phi_=atan2(WDaughter2Py_,WDaughter2Px_);
  WDaughter2Eta_=asinh(WDaughter2Pz_/sqrt(WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_));

  reconstructed_WDaughter2Px_=WDaughter2Px_;
  reconstructed_WDaughter2Py_=WDaughter2Py_;
  reconstructed_WDaughter2Pz_=WDaughter2Pz_;
  reconstructed_WDaughter2E_ =WDaughter2E_ ;
  reconstructed_WDaughter2Pt_=sqrt(WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_);
  reconstructed_WDaughter2Phi_=atan2(WDaughter2Py_,WDaughter2Px_);
  reconstructed_WDaughter2Eta_=asinh(WDaughter2Pz_/sqrt(WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_));
  reconstructed_WDaughter2Mass2_=max(WDaughter2E_*WDaughter2E_-WDaughter2Px_*WDaughter2Px_-WDaughter2Py_*WDaughter2Py_-WDaughter2Pz_*WDaughter2Pz_,0.);

  setBJetWidths(bJetPtWidth,bJetPhiWidth,bJetEtaWidth);
  setWDaughter1Widths(measJet1PtWidth,measJet1PhiWidth,measJet1EtaWidth);
  setTopMassWidth(sigmaMTop);
  setWMassWidth(sigmaMW);

  //printTopConstituents();
}

hadronicTopSystemChiSquare::~hadronicTopSystemChiSquare()
{
  //cout << "destructor" << endl;
}

void hadronicTopSystemChiSquare::printTopConstituents()
{
  cout << "Current hadronic top decay products:" << endl;
//  cout << "In the derived class:" << endl;
//  cout << "b-jet: "
//       << "\npx = " << bJetPx_
//       << "\npy = " << bJetPy_
//       << "\npz = " << bJetPz_
//       << "\ne  = " << bJetE_  
//       << endl;
//
//  cout << "first light quark: " 
//       << "\npx = " << WDaughter1Px_
//       << "\npy = " << WDaughter1Py_
//       << "\npz = " << WDaughter1Pz_
//       << "\ne  = " << WDaughter1E_
//       << endl;

  cout << "second light quark: "
       << "\npx = " << WDaughter2Px_
       << "\npy = " << WDaughter2Py_
       << "\npz = " << WDaughter2Pz_
       << "\ne  = " << WDaughter2E_
       << endl;

}

void hadronicTopSystemChiSquare::setEllipseAngle(double theta)
{
  theta_=theta;
  resetWDaughter2(theta_);
  calcWDaughter2Deltas();
  calcTopMomentum();
}

void hadronicTopSystemChiSquare::setWDaughter2(double px, double py, double pz)
{
  WDaughter2Px_=px;
  WDaughter2Py_=py;
  WDaughter2Pz_=pz;
  WDaughter2E_=sqrt(reconstructed_WDaughter2Mass2_+WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_+WDaughter2Pz_*WDaughter2Pz_);

  WDaughter2Pt_=sqrt(WDaughter2Px_*WDaughter2Px_+WDaughter2Py_*WDaughter2Py_);
  WDaughter2Phi_=atan2(WDaughter2Py_,WDaughter2Px_);
  WDaughter2Eta_=asinh(WDaughter2Pz_/WDaughter2Pt_);
}

//FIXME is this function necessary?
//typically, set ellipse angle -> calculate momentum -> recalculate deltas
//void hadronicTopSystemChiSquare::resetWDaughter2()
//{
//  double phiNew = reconstructed_WDaughter2Phi_+WDaughter2PhiDelta_*WDaughter2PhiWidth_;
//  double etaNew = reconstructed_WDaughter2Eta_+WDaughter2EtaDelta_*WDaughter2EtaWidth_;
//  double ptNew  = reconstructed_WDaughter2Pt_*exp(WDaughter2PtDelta_*WDaughter2PtWidth_) ;
//  double pxNew = ptNew*cos(phiNew);
//  double pyNew = ptNew*sin(phiNew);
//  double pzNew = ptNew*sinh(etaNew);
//  double pNew  = ptNew*cosh(etaNew);
//  double eNew  = sqrt(reconstructed_WDaughter2Mass2_ + pNew*pNew );
//
//  //reset values
//  WDaughter2Px_=pxNew;
//  WDaughter2Py_=pyNew;
//  WDaughter2Pz_=pzNew;
//  WDaughter2E_ =eNew ;
//
//  WDaughter2Pt_ =ptNew;
//  WDaughter2Phi_=phiNew;
//  WDaughter2Eta_=etaNew;
//
////  cout << "Second W daughter momentum:"
////       << "\npx = " << WDaughter2Px_
////       << "\npy = " << WDaughter2Py_
////       << "\npz = " << WDaughter2Pz_
////       << endl;
//}

void hadronicTopSystemChiSquare::getWDaughter2Deltas(double& ptDelta, double& phiDelta, double& etaDelta)
{
  //cout << "Current second W daughter pt  delta is " << ptDelta  << endl;
  //cout << "Current second W daughter phi delta is " << phiDelta << endl;
  //cout << "Current second W daughter eta delta is " << etaDelta << endl;
  ptDelta =WDaughter2PtDelta_;
  phiDelta=WDaughter2PhiDelta_;
  etaDelta=WDaughter2EtaDelta_;
}

void hadronicTopSystemChiSquare::setWDaughter2Deltas(double ptDelta, double phiDelta, double etaDelta)
{
//  cout << "Setting second W daughter deltas to" 
//       << "\ndelta pt  = " << ptDelta
//       << "\ndelta phi = " << phiDelta
//       << "\ndelta eta = " << etaDelta
//       << endl;
//  cout << "Second W daughter momentum (cartesian coordinates):"
//       << "\npx = " << WDaughter2Px_
//       << "\npy = " << WDaughter2Py_
//       << "\npz = " << WDaughter2Pz_
//       << endl;
//  cout << "Second W daughter momentum (polar coordinates)"
//       << "\npt  = " << WDaughter2Pt_
//       << "\nphi = " <<WDaughter2Phi_
//       << "\neta = " <<WDaughter2Eta_
//       << endl;
//  cout << "reconstructed second light quark: "
//       << "\npx = " << reconstructed_WDaughter2Px_
//       << "\npy = " << reconstructed_WDaughter2Py_
//       << "\npz = " << reconstructed_WDaughter2Pz_
//       << "\ne  = " << reconstructed_WDaughter2E_
//       << endl;

  WDaughter2PtDelta_ =ptDelta ;
  WDaughter2PhiDelta_=phiDelta;
  WDaughter2EtaDelta_=etaDelta;
}

void hadronicTopSystemChiSquare::calcWDaughter2Deltas()
{
//  cout << "Calculating second W daughter deltas" << endl;
//  cout << "Second W daughter momentum (cartesian coordinates):"
//       << "\npx = " << WDaughter2Px_
//       << "\npy = " << WDaughter2Py_
//       << "\npz = " << WDaughter2Pz_
//       << endl;
//  cout << "Second W daughter momentum (polar coordinates)"
//       << "\npt  = " << WDaughter2Pt_
//       << "\nphi = " << WDaughter2Phi_
//       << "\neta = " << WDaughter2Eta_
//       << endl;
//  cout << "Reconstructed second W daughter momentum (polar coordinates)"
//       << "\npt  = " << reconstructed_WDaughter2Pt_  << " with width " << WDaughter2PtWidth_
//       << "\nphi = " << reconstructed_WDaughter2Phi_ << " with width " << WDaughter2PhiWidth_
//       << "\neta = " << reconstructed_WDaughter2Eta_ << " with width " << WDaughter2EtaWidth_
//       << endl;

  //double deltaPt =abs(log(WDaughter2Pt_/reconstructed_WDaughter2Pt_)); //why is there an absolute value?
  double deltaPt =log(WDaughter2Pt_/reconstructed_WDaughter2Pt_);
  double deltaEta=WDaughter2Eta_-reconstructed_WDaughter2Eta_;         
  double deltaPhi=WDaughter2Phi_-reconstructed_WDaughter2Phi_;
  while(deltaPhi> 3.14159265359)
    {
      deltaPhi -= 2.*3.14159265359;
    }
  while(deltaPhi< -3.14159265359)
    {
      deltaPhi += 2.*3.14159265359;
    }
  //deltaPhi /= WDaughter2PhiWidth_;
  setWDaughter2Deltas(deltaPt,deltaPhi,deltaEta);
}

void hadronicTopSystemChiSquare::printWDaughter2()
{
  cout << "Current second light quark: "
       << "\npt  = " << WDaughter2Pt_
       << "\nphi = " << WDaughter2Phi_
       << "\neta = " << WDaughter2Eta_
       << endl;

  cout << "theta is " << theta_ << endl;

//  cout << "reconstructed second light quark: " 
//       << "\npt  = " << reconstructed_WDaughter2Pt_ 
//       << "\nphi = " << reconstructed_WDaughter2Phi_ 
//       << "\neta = " << reconstructed_WDaughter2Eta_ 
//       << endl;
//
  cout <<   "delta pt  = " << WDaughter2PtDelta_
       << "\ndelta phi = " << WDaughter2PhiDelta_
       << "\ndelta eta = " << WDaughter2EtaDelta_  
       << endl;

  //cout << "second W daughter ellipse in homogeneous coordinates:" << endl;
  //TMatrixD* a=getHomogeneousWDaughterEllipse();
  //a->Print();
}

void hadronicTopSystemChiSquare::calcChiSquare()
{
  //calcWDaughter2Deltas();
  chi2_=bJetPtDelta_*bJetPtDelta_
    +bJetPhiDelta_*bJetPhiDelta_
    +bJetEtaDelta_*bJetEtaDelta_
    +WDaughter1PtDelta_*WDaughter1PtDelta_
    +WDaughter1PhiDelta_*WDaughter1PhiDelta_
    +WDaughter1EtaDelta_*WDaughter1EtaDelta_
//    +WDaughter2PtDelta_*WDaughter2PtDelta_           //already counted
//    +WDaughter2PhiDelta_*WDaughter2PhiDelta_	       //in hadronic chi
//    +WDaughter2EtaDelta_*WDaughter2EtaDelta_         //square 
//    +breitWignerError(mTop_,sigmaMTop_,deltaMTop_)   //moved to inner minimization
    +breitWignerError(mW_,sigmaMW_,deltaMW_);
//    +deltaMW_*deltaMW_;
}

double hadronicTopSystemChiSquare::getChiSquare()
{
  calcChiSquare();
  return chi2_;
}

double hadronicTopSystemChiSquare::getHadronicChiSquare()
{
  //calcWDaughter2Deltas();
  //cout << "Current second W daughter pt delta is " << WDaughter2PtDelta_ << endl;
  //cout << "Current second W daughter phi delta is " << WDaughter2PhiDelta_ << endl;
  //cout << "Current second W daughter eta delta is " << WDaughter2EtaDelta_ << endl;
  double hadChi2=pow(WDaughter2PtDelta_,2)+pow(WDaughter2PhiDelta_,2)+pow(WDaughter2EtaDelta_,2);
  //cout << "Hadronic top decay: second W daughter chi square is " << hadChi2 << endl;
  return hadChi2;
}
