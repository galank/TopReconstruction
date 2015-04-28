#include "leptonicTopSystemChiSquare.h"

leptonicTopSystemChiSquare::leptonicTopSystemChiSquare(const double& bJetPx, const double& bJetPy, const double& bJetPz, const double& bJetE,
                                                       const double& bJetPtWidth, const double& bJetEtaWidth, const double& bJetPhiWidth,
                                                       const double& leptonPx, const double& leptonPy, const double& leptonPz, const double& leptonE,
                                                       const double& leptonPtWidth, const double& leptonEtaWidth, const double& leptonPhiWidth,
						       const double& mTop, const double& sigmaMTop,
                                                       const double& mW, const double& sigmaMW ) :
  topSystemChiSquare   (bJetPx, bJetPy, bJetPz, bJetE,
                        bJetPtWidth, bJetEtaWidth, bJetPhiWidth,
                        leptonPx, leptonPy, leptonPz, leptonE,
                        leptonPtWidth, leptonEtaWidth, leptonPhiWidth,
                        mTop, sigmaMTop,
                        mW, sigmaMW,
			0.
                        ),
  chi2_                 (0.)
{
  //cout << "constructor" << endl;

  setBJetWidths(bJetPtWidth,bJetPhiWidth,bJetEtaWidth);
  setWDaughter1Widths(leptonPtWidth,leptonPhiWidth,leptonEtaWidth);
  setTopMassWidth(sigmaMTop);
  setWMassWidth(sigmaMW);

  //printTopConstituents();
}

leptonicTopSystemChiSquare::~leptonicTopSystemChiSquare()
{
  //cout << "destructor" << endl;
}

void leptonicTopSystemChiSquare::printTopConstituents()
{
  cout << "Current leptonic top decay products:" << endl;
//  cout << "b-jet: "
//       << "\npx = " << bJetPx_
//       << "\npy = " << bJetPy_
//       << "\npz = " << bJetPz_
//       << "\ne  = " << bJetE_
//       << endl;
//
//  cout << "lepton: "
//       << "\npx = " << WDaughter1Px_
//       << "\npy = " << WDaughter1Py_
//       << "\npz = " << WDaughter1Pz_
//       << "\ne  = " << WDaughter1E_
//       << endl;
//
  cout << "neutrino: "
       << "\npx = " << WDaughter2Px_
       << "\npy = " << WDaughter2Py_
       << "\npz = " << WDaughter2Pz_
       << "\ne  = " << WDaughter2E_
       << endl;

  //cout << "neutrino mass squared is " << WDaughter2E_*WDaughter2E_-WDaughter2Px_*WDaughter2Px_-WDaughter2Py_*WDaughter2Py_-WDaughter2Pz_*WDaughter2Pz_ << endl;
}

void leptonicTopSystemChiSquare::setEllipseAngle(double theta)
{
  theta_=theta;
  resetWDaughter2(theta_);
  calcTopMomentum();
}

void leptonicTopSystemChiSquare::getWDaughter2Deltas(double& ptDelta, double& phiDelta, double& etaDelta)
{
  ptDelta =0.;
  phiDelta=0.;
  etaDelta=0.;
}

void leptonicTopSystemChiSquare::calcChiSquare()
{
  chi2_=bJetPtDelta_*bJetPtDelta_
    +bJetPhiDelta_*bJetPhiDelta_
    +bJetEtaDelta_*bJetEtaDelta_
    +WDaughter1PtDelta_*WDaughter1PtDelta_
    +WDaughter1PhiDelta_*WDaughter1PhiDelta_
    +WDaughter1EtaDelta_*WDaughter1EtaDelta_
//    +breitWignerError(mTop_,sigmaMTop_,deltaMTop_)   //moved to inner minimization
    +breitWignerError(mW_,sigmaMW_,deltaMW_);
//    +deltaMW_*deltaMW_;
}

double leptonicTopSystemChiSquare::getChiSquare()
{
  calcChiSquare();
  return chi2_;
}
