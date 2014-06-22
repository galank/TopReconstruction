#include "topSystemChiSquare.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1D.h"
#include "Math/QuantFuncMathCore.h"

class jetCompare
{
public:
  bool operator()(XYZTLorentzVector j1,XYZTLorentzVector j2){return j1.Pt()>j2.Pt();}
};


topSystemChiSquare::topSystemChiSquare(vector<XYZTLorentzVector> jets,
				       vector<double> jetPtWidths,
				       vector<double> jetPhiWidths,
				       vector<double> jetEtaWidths,
				       int bJet1, int lightJet1, XYZTLorentzVector lepton1, 
				       int bJet2, int lightJet2, XYZTLorentzVector lepton2,
				       XYZTLorentzVector METLorentzVector,
				       double mTop, double sigmaMTop,
				       double mW, double sigmaMW) :
  bJet1_    (bJet1),
  lightJet1_(lightJet1),
  bJet2_    (bJet2),
  lightJet2_(lightJet2),
  mTop_     (mTop),
  mW_       (mW),
  mNu_      (0.),
  sigmaMTop_(sigmaMTop),
  sigmaMW_  (sigmaMW),
  measuredMET_ (METLorentzVector),
  MET_         (METLorentzVector),
  measuredMETx_(METLorentzVector.Px()),
  METx_        (METLorentzVector.Px()),
  measuredMETy_(METLorentzVector.Py()),
  METy_        (METLorentzVector.Py()),
  lepton1LorentzVector_(lepton1),
  lepton2LorentzVector_(lepton2),
  jets_          (jets),
  //lightJets_     (jets),
  jetPtWidths_   (jetPtWidths),
  jetPhiWidths_  (jetPhiWidths),
  jetEtaWidths_  (jetEtaWidths),
//  lightJetPtWidths_ (jetPtWidths),
//  lightJetPhiWidths_(jetPhiWidths),
  nuEllipseOneOneCalc_(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,
		       bJet1PtWidth_,bJet1PhiWidth_,bJet1EtaWidth_,
		       lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
		       angleAdjusted_bJet1Px_, angleAdjusted_bJet1Py_, angleAdjusted_bJet1Pz_,
		       angleAdjusted_bJet1P_, angleAdjusted_bJet1E_),
  nuEllipseOneTwoCalc_(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,
		       bJet1PtWidth_,bJet1PhiWidth_,bJet1EtaWidth_,
		       lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
		       angleAdjusted_bJet1Px_, angleAdjusted_bJet1Py_, angleAdjusted_bJet1Pz_,
		       angleAdjusted_bJet1P_, angleAdjusted_bJet1E_),
  nuEllipseTwoOneCalc_(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,
		       bJet2PtWidth_,bJet2PhiWidth_,bJet2EtaWidth_,
		       lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
		       angleAdjusted_bJet2Px_, angleAdjusted_bJet2Py_, angleAdjusted_bJet2Pz_,
		       angleAdjusted_bJet2P_, angleAdjusted_bJet2E_),
  nuEllipseTwoTwoCalc_(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,
		       bJet2PtWidth_,bJet2PhiWidth_,bJet2EtaWidth_,
		       lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
		       angleAdjusted_bJet2Px_, angleAdjusted_bJet2Py_, angleAdjusted_bJet2Pz_,
		       angleAdjusted_bJet2P_, angleAdjusted_bJet2E_),
  lightJetChi2_(0.),
  bJetChi2_    (0.),
  chi2_        (0.),
  dx_          (0.),
  dy_          (0.),  
  lightJetChiSquare_(jets.size()-2,dx_,dy_),
  maxConsideredChiSquareRoot_(30.),
  intersectingEllipses_(false),
  intersectingEllipsesChi2_(-2)
{
  cout << "Setting up the top system" << endl;
  setBJets();
  setLightJets();
  setLightJetCollections();
  setLeptons();
  setupNeutrinoSolutions();
  setupNeutrinoEllipses();
  lightJetChiSquare_.setupEquations(lightJetPts_,lightJetPtWidths_,lightJetPhis_,lightJetPhiWidths_);
  setBJetAngleDeltas(0.,0.,0.,0.);
  setBJetPtDeltas(0.,0.);
}

topSystemChiSquare::topSystemChiSquare(const topSystemChiSquare& other) :
  bJet1_    (other.bJet1_),
  lightJet1_(other.lightJet1_),
  bJet2_    (other.bJet2_),
  lightJet2_(other.lightJet2_),
  mTop_     (other.mTop_),
  mW_       (other.mW_),
  mNu_      (other.mNu_),
  sigmaMTop_(other.sigmaMTop_),
  sigmaMW_  (other.sigmaMW_),
  measuredMET_ (other.measuredMET_),
  MET_         (other.MET_),
  measuredMETx_(other.measuredMETx_),
  METx_        (other.METx_),
  measuredMETy_(other.measuredMETy_),
  METy_        (other.METy_),
  lepton1LorentzVector_(other.lepton1LorentzVector_),
  lepton2LorentzVector_(other.lepton2LorentzVector_),
  jets_          (other.jets_),
  lightJets_     (other.jets_),
  jetPtWidths_   (other.jetPtWidths_),
  jetPhiWidths_  (other.jetPhiWidths_),
  jetEtaWidths_  (other.jetEtaWidths_),
  lightJetPtWidths_ (other.jetPtWidths_),
  lightJetPhiWidths_(other.jetPhiWidths_),
  nuEllipseOneOneCalc_(other.bJet1Px_,other.bJet1Py_,other.bJet1Pz_,other.bJet1E_,
		       other.bJet1PtWidth_,other.bJet1PhiWidth_,other.bJet1EtaWidth_,
		       other.lepton1Px_,other.lepton1Py_,other.lepton1Pz_,other.lepton1E_,
		       other.angleAdjusted_bJet1Px_, other.angleAdjusted_bJet1Py_, other.angleAdjusted_bJet1Pz_,
		       other.angleAdjusted_bJet1P_, other.angleAdjusted_bJet1E_),
  nuEllipseOneTwoCalc_(other.bJet1Px_,other.bJet1Py_,other.bJet1Pz_,other.bJet1E_,
		       other.bJet1PtWidth_,other.bJet1PhiWidth_,other.bJet1EtaWidth_,
		       other.lepton2Px_,other.lepton2Py_,other.lepton2Pz_,other.lepton2E_,
		       other.angleAdjusted_bJet1Px_, other.angleAdjusted_bJet1Py_, other.angleAdjusted_bJet1Pz_,
		       other.angleAdjusted_bJet1P_, other.angleAdjusted_bJet1E_),
  nuEllipseTwoOneCalc_(other.bJet2Px_,other.bJet2Py_,other.bJet2Pz_,other.bJet2E_,
		       other.bJet2PtWidth_,other.bJet2PhiWidth_,other.bJet2EtaWidth_,
		       other.lepton1Px_,other.lepton1Py_,other.lepton1Pz_,other.lepton1E_,
		       other.angleAdjusted_bJet2Px_, other.angleAdjusted_bJet2Py_, other.angleAdjusted_bJet2Pz_,
		       other.angleAdjusted_bJet2P_, other.angleAdjusted_bJet2E_),
  nuEllipseTwoTwoCalc_(other.bJet2Px_,other.bJet2Py_,other.bJet2Pz_,other.bJet2E_,
		       other.bJet2PtWidth_,other.bJet2PhiWidth_,other.bJet2EtaWidth_,
		       other.lepton2Px_,other.lepton2Py_,other.lepton2Pz_,other.lepton2E_,
		       other.angleAdjusted_bJet2Px_, other.angleAdjusted_bJet2Py_, other.angleAdjusted_bJet2Pz_,
		       other.angleAdjusted_bJet2P_, other.angleAdjusted_bJet2E_),
  lightJetChi2_(other.lightJetChi2_),
  bJetChi2_(other.bJetChi2_),
  chi2_(other.chi2_),
  dx_          (other.dx_),
  dy_          (other.dy_),
  lightJetChiSquare_(other.jets_.size()-2,dx_,dy_),
  maxConsideredChiSquareRoot_(other.maxConsideredChiSquareRoot_),
  intersectingEllipses_(other.intersectingEllipses_),
  intersectingEllipsesChi2_(other.intersectingEllipsesChi2_)
{
  setBJets();
  setLightJets();
  setLeptons();
  setupNeutrinoSolutions();
  setupNeutrinoEllipses();
  lightJetChiSquare_.setupEquations(lightJetPts_,lightJetPtWidths_,lightJetPhis_,lightJetPhiWidths_);
  setBJetAngleDeltas(0.,0.,0.,0.);
  setBJetPtDeltas(0.,0.);
}

topSystemChiSquare::~topSystemChiSquare()
{
}

void topSystemChiSquare::setBJets()
{
  //cout << "b-jet 1 has index " << bJet1_ << endl;

  bJet1LorentzVector_=jets_.at(bJet1_);
  bJet1PtWidth_ =log(1+jetPtWidths_.at(bJet1_)) ;
  bJet1PhiWidth_=jetPhiWidths_.at(bJet1_);
  bJet1EtaWidth_=jetEtaWidths_.at(bJet1_);

  bJet1Px_=bJet1LorentzVector_.Px();
  bJet1Py_=bJet1LorentzVector_.Py();
  bJet1Pz_=bJet1LorentzVector_.Pz();
  bJet1E_ =bJet1LorentzVector_.E();

  angleAdjusted_bJet1Px_=bJet1LorentzVector_.Px();
  angleAdjusted_bJet1Py_=bJet1LorentzVector_.Py();
  angleAdjusted_bJet1Pz_=bJet1LorentzVector_.Pz();
  angleAdjusted_bJet1P_ =bJet1LorentzVector_.P();
  angleAdjusted_bJet1E_=bJet1LorentzVector_.E();

  //cout << "bJet1 Unscaled Information:\n"
  //     << "Px: " << angleAdjusted_bJet1Px_
  //     << "\nPy: " << angleAdjusted_bJet1Py_
  //     << "\nPz: " << angleAdjusted_bJet1Pz_
  //     << "\nP : " << sqrt(angleAdjusted_bJet1Px_*angleAdjusted_bJet1Px_ + angleAdjusted_bJet1Py_*angleAdjusted_bJet1Py_ + angleAdjusted_bJet1Pz_*angleAdjusted_bJet1Pz_)
  //     << "\nE : " << angleAdjusted_bJet1E_ 
  //     << "\nM : " << sqrt(bJet1LorentzVector_.M2()) << endl;  


  reconstructed_bJet1Pt_ =bJet1LorentzVector_.Pt();
  reconstructed_bJet1Phi_=bJet1LorentzVector_.Phi();
  reconstructed_bJet1Eta_=bJet1LorentzVector_.Eta();
  reconstructed_bJet1Mass2_ = bJet1LorentzVector_.M2();

  //cout << "b-jet 2 has index " << bJet2_ << endl;

  bJet2LorentzVector_=jets_.at(bJet2_);
  bJet2PtWidth_ =log(1+jetPtWidths_.at(bJet2_)) ;
  bJet2PhiWidth_=jetPhiWidths_.at(bJet2_);
  bJet2EtaWidth_=jetEtaWidths_.at(bJet2_);

  bJet2Px_=bJet2LorentzVector_.Px();
  bJet2Py_=bJet2LorentzVector_.Py();
  bJet2Pz_=bJet2LorentzVector_.Pz();
  bJet2E_ =bJet2LorentzVector_.E();

  angleAdjusted_bJet2Px_=bJet2LorentzVector_.Px();
  angleAdjusted_bJet2Py_=bJet2LorentzVector_.Py();
  angleAdjusted_bJet2Pz_=bJet2LorentzVector_.Pz();
  angleAdjusted_bJet2P_ =bJet2LorentzVector_.P();
  angleAdjusted_bJet2E_=bJet2LorentzVector_.E();

  //cout << "bJet2 Unscaled Information:\n"
  //     << "Px: " << angleAdjusted_bJet2Px_
  //     << "\nPy: " << angleAdjusted_bJet2Py_
  //     << "\nPz: " << angleAdjusted_bJet2Pz_
  //     << "\nP : " << sqrt(angleAdjusted_bJet2Px_*angleAdjusted_bJet2Px_ + angleAdjusted_bJet2Py_*angleAdjusted_bJet2Py_ + angleAdjusted_bJet2Pz_*angleAdjusted_bJet2Pz_)
  //     << "\nE : " << angleAdjusted_bJet2E_
  //     << "\nM : " << sqrt(bJet2LorentzVector_.M2()) << endl;  

  reconstructed_bJet2Pt_ =bJet2LorentzVector_.Pt();
  reconstructed_bJet2Phi_=bJet2LorentzVector_.Phi();
  reconstructed_bJet2Eta_=bJet2LorentzVector_.Eta();
  reconstructed_bJet2Mass2_ = bJet2LorentzVector_.M2();
}

void topSystemChiSquare::setLightJets()
{
  //cout << "light jet 1 has index " << lightJet1_ << endl;
  lightJet1LorentzVector_=jets_.at(lightJet1_);
  //cout << "light jet 2 has index " << lightJet2_ << endl;
  lightJet2LorentzVector_=jets_.at(lightJet2_);
}

void topSystemChiSquare::setLightJetCollections()
{
  for(int iJet=0; iJet<int(jets_.size()); iJet++)
    {
      if( iJet==bJet1_ || iJet==bJet2_ )
	{
	  continue;
	}
      lightJets_.push_back(jets_.at(iJet));
      lightJetPts_ .push_back(jets_.at(iJet).Pt() );
      lightJetPhis_.push_back(jets_.at(iJet).Phi());
      lightJetPtWidths_ .push_back(jetPtWidths_.at(iJet) );
      lightJetPhiWidths_.push_back(jetPhiWidths_.at(iJet));
    }
}

void topSystemChiSquare::buildBestLightJets()
{
  vector<double>* minJetDeltasX = lightJetChiSquare_.getMinDeltasX();
  vector<double>* minJetDeltasY = lightJetChiSquare_.getMinDeltasY();
  if(minJetDeltasX->size() != lightJets_.size())
    { 
      cout << "bad size in light jet modification check!" << endl;
      chi2_ = -1.;
      bJetChi2_ = -1.;
      lightJetChi2_ = -1.;
      return;
    }
  vector<double>::iterator thisDeltaX = minJetDeltasX->begin();
  vector<double>::iterator thisDeltaY = minJetDeltasY->begin();
  lightJetsBest_.clear();
  int i(0.);
  for(vector<XYZTLorentzVector>::iterator thisJet=lightJets_.begin(); thisJet != lightJets_.end(); thisJet++, thisDeltaX++, thisDeltaY++)
    {
      i++;
      if((dxBest_ == 0. && dyBest_ == 0.) || lightJetChi2_ == 0.)
	{
	  lightJetsBest_.push_back(*thisJet);
	  //cout << "light jet " << i << " deltaPt /sigmaPt  = 0" << endl;
	  //cout << "light jet " << i << " deltaPhi/sigmaPhi = 0" << endl;
	}
      else
	{
	  double newPx = thisJet->Px()+*thisDeltaX;
	  double newPy = thisJet->Py()+*thisDeltaX;
	  double newPt = sqrt(newPx*newPx+newPy*newPy);
	  double ptScaleFactor = newPt/thisJet->Pt();

	  //cout << "light jet " << i << " Px       = " << thisJet->Px()  << endl;
	  //cout << "light jet " << i << " delta Px = " << *thisDeltaX    << endl;
	  //cout << "light jet " << i << " Py       = " << thisJet->Py()  << endl;
	  //cout << "light jet " << i << " delta Py = " << *thisDeltaY    << endl;

	  XYZTLorentzVector newJet(newPx,newPy,thisJet->Pz(),ptScaleFactor*thisJet->E()); 
	  lightJetsBest_.push_back(newJet);
	  //double deltaPt = abs(log(ptScaleFactor));
	  //cout << "old jet phi: " << thisJet->Phi() << " and new jet phi: " << newJet.Phi() << endl;
	  //double deltaPhi = thisJet->Phi() - newJet.Phi();
	  //while(deltaPhi> 3.14159265359)
	  //  {
	  //    deltaPhi-=2*3.14159265359;
	  //  }
	  //while(deltaPhi< -3.14159265359)
	  //  {
	  //    deltaPhi+=2*3.14159265359;
	  //  }
	  //cout << "delta phi is " << deltaPhi << endl;
	  //cout << "delta pt  is " << deltaPt  << endl;
	  //cout << "sigma phi is " << lightJetPhiWidths_[i-1] << endl;
	  //cout << "sigma pt  is " << lightJetPtWidths_[i-1] << endl;
	  //cout << "light jet " << i << " deltaPt /sigmaPt  = " << deltaPt /lightJetPtWidths_[i-1]  << endl;
	  //cout << "light jet " << i << " deltaPhi/sigmaPhi = " << deltaPhi/lightJetPhiWidths_[i-1] << endl;
	}
    }
  jetCompare compareJets;
  sort(lightJetsBest_.begin(),lightJetsBest_.end(),compareJets);
}

void topSystemChiSquare::setLeptons()
{
  lepton1Px_=lepton1LorentzVector_.Px();
  lepton1Py_=lepton1LorentzVector_.Py();
  lepton1Pz_=lepton1LorentzVector_.Pz();
  lepton1E_ =lepton1LorentzVector_.E() ;
  lepton2Px_=lepton2LorentzVector_.Px();
  lepton2Py_=lepton2LorentzVector_.Py();
  lepton2Pz_=lepton2LorentzVector_.Pz();
  lepton2E_ =lepton2LorentzVector_.E() ;
}

void topSystemChiSquare::resetAll()
{
  resetAngles();
  resetPt();
}

void topSystemChiSquare::resetAngles()
{
  double bJet1Phi = reconstructed_bJet1Phi_+bJet1PhiDelta_*bJet1PhiWidth_;
  double bJet1Eta = reconstructed_bJet1Eta_+bJet1EtaDelta_*bJet1EtaWidth_;
  double bJet2Phi = reconstructed_bJet2Phi_+bJet2PhiDelta_*bJet2PhiWidth_;
  double bJet2Eta = reconstructed_bJet2Eta_+bJet2EtaDelta_*bJet2EtaWidth_;
  angleAdjusted_bJet1Px_ = reconstructed_bJet1Pt_*cos(bJet1Phi);
  angleAdjusted_bJet1Py_ = reconstructed_bJet1Pt_*sin(bJet1Phi);
  angleAdjusted_bJet1Pz_ = reconstructed_bJet1Pt_*sinh(bJet1Eta);
  angleAdjusted_bJet1P_  = reconstructed_bJet1Pt_*cosh(bJet1Eta);
  angleAdjusted_bJet2Px_ = reconstructed_bJet2Pt_*cos(bJet2Phi);
  angleAdjusted_bJet2Py_ = reconstructed_bJet2Pt_*sin(bJet2Phi);
  angleAdjusted_bJet2Pz_ = reconstructed_bJet2Pt_*sinh(bJet2Eta);
  angleAdjusted_bJet2P_  = reconstructed_bJet2Pt_*cosh(bJet2Eta);
  angleAdjusted_bJet1E_  = sqrt(reconstructed_bJet1Mass2_ 
				+ angleAdjusted_bJet1P_*angleAdjusted_bJet1P_);
  angleAdjusted_bJet2E_  = sqrt(reconstructed_bJet2Mass2_ 
				+ angleAdjusted_bJet2P_*angleAdjusted_bJet2P_);
}

void topSystemChiSquare::resetPt()
{
  double bJet1PtScaleFactor = exp(bJet1PtDelta_*bJet1PtWidth_) ;
  double bJet2PtScaleFactor = exp(bJet2PtDelta_*bJet2PtWidth_) ;
  double bJet1PxNew = bJet1PtScaleFactor*angleAdjusted_bJet1Px_;
  double bJet1PyNew = bJet1PtScaleFactor*angleAdjusted_bJet1Py_;
  double bJet1PzNew = bJet1PtScaleFactor*angleAdjusted_bJet1Pz_;
  double bJet2PxNew = bJet2PtScaleFactor*angleAdjusted_bJet2Px_;
  double bJet2PyNew = bJet2PtScaleFactor*angleAdjusted_bJet2Py_;
  double bJet2PzNew = bJet2PtScaleFactor*angleAdjusted_bJet2Pz_;
  double bJet1ENew  = bJet1PtScaleFactor*angleAdjusted_bJet1E_ ;
  double bJet2ENew  = bJet2PtScaleFactor*angleAdjusted_bJet2E_ ;

  //recalculate MET
  METx_=measuredMETx_+bJet1LorentzVector_.Px()-bJet1PxNew+bJet2LorentzVector_.Px()-bJet2PxNew;
  METy_=measuredMETy_+bJet1LorentzVector_.Py()-bJet1PyNew+bJet2LorentzVector_.Py()-bJet2PyNew;

  //reset values
  bJet1Px_=bJet1PxNew;
  bJet1Py_=bJet1PyNew;
  bJet1Pz_=bJet1PzNew;
  bJet1E_ =bJet1ENew;
  bJet2Px_=bJet2PxNew;
  bJet2Py_=bJet2PyNew;
  bJet2Pz_=bJet2PzNew;
  bJet2E_ =bJet2ENew;

  //if(bJet1Pz_ > 1e4 || bJet2Pz_ > 1e4)
  //  {
  //    cout << "bJet1PtDelta_ : " << bJet1PtDelta_  << endl;
  //    cout << "bJet1PhiDelta_: " << bJet1PhiDelta_ << endl;
  //    cout << "bJet1EtaDelta_: " << bJet1EtaDelta_ << endl;
  //    cout << "bJet2PtDelta_ : " << bJet2PtDelta_  << endl;
  //    cout << "bJet2PhiDelta_: " << bJet2PhiDelta_ << endl;
  //    cout << "bJet2EtaDelta_: " << bJet2EtaDelta_ << endl;
  //  }

  //cout << "METx_    = " << METx_    << endl ;
  //cout << "METy_    = " << METy_    << endl ;
  //cout << "bJet1Px_ = " << bJet1Px_ << endl ;
  //cout << "bJet1Py_ = " << bJet1Py_ << endl ;
  //cout << "bJet1Pz_ = " << bJet1Pz_ << endl ;
  //cout << "bJet1E_  = " << bJet1E_  << endl ;
  //cout << "bJet2Px_ = " << bJet2Px_ << endl ;
  //cout << "bJet2Py_ = " << bJet2Py_ << endl ;
  //cout << "bJet2Pz_ = " << bJet2Pz_ << endl ;
  //cout << "bJet2E_  = " << bJet2E_  << endl ;
  //cout << "lepton1Px_ = " << lepton1Px_ << endl ;
  //cout << "lepton1Py_ = " << lepton1Py_ << endl ;
  //cout << "lepton1Pz_ = " << lepton1Pz_ << endl ;
  //cout << "lepton1E_  = " << lepton1E_  << endl ;
  //cout << "lepton2Px_ = " << lepton2Px_ << endl ;
  //cout << "lepton2Py_ = " << lepton2Py_ << endl ;
  //cout << "lepton2Pz_ = " << lepton2Pz_ << endl ;
  //cout << "lepton2E_  = " << lepton2E_  << endl ;

}

void topSystemChiSquare::setupNeutrinoSolutions()
{
  nuSolOne_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
			      bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
			      METx_,METy_,
			      mW_+sigmaMW_*deltaMW1_,mW_+sigmaMW_*deltaMW2_,
			      mTop_+sigmaMTop_*deltaMTop1_,mTop_+sigmaMTop_*deltaMTop2_);
  
  nuSolTwo_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
			      bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
			      METx_,METy_,
			      mW_+sigmaMW_*deltaMW1_,mW_+sigmaMW_*deltaMW2_,
			      mTop_+sigmaMTop_*deltaMTop1_,mTop_+sigmaMTop_*deltaMTop2_);
}

void topSystemChiSquare::setupNeutrinoEllipses()
{
  nuEllipseOneOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
  bJet1LogSFRangeOneOne_=nuEllipseOneOneCalc_.getBJetLogSFRange(nOneOneRanges_);

  nBJet1Ranges_ = nOneOneRanges_;
  bJet1LogSFRange_ = bJet1LogSFRangeOneOne_;

  nuEllipseTwoTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
  bJet2LogSFRangeTwoTwo_=nuEllipseTwoTwoCalc_.getBJetLogSFRange(nTwoTwoRanges_);

  nBJet2Ranges_ = nTwoTwoRanges_;
  bJet2LogSFRange_ = bJet2LogSFRangeTwoTwo_;

  nuEllipseOneOne_=nuEllipseOneOneCalc_.getHomogeneousNeutrinoEllipse();
  nuEllipseTwoTwo_=nuEllipseTwoTwoCalc_.getHomogeneousNeutrinoEllipse();

  nuEllipseOneTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
  bJet1LogSFRangeOneTwo_=nuEllipseOneTwoCalc_.getBJetLogSFRange(nOneTwoRanges_);

  nBJet1Ranges_ = nOneTwoRanges_;
  bJet1LogSFRange_ = bJet1LogSFRangeOneTwo_;

  nuEllipseTwoOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
  bJet2LogSFRangeTwoOne_=nuEllipseTwoOneCalc_.getBJetLogSFRange(nTwoOneRanges_);

  nBJet2Ranges_ = nTwoOneRanges_;
  bJet2LogSFRange_ = bJet2LogSFRangeTwoOne_;

  nuEllipseOneTwo_=nuEllipseOneTwoCalc_.getHomogeneousNeutrinoEllipse();
  nuEllipseTwoOne_=nuEllipseTwoOneCalc_.getHomogeneousNeutrinoEllipse();
}

void topSystemChiSquare::calcNeutrinoRanges()
{
  double bJet1PtDelta_temp(bJet1PtDelta_);
  double bJet2PtDelta_temp(bJet2PtDelta_);
  setBJetPtDeltas(0.,0.);

  if(pairingInfo_)
    {
      //cout << "OneOne" << endl;
      nuEllipseOneOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
      bJet1LogSFRangeOneOne_=nuEllipseOneOneCalc_.getBJetLogSFRange(nOneOneRanges_);
      //cout << "TwoTwo" << endl;
      nuEllipseTwoTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
      bJet2LogSFRangeTwoTwo_=nuEllipseTwoTwoCalc_.getBJetLogSFRange(nTwoTwoRanges_);
      if(nOneOneRanges_ > 0) nuEllipseOneOne_=nuEllipseOneOneCalc_.getHomogeneousNeutrinoEllipse();
      if(nTwoTwoRanges_ > 0) nuEllipseTwoTwo_=nuEllipseTwoTwoCalc_.getHomogeneousNeutrinoEllipse();
    }
  else
    {
      //cout << "OneTwo" << endl;
      nuEllipseOneTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
      bJet1LogSFRangeOneTwo_=nuEllipseOneTwoCalc_.getBJetLogSFRange(nOneTwoRanges_);
      //cout << "TwoOne" << endl;
      nuEllipseTwoOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
      bJet2LogSFRangeTwoOne_=nuEllipseTwoOneCalc_.getBJetLogSFRange(nTwoOneRanges_);
      if(nOneTwoRanges_ > 0) nuEllipseOneTwo_=nuEllipseOneTwoCalc_.getHomogeneousNeutrinoEllipse();
      if(nTwoOneRanges_ > 0) nuEllipseTwoOne_=nuEllipseTwoOneCalc_.getHomogeneousNeutrinoEllipse();
    }
  setBJetPtDeltas(bJet1PtDelta_temp,bJet2PtDelta_temp);
}

void topSystemChiSquare::calcNeutrinoEllipses()
{
  //cout << "calculating neutrino ellipses" << endl;
  if(pairingInfo_)
    {
      //cout <<"OneOne"<<endl;
      nuEllipseOneOneCalc_.calcNeutrinoEllipse();
      if(nuEllipseOneOneCalc_.badPoint()){
	cout << "OneOne bad point" << endl;
	//cout << "log bJet1 Pt Scaling factor is " << bJet1PtDelta_*bJet1PtWidth_ << endl;
	//cout << "bJet1 Pt Scaling factor is " << exp(bJet1PtDelta_*bJet1PtWidth_) << endl;
	//cout << "bJet1 Unscaled Information:\n"
	//     << "Px: " << angleAdjusted_bJet1Px_
	//     << "\nPy: " << angleAdjusted_bJet1Py_
	//     << "\nPz: " << angleAdjusted_bJet1Pz_
	//     << "\nP : " << sqrt(angleAdjusted_bJet1Px_*angleAdjusted_bJet1Px_ + angleAdjusted_bJet1Py_*angleAdjusted_bJet1Py_ + angleAdjusted_bJet1Pz_*angleAdjusted_bJet1Pz_)
	//     << "\nE : " << angleAdjusted_bJet1E_ << endl;
      }
      //cout <<"TwoTwo"<<endl;
      nuEllipseTwoTwoCalc_.calcNeutrinoEllipse();
      if(nuEllipseTwoTwoCalc_.badPoint())
	{
	  cout << "TwoTwo bad point" << endl;
	  //cout << "log bJet2 Pt Scaling factor is " << bJet2PtDelta_*bJet2PtWidth_ << endl;
	  //cout << "bJet2 Pt Scaling factor is " << exp(bJet2PtDelta_*bJet2PtWidth_) << endl;
	  //cout << "bJet2 Unscaled Information:\n"
	  //     << "Px: " << angleAdjusted_bJet2Px_
	  //     << "\nPy: " << angleAdjusted_bJet2Py_
	  //     << "\nPz: " << angleAdjusted_bJet2Pz_
	  //     << "\nP : " << sqrt(angleAdjusted_bJet2Px_*angleAdjusted_bJet2Px_ + angleAdjusted_bJet2Py_*angleAdjusted_bJet2Py_ + angleAdjusted_bJet2Pz_*angleAdjusted_bJet2Pz_)
	  //     << "\nE : " << angleAdjusted_bJet2E_ << endl;
	}
    }
  else
    {
      //cout <<"OneTwo"<<endl;
      nuEllipseOneTwoCalc_.calcNeutrinoEllipse();
      if(nuEllipseOneTwoCalc_.badPoint())
	{
	  cout << "OneTwo bad point" << endl;
	  //cout << "log bJet1 Pt Scaling factor is " << bJet1PtDelta_*bJet1PtWidth_ << endl;
	  //cout << "bJet1 Pt Scaling factor is " << exp(bJet1PtDelta_*bJet1PtWidth_) << endl;
	  //cout << "bJet1 Unscaled Information:\n"
	  //     << "Px: " << angleAdjusted_bJet1Px_
	  //     << "\nPy: " << angleAdjusted_bJet1Py_
	  //     << "\nPz: " << angleAdjusted_bJet1Pz_
	  //     << "\nP : " << sqrt(angleAdjusted_bJet1Px_*angleAdjusted_bJet1Px_ + angleAdjusted_bJet1Py_*angleAdjusted_bJet1Py_ + angleAdjusted_bJet1Pz_*angleAdjusted_bJet1Pz_)
	  //     << "\nE : " << angleAdjusted_bJet1E_ << endl;
	}
      //cout <<"TwoOne"<<endl;
      nuEllipseTwoOneCalc_.calcNeutrinoEllipse();
      if(nuEllipseTwoOneCalc_.badPoint())
	{
	  cout << "TwoOne bad point" << endl;
	  //cout << "log bJet2 Pt Scaling factor is " << bJet2PtDelta_*bJet2PtWidth_ << endl;
	  //cout << "bJet2 Pt Scaling factor is " << exp(bJet2PtDelta_*bJet2PtWidth_) << endl;
	  //cout << "bJet2 Unscaled Information:\n"
	  //     << "Px: " << angleAdjusted_bJet2Px_
	  //     << "\nPy: " << angleAdjusted_bJet2Py_
          //     << "\nPz: " << angleAdjusted_bJet2Pz_
	  //     << "\nP : " << sqrt(angleAdjusted_bJet2Px_*angleAdjusted_bJet2Px_ + angleAdjusted_bJet2Py_*angleAdjusted_bJet2Py_ + angleAdjusted_bJet2Pz_*angleAdjusted_bJet2Pz_)
	  //     << "\nE : " << angleAdjusted_bJet2E_ << endl;
	}
    }
}

bool topSystemChiSquare::calcNeutrinoSolutions()
{
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;
  if(pairingInfo_)
    {
      nuSolOne_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				  bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				  METx_,METy_,
				  mW_+sigmaMW_*deltaMW1_,mW_+sigmaMW_*deltaMW2_,
				  mTop_+sigmaMTop_*deltaMTop1_,mTop_+sigmaMTop_*deltaMTop2_);

      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);
      if(nu1E.size()>0) return false;

    }
  else
    {
      nuSolTwo_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				  bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				  METx_,METy_,
				  mW_+sigmaMW_*deltaMW1_,mW_+sigmaMW_*deltaMW2_,
				  mTop_+sigmaMTop_*deltaMTop1_,mTop_+sigmaMTop_*deltaMTop2_);

      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);
      if(nu1E.size()>0) return false;
    }
  return true;
}

void topSystemChiSquare::setBJetAngleDeltas(const double bJet1PhiDelta, const double bJet1EtaDelta, const double bJet2PhiDelta, const double bJet2EtaDelta)
{
  bJet1PhiDelta_ =bJet1PhiDelta;
  bJet1EtaDelta_ =bJet1EtaDelta;
  bJet2PhiDelta_ =bJet2PhiDelta;
  bJet2EtaDelta_ =bJet2EtaDelta;
  resetAngles();
}

void topSystemChiSquare::setStartingBJetAngleDeltas()
{
  if(whichEllipse_)
    {
      double bJet1Phi = reconstructed_bJet1Phi_+bJet1PhiDelta_*bJet1PhiWidth_;
      double bJet1Eta = reconstructed_bJet1Eta_+bJet1EtaDelta_*bJet1EtaWidth_;
      angleAdjusted_bJet1Px_ = reconstructed_bJet1Pt_*cos(bJet1Phi);
      angleAdjusted_bJet1Py_ = reconstructed_bJet1Pt_*sin(bJet1Phi);
      angleAdjusted_bJet1Pz_ = reconstructed_bJet1Pt_*sinh(bJet1Eta);
      angleAdjusted_bJet1P_  = reconstructed_bJet1Pt_*cosh(bJet1Eta);
      angleAdjusted_bJet1E_  = sqrt(reconstructed_bJet1Mass2_ 
				    + angleAdjusted_bJet1P_*angleAdjusted_bJet1P_);

      bJet1Px_ = angleAdjusted_bJet1Px_ ;
      bJet1Py_ = angleAdjusted_bJet1Py_ ;
      bJet1Pz_ = angleAdjusted_bJet1Pz_ ;
      bJet1E_  = angleAdjusted_bJet1E_  ;
    }
  else
    {
      double bJet2Phi = reconstructed_bJet2Phi_+bJet2PhiDelta_*bJet2PhiWidth_;
      double bJet2Eta = reconstructed_bJet2Eta_+bJet2EtaDelta_*bJet2EtaWidth_;
      angleAdjusted_bJet2Px_ = reconstructed_bJet2Pt_*cos(bJet2Phi);
      angleAdjusted_bJet2Py_ = reconstructed_bJet2Pt_*sin(bJet2Phi);
      angleAdjusted_bJet2Pz_ = reconstructed_bJet2Pt_*sinh(bJet2Eta);
      angleAdjusted_bJet2P_  = reconstructed_bJet2Pt_*cosh(bJet2Eta);
      angleAdjusted_bJet2E_  = sqrt(reconstructed_bJet2Mass2_ 
				    + angleAdjusted_bJet2P_*angleAdjusted_bJet2P_);

      bJet2Px_ = angleAdjusted_bJet2Px_ ;
      bJet2Py_ = angleAdjusted_bJet2Py_ ;
      bJet2Pz_ = angleAdjusted_bJet2Pz_ ;
      bJet2E_  = angleAdjusted_bJet2E_  ;
    }
}

void topSystemChiSquare::getStartingTopMassRange()
{
  double bJetPx (bJet2Px_) ;
  double bJetPy (bJet2Py_) ;
  double bJetPz (bJet2Pz_) ;
  double bJetE  (bJet2E_ ) ;
  double bJetP2 (bJet2Px_*bJet2Px_ + bJet2Py_*bJet2Py_ + bJet2Pz_*bJet2Pz_);
  double deltaMW = deltaMW2_;

  double leptonPx (lepton2Px_) ;
  double leptonPy (lepton2Py_) ;
  double leptonPz (lepton2Pz_) ;
  double leptonE  (lepton2E_ ) ;
  double leptonP2 (lepton2Px_*lepton2Px_ + lepton2Py_*lepton2Py_ + lepton2Pz_*lepton2Pz_);

  if(whichEllipse_)
    {
      bJetPx =bJet1Px_ ;
      bJetPy =bJet1Py_ ;
      bJetPz =bJet1Pz_ ;
      bJetE  =bJet1E_  ;
      bJetP2 =bJet2Px_*bJet2Px_ + bJet2Py_*bJet2Py_ + bJet2Pz_*bJet2Pz_;
      deltaMW = deltaMW1_;

      if(pairingInfo_)
	{
	  leptonPx = lepton1Px_ ;
	  leptonPy = lepton1Py_ ;
	  leptonPz = lepton1Pz_ ;
	  leptonE  = lepton1E_  ;
	  leptonP2 = lepton1Px_*lepton1Px_ + lepton1Py_*lepton1Py_ + lepton1Pz_*lepton1Pz_;
	}
    }
  else if(!pairingInfo_)
    {
      leptonPx = lepton1Px_ ;
      leptonPy = lepton1Py_ ;
      leptonPz = lepton1Pz_ ;
      leptonE  = lepton1E_  ;
      leptonP2 = lepton1Px_*lepton1Px_ + lepton1Py_*lepton1Py_ + lepton1Pz_*lepton1Pz_;
    }

  double bJetP = sqrt(bJetP2);
  double leptonP = sqrt(leptonP2);

  double leptonP4 = leptonP2*leptonP2;

  double bJetE2 = bJetE*bJetE;

  double leptonE2 = leptonE*leptonE;
  double leptonE3 = leptonE2*leptonE;
  double leptonE4 = leptonE2*leptonE2;

  double c = (leptonPx*bJetPx+leptonPy*bJetPy+leptonPz*bJetPz);
  double s2 = 1.-c*c/(bJetP2*leptonP2);
  c/=(bJetP*leptonP);
  double mnu2 = mNu_*mNu_;
  double mW = (mW_+deltaMW*sigmaMW_);
  double mW2 = mW*mW;
 
  double mTopEdge1 = (bJetE*leptonE3 - bJetE*leptonE*(mnu2 - mW2 + leptonP2) + leptonE2*(bJetE2 + mW2 - bJetP*(bJetP + c*leptonP)) + 
		      leptonP*((-bJetE2 - mW2 + bJetP2)*leptonP + c*bJetP*(mnu2 - mW2 + leptonP2)) - 
		      sqrt((leptonE4 + pow(mnu2 - mW2,2) + 2*(mnu2 + mW2)*leptonP2 + leptonP4 - 2*leptonE2*(mnu2 + mW2 + leptonP2))*
			   (pow(c*leptonE*bJetP - bJetE*leptonP,2) + bJetP2*(leptonE - leptonP)*(leptonE + leptonP)*s2)))/((leptonE - leptonP)*(leptonE + leptonP));
  
  double mTopEdge2 = (bJetE*leptonE3 - bJetE*leptonE*(mnu2 - mW2 + leptonP2) + leptonE2*(bJetE2 + mW2 - bJetP*(bJetP + c*leptonP)) + 
		      leptonP*((-bJetE2 - mW2 + bJetP2)*leptonP + c*bJetP*(mnu2 - mW2 + leptonP2)) + 
		      sqrt((leptonE4 + pow(mnu2 - mW2,2) + 2*(mnu2 + mW2)*leptonP2 + leptonP4 - 2*leptonE2*(mnu2 + mW2 + leptonP2))*
			   (pow(c*leptonE*bJetP - bJetE*leptonP,2) + bJetP2*(leptonE - leptonP)*(leptonE + leptonP)*s2)))/((leptonE - leptonP)*(leptonE + leptonP));

  double mTopEdgeLow = min(mTopEdge1,mTopEdge2);
  double mTopEdgeHigh = max(mTopEdge1,mTopEdge2);

  //cout << "checking top mass range:" << endl;

  //cout << "bJetPhiDelta_ " << inputDeltas[2] << "\n"
  //     << "bJetEtaDelta_ " << inputDeltas[3] << "\n"
  //     << "bJetPtDelta_  " << bJetPtDelta    << "\n"
  //     << "deltaMW_      " << inputDeltas[0] << "\n"
  //     << "deltaMTop_    " << inputDeltas[1] << "\n"
  //     << "W   mass is   " << mW_+sigmaMW_*inputDeltas[0] << "\n"
  //     << "top mass is   " << mTop_+sigmaMTop_*inputDeltas[1] << endl;

  //cout << "top mass squared edges are " << mTopEdgeLow << " and " << mTopEdgeHigh << endl;

  if(mTopEdgeLow < mW2 && mTopEdgeHigh > mW2)
    {
      deltaMTopRangeLow_ = (mW-mTop_)/sigmaMTop_;
      deltaMTopRangeHigh_ = max(0.,(sqrt(mTopEdgeHigh) - mTop_)/sigmaMTop_);
    }
  else if(mTopEdgeLow >= mTop_*mTop_)
    {
      deltaMTopRangeLow_ = (mW-mTop_)/sigmaMTop_;
      deltaMTopRangeHigh_ = (sqrt(mTopEdgeLow) - mTop_)/sigmaMTop_;
    }
  else if(mTopEdgeLow > mW2)
    {
      deltaMTopRangeLow_ = (sqrt(mTopEdgeLow) - mTop_)/sigmaMTop_;
      deltaMTopRangeHigh_ = max(0.,(sqrt(mTopEdgeHigh) - mTop_)/sigmaMTop_);
    }
  else
    {
      deltaMTopRangeLow_ = 1.;
      deltaMTopRangeHigh_ = -1.;
    }

  //cout << "setting low edge to " << mTop_ + deltaMTopRangeLow_*sigmaMTop_ << endl;
  //cout << "setting high edge to " << mTop_ + deltaMTopRangeHigh_*sigmaMTop_ << endl;
}

void topSystemChiSquare::setBJetPtDeltas(const double bJet1PtDelta, const double bJet2PtDelta)
{
  bJet1PtDelta_ =bJet1PtDelta ;
  bJet2PtDelta_ =bJet2PtDelta ;
  resetPt();
}

void topSystemChiSquare::guessBJetPtDeltas()
{
  whichEllipse_ = true;
  double bJet1PtDelta(getMinBJetPtDelta()) ;
  whichEllipse_ = false;
  double bJet2PtDelta(getMinBJetPtDelta()) ;
  setBJetPtDeltas(bJet1PtDelta,bJet2PtDelta);
}

void topSystemChiSquare::setupLightJetChiSquare()
{
  getDxDyFromEllipses();
  lightJetChi2_ =lightJetChiSquare_.getChiSquare();
}

void topSystemChiSquare::setupBJetChiSquare()
{ 
  calcBJetChiSquare();

  if(currentBestChi2_ > bJetChi2_ || (currentBestChi2_!=currentBestChi2_ && bJetChi2_==bJetChi2_))
    {
      currentBestChi2_ = bJetChi2_;
      currentBestbJet1PtDelta_ = bJet1PtDelta_;
      currentBestbJet2PtDelta_ = bJet2PtDelta_;
      currentBestTheta1_ = theta1_;
      currentBestTheta2_ = theta2_;
    }
}

void topSystemChiSquare::setupTotalChiSquare()
{
  calcTotalChiSquare();

  if(pairingInfo_)
    {
      if(chi2BestSamePairing_ > chi2_ || (chi2BestSamePairing_!=chi2BestSamePairing_ && chi2_==chi2_))
  	{
  	  chi2BestSamePairing_ = chi2_;
  	  theta1BestSamePairing_ = theta1_;
  	  theta2BestSamePairing_ = theta2_;
  	  bJet1PtDeltaBestSamePairing_ = bJet1PtDelta_;
  	  bJet2PtDeltaBestSamePairing_ = bJet2PtDelta_;
  	  bJet1PhiDeltaBestSamePairing_ = bJet1PhiDelta_;
  	  bJet1EtaDeltaBestSamePairing_= bJet1EtaDelta_;
  	  bJet2PhiDeltaBestSamePairing_ = bJet2PhiDelta_;
  	  bJet2EtaDeltaBestSamePairing_= bJet2EtaDelta_;
	  deltaMW1BestSamePairing_     = deltaMW1_     ;
	  deltaMW2BestSamePairing_     = deltaMW2_     ;
	  deltaMTop1BestSamePairing_   = deltaMTop1_   ;
	  deltaMTop2BestSamePairing_   = deltaMTop2_   ;
  	}
    }
  else
    {
      if(chi2BestOppositePairing_ > chi2_ || (chi2BestOppositePairing_!=chi2BestOppositePairing_ && chi2_==chi2_))
  	{
  	  chi2BestOppositePairing_ = bJetChi2_;
  	  theta1BestOppositePairing_ = theta1_;
  	  theta2BestOppositePairing_ = theta2_;
  	  bJet1PtDeltaBestOppositePairing_ = bJet1PtDelta_;
  	  bJet2PtDeltaBestOppositePairing_ = bJet2PtDelta_;
  	  bJet1PhiDeltaBestOppositePairing_ = bJet1PhiDelta_;
  	  bJet1EtaDeltaBestOppositePairing_= bJet1EtaDelta_;
  	  bJet2PhiDeltaBestOppositePairing_ = bJet2PhiDelta_;
  	  bJet2EtaDeltaBestOppositePairing_= bJet2EtaDelta_;
	  deltaMW1BestOppositePairing_     = deltaMW1_     ;
	  deltaMW2BestOppositePairing_     = deltaMW2_     ;
	  deltaMTop1BestOppositePairing_   = deltaMTop1_   ;
	  deltaMTop2BestOppositePairing_   = deltaMTop2_   ;
  	}
    }
}

void topSystemChiSquare::calcBJetChiSquare()
{
  bJetChi2_ = lightJetChi2_ 
    + bJet1PtDelta_*bJet1PtDelta_ 
    + bJet2PtDelta_*bJet2PtDelta_ ;
}

void topSystemChiSquare::calcTotalChiSquare()
{
  chi2_ = bJetChi2_
    + bJet1PhiDelta_*bJet1PhiDelta_ 
    + bJet2PhiDelta_*bJet2PhiDelta_
    + bJet1EtaDelta_*bJet1EtaDelta_ 
    + bJet2EtaDelta_*bJet2EtaDelta_ 
    + breitWignerError(mW_, sigmaMW_, deltaMW1_)
    + breitWignerError(mW_, sigmaMW_, deltaMW2_)
    + breitWignerError(mTop_, sigmaMTop_, deltaMTop1_)
    + breitWignerError(mTop_, sigmaMTop_, deltaMTop2_);
}

void topSystemChiSquare::getDxDyFromEllipses()
{
  double nu1Array[3] = {cos(theta1_),sin(theta1_),1.}; 
  double nu2Array[3] = {cos(theta2_),sin(theta2_),1.}; 

  TVectorD nu1Perp(3,nu1Array), nu2Perp(3,nu2Array);

  if(pairingInfo_)
    {
      nu1Perp*=*nuEllipseOneOne_;
      nu2Perp*=*nuEllipseTwoTwo_;
    }
  else
    {
      nu1Perp*=*nuEllipseOneTwo_;
      nu2Perp*=*nuEllipseTwoOne_;
    }

  //FIXME! should this be negative or positive?  I think this is right....
  dx_ = -nu1Perp[0] - nu2Perp[0] + METx_;
  dy_ = -nu1Perp[1] - nu2Perp[1] + METy_;
}

void topSystemChiSquare::buildBestNeutrinos()
{
  //calcNeutrinoEllipses();
  TVectorD nu1P(3), nu2P(3);
  if(pairingInfoBest_)
    {
      nu1P=*(nuEllipseOneOneCalc_.getNeutrinoMomentum(theta1Best_));
      nu2P=*(nuEllipseTwoTwoCalc_.getNeutrinoMomentum(theta2Best_));
    }
  else
    {
      nu1P=*(nuEllipseOneTwoCalc_.getNeutrinoMomentum(theta1Best_));
      nu2P=*(nuEllipseTwoOneCalc_.getNeutrinoMomentum(theta2Best_));
    }
  nu1Best_.SetPxPyPzE(nu1P[0],nu1P[1],nu1P[2],sqrt(pow(nu1P[0],2)+pow(nu1P[1],2)+pow(nu1P[2],2)));
  nu2Best_.SetPxPyPzE(nu2P[0],nu2P[1],nu2P[2],sqrt(pow(nu2P[0],2)+pow(nu2P[1],2)+pow(nu2P[2],2)));
}

bool topSystemChiSquare::getStartingValues()
{
  outerMin_->Clear();
  
  setBJetAngleDeltas(0.,0.,0.,0.);
  setBJetPtDeltas(0.,0.);

  //whichEllipse_ = true;
  //if(getMinBJetPtDelta()>0.)
  //  {
  //    for(int i=-5; i<6; i++)
  //	{
  //	  deltaMTop1_ = double(i)/sigmaMTop_;
  //	  for (int j =-5; j<6; j++)
  //	    {	
  //	      deltaMW1_ = double(i)/sigmaMW_;
  //	      cout << "for top mass = " << mTop_+sigmaMTop_*deltaMTop1_ << " and W mass = " << mW_+sigmaMW_*deltaMW1_ << endl;
  //	      cout << "pairing one two has Z2 = " << nuEllipseOneTwoCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_) << endl;
  //	      cout << "pairing one two has bJet Pt = " << getMinBJetPtDelta() << endl;
  //	    }
  //	}
  //  }
  //if(pairingInfo_)
  //  {
  //    whichEllipse_ = false;
  //    cout << "for top mass = " << mTop_ << " and W mass = " << mW_ << endl;
  //    cout << "pairing two two has Z2 = " << nuEllipseTwoTwoCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_) << endl;
  //    cout << "pairing two two has bJet Pt = " << getMinBJetPtDelta() << endl;
  //    if(abs(getMinBJetPtDelta())>0.)
  //	{
  //	  for(int i=-5; i<6; i++)
  //	    {
  //	      deltaMTop2_ = double(i)/sigmaMTop_;
  //	      for (int j =-5; j<6; j++)
  //		{
  //		  deltaMW2_ = double(j)/sigmaMW_;
  //		  cout << "for top mass = " << mTop_+sigmaMTop_*deltaMTop2_ << " and W mass = " << mW_+sigmaMW_*deltaMW2_ << endl;
  //		  //cout << "pairing one one has Z2 = " << nuEllipseOneOneCalc_.getZ2(mTop_+i,mW_+j,mNu_) << endl;
  //		  cout << "pairing two one has Z2 = " << nuEllipseTwoTwoCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_) << endl;
  //		  //cout << "pairing one two has Z2 = " << nuEllipseOneTwoCalc_.getZ2(mTop_+i,mW_+j,mNu_) << endl;
  //		  //cout << "pairing two one has Z2 = " << nuEllipseTwoOneCalc_.getZ2(mTop_+i,mW_+j,mNu_) << endl;
  //		  cout << "pairing two one has bJet Pt = " << getMinBJetPtDelta() << endl;
  //		}
  //	    }
  //	}
  //  }

  startingOuterMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  startingOuterMin_->SetMaxFunctionCalls(1000000);
  startingOuterMin_->SetTolerance(0.1);
  startingOuterMin_->SetPrintLevel(0);

  //Set up the functor
  ROOT::Math::Functor outerFunc(this,&topSystemChiSquare::startingValueOuterMinimizationOperator,3);

  startingOuterMin_->SetFunction(outerFunc);

  startingTopMassMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  startingTopMassMin_->SetMaxFunctionCalls(1000000);
  startingTopMassMin_->SetTolerance(0.1);
  startingTopMassMin_->SetPrintLevel(0);

  //Set up the functor
  ROOT::Math::Functor topMassFunc(this,&topSystemChiSquare::startingValueTopMassMinimizationOperator,1);

  startingTopMassMin_->SetFunction(topMassFunc);

  whichEllipse_ = true;

  startingValueChi2_ = getMinBJetPtDelta();
  deltaMTopBest_     = 0.;
  deltaMWBest_       = 0.;
  bJetPhiDeltaBest_  = 0.;
  bJetEtaDeltaBest_  = 0.;
  startingValueChi2_*=startingValueChi2_;

  bool startingValuesFound(true);

  if(startingValueChi2_ > 0.)
    {
      startingOuterMin_->SetVariable(0,"bJetPhiDelta",0.,0.1);
      startingOuterMin_->SetVariable(1,"bJetEtaDelta",0.,0.1);
      startingOuterMin_->SetLowerLimitedVariable(2,"deltaMW"     ,0.,0.1,-mW_/sigmaMW_);
      startingOuterMin_->Minimize();

      //startingOuterMin_->PrintResults();

      bJet1PhiDelta_ = bJetPhiDeltaBest_;  
      bJet1EtaDelta_ = bJetEtaDeltaBest_;  
      deltaMW1_      = deltaMWBest_;	   
      deltaMTop1_    = deltaMTopBest_;	   
      bJet1PtDelta_  = getMinBJetPtDelta();
    }
  else
    {
      cout << "exactly zero" << endl;
      bJet1PhiDelta_ = 0.;
      bJet1EtaDelta_ = 0.;
      deltaMW1_      = 0.;
      deltaMTop1_    = 0.;
      bJet1PtDelta_  = 0.;
    }

  cout << "initial values set in first top:"    << endl;
  cout << "bJet1PhiDelta_ : " << bJet1PhiDelta_ << endl;
  cout << "bJet1EtaDelta_ : " << bJet1EtaDelta_ << endl;
  cout << "deltaMW1_      : " << deltaMW1_      << endl;
  cout << "deltaMTop1_    : " << deltaMTop1_    << endl;
  cout << "bJet1PtDelta_  : " << bJet1PtDelta_  << endl;

  if(startingValueChi2_ > maxConsideredChiSquareRoot_*maxConsideredChiSquareRoot_) startingValuesFound = false;

  whichEllipse_ = false;
  startingOuterMin_->Clear();

  startingValueChi2_ = getMinBJetPtDelta();
  deltaMTopBest_     = 0.;
  deltaMWBest_       = 0.;
  bJetPhiDeltaBest_  = 0.;
  bJetEtaDeltaBest_  = 0.;
  startingValueChi2_*=startingValueChi2_;

  if(startingValueChi2_ > 0.)
    { 
      startingOuterMin_->SetVariable(0,"bJetPhiDelta",0.,0.1);
      startingOuterMin_->SetVariable(1,"bJetEtaDelta",0.,0.1);
      startingOuterMin_->SetLowerLimitedVariable(2,"deltaMW"     ,0.,0.1,-mW_/sigmaMW_);
      startingOuterMin_->Minimize();

      //startingOuterMin_->PrintResults();

      bJet2PhiDelta_ = bJetPhiDeltaBest_;
      bJet2EtaDelta_ = bJetEtaDeltaBest_;
      deltaMW2_      = deltaMWBest_;
      deltaMTop2_    = deltaMTopBest_;
      bJet2PtDelta_  = getMinBJetPtDelta();
    }
  else
    {
      cout << "exactly zero" << endl;
      bJet2PhiDelta_ = 0.;
      bJet2EtaDelta_ = 0.;
      deltaMW2_      = 0.;
      deltaMTop2_    = 0.;
      bJet2PtDelta_  = 0.;
    }

  cout << "initial values set in second top:"   << endl;
  cout << "bJet2PhiDelta_ : " << bJet2PhiDelta_ << endl;
  cout << "bJet2EtaDelta_ : " << bJet2EtaDelta_ << endl;
  cout << "deltaMW2_      : " << deltaMW2_      << endl;
  cout << "deltaMTop2_    : " << deltaMTop2_    << endl;
  cout << "bJet2PtDelta_  : " << bJet2PtDelta_  << endl;

  if(startingValueChi2_ > maxConsideredChiSquareRoot_*maxConsideredChiSquareRoot_) startingValuesFound = false;

  delete startingOuterMin_;
  delete startingTopMassMin_;

  setBJetAngleDeltas(bJet1PhiDelta_,bJet1EtaDelta_,bJet2PhiDelta_,bJet2EtaDelta_);
  setBJetPtDeltas(bJet1PtDelta_,bJet2PtDelta_);
  calcNeutrinoRanges();


  lightJetChi2_ = 0;
  calcBJetChiSquare();
  calcTotalChiSquare();
  
  bool noSolutions = calcNeutrinoSolutions();
  if(noSolutions && startingValuesFound)
    {
      minimizeLightJetChiSquare(20);
    }
  else
    { 
      lightJetChi2_ = 0;
      theta1_ = 0.;
      theta2_ = 0.;
    }
  calcBJetChiSquare();
  calcTotalChiSquare();

  cout << "Done finding starting values:\n"
       << "bJet1PhiDelta_ " << bJet1PhiDelta_ << "\n"
       << "bJet1EtaDelta_ " << bJet1EtaDelta_ << "\n"
       << "bJet2PhiDelta_ " << bJet2PhiDelta_ << "\n"
       << "bJet2EtaDelta_ " << bJet2EtaDelta_ << "\n"
       << "bJet1PtDelta_  " << bJet1PtDelta_  << "\n"
       << "bJet2PtDelta_  " << bJet2PtDelta_  << "\n"
       << "deltaMW1_      " << deltaMW1_      << "\n"
       << "deltaMW2_      " << deltaMW2_      << "\n"
       << "deltaMTop1_    " << deltaMTop1_    << "\n"
       << "deltaMTop2_    " << deltaMTop2_    << "\n"
       << "W 1   mass is " << mW_+sigmaMW_*deltaMW1_ << "\n"
       << "top 1 mass is " << mTop_+sigmaMTop_*deltaMTop1_ << "\n"
       << "W 2   mass is " << mW_+sigmaMW_*deltaMW2_ << "\n"
       << "top 2 mass is " << mTop_+sigmaMTop_*deltaMTop2_   << "\n"
       << "theta1_ is    " << theta1_ << "\n"
       << "theta2_ is    " << theta2_ << endl;

  cout << "Chi square starting values:\n"
       << "light jet component: " << lightJetChi2_ << "\n"
       << (noSolutions?"(non-intersecting ellipses)\n":"(intersecting ellipses)\n")
       << "b jet pt plus light jet chi square: " << bJetChi2_ << "\n"
       << "all component chi square " << chi2_ << endl;

  if(pairingInfo_)
    {
      chi2BestSamePairing_ = chi2_;
      theta1BestSamePairing_ = theta1_;
      theta2BestSamePairing_ = theta2_;
      bJet1PtDeltaBestSamePairing_ = bJet1PtDelta_;
      bJet2PtDeltaBestSamePairing_ = bJet2PtDelta_;
      bJet1PhiDeltaBestSamePairing_ = bJet1PhiDelta_;
      bJet1EtaDeltaBestSamePairing_= bJet1EtaDelta_;
      bJet2PhiDeltaBestSamePairing_ = bJet2PhiDelta_;
      bJet2EtaDeltaBestSamePairing_= bJet2EtaDelta_;
      deltaMW1BestSamePairing_     = deltaMW1_     ;
      deltaMW2BestSamePairing_     = deltaMW2_     ;
      deltaMTop1BestSamePairing_   = deltaMTop1_   ;
      deltaMTop2BestSamePairing_   = deltaMTop2_   ;
    }
  else
    {
      chi2BestOppositePairing_ = chi2_;
      theta1BestOppositePairing_ = theta1_;
      theta2BestOppositePairing_ = theta2_;
      bJet1PtDeltaBestOppositePairing_ = bJet1PtDelta_;
      bJet2PtDeltaBestOppositePairing_ = bJet2PtDelta_;
      bJet1PhiDeltaBestOppositePairing_ = bJet1PhiDelta_;
      bJet1EtaDeltaBestOppositePairing_= bJet1EtaDelta_;
      bJet2PhiDeltaBestOppositePairing_ = bJet2PhiDelta_;
      bJet2EtaDeltaBestOppositePairing_= bJet2EtaDelta_;
      deltaMW1BestOppositePairing_     = deltaMW1_     ;
      deltaMW2BestOppositePairing_     = deltaMW2_     ;
      deltaMTop1BestOppositePairing_   = deltaMTop1_   ;
      deltaMTop2BestOppositePairing_   = deltaMTop2_   ;
    }
  currentBestChi2_ = bJetChi2_;
  currentBestbJet1PtDelta_ = bJet1PtDelta_;
  currentBestbJet2PtDelta_ = bJet2PtDelta_;
  currentBestTheta1_ = theta1_;
  currentBestTheta2_ = theta2_;
  
  outerMin_->SetLimitedVariable(0,"bJetPhiDelta_1",bJet1PhiDelta_,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(1,"bJetEtaDelta_1",bJet1EtaDelta_,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(2,"bJetPhiDelta_2",bJet2PhiDelta_,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(3,"bJetEtaDelta_2",bJet2EtaDelta_,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(4,"deltaMW_1"     ,deltaMW1_     ,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(5,"deltaMW_2"     ,deltaMW2_     ,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(6,"deltaMTop_1"   ,deltaMTop1_   ,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);
  outerMin_->SetLimitedVariable(7,"deltaMTop_2"   ,deltaMTop2_   ,0.1,-maxConsideredChiSquareRoot_,maxConsideredChiSquareRoot_);

  if(chi2_ > maxConsideredChiSquareRoot_*maxConsideredChiSquareRoot_ || chi2_ != chi2_) return false;

  return true;
}

double topSystemChiSquare::startingValueOuterMinimizationOperator(const double* inputDeltas)
{				
  if(whichEllipse_)
    {
      bJet1PhiDelta_ = inputDeltas[0];
      bJet1EtaDelta_ = inputDeltas[1];
      deltaMW1_      = inputDeltas[2];
    }
  else
    {
      bJet2PhiDelta_ = inputDeltas[0];
      bJet2EtaDelta_ = inputDeltas[1];
      deltaMW2_      = inputDeltas[2];
    }

  setStartingBJetAngleDeltas();

  getStartingTopMassRange();

  //minimization of mass variables

  startingTopMassMin_->Clear();
  
  double topMassDelta(0.);
  if(deltaMTopRangeLow_ > 0.) topMassDelta = deltaMTopRangeLow_;
  if(deltaMTopRangeHigh_ < 0. ) topMassDelta = deltaMTopRangeHigh_;
  if(deltaMTopRangeLow_>=deltaMTopRangeHigh_) return 1e99;
  if(whichEllipse_)
    {
      deltaMTop1_    = topMassDelta;
    }
  else
    {
      deltaMTop2_    = topMassDelta;
    }

  double bJetPtDelta = getMinBJetPtDelta();

  //cout << "For W mass = " << mW_+inputDeltas[2]*sigmaMW_ << " i.e. " << inputDeltas[2] << " sigma\n"
  //     << " top mass set to " << mTop_+sigmaMTop_*topMassDelta << " in (" << mTop_+sigmaMTop_*deltaMTopRangeLow_ <<"," << mTop_+sigmaMTop_*deltaMTopRangeHigh_ << ")" << endl; 
  //
  //cout << "ellipse has Z2 = " << getZ2() << " and " 
  //     << "min bJet Pt = " << bJetPtDelta << endl;

  if(bJetPtDelta == 0.)
    {
      double tempChi2 = topMassDelta*topMassDelta
                        + inputDeltas[0]*inputDeltas[0] 
	                + inputDeltas[1]*inputDeltas[1] 
	                + breitWignerError(mW_, sigmaMW_, inputDeltas[2]); 

      //cout << "best top mass still " << mTop_+sigmaMTop_*topMassDelta << " where" << endl;
      //cout << "ellipse has Z2 = " << getZ2() << " and " 
      //	   << "min bJet Pt = " << bJetPtDelta << endl;
      //cout << "return chi^2 = " << tempChi2 << endl;


      if(tempChi2 < startingValueChi2_)
	{
	  startingValueChi2_ = tempChi2;
	  deltaMTopBest_     = 0.;
	  deltaMWBest_       = inputDeltas[2];
	  bJetPhiDeltaBest_  = inputDeltas[0];
	  bJetEtaDeltaBest_  = inputDeltas[1];
	}
      return tempChi2;
    }

  startingTopMassMin_->SetLimitedVariable(0,"topMassDelta",topMassDelta,0.1,deltaMTopRangeLow_,deltaMTopRangeHigh_);

  startingTopMassMin_->Minimize();

  double tempChi2 = startingTopMassMin_->MinValue()
                    + inputDeltas[0]*inputDeltas[0] 
                    + inputDeltas[1]*inputDeltas[1] 
                    + breitWignerError(mW_, sigmaMW_, inputDeltas[2]); 

  if(whichEllipse_)
    {
      deltaMTop1_    = startingTopMassMin_->X()[0];
    }
  else
    {
      deltaMTop2_    = startingTopMassMin_->X()[0];
    }

  bJetPtDelta = getMinBJetPtDelta();

  //cout << "best top mass found is " << mTop_+sigmaMTop_*startingTopMassMin_->X()[0] << " where" << endl;
  //cout << "ellipse has Z2 = " << getZ2() << " and " 
  //     << "min bJet Pt = " << bJetPtDelta << endl;

  if(tempChi2 < startingValueChi2_)
    {
      startingValueChi2_ = tempChi2;
      deltaMTopBest_     = startingTopMassMin_->X()[0];
      deltaMWBest_       = inputDeltas[2];
      bJetPhiDeltaBest_  = inputDeltas[0];
      bJetEtaDeltaBest_  = inputDeltas[1];
    }

  //cout << "return chi^2 = " << tempChi2 << endl;

  return tempChi2;
}

double topSystemChiSquare::getZ2()
{
  if(pairingInfo_)
    {
      if(whichEllipse_)
	{
	  return nuEllipseOneOneCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
  	}
      else
	{
	  return nuEllipseTwoTwoCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
	}
    }
  else
    {
      if(whichEllipse_)
	{
	  return nuEllipseOneTwoCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
  	}
      else
	{
	  return nuEllipseTwoOneCalc_.getZ2(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
	}
    }
}

double topSystemChiSquare::startingValueTopMassMinimizationOperator(const double* inputDeltas)
{				
  if(whichEllipse_)
    {
      deltaMTop1_    = inputDeltas[0];
    }
  else
    {
      deltaMTop2_    = inputDeltas[0];
    }

  //cout << "for delta m top = " << inputDeltas[0] << endl;
  double bJetPtDelta = getMinBJetPtDelta();
  //cout << "bjet pt delta is " << bJetPtDelta << endl;

  return bJetPtDelta*bJetPtDelta 
         + breitWignerError(mTop_, sigmaMTop_, inputDeltas[0]);  ;
}

double topSystemChiSquare::getMinBJetPtDelta()
{
  //cout << "getting min b jet pt delta" << endl;
  bool hasSolutions(false);
  double bJetLowerLimit(0);
  bool hasLowerLimit(false);
  double bJetUpperLimit(0);
  bool hasUpperLimit(false);

  if(getZ2()>0) return 0.;

  if(pairingInfo_)
    {
      if(whichEllipse_)
	{
	  nuEllipseOneOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
	  bJet1LogSFRangeOneOne_=nuEllipseOneOneCalc_.getBJetLogSFRange(nOneOneRanges_);
	  if(nOneOneRanges_ > 0)
	    {
	      hasSolutions = true;
	      bJetLowerLimit = bJet1LogSFRangeOneOne_.first.first;
	      hasLowerLimit = bJet1LogSFRangeOneOne_.first.second;
	      bJetUpperLimit = bJet1LogSFRangeOneOne_.second.first;
	      hasUpperLimit = bJet1LogSFRangeOneOne_.second.second;
	    }
	}
      else
	{
	  nuEllipseTwoTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
	  bJet2LogSFRangeTwoTwo_=nuEllipseTwoTwoCalc_.getBJetLogSFRange(nTwoTwoRanges_);
	  if(nTwoTwoRanges_ > 0)
	    {
	      hasSolutions = true;
	      bJetLowerLimit = bJet2LogSFRangeTwoTwo_.first.first;
	      hasLowerLimit = bJet2LogSFRangeTwoTwo_.first.second;
	      bJetUpperLimit = bJet2LogSFRangeTwoTwo_.second.first;
	      hasUpperLimit = bJet2LogSFRangeTwoTwo_.second.second;
	    }
	}
    }
  else
    {
      if(whichEllipse_)
	{
	  nuEllipseOneTwoCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop1_,mW_+sigmaMW_*deltaMW1_,mNu_);
	  bJet1LogSFRangeOneTwo_=nuEllipseOneTwoCalc_.getBJetLogSFRange(nOneTwoRanges_);
	  if(nOneTwoRanges_ > 0)
	    {
	      hasSolutions = true;
	      bJetLowerLimit = bJet1LogSFRangeOneTwo_.first.first;
	      hasLowerLimit = bJet1LogSFRangeOneTwo_.first.second;
	      bJetUpperLimit = bJet1LogSFRangeOneTwo_.second.first;
	      hasUpperLimit = bJet1LogSFRangeOneTwo_.second.second;
	    }
	}
      else
	{
	  nuEllipseTwoOneCalc_.setupEllipse(mTop_+sigmaMTop_*deltaMTop2_,mW_+sigmaMW_*deltaMW2_,mNu_);
	  bJet2LogSFRangeTwoOne_=nuEllipseTwoOneCalc_.getBJetLogSFRange(nTwoOneRanges_);
	  if(nTwoOneRanges_ > 0)
	    {
	      hasSolutions = true;
	      bJetLowerLimit = bJet2LogSFRangeTwoOne_.first.first;
	      hasLowerLimit = bJet2LogSFRangeTwoOne_.first.second;
	      bJetUpperLimit = bJet2LogSFRangeTwoOne_.second.first;
	      hasUpperLimit = bJet2LogSFRangeTwoOne_.second.second;
	    }
	}
    } 
  //cout << "bJet Range: (";
  //if(!hasLowerLimit)  cout << "-inf,";
  //else cout << bJetLowerLimit << ",";
  //if(!hasUpperLimit)  cout << "inf)" << endl;
  //else cout << bJetUpperLimit << ")" << endl;
  if(hasSolutions)
    {
      //cout << "found an edge";
      if(hasLowerLimit && hasUpperLimit)
	{
	  if(bJetLowerLimit < 0 && bJetUpperLimit > 0)
	    {
	      //cout << " around zero" << endl;
	      return 0;
	    }
	  if(bJetLowerLimit >= 0)
	    {
	      //cout << " at " << bJetLowerLimit << endl;
	      return bJetLowerLimit;
	    }
	  //cout <<" at " << bJetUpperLimit << endl;
	  return bJetUpperLimit;
	}
      if(hasLowerLimit)
	{
	  if(bJetLowerLimit > 0)
	    {
	      //cout << " at " << bJetLowerLimit << endl;
	      return bJetLowerLimit;
	    }
	  //cout << " around zero" << endl;
	  return 0;
	}
      if(hasUpperLimit)
	{
	  if(bJetUpperLimit < 0)
	    {
	      //cout << " at " << bJetUpperLimit << endl;
	      return bJetUpperLimit;
	    }
	  //cout << " around zero" << endl;
	  return 0;
	}
      //cout << " around zero" << endl;
      return 0;
    }
  //cout << "found no intervals";
  return 1e99;
}


double topSystemChiSquare::ellipseAngleMinimizationOperator(const double* inputAngles)
{
  theta1_ = inputAngles[0];
  theta2_ = inputAngles[1];
  setupLightJetChiSquare();
  return lightJetChi2_;
}

double topSystemChiSquare::bJetPtMinimizationOperator(const double* inputDeltas)
{
  bJet1PtDelta_ =inputDeltas[0];
  bJet2PtDelta_ =inputDeltas[1];

  resetPt();

  bool noSolutions = calcNeutrinoSolutions();
  if(noSolutions)
    {
      minimizeLightJetChiSquare(20);
    }
  else
    { 
      lightJetChi2_ = 0;
    }

  setupBJetChiSquare();
  return bJetChi2_;
}

double topSystemChiSquare::outerMinimizationOperator(const double* inputDeltas)
{
  bJet1PhiDelta_ =inputDeltas[0];
  bJet1EtaDelta_ =inputDeltas[1];
  bJet2PhiDelta_ =inputDeltas[2];
  bJet2EtaDelta_ =inputDeltas[3];
  deltaMW1_      =inputDeltas[4];
  deltaMW2_      =inputDeltas[5];
  deltaMTop1_    =inputDeltas[6];
  deltaMTop2_    =inputDeltas[7];

  resetAngles();
  //guessBJetPtDeltas();
  minimizeBJetPtChiSquare();
  return chi2_;
}

void topSystemChiSquare::minimizeLightJetChiSquareOnly()
{
  ellipseAngleMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  ellipseAngleMin_->SetMaxFunctionCalls(1000000);
  ellipseAngleMin_->SetTolerance(0.001);
  ellipseAngleMin_->SetPrintLevel(0);
      
  //Set up the functor
  ROOT::Math::Functor ellipseAngleFunc(this,&topSystemChiSquare::ellipseAngleMinimizationOperator,2);
      
  //Set up the minimization piece:
  ellipseAngleMin_->SetFunction(ellipseAngleFunc);
  pairingInfo_ = false;
  double tempChi2(0);
  for(unsigned int i = 0; i<2; i++)
    {
      if(i == 1) pairingInfo_ = true;
      minimizeLightJetChiSquare(50);
      if(i == 0)
	{
	  pairingInfoBest_ = pairingInfo_;
	  tempChi2 = ellipseAngleMin_->MinValue();
	  theta1Best_ = ellipseAngleMin_->X()[0];
	  theta2Best_ = ellipseAngleMin_->X()[1];
	}
      else 
	{
	  if((ellipseAngleMin_->MinValue() < tempChi2 && ellipseAngleMin_->MinValue() > 0) || (tempChi2 == 0 && ellipseAngleMin_->MinValue() > 0))
	    {
	      pairingInfoBest_ = pairingInfo_;
	      tempChi2 = ellipseAngleMin_->MinValue();
	      theta1Best_ = ellipseAngleMin_->X()[0];
	      theta2Best_ = ellipseAngleMin_->X()[1];
	    }
	}
    }
  delete ellipseAngleMin_;
  pairingInfo_ = pairingInfoBest_;
  theta1_ = theta1Best_;
  theta2_ = theta2Best_;
  lightJetChi2_ = tempChi2;
  bJetChi2_ = tempChi2;
  chi2_ = lightJetChi2_;
  getDxDyFromEllipses();
  dxBest_ = dx_;
  dyBest_ = dy_;
}

void topSystemChiSquare::minimizeLightJetChiSquare(int nPoints)
{
  ellipseAngleMin_->Clear();
  double theta1(0.),theta2(0.);

  calcNeutrinoEllipses();
  if(pairingInfo_)
    {
      if( nuEllipseOneOneCalc_.badPoint() || nuEllipseTwoTwoCalc_.badPoint() )
	{
	  cout << "Same Pairing topSystemChiSquare::minimizeLightJetPtChiSquare: ";
	  lightJetChi2_ = chi2BestSamePairing_ * 2.;
	  cout << "backing out" << endl;
	  return;
	}
    }
  else
    {
      if( nuEllipseOneTwoCalc_.badPoint() || nuEllipseTwoOneCalc_.badPoint() )
	{
	  cout << "Opposite Pairing topSystemChiSquare::minimizeLightJetPtChiSquare: ";
	  lightJetChi2_ = chi2BestOppositePairing_ * 2.;
	  cout << "backing out" << endl;
	  return;
	}
    }

  getDxDyFromEllipses();
  double startingChi2 =lightJetChiSquare_.getChiSquare();
  double twoPiOverN = 2.*3.14159265359/(double)nPoints;

  //cout << "starting distance is " << startingChi2 << endl;

  for(int iTheta1 = 0; iTheta1 < nPoints; iTheta1++)
    {
      for(int iTheta2 = 0; iTheta2 < nPoints; iTheta2++)
	{
	  if(iTheta1 == 0 && iTheta2 == 0) continue;
	  theta1_ = (double)iTheta1*twoPiOverN;
	  theta2_ = (double)iTheta2*twoPiOverN;
	  getDxDyFromEllipses();
	  double thisChi2 =lightJetChiSquare_.getChiSquare();
	  if(thisChi2 < startingChi2)
	    {
	      startingChi2 = thisChi2;
	      theta1 = theta1_;
	      theta2 = theta2_;
	    }
	}
    }

  //cout << "middle distance is " << startingChi2 << endl;  
  //cout << "at theta1 = " << theta1 << " and theta2 = " << theta2 << endl;

  //theta1_ = theta1;
  //theta2_ = theta2;
  //setupLightJetChiSquare();

  ellipseAngleMin_->SetVariable(0,"theta_1",theta1,0.02*3.14159265359);
  ellipseAngleMin_->SetVariable(1,"theta_2",theta2,0.02*3.14159265359);
  ellipseAngleMin_->Minimize();
  
  lightJetChi2_ = ellipseAngleMin_->MinValue();
  theta1_ = ellipseAngleMin_->X()[0];
  theta2_ = ellipseAngleMin_->X()[1];

  //cout << "ending distance is " << lightJetChi2_ << endl;  

  getDxDyFromEllipses();

  //cout <<"dx_ is " << dx_ << endl;
  //cout <<"dy_ is " << dy_ << endl;

  //cout << "theta1 ending is " << theta1_ << endl;
  //cout << "theta2 ending is " << theta2_ << endl;

}

void topSystemChiSquare::setBestCurrentChiSquare()
{
  bool noSolutions = calcNeutrinoSolutions();
  if(noSolutions)
    {
      minimizeLightJetChiSquare();
    }
  else
    { 
      lightJetChi2_ = 0;
    }
  setupBJetChiSquare();
}

bool topSystemChiSquare::setBJetPtRanges()
{
  calcNeutrinoRanges();
  if(pairingInfo_)
    {
      nBJet1Ranges_ = nOneOneRanges_;
      nBJet2Ranges_ = nTwoTwoRanges_;
      bJet1LogSFRange_ = bJet1LogSFRangeOneOne_;
      bJet2LogSFRange_ = bJet2LogSFRangeTwoTwo_;
      currentBestChi2_ =  2*chi2BestSamePairing_;
    }
  else
    {
      nBJet1Ranges_ = nOneTwoRanges_;
      nBJet2Ranges_ = nTwoOneRanges_;
      bJet1LogSFRange_ = bJet1LogSFRangeOneTwo_;
      bJet2LogSFRange_ = bJet2LogSFRangeTwoOne_;
      currentBestChi2_ =  2*chi2BestOppositePairing_;
    }

  if( nBJet1Ranges_ == 0 || nBJet2Ranges_ == 0 )
    {
      return false;
    }

  double bJet1LowerLimit(max(bJet1LogSFRange_.first.first,-maxConsideredChiSquareRoot_));

  if(!bJet1LogSFRange_.first.second) bJet1LowerLimit = -maxConsideredChiSquareRoot_;

  else if( bJet1LowerLimit > maxConsideredChiSquareRoot_) return false;

  double bJet1UpperLimit(min(bJet1LogSFRange_.second.first,maxConsideredChiSquareRoot_));

  if(!bJet1LogSFRange_.second.second) bJet1UpperLimit = maxConsideredChiSquareRoot_;

  else if(bJet1UpperLimit < -maxConsideredChiSquareRoot_) return false;

  double bJet2LowerLimit(max(bJet2LogSFRange_.first.first,-maxConsideredChiSquareRoot_));

  if(!bJet2LogSFRange_.first.second) bJet2LowerLimit = -maxConsideredChiSquareRoot_;

  else if(bJet2LowerLimit > maxConsideredChiSquareRoot_) return false;

  double bJet2UpperLimit(min(bJet2LogSFRange_.second.first,maxConsideredChiSquareRoot_));

  if(!bJet2LogSFRange_.second.second) bJet2UpperLimit = maxConsideredChiSquareRoot_;

  else if(bJet2UpperLimit < -maxConsideredChiSquareRoot_) return false;

  double bJet1PtDeltaStart(currentBestbJet1PtDelta_);

  if(bJet1LowerLimit > bJet1PtDeltaStart || bJet1UpperLimit < bJet1PtDeltaStart)
    {
      double width = abs(bJet1UpperLimit - bJet1LowerLimit);
      if(bJet1LowerLimit < 0 && bJet1UpperLimit > 0) bJet1PtDeltaStart = 0;
      else bJet1PtDeltaStart = (bJet1LowerLimit > 0)? 
	     bJet1LowerLimit + 0.01*width:
	     bJet1UpperLimit - 0.01*width; 
    }

  //if(bJet1LogSFRange_.first.second && bJet1LogSFRange_.second.second)
  //  {
  //    if(bJet1LowerLimit > bJet1PtDeltaStart || bJet1UpperLimit < bJet1PtDeltaStart)
  //	{
  //	  double width = abs(bJet1UpperLimit - bJet1LowerLimit);
  //	  bJet1PtDeltaStart = (bJet1LowerLimit > 0)? 
  //	    bJet1LowerLimit + 0.01*width:
  //	    bJet1UpperLimit - 0.01*width;
  //	}
  //  }
  //else if(bJet1LogSFRange_.first.second && bJet1LowerLimit > bJet1PtDeltaStart)
  //  {
  //    bJet1PtDeltaStart = bJet1LowerLimit + 0.01;
  //  }
  //else if(bJet1LogSFRange_.second.second && bJet1UpperLimit < bJet1PtDeltaStart)
  //  {
  //    bJet1PtDeltaStart = bJet1UpperLimit - 0.01;
  //  }

  double bJet2PtDeltaStart(currentBestbJet2PtDelta_);

  if(bJet2LowerLimit > bJet2PtDeltaStart || bJet2UpperLimit < bJet2PtDeltaStart)
    {
      double width = abs(bJet2UpperLimit - bJet2LowerLimit);
      if(bJet2LowerLimit < 0 && bJet2UpperLimit > 0) bJet2PtDeltaStart = 0;
      else bJet2PtDeltaStart = (bJet2LowerLimit > 0)? 
	     bJet2LowerLimit + 0.01*width:
	     bJet2UpperLimit - 0.01*width; 
    }

  //if(bJet2LogSFRange_.first.second && bJet2LogSFRange_.second.second)
  //  {
  //    if(bJet2LowerLimit > bJet2PtDeltaStart || bJet2UpperLimit < bJet2PtDeltaStart)
  //	{
  //	  double width = abs(bJet2UpperLimit - bJet2LowerLimit);
  //	  bJet2PtDeltaStart = (bJet2LowerLimit > 0)? 
  //	    bJet2LowerLimit + 0.01*width:
  //	    bJet2UpperLimit - 0.01*width;
  //	}
  //  }
  //else if(bJet2LogSFRange_.first.second && bJet2LowerLimit > bJet2PtDeltaStart)
  //  {
  //    bJet2PtDeltaStart = bJet2LowerLimit + 0.01;
  //  }
  //else if(bJet2LogSFRange_.second.second && bJet2UpperLimit < bJet2PtDeltaStart)
  //  {
  //    bJet2PtDeltaStart = bJet2UpperLimit - 0.01;
  //  }

  setBJetPtDeltas(bJet1PtDeltaStart,bJet2PtDeltaStart);

  ptMin_->Clear();

  ptMin_->SetLimitedVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart,0.1,bJet1LowerLimit,bJet1UpperLimit);
  ptMin_->SetLimitedVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart,0.1,bJet2LowerLimit,bJet2UpperLimit);

  //cout << "Setting b jet 1 range to be " << bJet1LowerLimit << " to " << bJet1UpperLimit << endl;
  //cout << "Setting b jet 2 range to be " << bJet2LowerLimit << " to " << bJet2UpperLimit << endl;

  return true;
  
  //if( bJet1LogSFRange_.first.second && bJet1LogSFRange_.second.second ) {
  //  //cout << "setting bjet1 limited" << endl;
  //  ptMin_->SetLimitedVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart,0.001,bJet1LowerLimit,bJet1UpperLimit);
  //}
  //else if( bJet1LogSFRange_.first.second ) {
  //  //cout << "setting bjet1 lower limited" << endl;
  //  ptMin_->SetLowerLimitedVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart,0.001,bJet1LowerLimit);
  //}
  //else if( bJet1LogSFRange_.second.second ) {
  //  //cout << "setting bjet1 upper limited" << endl;
  //  ptMin_->SetUpperLimitedVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart,0.001,bJet1UpperLimit);
  //}
  //else{
  //  //cout << "setting bjet1 unlimited" << endl;
  //  ptMin_->SetVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart,0.001);
  //}
  ////ptMin_->SetFixedVariable(0,"bJetPtDelta_1",bJet1PtDeltaStart);
  //
  //if( bJet2LogSFRange_.first.second && bJet2LogSFRange_.second.second ) {
  //  //cout << "setting bjet2 limited" << endl;
  //  ptMin_->SetLimitedVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart,0.001,bJet2LowerLimit,bJet2UpperLimit);
  //}
  //else if( bJet2LogSFRange_.first.second ){
  //  //cout << "setting bjet2 lower limited" << endl;
  //  ptMin_->SetLowerLimitedVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart,0.001,bJet2LowerLimit);
  //}
  //else if( bJet2LogSFRange_.second.second ){
  //  //cout << "setting bjet2 upper limited" << endl;
  //  ptMin_->SetUpperLimitedVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart,0.001,bJet2UpperLimit);
  //}
  //else{
  //  //cout << "setting bjet2 unlimited" << endl;
  //  ptMin_->SetVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart,0.001);
  //}
  ////ptMin_->SetFixedVariable(1,"bJetPtDelta_2",bJet2PtDeltaStart);
}

bool topSystemChiSquare::unstretchedBJetPtExists()
{
  calcNeutrinoRanges();

  pairingInfo_ = false;

  bool checkBit(true);

  for(unsigned int i = 0; i<2; i++)
    {
      if(i == 1) pairingInfo_ = true;

      if(pairingInfo_)
	{
	  nBJet1Ranges_ = nOneOneRanges_;
	  nBJet2Ranges_ = nTwoTwoRanges_;
	  bJet1LogSFRange_ = bJet1LogSFRangeOneOne_;
	  bJet2LogSFRange_ = bJet2LogSFRangeTwoTwo_;
	  currentBestChi2_ =  2*chi2BestSamePairing_;
	}
      else
	{
	  nBJet1Ranges_ = nOneTwoRanges_;
	  nBJet2Ranges_ = nTwoOneRanges_;
	  bJet1LogSFRange_ = bJet1LogSFRangeOneTwo_;
	  bJet2LogSFRange_ = bJet2LogSFRangeTwoOne_;
	  currentBestChi2_ =  2*chi2BestOppositePairing_;
	}

      if( nBJet1Ranges_ == 0 || nBJet2Ranges_ == 0 )
	{
	  checkBit = false;
	}

      double bJet1LowerLimit(bJet1LogSFRange_.first.first);

      if(bJet1LogSFRange_.first.second && bJet1LowerLimit > 0) checkBit = false;

      double bJet1UpperLimit(bJet1LogSFRange_.second.first);

      if(bJet1LogSFRange_.second.second && bJet1UpperLimit < 0) checkBit = false;

      double bJet2LowerLimit(bJet2LogSFRange_.first.first);

      if(bJet2LogSFRange_.first.second && bJet2LowerLimit > 0) checkBit = false;

      double bJet2UpperLimit(bJet2LogSFRange_.second.first);

      if(bJet2LogSFRange_.second.second && bJet2UpperLimit < 0) checkBit = false;

      if(checkBit) return true;
    }

  return false;

}

void topSystemChiSquare::minimizeBJetPtChiSquare()
{

  //cout << "running inner loop" << endl;

  if(!setBJetPtRanges()) return;

  setBestCurrentChiSquare();

  ptMin_->Minimize();

  double tempChi2 = ptMin_->MinValue();
  if(tempChi2 < currentBestChi2_)
    {
      bJet1PtDelta_ = ptMin_->X()[0];
      bJet2PtDelta_ = ptMin_->X()[1];
      setBestCurrentChiSquare();
    }
  setupTotalChiSquare();
}

void topSystemChiSquare::minimizeTotalChiSquare()
{
  pairingInfo_ = false;
  double tempChi2 = -1.;
  chi2BestSamePairing_ = 1e99;
  chi2BestOppositePairing_ = 1e99;
  //Run the minimization
  cout << "running the b-jet minimization" << endl;

  outerMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  outerMin_->SetMaxFunctionCalls(1000000);
  outerMin_->SetTolerance(0.001);
  outerMin_->SetPrintLevel(0);
  
  //Set up the functor
  ROOT::Math::Functor func(this,&topSystemChiSquare::outerMinimizationOperator,8);

  //Set up the minimization piece:
  outerMin_->SetFunction(func);

  ptMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  ptMin_->SetMaxFunctionCalls(1000000);
  ptMin_->SetTolerance(0.001);
  ptMin_->SetPrintLevel(0);
      
  //Set up the functor
  ROOT::Math::Functor ptFunc(this,&topSystemChiSquare::bJetPtMinimizationOperator,2);
      
  //Set up the minimization piece:
  ptMin_->SetFunction(ptFunc); 

  ellipseAngleMin_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  ellipseAngleMin_->SetMaxFunctionCalls(1000000);
  ellipseAngleMin_->SetTolerance(0.001);
  ellipseAngleMin_->SetPrintLevel(0);
      
  //Set up the functor
  ROOT::Math::Functor ellipseAngleFunc(this,&topSystemChiSquare::ellipseAngleMinimizationOperator,2);
      
  //Set up the minimization piece:
  ellipseAngleMin_->SetFunction(ellipseAngleFunc);

  for(unsigned int i = 0; i<2; i++)
    {
      cout << "i equals " << i << endl;

      //Set pairing of b-jets
      if(i == 1) pairingInfo_ = true;

      bool goOn = getStartingValues();
      if(!goOn)
	{
	  cout << "skipping this combination -- no good starting location" << endl;
	  continue;
	}
      cout << "doing this combination -- good starting location found" << endl;

      outerMin_->Minimize();
      double minValue = outerMin_->MinValue();
      //cout << "the minimum chi square is " << minValue << endl;
      //cout << "ending light jet chi square: " << lightJetChi2_ << endl;
      //outerMin_->SetPrintLevel(3);
      const double* values = outerMin_->X();
      if(tempChi2 < 0) 
       	{
       	  tempChi2=minValue;
       	  pairingInfoBest_   = pairingInfo_;
       	  bJet1PhiDeltaBest_ = values[0];
       	  bJet1EtaDeltaBest_ = values[1];
       	  bJet2PhiDeltaBest_ = values[2];
       	  bJet2EtaDeltaBest_ = values[3];
	  deltaMW1Best_      = values[4];
	  deltaMW2Best_      = values[5];
	  deltaMTop1Best_    = values[6];
	  deltaMTop2Best_    = values[7];
      	  theta1Best_        = currentBestTheta1_;
      	  theta2Best_        = currentBestTheta2_;
      	  bJet1PtDeltaBest_  = currentBestbJet1PtDelta_;
      	  bJet2PtDeltaBest_  = currentBestbJet2PtDelta_;
       	}
      else if(tempChi2 > minValue || (minValue==minValue && tempChi2!=tempChi2))
       	{
       	  tempChi2=minValue;
       	  pairingInfoBest_   = pairingInfo_;
       	  bJet1PhiDeltaBest_ = values[0];
       	  bJet1EtaDeltaBest_ = values[1];
       	  bJet2PhiDeltaBest_ = values[2];
       	  bJet2EtaDeltaBest_ = values[3];
	  deltaMW1Best_      = values[4];
	  deltaMW2Best_      = values[5];
	  deltaMTop1Best_    = values[6];
	  deltaMTop2Best_    = values[7];
      	  theta1Best_        = currentBestTheta1_;
      	  theta2Best_        = currentBestTheta2_;
      	  bJet1PtDeltaBest_  = currentBestbJet1PtDelta_;
      	  bJet2PtDeltaBest_  = currentBestbJet2PtDelta_;
      	}
      //double a[4] = {0,0,0,0};
      //outerMinimizationOperator(a);
    }
  cout << "done with minimization" << endl;
  if((chi2BestSamePairing_ < tempChi2 && chi2BestSamePairing_ < chi2BestOppositePairing_) 
     || (tempChi2!=tempChi2 && chi2BestSamePairing_ < chi2BestOppositePairing_)
     || (chi2BestSamePairing_ < tempChi2 && chi2BestOppositePairing_!=chi2BestOppositePairing_)
     || (tempChi2!=tempChi2 && chi2BestOppositePairing_!=chi2BestOppositePairing_ && chi2BestSamePairing_==chi2BestSamePairing_))
    {
      cout << "setting best values for same pairing" << endl;
      cout << "theta1Best_ = " << theta1BestSamePairing_ << endl;
      cout << "theta2Best_ = " << theta2BestSamePairing_ << endl;
      currentBestChi2_   = chi2BestSamePairing_;
      pairingInfoBest_   = true;
      theta1Best_        = theta1BestSamePairing_;
      theta2Best_        = theta2BestSamePairing_;
      bJet1PtDeltaBest_  = bJet1PtDeltaBestSamePairing_;
      bJet2PtDeltaBest_  = bJet2PtDeltaBestSamePairing_;
      bJet1PhiDeltaBest_ = bJet1PhiDeltaBestSamePairing_;
      bJet1EtaDeltaBest_ = bJet1EtaDeltaBestSamePairing_;
      bJet2PhiDeltaBest_ = bJet2PhiDeltaBestSamePairing_;
      bJet2EtaDeltaBest_ = bJet2EtaDeltaBestSamePairing_;
      deltaMW1Best_      = deltaMW1BestSamePairing_     ;
      deltaMW2Best_      = deltaMW2BestSamePairing_     ;
      deltaMTop1Best_    = deltaMTop1BestSamePairing_   ;
      deltaMTop2Best_    = deltaMTop2BestSamePairing_   ;
    }
  else if(chi2BestOppositePairing_ < tempChi2
     || tempChi2!=tempChi2
     || chi2BestSamePairing_!=chi2BestSamePairing_)
    {
      cout << "setting best values for opposite pairing" << endl;
      cout << "theta1Best_ = " << theta1BestOppositePairing_ << endl;
      cout << "theta2Best_ = " << theta2BestOppositePairing_ << endl;
      currentBestChi2_   = chi2BestOppositePairing_;
      pairingInfoBest_   = false;
      theta1Best_        = theta1BestOppositePairing_;
      theta2Best_        = theta2BestOppositePairing_;
      bJet1PtDeltaBest_  = bJet1PtDeltaBestOppositePairing_;
      bJet2PtDeltaBest_  = bJet2PtDeltaBestOppositePairing_;
      bJet1PhiDeltaBest_ = bJet1PhiDeltaBestOppositePairing_;
      bJet1EtaDeltaBest_ = bJet1EtaDeltaBestOppositePairing_;
      bJet2PhiDeltaBest_ = bJet2PhiDeltaBestOppositePairing_;
      bJet2EtaDeltaBest_ = bJet2EtaDeltaBestOppositePairing_;
      deltaMW1Best_      = deltaMW1BestOppositePairing_     ;
      deltaMW2Best_      = deltaMW2BestOppositePairing_     ;
      deltaMTop1Best_    = deltaMTop1BestOppositePairing_   ;
      deltaMTop2Best_    = deltaMTop2BestOppositePairing_   ;
    }
  else 
    {
      cout << "no best fit found" << endl;
      currentBestChi2_   = 1e99;
    }


  //pairingInfo_ = pairingInfoBest_;
  //setBJetAngleDeltas(bJet1PhiDeltaBest_,bJet1EtaDeltaBest_,bJet2PhiDeltaBest_,bJet2EtaDeltaBest_);
  //double angleDeltas[4] = {bJet1PhiDeltaBest_,bJet1EtaDeltaBest_,bJet2PhiDeltaBest_,bJet2EtaDeltaBest_};
  //outerMinimizationOperator(angleDeltas);
  //bJet1PtDeltaBest_ = currentBestbJet1PtDelta_;
  //bJet2PtDeltaBest_ = currentBestbJet2PtDelta_;
  //currentBestTheta1_ = currentBestTheta1_;
  //currentBestTheta2_ = currentBestTheta2_;
  delete ellipseAngleMin_;
  delete ptMin_;
  delete outerMin_; 

  setBestValues();

  //setBJetPtDeltas(0.,0.);
  //calcNeutrinoRanges();
  //setBJetPtDeltas(bJet1PtDeltaBest_,bJet2PtDeltaBest_);
  //theta1_ = theta1Best_;
  //theta2_ = theta2Best_;
  //calcNeutrinoEllipses();
  cout << "minimum chisquare for both pairings is " << currentBestChi2_ << endl;
  //bool noSolutions = calcNeutrinoSolutions();
  //if(noSolutions)
  //  {
  //    setupLightJetChiSquare();
  //    dxBest_ = dx_;
  //    dyBest_ = dy_;
  //  }
  //else
  //  { 
  //    lightJetChi2_ = 0;
  //    dxBest_ = 0.;
  //    dyBest_ = 0.;
  //  }
  cout << "Ending values for parameters:\n"
       << "bJet1PhiDeltaBest_ " << bJet1PhiDeltaBest_ << "\n"
       << "bJet1EtaDeltaBest_ " << bJet1EtaDeltaBest_ << "\n"
       << "bJet2PhiDeltaBest_ " << bJet2PhiDeltaBest_ << "\n"
       << "bJet2EtaDeltaBest_ " << bJet2EtaDeltaBest_ << "\n"
       << "bJet1PtDeltaBest_ "  << bJet1PtDeltaBest_  << "\n"
       << "bJet2PtDeltaBest_ "  << bJet2PtDeltaBest_  << "\n"
       << "deltaMW1Best_      " << deltaMW1Best_      << "\n"
       << "deltaMW2Best_      " << deltaMW2Best_      << "\n"
       << "deltaMTop1Best_    " << deltaMTop1Best_    << "\n"
      << "deltaMTop2Best_    " << deltaMTop2Best_    << "\n"
       << "W 1   mass is      " << mW_+sigmaMW_*deltaMW1Best_ << "\n"
       << "top 1 mass is      " << mTop_+sigmaMTop_*deltaMTop1Best_ << "\n"
       << "W 2   mass is      " << mW_+sigmaMW_*deltaMW2Best_ << "\n"
       << "top 2 mass is      " << mTop_+sigmaMTop_*deltaMTop2Best_ <<"\n"
       << "theta1Best_ is     " << theta1Best_      << "\n"
       << "theta2Best_ is     " << theta2Best_    << endl;

  cout << "Chi square starting values:\n"
       << "light jet component: " << lightJetChi2_ << "\n"
       << (lightJetChi2_!=0.?"(non-intersecting ellipses)\n":"(intersecting ellipses)\n")
       << "b jet pt plus light jet chi square: " << bJetChi2_ << "\n"
       << "all component chi square " << chi2_ << endl;

  buildBestLightJets();
  calcBJetChiSquare();
  calcTotalChiSquare();
  //getChiSquare();
  //cout << "the minimum chi square is " << chi2_ << endl;
  //minimizeLightJetChiSquareOnly();
  //setupBJetChiSquare();
  //cout << "ending light jet chi square: " << lightJetChi2_ << endl;
  //cout << "ending dx_: " << dx_ << " dy_: " << dy_ << endl;
  //lightJetChiSquare_.printResults();
}

void topSystemChiSquare::setBestValues()
{
      pairingInfo_   = pairingInfoBest_   ;
      setBJetAngleDeltas(bJet1PhiDeltaBest_,bJet1EtaDeltaBest_,bJet2PhiDeltaBest_,bJet2EtaDeltaBest_);
      deltaMW1_      = deltaMW1Best_      ;
      deltaMW2_      = deltaMW2Best_      ;
      deltaMTop1_    = deltaMTop1Best_    ;
      deltaMTop2_    = deltaMTop2Best_    ;
      setBJetPtDeltas(bJet1PtDeltaBest_,bJet2PtDeltaBest_);
      calcNeutrinoRanges();
      theta1_ = theta1Best_;
      theta2_ = theta2Best_;
      calcNeutrinoEllipses();
      bool noSolutions = calcNeutrinoSolutions();
      if(noSolutions)
	{
	  setupLightJetChiSquare();
	  dxBest_ = dx_;
	  dyBest_ = dy_;
	}
      else
	{ 
	  lightJetChi2_ = 0;
	  dxBest_ = 0.;
	  dyBest_ = 0.;
	}
      calcBJetChiSquare();
      calcTotalChiSquare();
      //cout << "Done setting best values:\n"
      // << "bJet1PhiDelta_ " << bJet1PhiDelta_ << "\n"
      // << "bJet1EtaDelta_ " << bJet1EtaDelta_ << "\n"
      // << "bJet2PhiDelta_ " << bJet2PhiDelta_ << "\n"
      // << "bJet2EtaDelta_ " << bJet2EtaDelta_ << "\n"
      // << "bJet1PtDelta_  " << bJet1PtDelta_  << "\n"
      // << "bJet2PtDelta_  " << bJet2PtDelta_  << "\n"
      // << "deltaMW1_      " << deltaMW1_      << "\n"
      // << "deltaMW2_      " << deltaMW2_      << "\n"
      // << "deltaMTop1_    " << deltaMTop1_    << "\n"
      // << "deltaMTop2_    " << deltaMTop2_    << "\n"
      // << "W 1   mass is " << mW_+sigmaMW_*deltaMW1_ << "\n"
      // << "top 1 mass is " << mTop_+sigmaMTop_*deltaMTop1_ << "\n"
      // << "W 2   mass is " << mW_+sigmaMW_*deltaMW2_ << "\n"
      // << "top 2 mass is " << mTop_+sigmaMTop_*deltaMTop2_   << "\n"
      // << "theta1_ is    " << theta1_ << "\n"
      // << "theta2_ is    " << theta2_ << endl;
}

double topSystemChiSquare::getChiSquare()
{
  calcTotalChiSquare();
  return chi2_;
}

double topSystemChiSquare::calcSolution()
{
  //First check whether there exist real solutions to the quartic equations
  //In this case, the neutrino ellipses intersect

  
  pairingInfo_ = true;
  bool solutionFlag1=!calcNeutrinoSolutions(); 
  pairingInfo_ = false;
  bool solutionFlag2=!calcNeutrinoSolutions();

  if(solutionFlag1||solutionFlag2) 
    {
      cout << "solution trivially attainable with intersecting ellipses!" << endl;
      intersectingEllipses_ = true;
    }
  if(solutionFlag1||solutionFlag2) return intersectingEllipsesChi2_;
  //If there are no real solutions, the ellipses don't touch
  //We need to minimize the global chi2
  minimizeTotalChiSquare();
  bool noSolutions = calcNeutrinoSolutions();
  if(!noSolutions)
    {
      cout << "created intersecting ellipses" << endl;
      //minimizeLightJetChiSquareOnly();
    }
  if(bJetChi2_<0) return -1.;
  //buildBestNeutrinos();
  return getChiSquare();
}

void topSystemChiSquare::fillBestMomenta(XYZTLorentzVector& bJet1, XYZTLorentzVector& bJet2,
					 XYZTLorentzVector& lightJet1, XYZTLorentzVector& lightJet2,
					 XYZTLorentzVector& neutrino1, XYZTLorentzVector& neutrino2,
					 XYZTLorentzVector& MET,
					 bool& pairingInfo,
					 vector<double>& minDeltas
					 )
{
  bJet1.SetPxPyPzE(0.,0.,0.,0.);
  bJet2.SetPxPyPzE(0.,0.,0.,0.);
  lightJet1.SetPxPyPzE(0.,0.,0.,0.);
  lightJet2.SetPxPyPzE(0.,0.,0.,0.);
  neutrino1.SetPxPyPzE(0.,0.,0.,0.);
  neutrino2.SetPxPyPzE(0.,0.,0.,0.);
  MET.SetPxPyPzE(0.,0.,0.,0.);
  pairingInfo = pairingInfoBest_;
  minDeltas.clear();
  minDeltas.assign(14,0.);

  pairingInfo_ = true;
  bool solutionFlag1=!calcNeutrinoSolutions();
  pairingInfo_ = false;
  bool solutionFlag2=!calcNeutrinoSolutions();
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;
  int minPzSol;
  double minPz;

  if(solutionFlag1&&solutionFlag2)
    {
      minPzSol = -1;
      minPz = 1.e99;
      bool bestSol1 = true;

      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);
      
      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          double thisPz = nu1pz.at(iSol)+nu2pz.at(iSol)+lepton1Pz_+lepton2Pz_+bJet1Pz_+bJet2Pz_ ;
          if(thisPz < minPz)
            {
              minPz = thisPz;
              minPzSol = iSol;
            }
        }

      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);

      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          double thisPz = nu1pz.at(iSol)+nu2pz.at(iSol)+lepton1Pz_+lepton2Pz_+bJet1Pz_+bJet2Pz_ ;
          if(thisPz < minPz)
            {
              minPz = thisPz;
              minPzSol = iSol;
	      bestSol1 = false;
            }
        }
      
      if(minPzSol > -1)
        {
	  if(bestSol1)
	    {
	      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
					       nu2E, nu2px, nu2py, nu2pz);
	    }
	  double px = reconstructed_bJet1Pt_*cos(reconstructed_bJet1Phi_);
	  double py = reconstructed_bJet1Pt_*sin(reconstructed_bJet1Phi_);
	  double pz = reconstructed_bJet1Pt_*sinh(reconstructed_bJet1Eta_);
	  double E = sqrt(reconstructed_bJet1Mass2_ + px*px + py*py + pz*pz) ;
          bJet1.SetPxPyPzE(px,py,pz,E);
	  px = reconstructed_bJet2Pt_*cos(reconstructed_bJet2Phi_);
          py = reconstructed_bJet2Pt_*sin(reconstructed_bJet2Phi_);
          pz = reconstructed_bJet2Pt_*sinh(reconstructed_bJet2Eta_);
          E = sqrt(reconstructed_bJet2Mass2_ + px*px + py*py + pz*pz) ;
          bJet2.SetPxPyPzE(px,py,pz,E);
          lightJet1 = lightJets_.at(0);
          lightJet2 = lightJets_.at(1);
          neutrino1.SetPxPyPzE(nu1px.at(minPzSol),nu1py.at(minPzSol),nu1pz.at(minPzSol),nu1E.at(minPzSol));
          neutrino2.SetPxPyPzE(nu2px.at(minPzSol),nu2py.at(minPzSol),nu2pz.at(minPzSol),nu2E.at(minPzSol));
        }

    }
  else if(solutionFlag1)
    {
      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);
      minPzSol = -1;
      minPz = 1.e99;
      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
	{
	  double thisPz = nu1pz.at(iSol)+nu2pz.at(iSol)+lepton1Pz_+lepton2Pz_+bJet1Pz_+bJet2Pz_ ;
	  if(thisPz < minPz)
	    {
	      minPz = thisPz;
	      minPzSol = iSol;
	    }
	}
      if(minPzSol > -1)
	{
	  //bJet1.SetPxPyPzE(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_);
	  //bJet2.SetPxPyPzE(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_);
	  double px = reconstructed_bJet1Pt_*cos(reconstructed_bJet1Phi_);
          double py = reconstructed_bJet1Pt_*sin(reconstructed_bJet1Phi_);
          double pz = reconstructed_bJet1Pt_*sinh(reconstructed_bJet1Eta_);
          double E = sqrt(reconstructed_bJet1Mass2_ + px*px + py*py + pz*pz) ;
          bJet1.SetPxPyPzE(px,py,pz,E);
          px = reconstructed_bJet2Pt_*cos(reconstructed_bJet2Phi_);
          py = reconstructed_bJet2Pt_*sin(reconstructed_bJet2Phi_);
          pz = reconstructed_bJet2Pt_*sinh(reconstructed_bJet2Eta_);
          E = sqrt(reconstructed_bJet2Mass2_ + px*px + py*py + pz*pz) ;
          bJet2.SetPxPyPzE(px,py,pz,E);
	  lightJet1 = lightJets_.at(0);
	  lightJet2 = lightJets_.at(1);
	  neutrino1.SetPxPyPzE(nu1px.at(minPzSol),nu1py.at(minPzSol),nu1pz.at(minPzSol),nu1E.at(minPzSol));
	  neutrino2.SetPxPyPzE(nu2px.at(minPzSol),nu2py.at(minPzSol),nu2pz.at(minPzSol),nu2E.at(minPzSol));
	}
    }
  else if(solutionFlag2)
    {
      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);
      minPzSol = -1;
      minPz = 1.e99;
      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          double thisPz = nu1pz.at(iSol)+nu2pz.at(iSol)+lepton1Pz_+lepton2Pz_+bJet1Pz_+bJet2Pz_ ;
          if(thisPz < minPz)
            {
              minPz = thisPz;
              minPzSol = iSol;
            }
        }
      if(minPzSol > -1)
        {
          //bJet1.SetPxPyPzE(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_);
          //bJet2.SetPxPyPzE(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_);
	  double px = reconstructed_bJet1Pt_*cos(reconstructed_bJet1Phi_);
          double py = reconstructed_bJet1Pt_*sin(reconstructed_bJet1Phi_);
          double pz = reconstructed_bJet1Pt_*sinh(reconstructed_bJet1Eta_);
          double E = sqrt(reconstructed_bJet1Mass2_ + px*px + py*py + pz*pz) ;
          bJet1.SetPxPyPzE(px,py,pz,E);
          px = reconstructed_bJet2Pt_*cos(reconstructed_bJet2Phi_);
          py = reconstructed_bJet2Pt_*sin(reconstructed_bJet2Phi_);
          pz = reconstructed_bJet2Pt_*sinh(reconstructed_bJet2Eta_);
          E = sqrt(reconstructed_bJet2Mass2_ + px*px + py*py + pz*pz) ;
          bJet2.SetPxPyPzE(px,py,pz,E);
          lightJet1 = lightJets_.at(0);
          lightJet2 = lightJets_.at(1);
          neutrino1.SetPxPyPzE(nu1px.at(minPzSol),nu1py.at(minPzSol),nu1pz.at(minPzSol),nu1E.at(minPzSol));
          neutrino2.SetPxPyPzE(nu2px.at(minPzSol),nu2py.at(minPzSol),nu2pz.at(minPzSol),nu2E.at(minPzSol));
	}
    }
  else if( chi2_ > 0. )
    {
      buildBestNeutrinos();
      buildBestLightJets();
      bJet1.SetPxPyPzE(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_);
      bJet2.SetPxPyPzE(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_);
      lightJet1 = lightJetsBest_.at(0);
      lightJet2 = lightJetsBest_.at(1);
      neutrino1 = nu1Best_;
      neutrino2 = nu2Best_;
      MET.SetPxPyPzE(METx_,METy_,0.,sqrt(METx_*METx_+METy_*METy_));
      pairingInfo = pairingInfoBest_;
      minDeltas.at(0) = bJet1PtDelta_;
      minDeltas.at(1) = bJet1PhiDelta_;
      minDeltas.at(2) = bJet1EtaDelta_;
      minDeltas.at(3) = bJet2PtDelta_;
      minDeltas.at(4) = bJet2PhiDelta_;
      minDeltas.at(5) = bJet2EtaDelta_;
      minDeltas.at(6) = deltaMTop1_;
      minDeltas.at(7) = deltaMTop2_;
      minDeltas.at(8) = deltaMW1_;
      minDeltas.at(9) = deltaMW2_;
      minDeltas.at(10) = lightJetChiSquare_.getMinDeltasX()->at(0);
      minDeltas.at(11) = lightJetChiSquare_.getMinDeltasY()->at(0);
      minDeltas.at(12) = lightJetChiSquare_.getMinDeltasX()->at(1);
      minDeltas.at(13) = lightJetChiSquare_.getMinDeltasY()->at(1);
    }
  //cout << "b-jet 1 px,py,pz,E " << bJet1.Px() << " " << bJet1.Py() << " " << bJet1.Pz() << " " << bJet1.E() << endl;
  //cout << "b-jet 2 px,py,pz,E " << bJet2.Px() << " " << bJet2.Py() << " " << bJet2.Pz() << " " << bJet2.E() << endl;
  //cout << "light jet 1 px,py,pz,E " << lightJet1.Px() << " " << lightJet1.Py() << " " << lightJet1.Pz() << " " << lightJet1.E() << endl;
  //cout << "light jet 2 px,py,pz,E " << lightJet2.Px() << " " << lightJet2.Py() << " " << lightJet2.Pz() << " " << lightJet2.E() << endl;
  //cout << "neutrino 1 px,py,pz,E " << neutrino1.Px() << " " << neutrino1.Py() << " " << neutrino1.Pz() << " " << neutrino1.E() << endl;
  //cout << "neutrino 2 px,py,pz,E " << neutrino2.Px() << " " << neutrino2.Py() << " " << neutrino2.Pz() << " " << neutrino2.E() << endl;
}

void topSystemChiSquare::fillBestMomenta(XYZTLorentzVector& bJet1, XYZTLorentzVector& bJet2,
					 XYZTLorentzVector& lightJet1, XYZTLorentzVector& lightJet2,
					 vector<XYZTLorentzVector>& neutrino1, vector<int>& pairing1,
					 vector<XYZTLorentzVector>& neutrino2, vector<int>& pairing2,
					 XYZTLorentzVector& MET,
					 bool& pairingInfo,
					 vector<double>& minDeltas
					 )
{
  bJet1.SetPxPyPzE(0.,0.,0.,0.);
  bJet2.SetPxPyPzE(0.,0.,0.,0.);
  lightJet1.SetPxPyPzE(0.,0.,0.,0.);
  lightJet2.SetPxPyPzE(0.,0.,0.,0.);
  neutrino1.clear();
  neutrino2.clear();
  pairing1.clear();
  pairing2.clear();
  MET.SetPxPyPzE(0.,0.,0.,0.);
  pairingInfo = pairingInfoBest_;
  minDeltas.clear();
  minDeltas.assign(14,0.);

  pairingInfo_ = true;
  bool solutionFlag1=!calcNeutrinoSolutions();
  pairingInfo_ = false;
  bool solutionFlag2=!calcNeutrinoSolutions();
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;

  if(solutionFlag1&&solutionFlag2)
    {
      //cout << "Intersecting ellipses for both pairings" << endl;
      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);
      
      //cout << "There are " << nu1E.size() << " solutions to the quartic equation for the correct pairing" << endl;
      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          //cout << "Adding a solution for the correct pairing" << endl;
          //cout << "Neutrino 1: " << nu1px.at(iSol) << " " << nu1py.at(iSol) << " " << nu1pz.at(iSol) << " " << nu1E.at(iSol) << endl;
          //cout << "Neutrino 2: " << nu2px.at(iSol) << " " << nu2py.at(iSol) << " " << nu2pz.at(iSol) << " " << nu2E.at(iSol) << endl;
	  neutrino1.push_back(XYZTLorentzVector(nu1px.at(iSol),nu1py.at(iSol),nu1pz.at(iSol),nu1E.at(iSol)));
	  pairing1.push_back(1);
	  neutrino2.push_back(XYZTLorentzVector(nu2px.at(iSol),nu2py.at(iSol),nu2pz.at(iSol),nu2E.at(iSol)));
          pairing2.push_back(1);
        }

      nu1px.clear();
      nu1py.clear();
      nu1pz.clear();
      nu1E.clear();
      nu2px.clear();
      nu2py.clear();
      nu2pz.clear();
      nu2E.clear();

      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);

      //cout << "There are " << nu1E.size() << " solutions to the quartic equation for the incorrect pairing" << endl;
      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          //cout << "Adding a solution for the incorrect pairing" << endl;
          //cout << "Neutrino 1: " << nu1px.at(iSol) << " " << nu1py.at(iSol) << " " << nu1pz.at(iSol) << " " << nu1E.at(iSol) << endl;
          //cout << "Neutrino 2: " << nu2px.at(iSol) << " " << nu2py.at(iSol) << " " << nu2pz.at(iSol) << " " << nu2E.at(iSol) << endl;
	  neutrino1.push_back(XYZTLorentzVector(nu1px.at(iSol),nu1py.at(iSol),nu1pz.at(iSol),nu1E.at(iSol)));
          pairing1.push_back(0);
          neutrino2.push_back(XYZTLorentzVector(nu2px.at(iSol),nu2py.at(iSol),nu2pz.at(iSol),nu2E.at(iSol)));
          pairing2.push_back(0);
        }
    }
  else if(solutionFlag1)
    {
      //cout << "Intersecting ellipses for the correct pairing" << endl;
      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);

      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
	  //cout << "Adding a solution for the correct pairing" << endl;
          neutrino1.push_back(XYZTLorentzVector(nu1px.at(iSol),nu1py.at(iSol),nu1pz.at(iSol),nu1E.at(iSol)));
          pairing1.push_back(1);
          neutrino2.push_back(XYZTLorentzVector(nu2px.at(iSol),nu2py.at(iSol),nu2pz.at(iSol),nu2E.at(iSol)));
          pairing2.push_back(1);
        }
    }
  else if(solutionFlag2)
    {
      //cout << "Intersecting ellipses for the incorrect pairing" << endl;
      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
                                       nu2E, nu2px, nu2py, nu2pz);

      for(unsigned int iSol=0; iSol<nu1E.size(); iSol++)
        {
          neutrino1.push_back(XYZTLorentzVector(nu1px.at(iSol),nu1py.at(iSol),nu1pz.at(iSol),nu1E.at(iSol)));
          pairing1.push_back(0);
          neutrino2.push_back(XYZTLorentzVector(nu2px.at(iSol),nu2py.at(iSol),nu2pz.at(iSol),nu2E.at(iSol)));
          pairing2.push_back(0);
        }
    }
  else if( chi2_ > 0. )
    {
      //cout << "No intersecting ellipses, doing a minimization" << endl;
      buildBestNeutrinos();
      buildBestLightJets();
      bJet1.SetPxPyPzE(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_);
      bJet2.SetPxPyPzE(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_);
      lightJet1 = lightJetsBest_.at(0);
      lightJet2 = lightJetsBest_.at(1);
      //cout << "Adding a solution for the minimization" << endl;
      neutrino1.push_back(nu1Best_);
      neutrino2.push_back(nu2Best_);
      pairing1.push_back(int(pairingInfoBest_));
      pairing2.push_back(int(pairingInfoBest_));
      MET.SetPxPyPzE(METx_,METy_,0.,sqrt(METx_*METx_+METy_*METy_));
      pairingInfo = pairingInfoBest_;
      minDeltas.at(0) = bJet1PtDelta_;
      minDeltas.at(1) = bJet1PhiDelta_;
      minDeltas.at(2) = bJet1EtaDelta_;
      minDeltas.at(3) = bJet2PtDelta_;
      minDeltas.at(4) = bJet2PhiDelta_;
      minDeltas.at(5) = bJet2EtaDelta_;
      minDeltas.at(6) = deltaMTop1_;
      minDeltas.at(7) = deltaMTop2_;
      minDeltas.at(8) = deltaMW1_;
      minDeltas.at(9) = deltaMW2_;
      minDeltas.at(10) = lightJetChiSquare_.getMinDeltasX()->at(0);
      minDeltas.at(11) = lightJetChiSquare_.getMinDeltasY()->at(0);
      minDeltas.at(12) = lightJetChiSquare_.getMinDeltasX()->at(1);
      minDeltas.at(13) = lightJetChiSquare_.getMinDeltasY()->at(1);
    }

  int maxSol = 4;
  for( int i = int(neutrino1.size()) ; i < 2*maxSol ; i++ )
    {
      neutrino1.push_back(XYZTLorentzVector(0.,0.,0.,0.));
      pairing1.push_back(-1);
    }
  for( int i = int(neutrino2.size()) ; i < 2*maxSol ; i++ )
    {
      neutrino2.push_back(XYZTLorentzVector(0.,0.,0.,0.));
      pairing2.push_back(-1);
    }
  if( int(neutrino1.size()) != 2*maxSol || neutrino1.size() != pairing1.size() ||
      int(neutrino2.size()) != 2*maxSol || neutrino2.size() != pairing2.size() )
    {
      cout << "Error occurred filling neutrino solution vectors!!!" << endl;
      cout << "There are " << neutrino1.size() << " neutrinos" << endl;
      cout << "There are " << neutrino2.size() << " anti-neutrinos" << endl;
      cout << "There are " << pairing1 .size() << " neutrino truth pairings" << endl;
      cout << "There are " << pairing2 .size() << " anti-neutrino truth pairings" << endl;
      neutrino1.clear(); 
      neutrino2.clear();
      pairing1.clear();
      pairing2.clear();
      neutrino1.assign(2*maxSol,XYZTLorentzVector(0.,0.,0.,0.));
      neutrino2.assign(2*maxSol,XYZTLorentzVector(0.,0.,0.,0.));
      pairing1.assign(2*maxSol,-1);
      pairing2.assign(2*maxSol,-1);
    }
}

void topSystemChiSquare::plotEllipses(TString plotName)
{
  int nPoints(2000);
  TGraph ellipse1(nPoints+1);
  TGraph ellipse2(nPoints+1);
  TMatrixD *nuEllipse1, *nuEllipse2;
  if(pairingInfoBest_)
    {
      nuEllipse1 = nuEllipseOneOne_;
      nuEllipse2 = nuEllipseTwoTwo_;
    }
  else
    {
      nuEllipse1 = nuEllipseOneTwo_;
      nuEllipse2 = nuEllipseTwoOne_;
    }

  double thisTheta = 0.;
  double theta1Array[3] = {1.,0.,1.}; 
  double theta2Array[3] = {1.,0.,1.}; 

  TVectorD nu1Perp(3,theta1Array);
  TVectorD nu2Perp(3,theta2Array);
  
  nu1Perp*=*nuEllipse1;
  nu2Perp*=*nuEllipse2;

  double nu1X(nu1Perp[0]), nu1Y(nu1Perp[1]);
  double nu2X(METx_ - nu2Perp[0]), nu2Y(METy_ - nu2Perp[1]);

  ellipse1.SetPoint(0,nu1X,nu1Y);
  ellipse2.SetPoint(0,nu2X,nu2Y);

  double thisMaxX(max(nu1X,nu2X)), thisMinX(min(nu1X,nu2X));
  double thisMaxY(max(nu1Y,nu2Y)), thisMinY(min(nu1Y,nu2Y));

  double maxX(thisMaxX),minX(thisMinX),maxY(thisMaxY),minY(thisMinY);

  double twoPiOverN = 2.*3.14159265359/(double)nPoints;
  for(int iPoint = 1; iPoint <= nPoints; iPoint++)
    {
      thisTheta = (double)iPoint * twoPiOverN;
      nu1Perp[0] = cos(thisTheta);
      nu1Perp[1] = sin(thisTheta);
      nu1Perp[2] = 1.;
      nu2Perp[0] = cos(thisTheta);
      nu2Perp[1] = sin(thisTheta);
      nu2Perp[2] = 1.;

      nu1Perp*=*nuEllipse1;
      nu2Perp*=*nuEllipse2;

      nu1X = nu1Perp[0];
      nu1Y = nu1Perp[1];
      nu2X = METx_ - nu2Perp[0];
      nu2Y = METy_ - nu2Perp[1];
      ellipse1.SetPoint(iPoint,nu1X,nu1Y);
      ellipse2.SetPoint(iPoint,nu2X,nu2Y);
      thisMaxX = max(nu1X,nu2X);
      thisMinX = min(nu1X,nu2X);
      thisMaxY = max(nu1Y,nu2Y);
      thisMinY = min(nu1Y,nu2Y);
      maxX = max(maxX,thisMaxX);
      maxY = max(maxY,thisMaxY);
      minX = min(minX,thisMinX);
      minY = min(minY,thisMinY);
    }

  //cout << "X is in ( " << minX << " , " << maxX << " )\n"
  //     << "Y is in ( " << minY << " , " << maxY << " )\n";
  //
  //cout << "setting X range to (" << minX - 0.05*(maxX-minX) << " , " << maxX + 0.05*(maxX-minX) << " )\n";
  //cout << "setting Y range to (" << minY - 0.05*(maxY-minY) << " , " << maxY + 0.05*(maxY-minY) << " )\n";

  //ellipse1.GetXaxis()->SetRangeUser(minX - 0.05*(maxX-minX),maxX + 0.05*(maxX-minX));
  //ellipse1.GetYaxis()->SetRangeUser(minY - 0.05*(maxY-minY),maxY + 0.05*(maxY-minY));  
  //ellipse2.GetXaxis()->SetRangeUser(minX - 0.05*(maxX-minX),maxX + 0.05*(maxX-minX));
  //ellipse2.GetYaxis()->SetRangeUser(minY - 0.05*(maxY-minY),maxY + 0.05*(maxY-minY));

  //ellipse1.Print("v");

  ellipse1.SetLineColor(kRed);
  ellipse2.SetLineColor(kBlue);

  nu1Perp[0] = cos(theta1Best_);
  nu1Perp[1] = sin(theta1Best_);
  nu1Perp[2] = 1.;
  nu2Perp[0] = cos(theta2Best_);
  nu2Perp[1] = sin(theta2Best_);
  nu2Perp[2] = 1.;

  nu1Perp*=*nuEllipse1;
  nu2Perp*=*nuEllipse2;

  TGraph points(2);
  points.SetPoint(0,nu1Perp[0], nu1Perp[1]);
  points.SetPoint(1,METx_ - nu2Perp[0],METy_ - nu2Perp[1]);
  points.SetMarkerStyle(24);
  points.SetMarkerSize(3);

  TH1D drawer("drawer","drawer",1,minX - 0.05*(maxX-minX),maxX + 0.05*(maxX-minX));

  TCanvas canv(plotName,plotName,800,800);
  drawer.Draw("");
  drawer.SetAxisRange(minY - 0.05*(maxY-minY),maxY + 0.05*(maxY-minY),"Y");
  drawer.Draw("AXIS");
  ellipse1.Draw("LSAME");
  ellipse2.Draw("LSAME");
  points.Draw("PSAME");

  canv.SaveAs(plotName+".pdf");

}

bool topSystemChiSquare::checkMasses()
{
  XYZTLorentzVector bJet1, bJet2, lightJet1, lightJet2, neutrino1, neutrino2, MET, W1, W2, top1, top2;
  bool pairingInfo;
  vector<double> minDeltas;
  fillBestMomenta(bJet1,bJet2,lightJet1,lightJet2,neutrino1,neutrino2,MET,pairingInfo,minDeltas);

  if( pairingInfo)
    {
      W1 = lepton1LorentzVector_ + neutrino1;
      W2 = lepton2LorentzVector_ + neutrino2;
      top1 = bJet1 + W1;
      top2 = bJet2 + W2;
    }
  else
    {
      W1 = lepton2LorentzVector_ + neutrino1;
      W2 = lepton1LorentzVector_ + neutrino2;
      top1 = bJet1 + W1;
      top2 = bJet2 + W2;
    }

  cout << "nu1 and W1 and top1 masses: " << neutrino1.M() << " " << W1.M() << " " << top1.M() << endl;
  cout << "nu1 and W2 and top2 masses: " << neutrino2.M() << " " << W2.M() << " " << top2.M() << endl;

  //if( double(W1.M()) != mW_ || double(W2.M()) != mW_ || double(top1.M()) != mTop_ || double(top2.M()) != mTop_ ) 
  //  {
  //    cout << "At least one of the W or top masses are wrong!" << endl;
  //    return false;
  //  }
  //
  //return true;

  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;

  setupNeutrinoSolutions();

  
  nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				   nu2E, nu2px, nu2py, nu2pz);

  unsigned int sizeFirstSolutions = nu1E.size();

  cout << "same pairing (" << sizeFirstSolutions << ") of them" << endl;

  for(unsigned int iNeutrino = 0; iNeutrino < sizeFirstSolutions; iNeutrino++)
    {
      neutrino1.SetPxPyPzE(nu1px[iNeutrino],nu1py[iNeutrino],nu1pz[iNeutrino],nu1E[iNeutrino]);
      neutrino2.SetPxPyPzE(nu2px[iNeutrino],nu2py[iNeutrino],nu2pz[iNeutrino],nu2E[iNeutrino]);
      W1 = lepton1LorentzVector_ + neutrino1;
      W2 = lepton2LorentzVector_ + neutrino2;
      top1 = bJet1 + W1;
      top2 = bJet2 + W2;
      cout << "for solution " << iNeutrino << endl;
      cout << "nu1 and W1 and top1 masses: " << neutrino1.M() << " " << W1.M() << " " << top1.M() << endl;
      cout << "nu1 and W2 and top2 masses: " << neutrino2.M() << " " << W2.M() << " " << top2.M() << endl;
    }

  nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				   nu2E, nu2px, nu2py, nu2pz);

  cout << "opposite pairing (" << nu1E.size() - sizeFirstSolutions << ") of them" << endl;
  
  for(unsigned int iNeutrino = sizeFirstSolutions; iNeutrino < nu1E.size(); iNeutrino++)
    {
      neutrino1.SetPxPyPzE(nu1px[iNeutrino],nu1py[iNeutrino],nu1pz[iNeutrino],nu1E[iNeutrino]);
      neutrino2.SetPxPyPzE(nu2px[iNeutrino],nu2py[iNeutrino],nu2pz[iNeutrino],nu2E[iNeutrino]);
      W1 = lepton2LorentzVector_ + neutrino1;
      W2 = lepton1LorentzVector_ + neutrino2;
      top1 = bJet1 + W1;
      top2 = bJet2 + W2;
      cout << "for solution " << iNeutrino << endl;
      cout << "nu1 and W1 and top1 masses: " << neutrino1.M() << " " << W1.M() << " " << top1.M() << endl;
      cout << "nu1 and W2 and top2 masses: " << neutrino2.M() << " " << W2.M() << " " << top2.M() << endl;
    }
  return true;
}

double topSystemChiSquare::breitWignerError(double& mass , double& width, const double& deltaMass)
{
  double scaledDeltaMass = deltaMass*width;
  double scaledDeltaMass2 = scaledDeltaMass*scaledDeltaMass;
  double Zscore = normal_quantile(0.31830988618379067154*atan2(scaledDeltaMass2 + 2.*scaledDeltaMass*mass,mass*width)+0.5,1.0);
  return Zscore*Zscore;
}

int topSystemChiSquare::getMinimizerStatusCode()
{
  //Minimizer status codes from http://root.cern.ch/root/html/ROOT__Minuit2__Minuit2Minimizer.html#ROOT__Minuit2__Minuit2Minimizer:Minimize
  int status = -1e9;
  if(intersectingEllipses_) status = -1;
  else
    {
      if(outerMin_) status = outerMin_->Status(); 
      {
	if(status>5) status = -2; 
	if(status<0) status = -3;
	if(status!=status) status = -4;
      }
    }
  //cout << "The minimizer status code is " << status << endl;
  return status; 
}
