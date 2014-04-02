#include "topSystemChiSquare.h"


topSystemChiSquare::topSystemChiSquare(vector<XYZTLorentzVector> jets,
				       vector<XYZTLorentzVector> jetWidths,
				       int bJet1, int lightJet1, XYZTLorentzVector lepton1, 
				       int bJet2, int lightJet2, XYZTLorentzVector lepton2,
				       XYZTLorentzVector METLorentzVector,
				       double mTop, double mW) :
  bJet1_    (bJet1),
  lightJet1_(lightJet1),
  bJet2_    (bJet2),
  lightJet2_(lightJet2),
  mTop_     (mTop),
  mW_       (mW),
  mNu_      (0.),
  measuredMET_ (METLorentzVector),
  MET_         (METLorentzVector),
  measuredMETx_(METLorentzVector.Px()),
  METx_        (METLorentzVector.Px()),
  measuredMETy_(METLorentzVector.Py()),
  METy_        (METLorentzVector.Py()),
  lepton1LorentzVector_(lepton1),
  lepton2LorentzVector_(lepton2),
  jets_          (jets),
  lightJets_     (jets),
  jetWidths_     (jetWidths),
  lightJetWidths_(jetWidths),
  lightJetChi2_(0.),
  bJetChi2_    (0.),
  chi2_        (0.),
  dx_          (0.),
  dy_          (0.),
  lightJetChiSquare_(jets.size()-2,dx_,dy_)
{
  setBJets();
  setLightJets();
  setLightJetCollections();
  setLeptons();
  setupNeutrinoSolutions();
  lightJetChiSquare_.setupEquations(lightJetPts_,lightJetPtWidths_,lightJetPhis_,lightJetPhiWidths_);
}

topSystemChiSquare::topSystemChiSquare(const topSystemChiSquare& other) :
  bJet1_    (other.bJet1_),
  lightJet1_(other.lightJet1_),
  bJet2_    (other.bJet2_),
  lightJet2_(other.lightJet2_),
  mTop_     (other.mTop_),
  mW_       (other.mW_),
  mNu_      (other.mNu_),
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
  jetWidths_     (other.jetWidths_),
  lightJetWidths_(other.jetWidths_),
  lightJetChi2_(other.lightJetChi2_),
  bJetChi2_(other.bJetChi2_),
  chi2_(other.chi2_),
  dx_          (other.dx_),
  dy_          (other.dy_),
  lightJetChiSquare_(other.jets_.size()-2,dx_,dy_)
{
  setBJets();
  setLightJets();
  setLeptons();
  setupNeutrinoSolutions();
  lightJetChiSquare_.setupEquations(lightJetPts_,lightJetPtWidths_,lightJetPhis_,lightJetPhiWidths_);
}

topSystemChiSquare::~topSystemChiSquare()
{
}

void topSystemChiSquare::setBJets()
{
  bJet1LorentzVector_=jets_.at(bJet1_);
  bJet1PtWidth_ =jetWidths_.at(bJet1_).Pt() ;
  bJet1PhiWidth_=jetWidths_.at(bJet1_).Phi();

  bJet1Px_=bJet1LorentzVector_.Px();
  bJet1Py_=bJet1LorentzVector_.Py();
  bJet1Pz_=bJet1LorentzVector_.Pz();
  bJet1E_ =bJet1LorentzVector_.E();

  reconstructed_bJet1Pt_ =bJet1LorentzVector_.Pt();
  reconstructed_bJet1Phi_=bJet1LorentzVector_.Phi();

  bJet2LorentzVector_=jets_.at(bJet2_);
  bJet2PtWidth_ =jetWidths_.at(bJet2_).Pt() ;
  bJet2PhiWidth_=jetWidths_.at(bJet2_).Phi();

  bJet2Px_=bJet2LorentzVector_.Px();
  bJet2Py_=bJet2LorentzVector_.Py();
  bJet2Pz_=bJet2LorentzVector_.Pz();
  bJet2E_ =bJet2LorentzVector_.E();

  reconstructed_bJet2Pt_ =bJet2LorentzVector_.Pt();
  reconstructed_bJet2Phi_=bJet2LorentzVector_.Phi();
}

void topSystemChiSquare::setLightJets()
{
  lightJet1LorentzVector_=jets_.at(lightJet1_);
  lightJet2LorentzVector_=jets_.at(lightJet2_);
}

void topSystemChiSquare::setLightJetCollections()
{
  int iJet=0;
  for(vector<XYZTLorentzVector>::iterator thisJet=lightJets_.begin(); thisJet!=lightJets_.end(); thisJet++)
    {
      if( iJet==bJet1_ || iJet==bJet2_ )
	{
	  lightJets_.erase(thisJet);
	  lightJetWidths_.erase(thisJet);
	  continue;
	}
      lightJetPts_ .push_back(lightJets_.at(iJet).Pt() );
      lightJetPhis_.push_back(lightJets_.at(iJet).Phi());
      lightJetPtWidths_ .push_back(lightJetWidths_.at(iJet).Pt() );
      lightJetPhiWidths_.push_back(lightJetWidths_.at(iJet).Phi());
      iJet++;
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
  for(vector<XYZTLorentzVector>::iterator thisJet=lightJets_.begin(); thisJet != lightJets_.end(); thisJet++, thisDeltaX++, thisDeltaY++)
    {
      lightJetsBest_.push_back(XYZTLorentzVector(thisJet->Px()+*thisDeltaX,thisJet->Py()+*thisDeltaY,thisJet->Pz(),thisJet->E()));
    }
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

void topSystemChiSquare::reset()
{
  double bJet1PtScaleFactor = exp(bJet1PtDelta_*bJet1PtWidth_);
  double bJet1Pt  = bJet1PtScaleFactor*reconstructed_bJet1Pt_;
  double bJet1Phi = reconstructed_bJet1Phi_+bJet1PhiDelta_*bJet1PhiWidth_;
  double bJet2PtScaleFactor = exp(bJet2PtDelta_*bJet2PtWidth_);
  double bJet2Pt  = bJet2PtScaleFactor*reconstructed_bJet2Pt_;
  double bJet2Phi = reconstructed_bJet2Phi_+bJet2PhiDelta_*bJet2PhiWidth_;
  double bJet1PxNew = bJet1Pt*sin(bJet1Phi);
  double bJet1PyNew = bJet1Pt*cos(bJet1Phi);
  double bJet2PxNew = bJet2Pt*sin(bJet2Phi);
  double bJet2PyNew = bJet2Pt*cos(bJet2Phi);
  double bJet1ENew  = bJet1PtScaleFactor*bJet1LorentzVector_.E();
  double bJet2ENew  = bJet2PtScaleFactor*bJet2LorentzVector_.E();

  //recalculate MET
  METx_=measuredMETx_+bJet1LorentzVector_.Px()-bJet1PxNew+bJet2LorentzVector_.Px()-bJet2PxNew;
  METy_=measuredMETy_+bJet1LorentzVector_.Py()-bJet1PyNew+bJet2LorentzVector_.Py()-bJet2PyNew;

  //reset values
  bJet1Px_=bJet1PxNew;
  bJet1Py_=bJet1PyNew;
  bJet1E_ =bJet1ENew;
  bJet2Px_=bJet2PxNew;
  bJet2Py_=bJet2PyNew;
  bJet2E_ =bJet2ENew;
}

void topSystemChiSquare::setupNeutrinoSolutions()
{
  nuSolOne_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
			      bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
			      METx_,METy_,mW_,mTop_);
  
  nuSolTwo_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
			      bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
			      METx_,METy_,mW_,mTop_);
}

void topSystemChiSquare::setupNeutrinoEllipses()
{
  nuEllipseOneOneCalc_.setupEllipse(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,
				    lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				    mTop_,mW_,mNu_);

  nuEllipseOneTwoCalc_.setupEllipse(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,
				    lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				    mTop_,mW_,mNu_);

  nuEllipseTwoOneCalc_.setupEllipse(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,
				    lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				    mTop_,mW_,mNu_);

  nuEllipseTwoTwoCalc_.setupEllipse(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,
				    lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				    mTop_,mW_,mNu_);
  nuEllipseOneOneCalc_.calcNeutrinoEllipse();
  nuEllipseOneTwoCalc_.calcNeutrinoEllipse();
  nuEllipseTwoOneCalc_.calcNeutrinoEllipse();
  nuEllipseTwoTwoCalc_.calcNeutrinoEllipse();
  nuEllipseOneOne_=nuEllipseOneOneCalc_.getHomogeneousNeutrinoEllipse();
  nuEllipseOneTwo_=nuEllipseOneTwoCalc_.getHomogeneousNeutrinoEllipse();
  nuEllipseTwoOne_=nuEllipseTwoOneCalc_.getHomogeneousNeutrinoEllipse();
  nuEllipseTwoTwo_=nuEllipseTwoTwoCalc_.getHomogeneousNeutrinoEllipse();
  bJet1LogSFRangeOneOne_=nuEllipseOneOneCalc_.getBJetLogSFRange(nOneOneRanges_);
  bJet1LogSFRangeOneTwo_=nuEllipseOneTwoCalc_.getBJetLogSFRange(nOneTwoRanges_);
  bJet2LogSFRangeTwoOne_=nuEllipseTwoOneCalc_.getBJetLogSFRange(nTwoOneRanges_);
  bJet2LogSFRangeTwoTwo_=nuEllipseTwoTwoCalc_.getBJetLogSFRange(nTwoTwoRanges_);
}

void topSystemChiSquare::calcNeutrinoEllipses()
{
  if(pairingInfo_)
    {
      nuEllipseOneOneCalc_.calcNeutrinoEllipse();
      nuEllipseOneTwoCalc_.calcNeutrinoEllipse();
    }
  else
    {
      nuEllipseTwoOneCalc_.calcNeutrinoEllipse();
      nuEllipseTwoTwoCalc_.calcNeutrinoEllipse();
    }
}

bool topSystemChiSquare::calcNeutrinoSolutions()
{
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;
  if(pairingInfo_)
    {
      nuSolOne_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				  bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				  METx_,METy_,mW_,mTop_);
      nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);
      if(nu1E.size()>0) return true;

    }
  else
    {
      nuSolTwo_.setupMeasurements(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
				  bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				  METx_,METy_,mW_,mTop_);

      nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				       nu2E, nu2px, nu2py, nu2pz);
      if(nu1E.size()>0) return true;
    }
  return false;
}

//void topSystemChiSquare::setupClosestApproaches()
//{
//  closestApproachOne_->setNeutrinoOneEllipse(nuEllipseOneOne_);
//  closestApproachOne_->setNeutrinoTwoEllipse(nuEllipseOneTwo_);
//  closestApproachOne_->getClosestApproach(d1_,theta1_);
//
//  closestApproachTwo_->setNeutrinoOneEllipse(nuEllipseTwoOne_);
//  closestApproachTwo_->setNeutrinoTwoEllipse(nuEllipseTwoTwo_);
//  closestApproachTwo_->getClosestApproach(d2_,theta2_);
//
//  //Keep the smallest as resulting from the most likely pairing of b-jets and leptons
//  d_=min(max(d1_,0.),max(d2_,0.));
//  theta_= (d_==d1_) ? theta1_ : (d_==d2_) ? theta2_ : 0.;
//}

void topSystemChiSquare::setBJetDeltas(double bJet1PtDelta, double bJet1PhiDelta, double bJet2PtDelta, double bJet2PhiDelta)
{
  bJet1PtDelta_ =bJet1PtDelta ;
  bJet1PhiDelta_=bJet1PhiDelta;
  bJet2PtDelta_ =bJet2PtDelta ;
  bJet2PhiDelta_=bJet2PhiDelta;
  reset();
}

void topSystemChiSquare::setupBJetChiSquare()
{

  //Set b-jet deltas
  //setBJetDeltas(bJet1PtDelta,bJet1PhiDelta,bJet2PtDelta,bJet2PhiDelta);
  
  bool noSolutions = calcNeutrinoSolutions();
  if(noSolutions)
    {
      //Recalculate the neutrino ellipses
      calcNeutrinoEllipses();

      //Minimize the light jet chi2 with the new MET
      minimizeLightJetChiSquare();
    }
  else{
    lightJetChi2_ = 0.;
    dxBest_=0;
    dyBest_=0;
  }
  calcBJetChiSquare();
  
}

void topSystemChiSquare::calcBJetChiSquare()
{
  bJetChi2_ = lightJetChi2_ + bJet1PtDelta_*bJet1PtDelta_ + bJet2PtDelta_*bJet2PtDelta_ + bJet1PhiDelta_*bJet1PhiDelta_ + bJet2PhiDelta_*bJet2PhiDelta_ ;
}

void topSystemChiSquare::getDxDyFromEllipses(double theta1, double theta2)
{
  double nu1Array[3] = {cos(theta1),sin(theta1),1.}; 
  double nu2Array[3] = {cos(theta2),sin(theta2),1.}; 

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
  dx_ = nu1Perp[0] + (METx_ - nu2Perp[0]);
  dy_ = nu1Perp[1] + (METy_ - nu2Perp[1]);
  theta1Best_ = theta1;
  theta2Best_ = theta2;
  dxBest_ = dx_;
  dyBest_ = dy_;
}

double topSystemChiSquare::lightJetMinimizationOperator(const double* inputThetas)
{
  double theta1 = inputThetas[0];
  double theta2 = inputThetas[1];

  getDxDyFromEllipses(theta1, theta2);

  return lightJetChiSquare_.getChiSquare();
}

double topSystemChiSquare::bJetMinimizationOperator(const double* inputDeltas)
{
  bJet1PtDelta_ =inputDeltas[0];
  bJet1PhiDelta_=inputDeltas[1];
  bJet2PtDelta_ =inputDeltas[2];
  bJet2PhiDelta_=inputDeltas[3];
  reset();
  setupBJetChiSquare();

  return bJetChi2_;
}

void topSystemChiSquare::minimizeLightJetChiSquare()
{
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.001);
  min->SetPrintLevel(0);
  
  //Set up the functor
  ROOT::Math::Functor func(this,&topSystemChiSquare::lightJetMinimizationOperator,2);

  //Set up the minimization piece:
  min->SetFunction(func);
  const TString ellipse1String ("theta1" );
  const TString ellipse2String ("theta2");
  TString parameterName;

  double theta1(0);
  double theta2(0);
  getDxDyFromEllipses(theta1,theta2);
  double startingChi2 = lightJetChiSquare_.getChiSquare();
  double piOver2 = 1.57079632679;
  for(int iTheta1 = 0; iTheta1 < 4; iTheta1++)
    {
      for(int iTheta2 = 0; iTheta2 < 4; iTheta2++)
	{
	  if(iTheta1 == 0 && iTheta2 == 0) continue;
	  getDxDyFromEllipses((double)iTheta1*piOver2,(double)iTheta2*piOver2);
	  double thisChi2 = lightJetChiSquare_.getChiSquare();
	  if(thisChi2 < startingChi2)
	    {
	      startingChi2 = thisChi2;
	      theta1 = (double)iTheta1*piOver2;
	      theta2 = (double)iTheta2*piOver2;
	    }
	}
    }

  
  min->SetVariable(0,ellipse1String.Data(),theta1,0.001);
  min->SetVariable(1,ellipse2String.Data(),theta2,0.001);

  //Run the minimization
  min->Minimize();
  lightJetChi2_=min->MinValue();
  getDxDyFromEllipses(min->X()[0],min->X()[1]);

  delete min; 
}

void topSystemChiSquare::minimizeBJetChiSquare()
{
  pairingInfo_ = false;
  bJetChi2_ = -1.;
  for(unsigned int i = 0; i<2; i++)
    {
      //Set pairing of b-jets
      if(i == 1) pairingInfo_ = true;

      if(pairingInfo_)
	{
	  nBJet1Ranges_ = nOneOneRanges_;
	  nBJet2Ranges_ = nTwoTwoRanges_;
	  bJet1LogSFRange_ = bJet1LogSFRangeOneOne_;
	  bJet2LogSFRange_ = bJet2LogSFRangeTwoTwo_;
	}
      else
	{
	  nBJet1Ranges_ = nOneTwoRanges_;
	  nBJet2Ranges_ = nTwoOneRanges_;
	  bJet1LogSFRange_ = bJet1LogSFRangeOneTwo_;
	  bJet2LogSFRange_ = bJet2LogSFRangeTwoOne_;
	}

      if( nBJet1Ranges_ == 0 || nBJet2Ranges_ == 0 )
	{
	  continue;
	}

      double bJet1PtDeltaStart(0.);
      if(bJet1LogSFRange_.first.second && bJet1LogSFRange_.second.second)
	{
	  if(bJet1LogSFRange_.first.first > 0 || bJet1LogSFRange_.second.first < 0.)
	    {
	      double width = abs(bJet1LogSFRange_.second.first - bJet1LogSFRange_.first.first);
	      bJet1PtDeltaStart = (abs(bJet1LogSFRange_.first.first) > abs(bJet1LogSFRange_.second.first))? 
		                  bJet1LogSFRange_.first.first  - (bJet1LogSFRange_.first.first /abs(bJet1LogSFRange_.first.first) )*0.01*width:
		                  bJet1LogSFRange_.second.first - (bJet1LogSFRange_.second.first/abs(bJet1LogSFRange_.second.first))*0.01*width;
	    }
	}
      double bJet2PtDeltaStart(0.);
      if(bJet2LogSFRange_.first.second && bJet2LogSFRange_.second.second)
	{
	  if(bJet2LogSFRange_.first.first > 0 || bJet2LogSFRange_.second.first < 0.)
	    {
	      double width = abs(bJet2LogSFRange_.second.first - bJet2LogSFRange_.first.first);
	      bJet2PtDeltaStart = (abs(bJet2LogSFRange_.first.first) > abs(bJet2LogSFRange_.second.first))? 
		                  bJet2LogSFRange_.first.first  - (bJet2LogSFRange_.first.first /abs(bJet2LogSFRange_.first.first) )*0.01*width:
		                  bJet2LogSFRange_.second.first - (bJet2LogSFRange_.second.first/abs(bJet2LogSFRange_.second.first))*0.01*width;
	    }
	}

      setBJetDeltas(bJet1PtDeltaStart,0.,bJet2PtDeltaStart,0.);

      ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
      min->SetMaxFunctionCalls(1000000);
      min->SetTolerance(0.001);
      min->SetPrintLevel(0);
  
      //Set up the functor
      ROOT::Math::Functor func(this,&topSystemChiSquare::bJetMinimizationOperator,4);

      //Set up the minimization piece:
      min->SetFunction(func);
      const TString bJetPtDeltaString ("bJetPtDelta_" );
      const TString bJetPhiDeltaString("bJetPhiDelta_");
      TString parameterName;

      parameterName=bJetPtDeltaString+"_0";
      min->SetVariable(0,parameterName.Data(),bJet1PtDelta_,0.001);
      if( bJet1LogSFRange_.first.second && bJet1LogSFRange_.first.second ) {
	min->SetLimitedVariable(0,parameterName.Data(),bJet1PtDelta_,0.001,bJet1LogSFRange_.first.first,bJet1LogSFRange_.second.first);
      }
      else if( bJet1LogSFRange_.first.second ) {
	min->SetLowerLimitedVariable(0,parameterName.Data(),bJet1PtDelta_,0.001,bJet1LogSFRange_.first.first);
      }
      else if( bJet1LogSFRange_.second.second ) {
	min->SetUpperLimitedVariable(0,parameterName.Data(),bJet1PtDelta_,0.001,bJet1LogSFRange_.second.first);
      }
      parameterName = bJetPhiDeltaString+"_0";
      min->SetVariable(1,parameterName.Data(),bJet1PhiDelta_,0.001);
      parameterName=bJetPtDeltaString+"_1";
      //min->SetVariable(2,parameterName.Data(),bJet2PtDelta_,0.001);
      if( bJet2LogSFRange_.first.second && bJet2LogSFRange_.second.second ) {
	min->SetLimitedVariable(2,parameterName.Data(),bJet2PtDelta_,0.001,bJet2LogSFRange_.first.first,bJet2LogSFRange_.second.first);
      }
      else if( bJet2LogSFRange_.first.second ){
	min->SetLowerLimitedVariable(2,parameterName.Data(),bJet2PtDelta_,0.001,bJet2LogSFRange_.first.first);
      }
      else if( bJet2LogSFRange_.second.second ){
	min->SetUpperLimitedVariable(2,parameterName.Data(),bJet2PtDelta_,0.001,bJet2LogSFRange_.second.first);
      }
      parameterName = bJetPhiDeltaString+"_1";
      min->SetVariable(3,parameterName.Data(),bJet2PhiDelta_,0.001);

      //Run the minimization
      min->Minimize();
      double minValue = min->MinValue();
      if(bJetChi2_ <0) 
	{
	  bJetChi2_=minValue;
	  pairingInfoBest_ = pairingInfo_;
	  bJet1PtDeltaBest_ = min->X()[0];
	  bJet1PhiDeltaBest_= min->X()[1];
	  bJet2PtDeltaBest_ = min->X()[2];
	  bJet2PhiDeltaBest_= min->X()[3];
	}
      else if(bJetChi2_ > minValue)
	{
	  bJetChi2_=minValue;
	  pairingInfoBest_ = pairingInfo_;
	  bJet1PtDeltaBest_ = min->X()[0];
	  bJet1PhiDeltaBest_= min->X()[1];
	  bJet2PtDeltaBest_ = min->X()[2];
	  bJet2PhiDeltaBest_= min->X()[3];
	}
      delete min; 
    }

  pairingInfo_ = pairingInfoBest_;
  bJet1PtDelta_ = bJet1PtDeltaBest_ ;
  bJet1PhiDelta_= bJet1PhiDeltaBest_;
  bJet2PtDelta_ = bJet2PtDeltaBest_ ;
  bJet2PhiDelta_= bJet2PhiDeltaBest_;
  reset();
  setupBJetChiSquare();
  buildBestLightJets();

}


double topSystemChiSquare::getChiSquare()
{
  chi2_=bJetChi2_;
  return chi2_;
}

double topSystemChiSquare::calcSolution()
{
  //First check whether there exist real solutions to the quartic equations
  //In this case, the neutrino ellipses intersect

  bool solutionFlag1=false, solutionFlag2=false;
  
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;

  nuSolOne_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				    nu2E, nu2px, nu2py, nu2pz);
  if(nu1E.size()>0) solutionFlag1=true;

  nuSolTwo_.getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				    nu2E, nu2px, nu2py, nu2pz);
  if(nu1E.size()>0) solutionFlag2=true;
  if(solutionFlag1||solutionFlag2) return chi2_;
  //If there are no real solutions, the ellipses don't touch
  //We need to minimize the global chi2
  minimizeBJetChiSquare();
  return getChiSquare();
}

//  LocalWords:  setupNeutrinoSolutions
