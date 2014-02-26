#include "topSystemChiSquare.h"


topSystemChiSquare::topSystemChiSquare(vector<math::XYZTLorentzVector> jets,
				     vector<math::XYZTLorentzVector> jetWidths,
				     int bJet1, int lightJet1, math::XYZTLorentzVector lepton1, 
				     int bJet2, int lightJet2, math::XYZTLorentzVector lepton2,
				     math::XYZTLorentzVector METLorentzVector,
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
  chi2_        (0.)
{
  setBJets();
  setLightJets();
  setLightJetCollections();
  setLeptons();
  setupNeutrinoSolutions();
  lightJetChiSquare_=new lightJetChiSquareMinimumSolver(&lightJetPts_,&lightJetPtWidths_,&lightJetPhis_,&lightJetPhiWidths_,&d_,&theta_);
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
  chi2_(other.chi2_)
{
  setBJets();
  setLightJets();
  setLeptons();
  setupNeutrinoSolutions();
}

topSystemChiSquare::~topSystemChiSquare()
{
  delete nuSolOne_;
  delete nuSolTwo_;
  if(nuEllipseOneOneCalc_!=NULL) delete nuEllipseOneOneCalc_;
  if(nuEllipseOneTwoCalc_!=NULL) delete nuEllipseOneTwoCalc_;
  if(nuEllipseTwoOneCalc_!=NULL) delete nuEllipseTwoOneCalc_;
  if(nuEllipseTwoTwoCalc_!=NULL) delete nuEllipseTwoTwoCalc_;
  if(closestApproachOne_!=NULL) delete closestApproachOne_;
  if(closestApproachTwo_!=NULL) delete closestApproachTwo_;
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

  bJet2LorentzVector_=jets_.at(bJet2_);
  bJet2PtWidth_ =jetWidths_.at(bJet2_).Pt() ;
  bJet2PhiWidth_=jetWidths_.at(bJet2_).Phi();

  bJet2Px_=bJet2LorentzVector_.Px();
  bJet2Py_=bJet2LorentzVector_.Py();
  bJet2Pz_=bJet2LorentzVector_.Pz();
  bJet2E_ =bJet2LorentzVector_.E();
}

void topSystemChiSquare::setLightJets()
{
  lightJet1LorentzVector_=jets_.at(lightJet1_);
  lightJet2LorentzVector_=jets_.at(lightJet2_);
}

void topSystemChiSquare::setLightJetCollections()
{
  int iJet=0;
  for(vector<math::XYZTLorentzVector>::iterator thisJet=lightJets_.begin(); thisJet!=lightJets_.end(); thisJet++)
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
  double bJet1PxNew=(1./bJet1PtDelta_)*bJet1LorentzVector_.Pt()*cos(bJet1LorentzVector_.Phi()+bJet1PhiDelta_);
  double bJet1PyNew=(1./bJet1PtDelta_)*bJet1LorentzVector_.Pt()*sin(bJet1LorentzVector_.Phi()+bJet1PhiDelta_);
  double bJet1PzNew=(1./bJet1PtDelta_)*bJet1LorentzVector_.Pt()*sinh(bJet1LorentzVector_.Eta());
  double bJet1ENew =(1./bJet1PtDelta_)*bJet1LorentzVector_.E();
  math::XYZTLorentzVector deltaBJet1;
  deltaBJet1.SetPxPyPzE(bJet1PxNew-bJet1Px_,bJet1PyNew-bJet1Py_,bJet1PzNew-bJet1Pz_,bJet1ENew-bJet1E_);

  double bJet2PxNew=(1./bJet2PtDelta_)*bJet2LorentzVector_.Pt()*cos(bJet2LorentzVector_.Phi()+bJet2PhiDelta_);
  double bJet2PyNew=(1./bJet2PtDelta_)*bJet2LorentzVector_.Pt()*sin(bJet2LorentzVector_.Phi()+bJet2PhiDelta_);
  double bJet2PzNew=(1./bJet2PtDelta_)*bJet2LorentzVector_.Pt()*sinh(bJet2LorentzVector_.Eta());
  double bJet2ENew =(1./bJet2PtDelta_)*bJet2LorentzVector_.E();
  math::XYZTLorentzVector deltaBJet2; 
  deltaBJet2.SetPxPyPzE(bJet2PxNew-bJet2Px_,bJet2PyNew-bJet2Py_,bJet2PzNew-bJet2Pz_,bJet2ENew-bJet2E_);

  //recalculate MET
  MET_+=deltaBJet1+deltaBJet2;
  METx_=MET_.Px();
  METy_=MET_.Py();

  //reset values
  bJet1Px_=bJet1PxNew;
  bJet1Py_=bJet1PyNew;
  bJet1Pz_=bJet1PzNew;
  bJet1E_ =bJet1ENew;
  bJet2Px_=bJet2PxNew;
  bJet2Py_=bJet2PyNew;
  bJet2Pz_=bJet2PzNew;
  bJet2E_ =bJet2ENew;
}

void topSystemChiSquare::setupNeutrinoSolutions()
{
  nuSolOne_ = new neutrinoSolutions(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				    bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
                                    METx_,METy_,mW_,mTop_);
  
  nuSolTwo_ = new neutrinoSolutions(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_,lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_,
                                    bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_,lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_,
				    METx_,METy_,mW_,mTop_);
}

void topSystemChiSquare::setupNeutrinoEllipses()
{
  nuEllipseOneOneCalc_->setBJet(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_); 
  nuEllipseOneTwoCalc_->setBJet(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_); 
  nuEllipseOneOneCalc_->setBJet(bJet1Px_,bJet1Py_,bJet1Pz_,bJet1E_); 
  nuEllipseOneTwoCalc_->setBJet(bJet2Px_,bJet2Py_,bJet2Pz_,bJet2E_); 

  nuEllipseOneOneCalc_->setLepton(lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_);
  nuEllipseOneOneCalc_->setLepton(lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_);
  nuEllipseOneOneCalc_->setLepton(lepton2Px_,lepton2Py_,lepton2Pz_,lepton2E_);
  nuEllipseOneOneCalc_->setLepton(lepton1Px_,lepton1Py_,lepton1Pz_,lepton1E_);

  nuEllipseOneOne_=nuEllipseOneOneCalc_->getHomogeneousNeutrinoEllipse();
  nuEllipseOneTwo_=nuEllipseOneTwoCalc_->getHomogeneousNeutrinoEllipse();
  nuEllipseTwoOne_=nuEllipseTwoOneCalc_->getHomogeneousNeutrinoEllipse();
  nuEllipseTwoTwo_=nuEllipseTwoTwoCalc_->getHomogeneousNeutrinoEllipse();
}

void topSystemChiSquare::setupClosestApproaches()
{
  closestApproachOne_->setNeutrinoOneEllipse(nuEllipseOneOne_);
  closestApproachOne_->setNeutrinoTwoEllipse(nuEllipseOneTwo_);
  closestApproachOne_->getClosestApproach(d1_,theta1_);

  closestApproachTwo_->setNeutrinoOneEllipse(nuEllipseTwoOne_);
  closestApproachTwo_->setNeutrinoTwoEllipse(nuEllipseTwoTwo_);
  closestApproachTwo_->getClosestApproach(d2_,theta2_);

  //Keep the smallest as resulting from the most likely pairing of b-jets and leptons
  d_=min(max(d1_,0.),max(d2_,0.));
  theta_= (d_==d1_) ? theta1_ : (d_==d2_) ? theta2_ : 0.;
}

void topSystemChiSquare::minimizeLightJetChiSquare()
{
  lightJetChiSquare_->minimize();
  lightJetChi2_=lightJetChiSquare_->getChiSquare();
}

void topSystemChiSquare::setBJetDeltas(double bJet1PtDelta, double bJet1PhiDelta, double bJet2PtDelta, double bJet2PhiDelta)
{
  bJet1PtDelta_ =bJet1PtDelta ;
  bJet1PhiDelta_=bJet1PhiDelta;
  bJet2PtDelta_ =bJet2PtDelta ;
  bJet2PhiDelta_=bJet2PhiDelta;
  reset();
}

void topSystemChiSquare::setupBJetChiSquare(double bJet1PtDelta, double bJet1PhiDelta, double bJet2PtDelta, double bJet2PhiDelta)
{
  //Set b-jet deltas
  setBJetDeltas(bJet1PtDelta,bJet1PhiDelta,bJet2PtDelta,bJet2PhiDelta);

  //Recalculate the neutrino ellipses
  setupNeutrinoEllipses();
  setupClosestApproaches();

  //Minimize the light jet chi2 with the new MET
  //Since lightJetChiSquare_ is constructed with pointers, the values are already updated
  minimizeLightJetChiSquare();
}

double topSystemChiSquare::calcBJetChiSquare()
{
  double bJetChi2=
    pow((bJet1PtDelta_-1.)/bJet1PtWidth_,2)+pow((bJet1PhiDelta_-1.)/bJet1PhiWidth_,2)+
    pow((bJet2PtDelta_-1.)/bJet2PtWidth_,2)+pow((bJet2PhiDelta_-1.)/bJet2PhiWidth_,2);
  return bJetChi2;
}

double topSystemChiSquare::operator()(const double* inputDeltas)
{
  double bJet1PtDelta =inputDeltas[0];
  double bJet1PhiDelta=inputDeltas[1];
  double bJet2PtDelta =inputDeltas[2];
  double bJet2PhiDelta=inputDeltas[3];
  setupBJetChiSquare(bJet1PtDelta,bJet1PhiDelta,bJet2PtDelta,bJet2PhiDelta);

  return calcBJetChiSquare();
}

void topSystemChiSquare::minimizeBJetChiSquare()
{
  //Initial choice?
  setBJetDeltas(0,0,0,0);

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.001);
  min->SetPrintLevel(0);
  
  //Set up the functor
  ROOT::Math::Functor func(this,&topSystemChiSquare::operator(),4);

  //Set up the minimization piece:
  min->SetFunction(func);
  const TString bJetPtDeltaString ("bJetPtDelta_" );
  const TString bJetPhiDeltaString("bJetPhiDelta_");
  TString parameterName;
  
  parameterName=bJetPtDeltaString+"_0";;
  min->SetVariable(0,parameterName.Data(),bJet1PtDelta_,0.001);
  parameterName = bJetPhiDeltaString+"_0";
  min->SetVariable(1,parameterName.Data(),bJet1PhiDelta_,0.001);
  parameterName=bJetPtDeltaString+"_1";;
  min->SetVariable(2,parameterName.Data(),bJet2PtDelta_,0.001);
  parameterName = bJetPhiDeltaString+"_1";
  min->SetVariable(3,parameterName.Data(),bJet2PhiDelta_,0.001);

  //Run the minimization
  min->Minimize();
  bJetChi2_=min->MinValue();
  
  delete min; 
  
  //return bJetChi2_;

}

double topSystemChiSquare::getChiSquare()
{
  chi2_=bJetChi2_+lightJetChi2_;
  return chi2_;
}

double topSystemChiSquare::calcSolution()
{
  //First check whether there exist real solutions to the quartic equations
  //In this case, the neutrino ellipses intersect

  bool solutionFlag1=false, solutionFlag2=false;
  
  vector<double> nu1E,  nu1px,  nu1py,  nu1pz,  nu2E,  nu2px,  nu2py,  nu2pz;

  nuSolOne_->getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				    nu2E, nu2px, nu2py, nu2pz);
  if(nu1E.size()>0) solutionFlag1=true;

  nuSolTwo_->getRealNeutrinoVectors(nu1E, nu1px, nu1py, nu1pz,
				    nu2E, nu2px, nu2py, nu2pz);
  if(nu1E.size()>0) solutionFlag2=true;

  if(solutionFlag1||solutionFlag2) return chi2_;


  //If there are no real solutions, the ellipses don't touch
  //We need to minimize the global chi2

  minimizeBJetChiSquare();


  return getChiSquare();
}

//  LocalWords:  setupNeutrinoSolutions
