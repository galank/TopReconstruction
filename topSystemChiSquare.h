#ifndef TOPSYSTEMCHISQUARE
#define TOPSYSTEMCHISQUARE

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/GenVector/LorentzVector.h"

#include "neutrinoSolutions.h"
#include "NeutrinoEllipseCalculator.h"
#include "lightJetChiSquareMinimumSolver.h"

using namespace std;
using namespace ROOT::Math;
typedef LorentzVector<PxPyPzE4D<double> > XYZTLorentzVector;


class topSystemChiSquare
{

 private:

  int bJet1_, lightJet1_;
  int bJet2_, lightJet2_;

  double mTop_, mW_, mNu_;

  double sigmaMTop_, sigmaMW_;

  double deltaMTop1_, deltaMTop2_;
  double deltaMW1_, deltaMW2_;

  double deltaMTop1Best_, deltaMTop2Best_;
  double deltaMW1Best_, deltaMW2Best_;

  double deltaMTop1BestSamePairing_, deltaMTop2BestSamePairing_;
  double deltaMW1BestSamePairing_, deltaMW2BestSamePairing_;

  double deltaMTop1BestOppositePairing_, deltaMTop2BestOppositePairing_;
  double deltaMW1BestOppositePairing_, deltaMW2BestOppositePairing_;

  XYZTLorentzVector measuredMET_, MET_;

  double measuredMETx_, METx_;
  double measuredMETy_, METy_;
  double measuredMETphi_, METphi_;

  XYZTLorentzVector bJet1LorentzVector_, lightJet1LorentzVector_, lepton1LorentzVector_;
  XYZTLorentzVector bJet2LorentzVector_, lightJet2LorentzVector_, lepton2LorentzVector_;

  double bJet1Px_, bJet1Py_, bJet1Pz_, bJet1E_, lepton1Px_, lepton1Py_, lepton1Pz_, lepton1E_;
  double bJet2Px_, bJet2Py_, bJet2Pz_, bJet2E_, lepton2Px_, lepton2Py_, lepton2Pz_, lepton2E_;

  double reconstructed_bJet1Pt_, reconstructed_bJet1Phi_,reconstructed_bJet1Eta_, reconstructed_bJet1Mass2_;
  double reconstructed_bJet2Pt_, reconstructed_bJet2Phi_,reconstructed_bJet2Eta_, reconstructed_bJet2Mass2_;

  double angleAdjusted_bJet1Px_, angleAdjusted_bJet1Py_, angleAdjusted_bJet1Pz_, angleAdjusted_bJet1P_, angleAdjusted_bJet1E_;
  double angleAdjusted_bJet2Px_, angleAdjusted_bJet2Py_, angleAdjusted_bJet2Pz_, angleAdjusted_bJet2P_, angleAdjusted_bJet2E_;

  vector<XYZTLorentzVector> jets_, lightJets_;
  //vector<XYZTLorentzVector> jetWidths_, lightJetWidths_;
  vector<double> jetPtWidths_;
  vector<double> jetPhiWidths_;
  vector<double> jetEtaWidths_;

  vector<XYZTLorentzVector> lightJetsBest_;

  XYZTLorentzVector nu1Best_, nu2Best_;

  double bJet1PtWidth_, bJet1PhiWidth_, bJet1EtaWidth_;
  double bJet2PtWidth_, bJet2PhiWidth_, bJet2EtaWidth_;
  vector<double> lightJetPts_ , lightJetPtWidths_ , lightJetPtDeltas_ ;
  vector<double> lightJetPhis_, lightJetPhiWidths_, lightJetPhiDeltas_;

  double bJet1PtDelta_, bJet1PhiDelta_, bJet1EtaDelta_;
  double bJet2PtDelta_, bJet2PhiDelta_, bJet2EtaDelta_;

  double bJet1PtDeltaBest_, bJet1PhiDeltaBest_, bJet1EtaDeltaBest_;
  double bJet2PtDeltaBest_, bJet2PhiDeltaBest_, bJet2EtaDeltaBest_;

  double bJet1PtDeltaBestSamePairing_, bJet1PhiDeltaBestSamePairing_, bJet1EtaDeltaBestSamePairing_;
  double bJet2PtDeltaBestSamePairing_, bJet2PhiDeltaBestSamePairing_, bJet2EtaDeltaBestSamePairing_;

  double bJet1PtDeltaBestOppositePairing_, bJet1PhiDeltaBestOppositePairing_, bJet1EtaDeltaBestOppositePairing_;
  double bJet2PtDeltaBestOppositePairing_, bJet2PhiDeltaBestOppositePairing_, bJet2EtaDeltaBestOppositePairing_;

  neutrinoSolutions nuSolOne_, nuSolTwo_;
  
  NeutrinoEllipseCalculator nuEllipseOneOneCalc_, nuEllipseOneTwoCalc_;
  NeutrinoEllipseCalculator nuEllipseTwoOneCalc_, nuEllipseTwoTwoCalc_;

  TMatrixD *nuEllipseOneOne_, *nuEllipseOneTwo_;
  TMatrixD *nuEllipseTwoOne_, *nuEllipseTwoTwo_;

  int nOneOneRanges_, nOneTwoRanges_, nTwoOneRanges_, nTwoTwoRanges_, nBJet1Ranges_, nBJet2Ranges_;

  pair<pair<double, bool>,pair<double, bool> > bJet1LogSFRangeOneOne_;
  pair<pair<double, bool>,pair<double, bool> > bJet1LogSFRangeOneTwo_;
  pair<pair<double, bool>,pair<double, bool> > bJet2LogSFRangeTwoOne_;
  pair<pair<double, bool>,pair<double, bool> > bJet2LogSFRangeTwoTwo_;
  pair<pair<double, bool>,pair<double, bool> > bJet1LogSFRange_;
  pair<pair<double, bool>,pair<double, bool> > bJet2LogSFRange_;

  double lightJetChi2_, bJetChi2_, chi2_;
  double chi2BestSamePairing_, chi2BestOppositePairing_;
  double currentBestChi2_, currentBestbJet1PtDelta_, currentBestbJet2PtDelta_, currentBestTheta1_, currentBestTheta2_;
  double startingValueChi2_;
  double deltaMTopBest_,deltaMWBest_,bJetPhiDeltaBest_,bJetEtaDeltaBest_;

  bool pairingInfo_,pairingInfoBest_;
  bool whichEllipse_;

  //ClosestApproach *closestApproachOne_, *closestApproachTwo_;
  double d1_, theta1_, d2_, theta2_, d_, theta_,dx_,dy_;
  double dxBest_,dyBest_,theta1Best_,theta2Best_;

  double theta1BestSamePairing_,theta2BestSamePairing_;
  double theta1BestOppositePairing_,theta2BestOppositePairing_;

  lightJetChiSquareMinimumSolver lightJetChiSquare_;

  double maxConsideredChiSquareRoot_;

  bool intersectingEllipses_;
  int intersectingEllipsesChi2_;

  ROOT::Math::Minimizer* outerMin_;
  ROOT::Math::Minimizer* ptMin_;
  ROOT::Math::Minimizer* ellipseAngleMin_;

  ROOT::Math::Minimizer* startingOuterMin_;
  ROOT::Math::Minimizer* startingTopMassMin_;

  double deltaMTopRangeLow_, deltaMTopRangeHigh_;

  void setBJets();
  void setLightJets();
  void setLightJetCollections();
  void setLeptons();

  //bool sortJets(XYZTLorentzVector , XYZTLorentzVector );

  void resetAll();
  void resetAngles();
  void resetPt();

  void setupNeutrinoSolutions();
  void setupNeutrinoEllipses();
  //void setupClosestApproaches();

  void calcInitialLightJetDeltas(vector<double>& , vector<double>& );
  void buildBestLightJets();
  void calcNeutrinoEllipses();
  void calcNeutrinoRanges();
  bool calcNeutrinoSolutions();
  void getDxDyFromEllipses();
  void buildBestNeutrinos();
  void setBJetAngleDeltas(const double , const double , const double , const double );
  void setBJetPtDeltas(const double , const double );
  void guessBJetPtDeltas();
  void setBestValues();
  void setBestCurrentChiSquare();
  bool setBJetPtRanges();
  void setupTotalChiSquare();
  void setupBJetChiSquare();
  void setupLightJetChiSquare();
  void setStartingBJetAngleDeltas();
  void calcBJetChiSquare();
  void calcTotalChiSquare();
  double getMinBJetPtDelta();
  bool getStartingValues();
  void getStartingTopMassRange();
  double getZ2();
  double breitWignerError(double&,double&,const double&);

 public:

  topSystemChiSquare(vector<XYZTLorentzVector> , vector<double> , vector<double> ,  vector<double> ,
		     int, int, XYZTLorentzVector, 
		     int, int, XYZTLorentzVector,
		     XYZTLorentzVector , double, double, double, double );

  topSystemChiSquare(const topSystemChiSquare& other);

  ~topSystemChiSquare();

  double startingValueOuterMinimizationOperator(const double* );
  double startingValueTopMassMinimizationOperator(const double* );
  double ellipseAngleMinimizationOperator(const double* );
  double bJetPtMinimizationOperator(const double* );
  double outerMinimizationOperator(const double* );

  bool unstretchedBJetPtExists();

  void plotEllipses(TString);
  
  //void minimizeLightJetChiSquare();
  void minimizeBJetPtChiSquare();
  void minimizeLightJetChiSquareOnly();
  void minimizeLightJetChiSquare(int nPoints=10);
  void minimizeTotalChiSquare();
  double getChiSquare();
  int getMinimizerStatusCode(); 

  double calcSolution();

  void fillBestMomenta(XYZTLorentzVector& , XYZTLorentzVector& ,
		       XYZTLorentzVector& , XYZTLorentzVector& ,
		       XYZTLorentzVector& , XYZTLorentzVector& ,
		       XYZTLorentzVector& ,
		       bool& , vector<double>& );
  void fillBestMomenta(XYZTLorentzVector& , XYZTLorentzVector& ,
                       XYZTLorentzVector& , XYZTLorentzVector& ,
                       vector<XYZTLorentzVector>& , vector<int>& ,
		       vector<XYZTLorentzVector>& , vector<int>& ,
		       XYZTLorentzVector& ,
                       bool& , vector<double>& );

  bool checkMasses();
};

#endif
