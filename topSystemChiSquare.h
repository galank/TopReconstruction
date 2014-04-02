#ifndef TOPSYSTEMCHISQUARE
#define TOPSYSTEMCHISQUARE

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

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

  XYZTLorentzVector measuredMET_, MET_;

  double measuredMETx_, METx_;
  double measuredMETy_, METy_;
  double measuredMETphi_, METphi_;

  XYZTLorentzVector bJet1LorentzVector_, lightJet1LorentzVector_, lepton1LorentzVector_;
  XYZTLorentzVector bJet2LorentzVector_, lightJet2LorentzVector_, lepton2LorentzVector_;

  double bJet1Px_, bJet1Py_, bJet1Pz_, bJet1E_, lepton1Px_, lepton1Py_, lepton1Pz_, lepton1E_;
  double bJet2Px_, bJet2Py_, bJet2Pz_, bJet2E_, lepton2Px_, lepton2Py_, lepton2Pz_, lepton2E_;

  double reconstructed_bJet1Pt_, reconstructed_bJet1Phi_;
  double reconstructed_bJet2Pt_, reconstructed_bJet2Phi_;

  vector<XYZTLorentzVector> jets_, lightJets_;
  vector<XYZTLorentzVector> jetWidths_, lightJetWidths_;

  vector<XYZTLorentzVector> lightJetsBest_;

  double bJet1PtWidth_, bJet1PhiWidth_;
  double bJet2PtWidth_, bJet2PhiWidth_;
  vector<double> lightJetPts_ , lightJetPtWidths_ , lightJetPtDeltas_ ;
  vector<double> lightJetPhis_, lightJetPhiWidths_, lightJetPhiDeltas_;

  double bJet1PtDelta_, bJet1PhiDelta_;
  double bJet2PtDelta_, bJet2PhiDelta_;

  double bJet1PtDeltaBest_, bJet1PhiDeltaBest_;
  double bJet2PtDeltaBest_, bJet2PhiDeltaBest_;

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

  bool pairingInfo_,pairingInfoBest_;

  //ClosestApproach *closestApproachOne_, *closestApproachTwo_;
  double d1_, theta1_, d2_, theta2_, d_, theta_,dx_,dy_;
  double dxBest_,dyBest_,theta1Best_,theta2Best_;

  lightJetChiSquareMinimumSolver lightJetChiSquare_;

  void setBJets();
  void setLightJets();
  void setLightJetCollections();
  void setLeptons();

  void reset();

  void setupNeutrinoSolutions();
  void setupNeutrinoEllipses();
  //void setupClosestApproaches();

  void calcInitialLightJetDeltas(vector<double>& , vector<double>& );
  void buildBestLightJets();
  void calcNeutrinoEllipses();
  bool calcNeutrinoSolutions();
  void getDxDyFromEllipses(double, double);

 public:

  topSystemChiSquare(vector<XYZTLorentzVector> , vector<XYZTLorentzVector> ,
		    int , int, XYZTLorentzVector , 
		    int , int, XYZTLorentzVector ,
		    XYZTLorentzVector , double , double );

  topSystemChiSquare(const topSystemChiSquare& other);

  ~topSystemChiSquare();

  //double operator()(const double* );
  double lightJetMinimizationOperator(const double* );
  double bJetMinimizationOperator(const double* );

  void setBJetDeltas(double , double , double , double );
  void setupBJetChiSquare();

  void calcBJetChiSquare();
  
  void minimizeLightJetChiSquare();
  void minimizeBJetChiSquare();
  double getChiSquare();

  double calcSolution();
};

#endif
