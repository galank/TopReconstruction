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
#include "DataFormats/Math/interface/LorentzVector.h"

#include "neutrinoSolutions.h"
#include "NeutrinoEllipseCalculator.cxx"
#include "ClosestApproach.cxx"
//#include "lightJetChiSquare.cxx"
#include "lightJetChiSquareMinimumSolver.cxx"

using namespace std;

class topSystemChiSquare
{

 private:

  double bJet1_, lightJet1_;
  double bJet2_, lightJet2_;

  double mTop_, mW_, mNu_;

  math::XYZTLorentzVector measuredMET_, MET_;

  double measuredMETx_, METx_;
  double measuredMETy_, METy_;
  double measuredMETphi_, METphi_;

  math::XYZTLorentzVector bJet1LorentzVector_, lightJet1LorentzVector_, lepton1LorentzVector_;
  math::XYZTLorentzVector bJet2LorentzVector_, lightJet2LorentzVector_, lepton2LorentzVector_;

  double bJet1Px_, bJet1Py_, bJet1Pz_, bJet1E_, lepton1Px_, lepton1Py_, lepton1Pz_, lepton1E_;
  double bJet2Px_, bJet2Py_, bJet2Pz_, bJet2E_, lepton2Px_, lepton2Py_, lepton2Pz_, lepton2E_;

  vector<math::XYZTLorentzVector> jets_, lightJets_;
  vector<math::XYZTLorentzVector> jetWidths_, lightJetWidths_;

  double bJet1PtWidth_, bJet1PhiWidth_;
  double bJet2PtWidth_, bJet2PhiWidth_;
  vector<double> lightJetPts_ , lightJetPtWidths_ , lightJetPtDeltas_ ;
  vector<double> lightJetPhis_, lightJetPhiWidths_, lightJetPhiDeltas_;

  double bJet1PtDelta_, bJet1PhiDelta_;
  double bJet2PtDelta_, bJet2PhiDelta_;

  neutrinoSolutions *nuSolOne_, *nuSolTwo_;
  
  NeutrinoEllipseCalculator *nuEllipseOneOneCalc_, *nuEllipseOneTwoCalc_;
  NeutrinoEllipseCalculator *nuEllipseTwoOneCalc_, *nuEllipseTwoTwoCalc_;

  TMatrixD nuEllipseOneOne_, nuEllipseOneTwo_;
  TMatrixD nuEllipseTwoOne_, nuEllipseTwoTwo_;

  ClosestApproach *closestApproachOne_, *closestApproachTwo_;
  double d1_, theta1_, d2_, theta2_, d_, theta_;

  lightJetChiSquareMinimumSolver *lightJetChiSquare_;
  double lightJetChi2_, bJetChi2_, chi2_;

  void setBJets();
  void setLightJets();
  void setLightJetCollections();
  void setLeptons();

  void reset();

  void setupNeutrinoSolutions();
  void setupNeutrinoEllipses();
  void setupClosestApproaches();

  void calcInitialLightJetDeltas(vector<double>& , vector<double>& );

 public:

  topSystemChiSquare(vector<math::XYZTLorentzVector> , vector<math::XYZTLorentzVector> ,
		    int , int, math::XYZTLorentzVector , 
		    int , int, math::XYZTLorentzVector ,
		    math::XYZTLorentzVector , double , double );

  topSystemChiSquare(const topSystemChiSquare& other);

  ~topSystemChiSquare();

  double operator()(const double* );

  void setBJetDeltas(double , double , double , double );
  void setupBJetChiSquare(double , double , double , double );

  double calcBJetChiSquare();
  
  void minimizeLightJetChiSquare();
  void minimizeBJetChiSquare();
  double getChiSquare();

  double calcSolution();
};

#endif
