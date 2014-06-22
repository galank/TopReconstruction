// C++ implementation of arXiv:1305.1872
// Analytic solutions for neutrino momenta in decay of top quarks
// Gala Nicolas Kaufman (gnn4@cornell.edu)

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <list>
#include <utility>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TMath.h>
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

#include "Math/RootFinderAlgorithms.h"
#include "Math/Polynomial.h"

using namespace std;
using namespace ROOT::Math;

class NeutrinoEllipseCalculator{

 private:
  //inputs: lepton and b-jet momenta
  //TLorentzVector bJet_;
  //TLorentzVector lepton_;

  const double& bJetPx_;
  const double& bJetPy_;
  const double& bJetPz_;
  const double& bJetE_ ;

  const double& bJetPtWidth_;
  const double& bJetPhiWidth_;
  const double& bJetEtaWidth_;

  const double& leptonPx_;
  const double& leptonPy_;
  const double& leptonPz_;
  const double& leptonE_ ;

  const double& unscaled_bJetPx_;
  const double& unscaled_bJetPy_;
  const double& unscaled_bJetPz_;
  const double& unscaled_bJetE_ ;
  const double& unscaled_bJetP_ ;

  double bJetP2_, bJetP_, bJetMass2_, bJetE_nonConst_;
  double bJetBeta_,   bJetBeta2_,   bJetGamma_,   bJetGamma2_;
  double leptonBeta_, leptonBeta2_, leptonGamma_, leptonGamma2_;
  double leptonP2_, leptonP_, leptonMass2_;
  double leptonPhi_, leptonTheta_;

  //particle masses
  double mW_, mt_, mnu_;
  double mW2_, mt2_, mnu2_;

  //parameters
  double x0_, x0p_;
  double Sx_, Sy_;
  double epsilon2_;

  double c_, s_,c2_,s2_; //cosine and sine of theta_{b,l}

  double omega_;
  double Omega_;
  double x1_, y1_;
  double Z2_;

  //matrices
  TMatrixD Ab_;
  TMatrixD Al_;

  TMatrixD Htilde_;
  TMatrixD H_;
  TMatrixD Hperp_;
  TMatrixD HperpInv_;

  TMatrixD Nperp_;

  TVectorD nuPerp_;
  TVectorD pNu_;

  int nRanges_;
  pair<pair<double,bool>,pair<double,bool> > bJetLogSFRange_;

 public:

  void setBJetFactors();
  void setTempBJetFactors(double , double , double , double);
  void setLeptonFactors();

  double getZ2(double, double, double);

  void setAngles();
  void setTempAngles(double,double,double);

  void initializeMatrices();

  TMatrixD rotationMatrix(int, double);

  void Wsurface();
  void bJetEllipsoid();
  void leptonEllipsoid();
  void calcZ2();
  void neutrinoSolution();
  void labSystemTransform();

  bool errorFlag_;

  NeutrinoEllipseCalculator(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
  NeutrinoEllipseCalculator(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, double , double , double );
  ~NeutrinoEllipseCalculator();

  bool badPoint(){return errorFlag_;};

  void setupEllipse(double , double , double );

  //void setBJet(const double , const double, const double , const double );
  //void setLepton(const double , const double, const double , const double );

  void setTopMass(double& mTop){ mt_=mTop; };
  void setWBosonMass(double& mW){ mW_=mW; };
  void setNeutrinoMass(double& mNu){ mnu_=mNu; };
  void setMasses(double& , double& , double& );
  
  TMatrixD* getHomogeneousNeutrinoEllipse();
  TMatrixD* getExtendedNeutrinoEllipse();

  TVectorD* getNeutrinoMomentum(double theta);

  void calcNeutrinoEllipse();
  void calcExtendedNeutrinoEllipse();
  
  void calcBJetCorrection();
  pair<pair<double,bool>,pair<double,bool> > getBJetLogSFRange(int& nRanges);

  void print3By3Matrix(const TMatrixD& m);
  void print4By4Matrix(const TMatrixD& m);

};
