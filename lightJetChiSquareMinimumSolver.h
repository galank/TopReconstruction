#ifndef LIGHTJETCHISQUAREMINIMUMSOLVER
#define LIGHTJETCHISQUAREMINIMUMSOLVER

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBase.h"

#include "Math/GenVector/LorentzVector.h"

using namespace std;
using namespace ROOT::Math;
typedef LorentzVector<PxPyPzE4D<double> > XYZTLorentzVector;

class lightJetChiSquareMinimumSolver
{

private:

  bool do3D_;

  vector<double> jetPxWidths2_ , jetPyWidths2_ , jetPzWidths2_ , jetPxPyWidths_ , jetPxPzWidths_ , jetPyPzWidths_;

  const double &dx_, &dy_ , &dz_;
  double dxCheck_, dyCheck_, dzCheck_;

  vector<TMatrixD> jetSigmas2D_;
  vector<TMatrixD> jetSigmas3D_;

  TDecompBase* inverter2D_;
  TDecompBase* inverter3D_;

  TMatrixD inverseSumSigmas2D_;
  TMatrixD inverseSumSigmas3D_;

  vector<double> minDeltasX_;
  vector<double> minDeltasY_;
  vector<double> minDeltasZ_;

  double chi2_;

  unsigned int nJets_;

  void setCartesianWidths(vector<XYZTLorentzVector>& , vector<double>& , vector<double>& , vector<double>& );

  void calcMin();

  void calcSigmas();

  void checkSize(vector<XYZTLorentzVector>& , vector<double>& , vector<double>& , vector<double>& );

 public:

  lightJetChiSquareMinimumSolver(vector<XYZTLorentzVector>& , vector<double>& , vector<double>& , vector<double>& , double& , double& , double& , bool );
  lightJetChiSquareMinimumSolver(int, double& , double& , double& , bool );

  lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver& other);

  void setupEquations(vector<XYZTLorentzVector>& , vector<double>& , vector<double>&, vector<double>&);
  void printResults();

  int nJets(){return nJets_;};

  double jetPxWidth2(int i){return jetPxWidths2_.at(i);};
  double jetPyWidth2(int i){return jetPyWidths2_.at(i);};
  double jetPzWidth2(int i){return jetPzWidths2_.at(i);};
  double jetPxPyWidth(int i){return jetPxPyWidths_.at(i);};
  double jetPxPzWidth(int i){return jetPxPzWidths_.at(i);};
  double jetPyPzWidth(int i){return jetPyPzWidths_.at(i);};

  //void setRecoil(double , double, double );

  vector<double>* getMinDeltasX(){return &minDeltasX_;};
  vector<double>* getMinDeltasY(){return &minDeltasY_;};
  vector<double>* getMinDeltasZ(){return &minDeltasZ_;};

  ~lightJetChiSquareMinimumSolver();

  double getChiSquare();
  
};

#endif
