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

using namespace std;

class lightJetChiSquareMinimumSolver
{

private:

  vector<double> jetPxWidths2_ , jetPyWidths2_ , jetPxPyWidths_;

  const double &dx_, &dy_;
  double dxCheck_, dyCheck_;

  TMatrixD A_;
  TDecompBase* solver_;
  TVectorD b_;

  vector<double> minDeltasX_;
  vector<double> minDeltasY_;

  double chi2_;

  unsigned int nJets_;

  void setCartesianWidths(vector<double>& jetPts , vector<double>& jetPtWidths , 
			  vector<double>& jetPhis, vector<double>& jetPhiWidths);

  void calcLinearCoefficients();

  void calcVector();

  void checkSize(vector<double>& jetPts , vector<double>& jetPtWidths , 
		 vector<double>& jetPhis, vector<double>& jetPhiWidths);

 public:

  lightJetChiSquareMinimumSolver(vector<double>& , vector<double>& , vector<double>& , vector<double>& , double& , double& );
  lightJetChiSquareMinimumSolver(int, double& , double& );

  lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver& other);

  void setupEquations(vector<double>& , vector<double>& , vector<double>&, vector<double>&);

  int nJets(){return nJets_;};
  double jetPxWidth2(int i){return jetPxWidths2_.at(i);};
  double jetPyWidth2(int i){return jetPyWidths2_.at(i);};
  double jetPxPyWidth(int i){return jetPxPyWidths_.at(i);};
  vector<double>* getMinDeltasX(){return &minDeltasX_;};
  vector<double>* getMinDeltasY(){return &minDeltasY_;};

  ~lightJetChiSquareMinimumSolver();

  double getChiSquare();
  
};

#endif
