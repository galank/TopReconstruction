#include "lightJetChiSquareMinimumSolver.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"


lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(vector<double>& jetPts , vector<double>& jetPtWidths , 
							       vector<double>& jetPhis, vector<double>& jetPhiWidths,
							       double& dx, double& dy) :
  jetPxWidths2_ (vector<double>(jetPts.size(),0.)),
  jetPyWidths2_ (vector<double>(jetPts.size(),0.)),
  jetPxPyWidths_(vector<double>(jetPts.size(),0.)),
  dx_           (dx),
  dy_           (dy),
  dxCheck_      (0.),
  dyCheck_      (0.),
  A_            (2*(jetPts.size()-1),2*(jetPts.size()-1)),
  solver_       (new TDecompLU(2*(jetPts.size()-1))),
  b_            (2*(jetPts.size()-1)),
  minDeltasX_   (jetPts.size(),0.),
  minDeltasY_   (jetPts.size(),0.),
  chi2_         (0.),
  nJets_        (jetPts.size())
{
  checkSize(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  setCartesianWidths(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  calcLinearCoefficients();
  calcVector();
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver& other) :
  jetPxWidths2_ (other.jetPxWidths2_),
  jetPyWidths2_ (other.jetPyWidths2_),
  jetPxPyWidths_(other.jetPxPyWidths_),
  dx_           (other.dx_),
  dy_           (other.dy_),
  dxCheck_      (other.dxCheck_),
  dyCheck_      (other.dyCheck_),
  A_            (other.A_),
  solver_       (dynamic_cast<TDecompBase*>(other.solver_->Clone())),
  b_            (other.b_),
  minDeltasX_   (other.minDeltasX_),
  minDeltasY_   (other.minDeltasY_),
  chi2_         (other.chi2_),
  nJets_        (other.nJets_)
{
}

lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
  //delete all the pointers
  delete solver_;
}

void lightJetChiSquareMinimumSolver::setCartesianWidths(vector<double>& jetPts , vector<double>& jetPtWidths , 
							vector<double>& jetPhis, vector<double>& jetPhiWidths)
{
  for( unsigned int i = 0 ; i < nJets_;  i++)
    {
      double expTwoSigmaPtMinusTwoSigmaPhi  = exp(2.*log(1+jetPtWidths .at(i))-2.*jetPhiWidths.at(i));
      double p = 0.5*pow(jetPts.at(i),2)*expTwoSigmaPtMinusTwoSigmaPhi;
      double expTwoSigmaPhi = exp(2.*jetPhiWidths.at(i));
      double cosTwoPhi = cos(2.*jetPhis.at(i));
      double sinTwoPhi = sin(2.*jetPhis.at(i));
      jetPxWidths2_.at(i)   =  p * (expTwoSigmaPhi + cosTwoPhi) ;
      jetPyWidths2_.at(i)   =  p * (expTwoSigmaPhi - cosTwoPhi) ;
      jetPxPyWidths_.at(i) = p * sinTwoPhi;      
    }
}

void lightJetChiSquareMinimumSolver::calcLinearCoefficients()
{
  double nthSigmaX2 = jetPxWidths2_.at(nJets_);
  double nthSigmaY2 = jetPyWidths2_.at(nJets_);
  double nthSigmaXY = jetPxPyWidths_.at(nJets_);
  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      for( unsigned int j = 0 ; j < nJets_-1 ; j++ )
        {
	  A_[i][j] = nthSigmaX2;
          A_[nJets_-1 + i][j] = -1*nthSigmaXY;
          A_[nJets_-1 + i][nJets_-1 + j] = nthSigmaY2;
	  if(i == j)
	    {
	      A_[i][j] += jetPxWidths2_.at(i);
	      A_[nJets_-1 + i][j] -= jetPxPyWidths_.at(i);
	      A_[nJets_-1 + i][nJets_-1 + j] += jetPyWidths2_.at(i);
	    }
          A_[i][nJets_-1 + j] = A_[nJets_-1 + i][j];
        }
    }
  dynamic_cast<TDecompLU*>(solver_)->SetMatrix(A_);
  bool checkDecomposition = solver_->Decompose();
  if(!checkDecomposition)
    {
      cout <<"singular matrix -- using SVD decomposition -- experimental!"<<endl;
      delete solver_;
      solver_ = new TDecompSVD(2*nJets_-1,2*nJets_-1);
      dynamic_cast<TDecompSVD*>(solver_)->SetMatrix(A_);
      solver_->Decompose();
    }
}

void lightJetChiSquareMinimumSolver::calcVector()
{
  if(dxCheck_ == dx_ && dyCheck_ == dy_) return;
  dxCheck_ = dx_;
  dyCheck_ = dy_;
  double nthSigmaX2 = jetPxWidths2_.at(nJets_);
  double nthSigmaY2 = jetPyWidths2_.at(nJets_);
  double nthSigmaXY = jetPxPyWidths_.at(nJets_);
  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      b_[i]            = dx_ * nthSigmaX2 - dy_ * nthSigmaXY;
      b_[nJets_-1 + i] = dy_ * nthSigmaY2 - dx_ * nthSigmaXY;
    }
  
  solver_->Solve(b_);

  minDeltasX_.at(nJets_-1) = dx_;
  minDeltasY_.at(nJets_-1) = dy_;
  chi2_ = 0;
  
  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      minDeltasX_.at(i) = b_[i];
      minDeltasX_.at(nJets_-1) -= b_[i];
      minDeltasY_.at(i) = b_[nJets_-1 + i];
      minDeltasX_.at(nJets_-1) -= b_[nJets_-1 + i];
      chi2_ += ((minDeltasX_.at(i)*minDeltasX_.at(i))*jetPxWidths2_.at(i) +
		(minDeltasY_.at(i)*minDeltasY_.at(i))*jetPyWidths2_.at(i) -
		2*(minDeltasX_.at(i)*minDeltasY_.at(i))*jetPxPyWidths_.at(i))/
	       (jetPxWidths2_.at(i)*jetPyWidths2_.at(i) - jetPxPyWidths_.at(i));
    }

  chi2_ += ((minDeltasX_.at(nJets_-1)*minDeltasX_.at(nJets_-1))*jetPxWidths2_.at(nJets_-1) +
	    (minDeltasY_.at(nJets_-1)*minDeltasY_.at(nJets_-1))*jetPyWidths2_.at(nJets_-1) -
	    2*(minDeltasX_.at(nJets_-1)*minDeltasY_.at(nJets_-1))*jetPxPyWidths_.at(nJets_-1))/
           (jetPxWidths2_.at(nJets_-1)*jetPyWidths2_.at(nJets_-1) - jetPxPyWidths_.at(nJets_-1));

}

double lightJetChiSquareMinimumSolver::getChiSquare()
{
  calcVector();
  return chi2_;
}

void lightJetChiSquareMinimumSolver::checkSize(vector<double>& jetPts , vector<double>& jetPtWidths , 
					       vector<double>& jetPhis, vector<double>& jetPhiWidths)
{
  if( jetPts.size() != jetPtWidths.size() )
    {
      cout << "Unequal number of jet pTs and widths!" << endl;
      return;
    }
  if( jetPhis.size() != jetPhiWidths.size() )
    {
      cout << "Unequal number of jet phis and widths!" << endl;
      return;
    }
  if( jetPts.size() != jetPhis.size() )
    {
      cout << "Unequal number of jet pTs and phis!" << endl;
      return;
    }
  return;
}
