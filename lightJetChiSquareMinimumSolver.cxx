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
  //calcVector();
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(int nObjects,
							       double& dx, double& dy) :
  jetPxWidths2_ (vector<double>(nObjects,0.)),
  jetPyWidths2_ (vector<double>(nObjects,0.)),
  jetPxPyWidths_(vector<double>(nObjects,0.)),
  dx_           (dx),
  dy_           (dy),
  dxCheck_      (0.),
  dyCheck_      (0.),
  A_            (2*(nObjects-1),2*(nObjects-1)),
  solver_       (new TDecompLU(2*(nObjects-1))),
  b_            (2*(nObjects-1)),
  minDeltasX_   (nObjects,0.),
  minDeltasY_   (nObjects,0.),
  chi2_         (0.),
  nJets_        (nObjects)
{
  //checkSize(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  //setCartesianWidths(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  //calcLinearCoefficients();
  //calcVector();
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
  calcLinearCoefficients();
}

lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
  //delete all the pointers
  delete solver_;
}

void lightJetChiSquareMinimumSolver::setupEquations(vector<double>& jetPts , vector<double>& jetPtWidths , 
						    vector<double>& jetPhis, vector<double>& jetPhiWidths)
{
  checkSize(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  if(jetPts.size() != jetPxWidths2_.size())
    {
      cout << "Unequal number of cartesian and radial jets!" << endl;
      return;
    }
  setCartesianWidths(jetPts, jetPtWidths, jetPhis, jetPhiWidths);
  calcLinearCoefficients();
}

void lightJetChiSquareMinimumSolver::setCartesianWidths(vector<double>& jetPts , vector<double>& jetPtWidths , 
							vector<double>& jetPhis, vector<double>& jetPhiWidths)
{
  for( unsigned int i = 0 ; i < nJets_;  i++)
    {

      double halfPt2 = 0.5*jetPts.at(i)*jetPts.at(i);
      double sigmaPt2 = log(1+jetPtWidths.at(i));
      sigmaPt2*=sigmaPt2;
      double sigmaPhi2 = jetPhiWidths.at(i)*jetPhiWidths.at(i);

      double expSigmaPt2 = exp(sigmaPt2);
      double expTwoSigmaPt2 = expSigmaPt2*expSigmaPt2;
      double expMinusSigmaPhi2 = exp(-sigmaPhi2);
      double expMinusTwoSigmaPhi2 = expMinusSigmaPhi2*expMinusSigmaPhi2;
      double cosTwoPhi = cos(2.*jetPhis.at(i));
      double sinTwoPhi = sin(2.*jetPhis.at(i));
      jetPxWidths2_.at(i)   =  halfPt2* ( expTwoSigmaPt2 * (1 + cosTwoPhi*expMinusTwoSigmaPhi2)
				      - expSigmaPt2 * (1 + cosTwoPhi) * expMinusSigmaPhi2) ;
      jetPyWidths2_.at(i)   =  halfPt2* ( expTwoSigmaPt2 * (1 - cosTwoPhi*expMinusTwoSigmaPhi2)
				      - expSigmaPt2 * (1 - cosTwoPhi) * expMinusSigmaPhi2) ;
      jetPxPyWidths_.at(i) = halfPt2 * sinTwoPhi *( expTwoSigmaPt2 *expMinusTwoSigmaPhi2
						- expSigmaPt2 * expMinusSigmaPhi2);
      //cout << "calculating widths:\n"
      //     << "pt  is " << jetPts.at(i)  << " with width of " << log(1+jetPtWidths .at(i)) << "\n"
      //     << "phi is " << jetPhis.at(i) << " with width of " << log(1+jetPhiWidths .at(i)) << "\n"
      //     << "px  is " << jetPts.at(i)*cos(jetPhis.at(i)) << " with width of " << sqrt(jetPxWidths2_.at(i)) << "\n"
      //     << "py  is " << jetPts.at(i)*sin(jetPhis.at(i)) << " with width of " << sqrt(jetPyWidths2_.at(i)) << "\n"
      //     << "correlation coefficient is " << jetPxPyWidths_.at(i)/(sqrt(jetPxWidths2_.at(i))*sqrt(jetPyWidths2_.at(i))) << endl;
    }
}

void lightJetChiSquareMinimumSolver::calcLinearCoefficients()
{
  //cout << "calculating linear coefficients" << endl;
  double nthSigmaX2 = jetPxWidths2_.at(nJets_-1);
  double nthSigmaY2 = jetPyWidths2_.at(nJets_-1);
  double nthSigmaXY = jetPxPyWidths_.at(nJets_-1);
  double nthSigmaDenominator = (nthSigmaX2*nthSigmaY2-nthSigmaXY*nthSigmaXY);
  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      double ithSigmaX2 = jetPxWidths2_.at(i);
      double ithSigmaY2 = jetPyWidths2_.at(i);
      double ithSigmaXY = jetPxPyWidths_.at(i);
      double ithSigmaRatio = nthSigmaDenominator/(ithSigmaX2*ithSigmaY2-ithSigmaXY*ithSigmaXY);
      for( unsigned int j = 0 ; j < nJets_-1 ; j++ )
        {
	  //cout << "Index i=" << i << " , j=" << j << endl;
	  A_[i][j] = nthSigmaY2;
          A_[nJets_-1 + i][j] = -1*nthSigmaXY;
          A_[nJets_-1 + i][nJets_-1 + j] = nthSigmaX2;
	  if(i == j)
	    {
	      A_[i][j] += jetPyWidths2_.at(i)*ithSigmaRatio;
	      A_[nJets_-1 + i][j] -= jetPxPyWidths_.at(i)*ithSigmaRatio;
	      A_[nJets_-1 + i][nJets_-1 + j] += jetPxWidths2_.at(i)*ithSigmaRatio;
	    }
          A_[i][nJets_-1 + j] = A_[nJets_-1 + i][j];
        }
    }

  //A_.Print();

  //cout << "setting matrix" << endl;
  dynamic_cast<TDecompLU*>(solver_)->SetMatrix(A_);
  bool checkDecomposition = solver_->Decompose();
  if(!checkDecomposition)
    {
      cout <<"singular matrix -- using SVD decomposition -- experimental!"<<endl;
      delete solver_;
      solver_ = new TDecompSVD(2*(nJets_-1),2*(nJets_-1));
      dynamic_cast<TDecompSVD*>(solver_)->SetMatrix(A_);
      dynamic_cast<TDecompSVD*>(solver_)->Decompose();
      //dynamic_cast<TDecompSVD*>(solver_)->Print();
    }
}

void lightJetChiSquareMinimumSolver::calcVector()
{
  if(dxCheck_ == dx_ && dyCheck_ == dy_) return;
  //cout << "calculating minimum chi^2" << endl;
  dxCheck_ = dx_;
  dyCheck_ = dy_; 
  //chi2_ = dx_*dx_ + dy_*dy_;
  //return;
  double nthSigmaX2 = jetPxWidths2_.at(nJets_-1);
  double nthSigmaY2 = jetPyWidths2_.at(nJets_-1);
  double nthSigmaXY = jetPxPyWidths_.at(nJets_-1);
  //cout << "dx_ = " << dx_ << endl;
  //cout << "dy_ = " << dy_ << endl;
  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      b_[i]            = dx_ * nthSigmaY2 - dy_ * nthSigmaXY;
      b_[nJets_-1 + i] = dy_ * nthSigmaX2 - dx_ * nthSigmaXY;
    }

  //b_.Print();
  
  solver_->Solve(b_);

  minDeltasX_.at(nJets_-1) = dx_;
  minDeltasY_.at(nJets_-1) = dy_;
  chi2_ = 0;

  for( unsigned int i = 0 ; i < nJets_-1 ; i++ )
    {
      minDeltasX_.at(i) = b_[i];
      minDeltasX_.at(nJets_-1) -= b_[i];
      minDeltasY_.at(i) = b_[nJets_-1 + i];
      minDeltasY_.at(nJets_-1) -= b_[nJets_-1 + i];
      chi2_ += ((minDeltasX_.at(i)*minDeltasX_.at(i))*jetPyWidths2_.at(i) +
		(minDeltasY_.at(i)*minDeltasY_.at(i))*jetPxWidths2_.at(i) -
		2*(minDeltasX_.at(i)*minDeltasY_.at(i))*jetPxPyWidths_.at(i))/
	       (jetPxWidths2_.at(i)*jetPyWidths2_.at(i) - jetPxPyWidths_.at(i)*jetPxPyWidths_.at(i));
    }  

  chi2_ += ((minDeltasX_.at(nJets_-1)*minDeltasX_.at(nJets_-1))*jetPyWidths2_.at(nJets_-1) +
	    (minDeltasY_.at(nJets_-1)*minDeltasY_.at(nJets_-1))*jetPxWidths2_.at(nJets_-1) -
	    2*(minDeltasX_.at(nJets_-1)*minDeltasY_.at(nJets_-1))*jetPxPyWidths_.at(nJets_-1))/
           (jetPxWidths2_.at(nJets_-1)*jetPyWidths2_.at(nJets_-1) - jetPxPyWidths_.at(nJets_-1)*jetPxPyWidths_.at(nJets_-1));

}

void lightJetChiSquareMinimumSolver::printResults()
{

  for( unsigned int i = 0 ; i < nJets_ ; i++ )
    {
      cout << "delta px " << i+1 << " = " << minDeltasX_.at(i) << endl;
      cout << "delta py " << i+1 << " = " << minDeltasY_.at(i) << endl;
    }
  
}

double lightJetChiSquareMinimumSolver::getChiSquare()
{
  calcVector();
  //vector<double>::iterator thisDeltaX = minDeltasX_.begin();
  //double deltaXCheck(0.);
  //double deltaYCheck(0.);
  //for(vector<double>::iterator thisDeltaY = minDeltasY_.begin(); thisDeltaY != minDeltasY_.end(); thisDeltaX++, thisDeltaY++)
  //  {
  //    deltaXCheck+=*thisDeltaX;
  //    deltaYCheck+=*thisDeltaY;
  //  }
  //cout << "delta x = " << dx_ << " and delta x check = " << deltaXCheck << endl;
  //cout << "delta y = " << dy_ << " and delta y check = " << deltaYCheck << endl;
  return chi2_;
}

void lightJetChiSquareMinimumSolver::checkSize(vector<double>& jetPts , vector<double>& jetPtWidths , 
					       vector<double>& jetPhis, vector<double>& jetPhiWidths)
{
  if( jetPts.size() != jetPtWidths.size() )
    {
      cout << "Unequal number of jet pTs and widths!" << endl;
      cout << "there are " << jetPts.size() << " jet pts and " << jetPtWidths.size() << " jet pt widths" << endl;
      return;
    }
  if( jetPhis.size() != jetPhiWidths.size() )
    {
      cout << "Unequal number of jet phis and widths!" << endl;
      cout << "there are " << jetPhis.size() << " jet phis and " << jetPhiWidths.size() << " jet phi widths" << endl;
      return;
    }
  if( jetPts.size() != jetPhis.size() )
    {
      cout << "Unequal number of jet pTs and phis!" << endl;
      return;
    }
  return;
}
