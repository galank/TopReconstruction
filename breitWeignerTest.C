#include "Math/QuantFuncMathCore.h"
#include <iostream>

using namespace std;
using namespace ROOT::Math;

double breitWignerError(double& mass ,double& width, const double& deltaMass)
{
  double scaledDeltaMass = deltaMass*width;
  double scaledDeltaMass2 = scaledDeltaMass*scaledDeltaMass;
  double Zscore = 1.4142135623730950488*normal_quantile(0.63661977236758134308*atan2(scaledDeltaMass2 + 2.*scaledDeltaMass,mass*width),1.0);
  return Zscore*Zscore;
}

void breitWignerTest()
{
  double mtop = 173.;
  double wtop = 2.0;
  double mW = 80.4;
  double wW = 2.09;

  for(int i = 0; i<20; i++)
    {
      double sigma = (i - 10.)/5.;
      cout << "for input value of " << sigma << " which has nominal chi-square of " << sigma*sigma << "\n"
	   << "we find that for a top of mass 173 GeV and width 2.00 GeV the adjusted chi-square is " << breitWignerError(mtop,wtop,sigma) << "\n"
	   << "we find that for a W   of mass 80.4GeV and width 2.09 GeV the adjusted chi-square is " << breitWignerError(mW,wW,sigma) << endl;
    }
}

