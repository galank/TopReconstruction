#ifndef HADRONICTOPSYSTEMCHISQUARE
#define HADRONICTOPSYSTEMCHISQUARE 

#include "topSystemChiSquare.h"

class hadronicTopSystemChiSquare : public topSystemChiSquare
{
 private:
//  double WDaughter2Px_;
//  double WDaughter2Py_;
//  double WDaughter2Pz_;
//  double WDaughter2E_ ;
//
//  double WDaughter2Pt_ ;
//  double WDaughter2Phi_;
//  double WDaughter2Eta_;

  double WDaughter2PtWidth_ ;
  double WDaughter2PhiWidth_;
  double WDaughter2EtaWidth_;

  double WDaughter2PtDelta_ ;
  double WDaughter2PhiDelta_;
  double WDaughter2EtaDelta_;

  double reconstructed_WDaughter2Px_;
  double reconstructed_WDaughter2Py_;
  double reconstructed_WDaughter2Pz_;
  double reconstructed_WDaughter2E_ ;
	 
  double reconstructed_WDaughter2Pt_ ;
  double reconstructed_WDaughter2Phi_;
  double reconstructed_WDaughter2Eta_;
  double reconstructed_WDaughter2Mass2_;

  double chi2_;
//
//  double topPx_;
//  double topPy_;
//  double topPz_;
//  double topE_ ;


 public:
  
  hadronicTopSystemChiSquare(const double&, const double&, const double&, const double&,
			     const double&, const double&, const double&,
			     const double&, const double&, const double&, const double&,
			     const double&, const double&, const double&,
			     const double&, const double&, const double&, const double&,
			     const double&, const double&, const double&,
			     const double&, const double&,
			     const double&, const double& ) ;

  ~hadronicTopSystemChiSquare();

  virtual void printTopConstituents();

  virtual void setEllipseAngle(double);

  void setWDaughter2(double, double, double);
  //void resetWDaughter2();
  void calcWDaughter2Deltas();
  void setWDaughter2Deltas(double , double , double );
  virtual void getWDaughter2Deltas(double& , double& , double& );
  virtual void printWDaughter2();

  void calcChiSquare();
  double getChiSquare();
  virtual double getHadronicChiSquare();

};

#endif
