#define topReconstructionFromLHE_cxx
#include "topReconstructionFromLHE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include "topEventMinimizer.h"
#include "neutrinoSolutions.h"

using namespace std;
using namespace ROOT::Math;


void topReconstructionFromLHE::printVector(XYZTLorentzVector& v)
{
  cout << "px = " << v.Px() << endl;
  cout << "py = " << v.Py() << endl;
  cout << "pz = " << v.Pz() << endl;
  cout << "E  = " << v.E () << endl;
}

void topReconstructionFromLHE::Loop(TString dir, int whichLoop, int maxLoops)
{
//   In a ROOT session, you can do:
//      Root > .L topReconstructionFromLHE.C
//      Root > topReconstructionFromLHE t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TRandom3 rand;
   rand.SetSeed(10);

   initOutput(dir,whichLoop);

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;

   vector<XYZTLorentzVector> nonTopObjects;
   vector<double> nonTopObjectPtWidths;
   vector<double> nonTopObjectEtaWidths;
   vector<double> nonTopObjectPhiWidths;


   int jStart = whichLoop*(nentries/maxLoops) + ((whichLoop>(maxLoops-nentries%maxLoops))?(whichLoop+nentries%maxLoops-maxLoops):0);

   int jFinish = jStart + (nentries+whichLoop)/maxLoops;

   for (Long64_t jentry=jStart; jentry<jFinish; jentry++) {
     cout << "BEGINNING BRANCH NUMBER " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if (jentry > 5) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      XYZTLorentzVector lepton, antiLepton;
      XYZTLorentzVector neutrino, antiNeutrino;
      XYZTLorentzVector bottomQuark, antiBottomQuark;
      XYZTLorentzVector lightParton1, lightParton2;
      XYZTLorentzVector topQuark, antiTopQuark;
      XYZTLorentzVector recoTopQuark(0,0,0,0), recoAntiTopQuark(0,0,0,0);
      XYZTLorentzVector WPlus, WMinus;
      XYZTLorentzVector recoWPlus(0,0,0,0), recoWMinus(0,0,0,0);
      double METx(0.), METy(0.);
      int iJet = 0;
      for(int iParticle = 0 ; iParticle < n_particles; iParticle++)
	{
	  if( PID->at(iParticle) == 6 ) topQuark.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	  if( PID->at(iParticle) == -6 ) antiTopQuark.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	  if( PID->at(iParticle) == 24 ) WPlus.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	  if( PID->at(iParticle) == -24 ) WMinus.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	  if( PID->at(iParticle) == 11 || PID->at(iParticle) == 13 || PID->at(iParticle) == 15 ) 
	    {
	      lepton.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoWMinus += lepton;
	      //recoAntiTopQuark += lepton;
	    }
	  if( PID->at(iParticle) == -11 || PID->at(iParticle) == -13 || PID->at(iParticle) == -15 )
	    {
	      antiLepton.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoWPlus += antiLepton;
	      //recoTopQuark += antiLepton;
	    }
	  if( PID->at(iParticle) == 12 || PID->at(iParticle) == 14 || PID->at(iParticle) == 16 )
	    {
	      neutrino.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoWPlus += neutrino;
	      //recoTopQuark += neutrino;
	    }
	  if( PID->at(iParticle) == -12 || PID->at(iParticle) == -14 || PID->at(iParticle) == -16 ) 
	    {
	      antiNeutrino.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoWMinus += antiNeutrino;
	      //recoAntiTopQuark += antiNeutrino;
	    }
	  if( PID->at(iParticle) == 5 ) 
	    {
	      bottomQuark.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoTopQuark += bottomQuark;
	    }
	  if( PID->at(iParticle) == -5 )
	    {
	      antiBottomQuark.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
	      //recoAntiTopQuark += antiBottomQuark;
	    }
	  if( status->at(iParticle) == 1 )
	    {
	      if(abs(PID->at(iParticle)) == 1 || 
		 abs(PID->at(iParticle)) == 2 ||  
		 abs(PID->at(iParticle)) == 3 ||  
		 abs(PID->at(iParticle)) == 4 ||  
		 abs(PID->at(iParticle)) == 21 )
		{
		  if (iJet == 0)
		    {
		      lightParton1.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
		      iJet++;
		    }
		  else if (iJet == 1)
		    {
		      lightParton2.SetPxPyPzE(P_X->at(iParticle),P_Y->at(iParticle),P_Z->at(iParticle),E->at(iParticle));
		      iJet++;
		    }
		  else cout <<"too many light partons"<<endl;
		}
	    }
	} //end loop over particles

      //cout << "true W+             mass is " << WPlus.M()        << endl;
      //cout << "true top quark      mass is " << topQuark.M()     << endl;
      //cout << "true W-             mass is " << WMinus.M()       << endl;
      //cout << "true anti top quark mass is " << antiTopQuark.M() << endl;


      //cout << "true neutrino momentum:" << endl; 
      //printVector(neutrino);
      //cout << "true anti-neutrino momentum:" << endl;
      //printVector(antiNeutrino);


      //RESET EVERYTHING
      recoWMinus.SetPxPyPzE(0,0,0,0);
      recoAntiTopQuark.SetPxPyPzE(0,0,0,0);
      recoWPlus.SetPxPyPzE(0,0,0,0);
      recoTopQuark.SetPxPyPzE(0,0,0,0);
      METx = 0.;
      METy = 0.;
      nonTopObjects.clear();
      nonTopObjectPtWidths.clear();
      nonTopObjectEtaWidths.clear();
      nonTopObjectPhiWidths.clear();
      
      
      //for fully leptonic decays
      //METx -= lepton.Px() + antiLepton.Px();
      //METy -= lepton.Py() + antiLepton.Py();
      //simulate hadronically decaying \bar{t} quark, so smear the lepton and subtract it from the MET separately
      METx -= antiLepton.Px();
      METy -= antiLepton.Py();


      XYZTLorentzVector smearedBottomQuark;
      XYZTLorentzVector smearedAntiBottomQuark;
      XYZTLorentzVector smearedLightQuark1;
      XYZTLorentzVector smearedLightQuark2;
      XYZTLorentzVector smearedLightParton1;
      XYZTLorentzVector smearedLightParton2;
      XYZTLorentzVector tempP4;

      double bJet1PtSmear = log(1.+0.1)*rand.Gaus();
      double bJet2PtSmear = log(1.+0.1)*rand.Gaus();
      double lightJet1PtSmear = log(1.+0.1)*rand.Gaus();
      double lightJet2PtSmear = log(1.+0.1)*rand.Gaus();
      double bJet1PhiSmear = 0.01*rand.Gaus();
      double bJet2PhiSmear = 0.01*rand.Gaus();
      double lightJet1PhiSmear = 0.01*rand.Gaus();
      double lightJet2PhiSmear = 0.01*rand.Gaus();
      double bJet1EtaSmear = 0.01*rand.Gaus();
      double bJet2EtaSmear = 0.01*rand.Gaus();
      //smear the lepton and anti-neutrino as if they were jets, to simulate the \bar{t} quark decaying hadronically
      double leptonPtSmear = log(1.+0.1)*rand.Gaus();
      double leptonPhiSmear = 0.01*rand.Gaus();
      double leptonEtaSmear = 0.01*rand.Gaus();
      double antiNeutrinoPtSmear = log(1.+0.1)*rand.Gaus();
      double antiNeutrinoPhiSmear = 0.01*rand.Gaus();
      double antiNeutrinoEtaSmear = 0.01*rand.Gaus();


      
      double pt(bottomQuark.Pt()*exp(bJet1PtSmear));
      double phi(bottomQuark.Phi()+bJet1PhiSmear);
      double eta(bottomQuark.Eta()+bJet1EtaSmear);
      double energy(sqrt((bottomQuark.M2()+bottomQuark.Pt()*bottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet1PtSmear));
      
      smearedBottomQuark.SetPx(pt*cos(phi));
      smearedBottomQuark.SetPy(pt*sin(phi));
      smearedBottomQuark.SetPz(pt*sinh(eta));
      smearedBottomQuark.SetE(energy);
      
      recoTopQuark+=smearedBottomQuark;
	  
      METx -= smearedBottomQuark.Px();
      METy -= smearedBottomQuark.Py();
	  
      pt = antiBottomQuark.Pt()*exp(bJet2PtSmear);
      phi = antiBottomQuark.Phi()+bJet2PhiSmear;
      eta = antiBottomQuark.Eta()+bJet2EtaSmear;
      energy = sqrt((antiBottomQuark.M2()+antiBottomQuark.Pt()*antiBottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet2PtSmear);
      
      smearedAntiBottomQuark.SetPx(pt*cos(phi));
      smearedAntiBottomQuark.SetPy(pt*sin(phi));
      smearedAntiBottomQuark.SetPz(pt*sinh(eta));
      smearedAntiBottomQuark.SetE(energy);
      
      recoAntiTopQuark+=smearedAntiBottomQuark;
	  
      METx -= smearedAntiBottomQuark.Px();
      METy -= smearedAntiBottomQuark.Py();

      //smear the lepton to simulate a hadronic \bar{t} decay
      pt = lepton.Pt()*exp(leptonPtSmear);
      phi = lepton.Phi()+leptonPhiSmear;
      eta = lepton.Eta()+leptonEtaSmear;
      energy = sqrt((lepton.M2()+lepton.Pt()*lepton.Pt()*cosh(eta)*cosh(eta)))*exp(leptonPtSmear);
      
      smearedLightQuark1.SetPx(pt*cos(phi));
      smearedLightQuark1.SetPy(pt*sin(phi));
      smearedLightQuark1.SetPz(pt*sinh(eta));
      smearedLightQuark1.SetE(energy);
      
      METx -= smearedLightQuark1.Px();
      METy -= smearedLightQuark1.Py();

      //also smear the anti-neutrino
      pt = antiNeutrino.Pt()*exp(antiNeutrinoPtSmear);
      phi = antiNeutrino.Phi()+antiNeutrinoPhiSmear;
      eta = antiNeutrino.Eta()+antiNeutrinoEtaSmear;
      energy = sqrt((antiNeutrino.M2()+antiNeutrino.Pt()*antiNeutrino.Pt()*cosh(eta)*cosh(eta)))*exp(antiNeutrinoPtSmear);

      smearedLightQuark2.SetPx(pt*cos(phi));
      smearedLightQuark2.SetPy(pt*sin(phi));
      smearedLightQuark2.SetPz(pt*sinh(eta));
      smearedLightQuark2.SetE(energy);

      METx -= smearedLightQuark2.Px();
      METy -= smearedLightQuark2.Py();

	  
      pt = lightParton1.Pt()*exp(lightJet1PtSmear);
      phi = lightParton1.Phi()+lightJet1PhiSmear;
      eta = lightParton1.Eta();
      energy = lightParton1.E()*exp(lightJet1PtSmear);
      
      smearedLightParton1.SetPx(pt*cos(phi));
      smearedLightParton1.SetPy(pt*sin(phi));
      smearedLightParton1.SetPz(pt*sinh(eta));
      smearedLightParton1.SetE(energy);
      METx -= smearedLightParton1.Px();
      METy -= smearedLightParton1.Py();

	  
      pt = lightParton2.Pt()*exp(lightJet2PtSmear);
      phi = lightParton2.Phi()+lightJet2PhiSmear;
      eta = lightParton2.Eta();
      energy = lightParton2.E()*exp(lightJet2PtSmear);
      
      smearedLightParton2.SetPx(pt*cos(phi));
      smearedLightParton2.SetPy(pt*sin(phi));
      smearedLightParton2.SetPz(pt*sinh(eta));
      smearedLightParton2.SetE(energy);
      METx -= smearedLightParton2.Px();
      METy -= smearedLightParton2.Py();





      recoWMinus += smearedLightQuark1;
      recoWMinus += smearedLightQuark2;
      recoAntiTopQuark+=recoWMinus;

      recoWPlus += antiLepton;
      recoWPlus += neutrino;
      recoTopQuark+=recoWPlus;



      int N=2;
      
      nonTopObjectPtWidths =vector<double>(N,0.1);
      nonTopObjectPhiWidths=vector<double>(N,0.01);
      nonTopObjectEtaWidths=vector<double>(N,0.01);

      nonTopObjects.push_back(smearedLightParton1);
      nonTopObjects.push_back(smearedLightParton2);


      //cout << "METx is " << METx << endl;				     
      //cout << "real METx is " << neutrino.Px() + antiNeutrino.Px() << endl;
      //cout << "METy is " << METy << endl;				     
      //cout << "real METy is " << neutrino.Py() + antiNeutrino.Py() << endl;
      //cout << "neutrino       mass is      " << neutrino.M() << endl;
      //printVector(neutrino);
      //cout << "anti-lepton    mass is      " << antiLepton.M() << endl;
      //printVector(antiLepton);
      //cout << "reco W+        mass is      " << recoWPlus.M() << endl;
      //printVector(recoWPlus);
      //cout << "b quark        mass is      " << smearedBottomQuark.M() << endl;
      //printVector(smearedBottomQuark);
      //cout << "reco top quark mass is      " << recoTopQuark.M() << endl;
      //printVector(recoTopQuark);
      //cout << "anti-neutrino  mass is      " << antiNeutrino.M() << endl;
      //printVector(antiNeutrino);
      //cout << "lepton         mass is      " << lepton.M() << endl;
      //printVector(lepton);
      //cout << "smeared   1st light quark mass is " << smearedLightQuark1.M() << endl;
      //printVector(smearedLightQuark1);
      //cout << "smeared   2nd light quark mass is " << smearedLightQuark2.M() << endl;
      //printVector(smearedLightQuark2);
      //cout << "reco W-        mass is      " << recoWMinus.M() << endl;
      //printVector(recoWMinus);
      //cout << "b-bar quark    mass is      " << smearedAntiBottomQuark.M() << endl;
      //printVector(smearedAntiBottomQuark);
      //cout << "reco anti top quark mass is " << recoAntiTopQuark.M() << endl;
      //printVector(recoAntiTopQuark);
      //cout << "first light parton  mass is " << smearedLightParton1.M() << endl;
      //printVector(smearedLightParton1);
      //cout << "second light parton mass is " << smearedLightParton2.M() << endl;
      //printVector(smearedLightParton2);


	  


      double bJet1Px = smearedAntiBottomQuark.Px();
      double bJet1Py = smearedAntiBottomQuark.Py();
      double bJet1Pz = smearedAntiBottomQuark.Pz();
      double bJet1E  = smearedAntiBottomQuark.E() ;

      double bJet2Px = smearedBottomQuark.Px();
      double bJet2Py = smearedBottomQuark.Py();
      double bJet2Pz = smearedBottomQuark.Pz();
      double bJet2E  = smearedBottomQuark.E() ;

      double lightJet1Px = smearedLightQuark1.Px();
      double lightJet1Py = smearedLightQuark1.Py();
      double lightJet1Pz = smearedLightQuark1.Pz();
      double lightJet1E  = smearedLightQuark1.E() ;

      double lightJet2Px = smearedLightQuark2.Px();
      double lightJet2Py = smearedLightQuark2.Py();
      double lightJet2Pz = smearedLightQuark2.Pz();
      double lightJet2E  = smearedLightQuark2.E() ;

      double antiLeptonPx = antiLepton.Px();
      double antiLeptonPy = antiLepton.Py();
      double antiLeptonPz = antiLepton.Pz();
      double antiLeptonE  = antiLepton.E() ;

      //double leptonPx = lepton.Px();
      //double leptonPy = lepton.Py();
      //double leptonPz = lepton.Pz();
      //double leptonE  = lepton.E() ;




      double mTop=173., mW=80.4;
      double sigmaMTop=2.0, sigmaMW=2.085;
      double sigmaPtJet=0.1, sigmaPhiJet=0.01, sigmaEtaJet=0.01;
      double sigmaPtLep=0., sigmaPhiLep=0., sigmaEtaLep=0.;

      topEventMinimizer ev(nonTopObjects,
			   nonTopObjectPtWidths,
			   nonTopObjectEtaWidths,
			   nonTopObjectPhiWidths,
			   mTop,sigmaMTop,
			   mW,sigmaMW
			   );

      ev.addLeptonicTop(bJet2Px,bJet2Py,bJet2Pz,bJet2E,
                        sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
                        antiLeptonPx,antiLeptonPy,antiLeptonPz,antiLeptonE,
                        sigmaPtLep,sigmaEtaLep,sigmaPhiLep,
                        mTop,sigmaMTop,
                        mW,sigmaMW);
      
      ev.addHadronicTop(bJet1Px,bJet1Py,bJet1Pz,bJet1E,
			sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			lightJet1Px,lightJet1Py,lightJet1Pz,lightJet1E,
			sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			lightJet2Px,lightJet2Py,lightJet2Pz,lightJet2E,
			sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			mTop,sigmaMTop,
			mW,sigmaMW);

      
      //ev.printTopConstituents();
      ev.initializeDeltas();
      //ev.calcTopMassRanges();

      //ev.findStartingValues(50);
      //TString thisPlotName="ellipses_";
      //thisPlotName+=jentry;
      //ev.plotEllipses(thisPlotName);

      //ev.minimizeNonTopChiSquare();

      ev.minimizeTotalChiSquare();
      

      //ev.printTopConstituents();






      //Fill branches

      eventNumber=jentry;

      //Minimizer stuff
      nonTopChi2=ev.getNonTopChiSquare();
      hadronicChi2=ev.getHadronicChiSquare();
      topMassChi2=ev.getTopMassChiSquare();
      topSystemChi2=ev.getTopChiSquare();
      totalChi2=ev.getChiSquare();

      //cout << "Total chi^2 is " << totalChi2_ << endl;
      //cout << "Hadronic chi^2 is " << hadronicChi2_ << endl;
      //cout << "Top mass chi^2 is " << topMassChi2_ << endl;
      //cout << "Non-top chi^2 is " << nonTopChi2_ << endl;
      //cout << "Top system chi^2 is " << topSystemChi2_ << endl;

      innerMinStatus=ev.getInnerMinimizerStatus();
      outerMinStatus=ev.getOuterMinimizerStatus();
      outerMinEdm=ev.getOuterMinimizerEdm();

      //True values
      top1TruePt        =topQuark.Pt()    ; 
      top1TrueEta       =topQuark.Eta()   ; 
      top1TruePhi       =topQuark.Phi()   ;
      top1TrueE         =topQuark.E()     ;
      top1TrueMass      =sqrt(topQuark.M2()); 
      W1TruePt          =WPlus.Pt()       ;
      W1TrueEta         =WPlus.Eta()      ;
      W1TruePhi         =WPlus.Phi()      ;
      W1TrueE           =WPlus.E()        ;
      W1TrueMass        =sqrt(WPlus.M2()) ;
      bJet1TruePt       =bottomQuark.Pt() ;
      bJet1TrueEta      =bottomQuark.Eta();
      bJet1TruePhi      =bottomQuark.Phi();
      bJet1TrueE        =bottomQuark.E()  ;
      W1Daughter1TruePt =antiLepton.Pt()  ;
      W1Daughter1TrueEta=antiLepton.Eta() ;
      W1Daughter1TruePhi=antiLepton.Phi() ;
      W1Daughter1TrueE  =antiLepton.E()   ;
      W1Daughter2TruePt =neutrino.Pt()    ;
      W1Daughter2TrueEta=neutrino.Eta()   ;
      W1Daughter2TruePhi=neutrino.Phi()   ;
      W1Daughter2TrueE  =neutrino.E()     ;
      lightJet1TruePx   =lightParton1.Px(); 
      lightJet1TruePy   =lightParton1.Py();
      lightJet1TruePz   =lightParton1.Pz();
      lightJet1TrueE    =lightParton1.E() ;
      lightJet1TruePt   =lightParton1.Pt();


      top2TruePt        =antiTopQuark.Pt()    ;
      top2TrueEta       =antiTopQuark.Eta()   ;
      top2TruePhi       =antiTopQuark.Phi()   ;
      top2TrueE         =antiTopQuark.E()     ;
      top2TrueMass      =sqrt(antiTopQuark.M2());
      W2TruePt          =WMinus.Pt()          ;
      W2TrueEta         =WMinus.Eta()         ;
      W2TruePhi         =WMinus.Phi()         ;
      W2TrueE           =WMinus.E()           ;
      W2TrueMass        =sqrt(WMinus.M2())    ;
      bJet2TruePt       =antiBottomQuark.Pt() ;
      bJet2TrueEta      =antiBottomQuark.Eta();
      bJet2TruePhi      =antiBottomQuark.Phi();
      bJet2TrueE        =antiBottomQuark.E()  ;
      W2Daughter1TruePt =antiNeutrino.Pt()    ;
      W2Daughter1TrueEta=antiNeutrino.Eta()   ;
      W2Daughter1TruePhi=antiNeutrino.Phi()   ;
      W2Daughter1TrueE  =antiNeutrino.E()     ;
      W2Daughter2TruePt =lepton.Pt()          ;
      W2Daughter2TrueEta=lepton.Eta()         ;
      W2Daughter2TruePhi=lepton.Phi()         ;
      W2Daughter2TrueE  =lepton.E()           ;
      lightJet2TruePx   =lightParton2.Px()    ;
      lightJet2TruePy   =lightParton2.Py()    ;
      lightJet2TruePz   =lightParton2.Pz()    ;
      lightJet2TrueE    =lightParton2.E()     ;		
      lightJet2TruePt   =lightParton2.Pt()    ;



      //Smeared values
      top1SmearedPt        =recoTopQuark.Pt()       ;
      top1SmearedEta       =recoTopQuark.Eta()      ;
      top1SmearedPhi       =recoTopQuark.Phi()      ;
      top1SmearedE         =recoTopQuark.E()        ;
      top1SmearedMass      =sqrt(recoTopQuark.M2()) ;
      W1SmearedPt          =recoWPlus.Pt()          ;
      W1SmearedEta         =recoWPlus.Eta()         ;
      W1SmearedPhi         =recoWPlus.Phi()         ;
      W1SmearedE           =recoWPlus.E()           ;
      W1SmearedMass        =sqrt(recoWPlus.M2())    ;
      bJet1SmearedPt       =smearedBottomQuark.Pt() ;
      bJet1SmearedEta      =smearedBottomQuark.Eta();
      bJet1SmearedPhi      =smearedBottomQuark.Phi();
      bJet1SmearedE        =smearedBottomQuark.E()  ;
      W1Daughter1SmearedPt =antiLepton.Pt()         ;
      W1Daughter1SmearedEta=antiLepton.Eta()        ;
      W1Daughter1SmearedPhi=antiLepton.Phi()        ;
      W1Daughter1SmearedE  =antiLepton.E()          ;
      W1Daughter2SmearedPt =-1.                     ;
      W1Daughter2SmearedEta=-1.                     ;
      W1Daughter2SmearedPhi=-1.                     ;
      W1Daughter2SmearedE  =-1.                     ;
      lightJet1SmearedPx   =smearedLightParton1.Px();
      lightJet1SmearedPy   =smearedLightParton1.Py();
      lightJet1SmearedPz   =smearedLightParton1.Pz();
      lightJet1SmearedE    =smearedLightParton1.E() ;
      lightJet1SmearedPt   =smearedLightParton1.Pt();

      top2SmearedPt        =recoAntiTopQuark.Pt()       ;
      top2SmearedEta       =recoAntiTopQuark.Eta()      ;
      top2SmearedPhi       =recoAntiTopQuark.Phi()      ;
      top2SmearedE         =recoAntiTopQuark.E()        ;
      top2SmearedMass      =sqrt(recoAntiTopQuark.M2()) ;
      W2SmearedPt          =recoWMinus.Pt()             ;
      W2SmearedEta         =recoWMinus.Eta()            ;
      W2SmearedPhi         =recoWMinus.Phi()            ;
      W2SmearedE           =recoWMinus.E()              ;
      W2SmearedMass        =sqrt(recoWMinus.M2())       ;
      bJet2SmearedPt       =smearedAntiBottomQuark.Pt() ;
      bJet2SmearedEta      =smearedAntiBottomQuark.Eta();
      bJet2SmearedPhi      =smearedAntiBottomQuark.Phi();
      bJet2SmearedE        =smearedAntiBottomQuark.E()  ;
      W2Daughter1SmearedPt =smearedLightQuark2.Pt()     ;
      W2Daughter1SmearedEta=smearedLightQuark2.Eta()    ;
      W2Daughter1SmearedPhi=smearedLightQuark2.Phi()    ;
      W2Daughter1SmearedE  =smearedLightQuark2.E()      ;
      W2Daughter2SmearedPt =smearedLightQuark1.Pt()     ;
      W2Daughter2SmearedEta=smearedLightQuark1.Eta()    ;
      W2Daughter2SmearedPhi=smearedLightQuark1.Phi()    ;
      W2Daughter2SmearedE  =smearedLightQuark1.E()      ;
      lightJet2SmearedPx   =smearedLightParton2.Px()    ;
      lightJet2SmearedPy   =smearedLightParton2.Py()    ;
      lightJet2SmearedPz   =smearedLightParton2.Pz()    ;
      lightJet2SmearedE    =smearedLightParton2.E()     ;
      lightJet2SmearedPt   =smearedLightParton2.Pt()    ;



      //Best values


      //Minimizer deltas

      vector<double> bJetPtDeltas; 
      vector<double> bJetPhiDeltas;
      vector<double> bJetEtaDeltas;
      vector<double> firstWDaughterPtDeltas;
      vector<double> firstWDaughterPhiDeltas;
      vector<double> firstWDaughterEtaDeltas;
      vector<double> secondWDaughterPtDeltas;
      vector<double> secondWDaughterPhiDeltas;
      vector<double> secondWDaughterEtaDeltas;
      vector<double> WMassDeltas;
      vector<double> topMassDeltas;
      vector<double> lightJetPxDeltas;
      vector<double> lightJetPyDeltas;

      ev.getBestDeltas(bJetPtDeltas, bJetPhiDeltas, bJetEtaDeltas,
		       firstWDaughterPtDeltas, firstWDaughterPhiDeltas, firstWDaughterEtaDeltas,
		       secondWDaughterPtDeltas, secondWDaughterPhiDeltas, secondWDaughterEtaDeltas,
		       WMassDeltas, topMassDeltas,
		       lightJetPxDeltas, lightJetPyDeltas);


      bJet1PtDelta =bJetPtDeltas .at(0); 
      bJet1PhiDelta=bJetPhiDeltas.at(0); 
      bJet1EtaDelta=bJetEtaDeltas.at(0);

      bJet2PtDelta =bJetPtDeltas .at(1);
      bJet2PhiDelta=bJetPhiDeltas.at(1);
      bJet2EtaDelta=bJetEtaDeltas.at(1);

      W1Daughter1PtDelta =firstWDaughterPtDeltas .at(0); 
      W1Daughter1PhiDelta=firstWDaughterPhiDeltas.at(0);
      W1Daughter1EtaDelta=firstWDaughterEtaDeltas.at(0);

      W2Daughter1PtDelta =firstWDaughterPtDeltas .at(1);
      W2Daughter1PhiDelta=firstWDaughterPhiDeltas.at(1);
      W2Daughter1EtaDelta=firstWDaughterEtaDeltas.at(1);

      W1Daughter2PtDelta =secondWDaughterPtDeltas .at(0);
      W1Daughter2PhiDelta=secondWDaughterPhiDeltas.at(0);
      W1Daughter2EtaDelta=secondWDaughterEtaDeltas.at(0);

      W2Daughter2PtDelta =secondWDaughterPtDeltas .at(1);
      W2Daughter2PhiDelta=secondWDaughterPhiDeltas.at(1);
      W2Daughter2EtaDelta=secondWDaughterEtaDeltas.at(1);

      W1MassDelta=WMassDeltas.at(0);
      W2MassDelta=WMassDeltas.at(1);

      top1MassDelta=topMassDeltas.at(0);
      top2MassDelta=topMassDeltas.at(1);

      lightJet1PxDelta=lightJetPxDeltas.at(0);
      lightJet2PxDelta=lightJetPxDeltas.at(1);

      lightJet1PyDelta=lightJetPyDeltas.at(0);
      lightJet2PyDelta=lightJetPyDeltas.at(1);


      //Object momenta

      ev.getTop         (0,top1BestPx       ,top1BestPy       ,top1BestPz       ,top1BestE       );
      ev.getBJet        (0,bJet1BestPx      ,bJet1BestPy      ,bJet1BestPz      ,bJet1BestE      );
      ev.getW           (0,W1BestPx         ,W1BestPy         ,W1BestPz         ,W1BestE         );
      ev.getWDaughter1  (0,W1Daughter1BestPx,W1Daughter1BestPy,W1Daughter1BestPz,W1Daughter1BestE);
      ev.getWDaughter2  (0,W1Daughter2BestPx,W1Daughter2BestPy,W1Daughter2BestPz,W1Daughter2BestE);
      ev.getNonTopObject(0,lightJet1BestPx,lightJet1BestPy);

      XYZTLorentzVector top1Best;
      top1Best.SetPxPyPzE(top1BestPx,top1BestPy,top1BestPz,top1BestE);
      XYZTLorentzVector bJet1Best;
      bJet1Best.SetPxPyPzE(bJet1BestPx,bJet1BestPy,bJet1BestPz,bJet1BestE);
      XYZTLorentzVector W1Best;
      W1Best.SetPxPyPzE(W1BestPx,W1BestPy,W1BestPz,W1BestE);
      XYZTLorentzVector W1Daughter1Best;
      W1Daughter1Best.SetPxPyPzE(W1Daughter1BestPx,W1Daughter1BestPy,W1Daughter1BestPz,W1Daughter1BestE);
      XYZTLorentzVector W1Daughter2Best;
      W1Daughter2Best.SetPxPyPzE(W1Daughter2BestPx,W1Daughter2BestPy,W1Daughter2BestPz,W1Daughter2BestE);
      XYZTLorentzVector lightJet1Best;
      lightJet1BestPz=smearedLightParton1.Pz();
      lightJet1BestE =sqrt(smearedLightParton1.M2()+pow(lightJet1BestPx,2)+pow(lightJet1BestPy,2)+pow(smearedLightParton1.Pz(),2));
      lightJet1Best.SetPxPyPzE(lightJet1BestPx,lightJet1BestPy,lightJet1BestPz,lightJet1BestE);

      top1BestPt        =top1Best.Pt() ;
      top1BestEta       =top1Best.Eta();
      top1BestPhi       =top1Best.Phi();
      top1BestMass      =max(0.,sqrt(top1Best.M2()));
      bJet1BestPt       =bJet1Best.Pt() ;
      bJet1BestEta      =bJet1Best.Eta();
      bJet1BestPhi      =bJet1Best.Phi();
      W1Daughter1BestPt =W1Daughter1Best.Pt() ;
      W1Daughter1BestEta=W1Daughter1Best.Eta();
      W1Daughter1BestPhi=W1Daughter1Best.Phi();
      W1Daughter2BestPt =W1Daughter2Best.Pt() ;
      W1Daughter2BestEta=W1Daughter2Best.Eta();
      W1Daughter2BestPhi=W1Daughter2Best.Phi();
      W1BestPt          =W1Best.Pt() ;
      W1BestEta         =W1Best.Eta();
      W1BestPhi         =W1Best.Phi();
      W1BestMass        =max(0.,sqrt(W1Best.M2()));
      lightJet1BestPt   =lightJet1Best.Pt();


      ev.getTop         (1,top2BestPx       ,top2BestPy       ,top2BestPz       ,top2BestE       );
      ev.getBJet        (1,bJet2BestPx      ,bJet2BestPy      ,bJet2BestPz      ,bJet2BestE      );
      ev.getWDaughter1  (1,W2Daughter1BestPx,W2Daughter1BestPy,W2Daughter1BestPz,W2Daughter1BestE);
      ev.getWDaughter2  (1,W2Daughter2BestPx,W2Daughter2BestPy,W2Daughter2BestPz,W2Daughter2BestE);
      ev.getW           (1,W2BestPx         ,W2BestPy         ,W2BestPz         ,W2BestE         );
      ev.getNonTopObject(1,lightJet2BestPx,lightJet2BestPy);


      XYZTLorentzVector top2Best;
      top2Best.SetPxPyPzE(top2BestPx,top2BestPy,top2BestPz,top2BestE);
      XYZTLorentzVector bJet2Best;
      bJet2Best.SetPxPyPzE(bJet2BestPx,bJet2BestPy,bJet2BestPz,bJet2BestE);
      XYZTLorentzVector W2Best;
      W2Best.SetPxPyPzE(W2BestPx,W2BestPy,W2BestPz,W2BestE);
      XYZTLorentzVector W2Daughter1Best;
      W2Daughter1Best.SetPxPyPzE(W2Daughter1BestPx,W2Daughter1BestPy,W2Daughter1BestPz,W2Daughter1BestE);
      XYZTLorentzVector W2Daughter2Best;
      W2Daughter2Best.SetPxPyPzE(W2Daughter2BestPx,W2Daughter2BestPy,W2Daughter2BestPz,W2Daughter2BestE);
      XYZTLorentzVector lightJet2Best;
      lightJet2BestPz=smearedLightParton2.Pz();
      lightJet2BestE =sqrt(smearedLightParton2.M2()+pow(lightJet2BestPx,2)+pow(lightJet2BestPy,2)+pow(smearedLightParton2.Pz(),2));
      lightJet2Best.SetPxPyPzE(lightJet2BestPx,lightJet2BestPy,lightJet2BestPz,lightJet2BestE);

      top2BestPt        =top2Best.Pt() ;
      top2BestEta       =top2Best.Eta();
      top2BestPhi       =top2Best.Phi();
      top2BestMass      =max(0.,sqrt(top2Best.M2()));
      bJet2BestPt       =bJet2Best.Pt() ;
      bJet2BestEta      =bJet2Best.Eta();
      bJet2BestPhi      =bJet2Best.Phi();
      W2Daughter1BestPt =W2Daughter1Best.Pt() ;
      W2Daughter1BestEta=W2Daughter1Best.Eta();
      W2Daughter1BestPhi=W2Daughter1Best.Phi();
      W2Daughter2BestPt =W2Daughter2Best.Pt() ;
      W2Daughter2BestEta=W2Daughter2Best.Eta();
      W2Daughter2BestPhi=W2Daughter2Best.Phi();
      W2BestPt          =W2Best.Pt() ;
      W2BestEta         =W2Best.Eta();
      W2BestPhi         =W2Best.Phi();
      W2BestMass        =max(0.,sqrt(W2Best.M2()));
      lightJet2BestPt   =lightJet2Best.Pt();



      //Compare true and smeared
      top1DeltaPtTrueSmeared       =top1TruePt-top1SmearedPt; 
      top1DeltaRTrueSmeared        =deltaR(topQuark,recoTopQuark); 
      top1DeltaMTrueSmeared        =top1TrueMass-top1SmearedMass;
      W1DeltaPtTrueSmeared         =W1TruePt-W1SmearedPt; 
      W1DeltaRTrueSmeared          =deltaR(WPlus,recoWPlus); 
      W1DeltaMTrueSmeared          =W1TrueMass-W1SmearedMass;
      bJet1DeltaPtTrueSmeared      =bJet1TruePt-bJet1SmearedPt; 
      bJet1DeltaRTrueSmeared       =deltaR(bottomQuark,smearedBottomQuark); 
      W1Daughter1DeltaPtTrueSmeared=0;
      W1Daughter1DeltaRTrueSmeared =0;
      W1Daughter2DeltaPtTrueSmeared=-1;
      W1Daughter2DeltaRTrueSmeared =-1;
      lightJet1DeltaPxTrueSmeared  =lightJet1TruePx-lightJet1SmearedPx;
      lightJet1DeltaPyTrueSmeared  =lightJet1TruePy-lightJet1SmearedPy;
      lightJet1DeltaRTrueSmeared   =deltaR(lightParton1,smearedLightParton1);

      top2DeltaPtTrueSmeared       =top2TruePt-top2SmearedPt;
      top2DeltaRTrueSmeared        =deltaR(antiTopQuark,recoAntiTopQuark);
      top2DeltaMTrueSmeared        =top2TrueMass-top2SmearedMass;
      W2DeltaPtTrueSmeared         =W2TruePt-W2SmearedPt;
      W2DeltaRTrueSmeared          =deltaR(WMinus,recoWMinus);
      W2DeltaMTrueSmeared          =W2TrueMass-W2SmearedMass;
      bJet2DeltaPtTrueSmeared      =bJet2TruePt-bJet2SmearedPt;
      bJet2DeltaRTrueSmeared       =deltaR(antiBottomQuark,smearedAntiBottomQuark);
      W2Daughter1DeltaPtTrueSmeared=W2Daughter1TruePt-W2Daughter1SmearedPt;
      W2Daughter1DeltaRTrueSmeared =deltaR(antiNeutrino,smearedLightQuark2);
      W2Daughter2DeltaPtTrueSmeared=W2Daughter2TruePt-W2Daughter2SmearedPt;
      W2Daughter2DeltaRTrueSmeared =deltaR(lepton,smearedLightQuark1);
      lightJet2DeltaPxTrueSmeared  =lightJet2TruePx-lightJet2SmearedPx;
      lightJet2DeltaPyTrueSmeared  =lightJet2TruePy-lightJet2SmearedPy;
      lightJet2DeltaRTrueSmeared   =deltaR(lightParton2,smearedLightParton2);



      //Compare true and corrected
      top1DeltaPtTrueBest       =top1TruePt-top1BestPt;
      top1DeltaRTrueBest        =deltaR(topQuark,top1Best);
      top1DeltaMTrueBest        =top1TrueMass-top1BestMass;
      W1DeltaPtTrueBest         =W1TruePt-W1BestPt;
      W1DeltaRTrueBest          =deltaR(WPlus,W1Best);
      W1DeltaMTrueBest          =W1TrueMass-W1BestMass;
      bJet1DeltaPtTrueBest      =bJet1TruePt-bJet1BestPt;
      bJet1DeltaRTrueBest       =deltaR(bottomQuark,bJet1Best);
      W1Daughter1DeltaPtTrueBest=W1Daughter1TruePt-W1Daughter1BestPt;
      W1Daughter1DeltaRTrueBest =deltaR(antiLepton,W1Daughter1Best);
      W1Daughter2DeltaPtTrueBest=W1Daughter2TruePt-W1Daughter2BestPt;
      W1Daughter2DeltaRTrueBest =deltaR(neutrino,W1Daughter2Best);
      lightJet1DeltaPxTrueBest  =lightJet1TruePx-lightJet1BestPx;
      lightJet1DeltaPyTrueBest  =lightJet1TruePy-lightJet1BestPy;
      lightJet1DeltaRTrueBest   =deltaR(lightParton1,lightJet1Best);

      top2DeltaPtTrueBest       =top2TruePt-top2BestPt;
      top2DeltaRTrueBest        =deltaR(antiTopQuark,top2Best);
      top2DeltaMTrueBest        =top2TrueMass-top2BestMass;
      W2DeltaPtTrueBest         =W2TruePt-W2BestPt;
      W2DeltaRTrueBest          =deltaR(WMinus,W2Best);
      W2DeltaMTrueBest          =W2TrueMass-W2BestMass;
      bJet2DeltaPtTrueBest      =bJet2TruePt-bJet2BestPt;
      bJet2DeltaRTrueBest       =deltaR(antiBottomQuark,bJet2Best);
      W2Daughter1DeltaPtTrueBest=W2Daughter1TruePt-W2Daughter1BestPt;
      W2Daughter1DeltaRTrueBest =deltaR(antiNeutrino,W2Daughter1Best);
      W2Daughter2DeltaPtTrueBest=W2Daughter2TruePt-W2Daughter2BestPt;
      W2Daughter2DeltaRTrueBest =deltaR(lepton,W2Daughter2Best);
      lightJet2DeltaPxTrueBest  =lightJet2TruePx-lightJet2BestPx;
      lightJet2DeltaPyTrueBest  =lightJet2TruePy-lightJet2BestPy;
      lightJet2DeltaRTrueBest   =deltaR(lightParton2,lightJet2Best);


      //Compare smeared and corrected
      top1DeltaPtSmearedBest       =top1SmearedPt-top1BestPt;
      top1DeltaRSmearedBest        =deltaR(recoTopQuark,top1Best);
      top1DeltaMSmearedBest        =top1SmearedMass-top1BestMass;
      W1DeltaPtSmearedBest         =W1SmearedPt-W1BestPt;
      W1DeltaRSmearedBest          =deltaR(recoWPlus,W1Best);
      W1DeltaMSmearedBest          =W1SmearedMass-W1BestMass;
      bJet1DeltaPtSmearedBest      =bJet1SmearedPt-bJet1BestPt;
      bJet1DeltaRSmearedBest       =deltaR(smearedBottomQuark,bJet1Best);
      W1Daughter1DeltaPtSmearedBest=W1Daughter1SmearedPt-W1Daughter1BestPt;
      W1Daughter1DeltaRSmearedBest =deltaR(antiLepton,W1Daughter1Best);
      W1Daughter2DeltaPtSmearedBest=-1;
      W1Daughter2DeltaRSmearedBest =-1;
      lightJet1DeltaPxSmearedBest  =lightJet1SmearedPx-lightJet1BestPx;
      lightJet1DeltaPySmearedBest  =lightJet1SmearedPy-lightJet1BestPy;
      lightJet1DeltaRSmearedBest   =deltaR(smearedLightParton1,lightJet1Best);

      top2DeltaPtSmearedBest       =top2SmearedPt-top2BestPt;
      top2DeltaRSmearedBest        =deltaR(recoAntiTopQuark,top2Best);
      top2DeltaMSmearedBest        =top2SmearedMass-top2BestMass;
      W2DeltaPtSmearedBest         =W2SmearedPt-W2BestPt;
      W2DeltaRSmearedBest          =deltaR(recoWMinus,W2Best);
      W2DeltaMSmearedBest          =W2SmearedMass-W2BestMass;
      bJet2DeltaPtSmearedBest      =bJet2SmearedPt-bJet2BestPt;
      bJet2DeltaRSmearedBest       =deltaR(smearedAntiBottomQuark,bJet2Best);
      W2Daughter1DeltaPtSmearedBest=W2Daughter1SmearedPt-W2Daughter1BestPt;
      W2Daughter1DeltaRSmearedBest =deltaR(smearedLightQuark2,W2Daughter1Best);
      W2Daughter2DeltaPtSmearedBest=W2Daughter2SmearedPt-W2Daughter2BestPt;
      W2Daughter2DeltaRSmearedBest =deltaR(smearedLightQuark1,W2Daughter2Best);
      lightJet2DeltaPxSmearedBest  =lightJet2SmearedPx-lightJet2BestPx;
      lightJet2DeltaPySmearedBest  =lightJet2SmearedPy-lightJet2BestPy;
      lightJet2DeltaRSmearedBest   =deltaR(smearedLightParton2,lightJet2Best);

      
      outTree->Fill();

   } //end loop over events

   outFile->cd();
   outFile->Write();
   outFile->Close();
   
}
