#define ttbarReconstructionFromLHE_cxx
#include "ttbarReconstructionFromLHE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include "topSystemChiSquare.h"

using namespace std;
using namespace ROOT::Math;


void ttbarReconstructionFromLHE::printVector(XYZTLorentzVector& v)
{
  cout << "px = " << v.Px() << endl;
  cout << "py = " << v.Py() << endl;
  cout << "pz = " << v.Pz() << endl;
  cout << "E  = " << v.E () << endl;
}

void ttbarReconstructionFromLHE::Loop(int nSmearings, int whichLoop, int maxLoops)
{
//   In a ROOT session, you can do:
//      Root > .L ttbarReconstructionFromLHE.C
//      Root > ttbarReconstructionFromLHE t
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
   rand.SetSeed(whichLoop+40);

   initOutput(whichLoop);

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;

   int jStart = whichLoop*(nentries/maxLoops) + ((whichLoop>(maxLoops-nentries%maxLoops))?(whichLoop+nentries%maxLoops-maxLoops):0);

   //for(int iStart = 0; iStart < whichLoop; iStart++)
   //  {
   //    jStart += (nentries+iStart)/maxLoops;
   //  }

   int jFinish = jStart + (nentries+whichLoop)/maxLoops;

   for (Long64_t jentry=jStart; jentry<jFinish;jentry++) {
     cout << "BEGINNING BRANCH NUMBER " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
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
	}

      recoWMinus += lepton;
      recoWMinus += antiNeutrino;
      recoAntiTopQuark+=recoWMinus;

      recoWPlus += antiLepton;
      recoWPlus += neutrino;
      recoTopQuark+=recoWPlus;

      //Increase statistics by smearing one event multiple times
      for( int is = 0; is < nSmearings; is++ )
	{
	  eventNumber = jentry;
          smearingIndex = is; 
	  
	  METx -= lepton.Px() + antiLepton.Px();
	  METy -= lepton.Py() + antiLepton.Py();
	  vector<XYZTLorentzVector> recoJets;
	  int i_bottomQuark(0), i_antiBottomQuark(1);
	  int i_lightParton1(2), i_lightParton2(3);
	  XYZTLorentzVector smearedBottomQuark;
	  XYZTLorentzVector smearedAntiBottomQuark;
	  XYZTLorentzVector smearedLightParton1;
	  XYZTLorentzVector smearedLightParton2;
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

	  //cout << "bJet1PtSmear      :" << bJet1PtSmear      << endl;
	  //cout << "bJet2PtSmear      :" << bJet2PtSmear      << endl;
	  //cout << "lightJet1PtSmear  :" << lightJet1PtSmear  << endl;
	  //cout << "lightJet2PtSmear  :" << lightJet2PtSmear  << endl;
	  //cout << "bJet1PhiSmear     :" << bJet1PhiSmear     << endl;
	  //cout << "bJet2PhiSmear     :" << bJet2PhiSmear     << endl;
	  //cout << "lightJet1PhiSmear :" << lightJet1PhiSmear << endl;
	  //cout << "lightJet2PhiSmear :" << lightJet2PhiSmear << endl;
	  //cout << "bJet1EtaSmear     :" << bJet1EtaSmear     << endl;
	  //cout << "bJet2EtaSmear     :" << bJet2EtaSmear     << endl;
	  
	  vector<double> jetPtResolutions(4,0.1);
	  vector<double> jetPhiResolutions(4,0.01);
	  vector<double> jetEtaResolutions(4,0.01);
	  
	  double pt(bottomQuark.Pt()*exp(bJet1PtSmear));
	  double phi(bottomQuark.Phi()+bJet1PhiSmear);
	  double eta(bottomQuark.Eta()+bJet1EtaSmear);
	  double energy(sqrt((bottomQuark.M2()+bottomQuark.Pt()*bottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet1PtSmear));
	  
	  smearedBottomQuark.SetPx(pt*cos(phi));
	  smearedBottomQuark.SetPy(pt*sin(phi));
	  smearedBottomQuark.SetPz(pt*sinh(eta));
	  smearedBottomQuark.SetE(energy);
	  
	  recoTopQuark+=smearedBottomQuark;
	  
	  //cout << "Bottom Quark:"
	  //	   << "\nPx: " << bottomQuark.Px()
	  //	   << "\nPy: " << bottomQuark.Py()
	  //	   << "\nPz: " << bottomQuark.Pz()
	  //	   << "\nE : " << bottomQuark.E() << endl;
	  //
	  //cout << "Smeared Bottom Quark:"
	  //	   << "\nPx: " << smearedBottomQuark.Px()
	  //	   << "\nPy: " << smearedBottomQuark.Py()
	  //	   << "\nPz: " << smearedBottomQuark.Pz()
	  //	   << "\nE : " << smearedBottomQuark.E() << endl;
	  
	  METx -= smearedBottomQuark.Px();
	  METy -= smearedBottomQuark.Py();
	  //METx += bottomQuark.Px();
	  //METy += bottomQuark.Py();
	  
	  pt = antiBottomQuark.Pt()*exp(bJet2PtSmear);
	  phi = antiBottomQuark.Phi()+bJet2PhiSmear;
	  eta = antiBottomQuark.Eta()+bJet2EtaSmear;
	  energy = sqrt((antiBottomQuark.M2()+antiBottomQuark.Pt()*antiBottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet2PtSmear);
	  
	  smearedAntiBottomQuark.SetPx(pt*cos(phi));
	  smearedAntiBottomQuark.SetPy(pt*sin(phi));
	  smearedAntiBottomQuark.SetPz(pt*sinh(eta));
	  smearedAntiBottomQuark.SetE(energy);
	  
	  recoAntiTopQuark+=smearedAntiBottomQuark;
	  
	  //cout << "Anti Bottom Quark:"
	  //	   << "\nPx: " << antiBottomQuark.Px()
	  //	   << "\nPy: " << antiBottomQuark.Py()
	  //	   << "\nPz: " << antiBottomQuark.Pz()
	  //	   << "\nE : " << antiBottomQuark.E() << endl;
	  //
	  //cout << "Smeared Anti Bottom Quark:"
	  //	   << "\nPx: " << smearedAntiBottomQuark.Px()
	  //	   << "\nPy: " << smearedAntiBottomQuark.Py()
	  //	   << "\nPz: " << smearedAntiBottomQuark.Pz()
	  //	   << "\nE : " << smearedAntiBottomQuark.E() << endl;
	  
	  METx -= smearedAntiBottomQuark.Px();
	  METy -= smearedAntiBottomQuark.Py();
	  //METx += antiBottomQuark.Px();
	  //METy += antiBottomQuark.Py();
	  
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
	  //METx += lightParton1.Px();
	  //METy += lightParton1.Py();
	  
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
	  //METx += lightParton2.Px();
	  //METy += lightParton2.Py();
	  
	  recoJets.push_back(smearedBottomQuark    );
	  recoJets.push_back(smearedAntiBottomQuark);
	  recoJets.push_back(smearedLightParton1   );
	  recoJets.push_back(smearedLightParton2   );

	  //generator-level 4-momentum
	  genTopQuarkPt  = topQuark.Pt();
	  genTopQuarkEta = topQuark.Eta();
	  genTopQuarkPhi = topQuark.Phi();
	  genTopQuarkE   = topQuark.E();
	  genTopQuarkM   = sqrt(topQuark.M2());
	  genAntiTopQuarkPt  = antiTopQuark.Pt();
	  genAntiTopQuarkEta = antiTopQuark.Eta();
	  genAntiTopQuarkPhi = antiTopQuark.Phi();
	  genAntiTopQuarkE   = antiTopQuark.E();
	  genAntiTopQuarkM   = sqrt(antiTopQuark.M2());
	  
	  genWPlusPt  = WPlus.Pt();
	  genWPlusEta = WPlus.Eta();
	  genWPlusPhi = WPlus.Phi();
	  genWPlusE   = WPlus.E();
	  genWPlusM   = sqrt(WPlus.M2());
	  genWMinusPt  = WMinus.Pt();
	  genWMinusEta = WMinus.Eta();
	  genWMinusPhi = WMinus.Phi();
	  genWMinusE   = WMinus.E();
	  genWMinusM   = sqrt(WMinus.M2());
	  
	  genBottomQuarkPt  = bottomQuark.Pt();
	  genBottomQuarkEta = bottomQuark.Eta();
	  genBottomQuarkPhi = bottomQuark.Phi();
	  genBottomQuarkE   = bottomQuark.E();
	  genAntiBottomQuarkPt  = antiBottomQuark.Pt();
	  genAntiBottomQuarkEta = antiBottomQuark.Eta();
	  genAntiBottomQuarkPhi = antiBottomQuark.Phi();
	  genAntiBottomQuarkE   = antiBottomQuark.E();
	  
	  genNeutrinoPt  = neutrino.Pt();
	  genNeutrinoEta = neutrino.Eta();
	  genNeutrinoPhi = neutrino.Phi();
	  genNeutrinoE   = neutrino.E();
	  genAntiNeutrinoPt  = antiNeutrino.Pt();
	  genAntiNeutrinoEta = antiNeutrino.Eta();
	  genAntiNeutrinoPhi = antiNeutrino.Phi();
	  genAntiNeutrinoE   = antiNeutrino.E();

	  genLeptonPt  = lepton.Pt();
	  genLeptonEta = lepton.Eta();
	  genLeptonPhi = lepton.Phi();
	  genLeptonE   = lepton.E();
	  genAntiLeptonPt  = antiLepton.Pt();
	  genAntiLeptonEta = antiLepton.Eta();
	  genAntiLeptonPhi = antiLepton.Phi();
	  genAntiLeptonE   = antiLepton.E();
	  
	  genLightParton1Pt  = lightParton1.Pt();
	  genLightParton1Eta = lightParton1.Eta();
	  genLightParton1Phi = lightParton1.Phi();
	  genLightParton1E   = lightParton1.E();
	  genLightParton2Pt  = lightParton2.Pt();
	  genLightParton2Eta = lightParton2.Eta();
	  genLightParton2Phi = lightParton2.Phi();
	  genLightParton2E   = lightParton2.E();


	  //smeared 4-momentum
	  smearedTopQuarkPt  = recoTopQuark.Pt();
	  smearedTopQuarkEta = recoTopQuark.Eta();
	  smearedTopQuarkPhi = recoTopQuark.Phi();
	  smearedTopQuarkE   = recoTopQuark.E();
	  smearedTopQuarkM   = sqrt(recoTopQuark.M2());
	  smearedAntiTopQuarkPt  = recoAntiTopQuark.Pt();
	  smearedAntiTopQuarkEta = recoAntiTopQuark.Eta();
	  smearedAntiTopQuarkPhi = recoAntiTopQuark.Phi();
	  smearedAntiTopQuarkE   = recoAntiTopQuark.E();
	  smearedAntiTopQuarkM   = sqrt(recoAntiTopQuark.M2());
	  
	  smearedWPlusPt  = recoWPlus.Pt();
	  smearedWPlusEta = recoWPlus.Eta();
	  smearedWPlusPhi = recoWPlus.Phi();
	  smearedWPlusE   = recoWPlus.E();
	  smearedWPlusM   = sqrt(recoWPlus.M2());
	  smearedWMinusPt  = recoWMinus.Pt();
	  smearedWMinusEta = recoWMinus.Eta();
	  smearedWMinusPhi = recoWMinus.Phi();
	  smearedWMinusE   = recoWMinus.E();
	  smearedWMinusM   = sqrt(recoWMinus.M2());
	  
	  smearedBottomQuarkPt  = smearedBottomQuark.Pt();
	  smearedBottomQuarkEta = smearedBottomQuark.Eta();
	  smearedBottomQuarkPhi = smearedBottomQuark.Phi();
	  smearedBottomQuarkE   = smearedBottomQuark.E();
	  smearedAntiBottomQuarkPt  = smearedAntiBottomQuark.Pt();
	  smearedAntiBottomQuarkEta = smearedAntiBottomQuark.Eta();
	  smearedAntiBottomQuarkPhi = smearedAntiBottomQuark.Phi();
	  smearedAntiBottomQuarkE   = smearedAntiBottomQuark.E();
	  //deltaRBottomQuarkGenSmeared     = deltaR(    bottomQuark,smearedBottomQuark);
	  //deltaRAntiBottomQuarkGenSmeared = deltaR(antiBottomQuark,smearedAntiBottomQuark);
	  
	  smearedLightParton1Pt  = smearedLightParton1.Pt();
	  smearedLightParton1Eta = smearedLightParton1.Eta();
	  smearedLightParton1Phi = smearedLightParton1.Phi();
	  smearedLightParton1E   = smearedLightParton1.E();
	  smearedLightParton2Pt  = smearedLightParton2.Pt();
	  smearedLightParton2Eta = smearedLightParton2.Eta();
	  smearedLightParton2Phi = smearedLightParton2.Phi();
	  smearedLightParton2E   = smearedLightParton2.E();
	  //deltaRLightParton1GenSmeared = deltaR(lightParton1,smearedLightParton1);
	  //deltaRLightParton2GenSmeared = deltaR(lightParton2,smearedLightParton2);

	  deltaPtGenSmearedTopQuark = recoTopQuark.Pt() - topQuark.Pt();
	  deltaRGenSmearedTopQuark  = deltaR(recoTopQuark,topQuark);
	  deltaMGenSmearedTopQuark  = sqrt(recoTopQuark.M2()) - sqrt(topQuark.M2());
	  deltaPtGenSmearedAntiTopQuark = recoAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenSmearedAntiTopQuark  = deltaR(recoAntiTopQuark,antiTopQuark);
          deltaMGenSmearedAntiTopQuark  = sqrt(recoAntiTopQuark.M2()) -sqrt(antiTopQuark.M2());

	  deltaPtGenSmearedWPlus = recoWPlus.Pt() - WPlus.Pt();
          deltaRGenSmearedWPlus  = deltaR(recoWPlus,WPlus);
          deltaMGenSmearedWPlus  = sqrt(recoWPlus.M2()) - sqrt(WPlus.M2());
	  deltaPtGenSmearedWMinus = recoWMinus.Pt() - WMinus.Pt();
          deltaRGenSmearedWMinus  = deltaR(recoWMinus,WMinus);
          deltaMGenSmearedWMinus  = sqrt(recoWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenSmearedBottomQuark = smearedBottomQuark.Pt() - bottomQuark.Pt();
          deltaRGenSmearedBottomQuark  = deltaR(smearedBottomQuark,bottomQuark);
	  deltaPtGenSmearedAntiBottomQuark = smearedAntiBottomQuark.Pt() - antiBottomQuark.Pt();
          deltaRGenSmearedAntiBottomQuark  = deltaR(smearedAntiBottomQuark,antiBottomQuark);

	  deltaPtGenSmearedLightParton1 = smearedLightParton1.Pt() - lightParton1.Pt();
	  deltaPxGenSmearedLightParton1 = smearedLightParton1.Px() - lightParton1.Px();
	  deltaPyGenSmearedLightParton1 = smearedLightParton1.Py() - lightParton1.Py();
          deltaRGenSmearedLightParton1  = deltaR(smearedLightParton1,lightParton1);
	  deltaPtGenSmearedLightParton2 = smearedLightParton2.Pt() - lightParton2.Pt();
	  deltaPxGenSmearedLightParton2 = smearedLightParton2.Px() - lightParton2.Px();
          deltaPyGenSmearedLightParton2 = smearedLightParton2.Py() - lightParton2.Py();
          deltaRGenSmearedLightParton2  = deltaR(smearedLightParton2,lightParton2);

	  deltaRGenBottomQuarks = deltaR(bottomQuark,antiBottomQuark);
	  deltaRSmearedBottomQuarks = deltaR(smearedBottomQuark,smearedAntiBottomQuark);

	  deltaRGenLightPartons = deltaR(lightParton1,lightParton2);
	  deltaRSmearedLightPartons = deltaR(smearedLightParton1,smearedLightParton2);

	  deltaRGenTopQuarks = deltaR(topQuark,antiTopQuark);
          deltaRSmearedTopQuarks = deltaR(recoTopQuark,recoAntiTopQuark);

	  deltaRGenWBosons = deltaR(WPlus,WMinus);
	  deltaRSmearedWBosons = deltaR(recoWPlus,recoWMinus);

	  deltaRGenNeutrinos = deltaR(neutrino,antiNeutrino);
	  	  
	  //if(jentry !=9) continue;
	  
	  cout << "METx is " << METx << endl;
	  cout << "real METx is " << neutrino.Px() + antiNeutrino.Px() << endl;
	  cout << "METy is " << METy << endl;
	  //cout << "real METy is " << neutrino.Py() + antiNeutrino.Py() << endl;
	  //cout << "neutrino       mass is      " << neutrino.M() << endl;
	  //printVector(neutrino);
	  //cout << "anti-lepton    mass is      " << antiLepton.M() << endl;
	  //printVector(antiLepton);
	  cout << "reco W+        mass is      " << recoWPlus.M() << endl;
	  //printVector(recoWPlus);
	  //cout << "b quark        mass is      " << smearedBottomQuark.M() << endl;
	  //printVector(smearedBottomQuark);
	  cout << "reco top quark mass is      " << recoTopQuark.M() << endl;
	  //printVector(recoTopQuark);
	  //cout << "anti-neutrino  mass is      " << antiNeutrino.M() << endl;
	  //printVector(antiNeutrino);
	  //cout << "lepton         mass is      " << lepton.M() << endl;
	  //printVector(antiLepton);
	  cout << "reco W-        mass is      " << recoWMinus.M() << endl;
	  //printVector(recoWMinus);
	  //cout << "b-bar quark    mass is      " << smearedAntiBottomQuark.M() << endl;
	  //printVector(smearedAntiBottomQuark);
	  cout << "reco anti top quark mass is " << recoAntiTopQuark.M() << endl;
	  //printVector(recoAntiTopQuark);
	  
	  //double bQuarkPx(smearedBottomQuark.Px());
	  //double bQuarkPy(smearedBottomQuark.Py());
	  //double bQuarkPz(smearedBottomQuark.Pz());
	  //double bQuarkP (smearedBottomQuark.P ());
	  //double bQuarkE (smearedBottomQuark.E ());
	  //
	  //double antiLeptonPx(antiLepton.Px());
	  //double antiLeptonPy(antiLepton.Py());
	  //double antiLeptonPz(antiLepton.Pz());
	  //double antiLeptonE (antiLepton.E ());
	  //
	  //double bBarQuarkPx(smearedAntiBottomQuark.Px());
	  //double bBarQuarkPy(smearedAntiBottomQuark.Py());
	  //double bBarQuarkPz(smearedAntiBottomQuark.Pz());
	  //double bBarQuarkP (smearedAntiBottomQuark.P ());
	  //double bBarQuarkE (smearedAntiBottomQuark.E ());
	  //
	  //double leptonPx(lepton.Px());
	  //double leptonPy(lepton.Py());
	  //double leptonPz(lepton.Pz());
	  //double leptonE (lepton.E ());
	  //
	  //double ptWidth(log(1+0.1)),phiWidth(0.01),etaWidth(0.01);
	  //
	  //NeutrinoEllipseCalculator nu(bQuarkPx,bQuarkPy,bQuarkPz,bQuarkE,
	  //				   ptWidth,phiWidth,etaWidth,
	  //				   antiLeptonPx,antiLeptonPy,antiLeptonPz,antiLeptonE,
	  //				   bQuarkPx,bQuarkPy,bQuarkPz,
	  //				   bQuarkP, bQuarkE);
	  //
	  //cout << "for top quark to neutrino we have:\n"
	  //	   << "for actual masses Z2 = " << nu.getZ2(recoTopQuark.M(),recoWPlus.M(),0.)
	  //	   << "\nfor nominal masses Z2 = " << nu.getZ2(173.,80.4,0.) << endl;
	  //	
	  //
	  //NeutrinoEllipseCalculator nuBar(bBarQuarkPx,bBarQuarkPy,bBarQuarkPz,bBarQuarkE,
	  //				      ptWidth,phiWidth,etaWidth,
	  //				      leptonPx,leptonPy,leptonPz,leptonE,
	  //				      bBarQuarkPx,bBarQuarkPy,bBarQuarkPz,
	  //				      bBarQuarkP, bBarQuarkE);
	  //
	  //cout << "for anti top quark to anti neutrino we have:\n"
	  //	   << "for actual masses Z2 = " << nuBar.getZ2(recoAntiTopQuark.M(),recoWMinus.M(),0.)
	  //	   << "\nfor nominal masses Z2 = " << nuBar.getZ2(173.,80.4,0.) << endl;
	  //
	  //return;
	  
	  //continue;

	  XYZTLorentzVector MET(METx,METy,0,0);
	  topSystemChiSquare system(recoJets,
				    jetPtResolutions,
				    jetPhiResolutions,
				    jetEtaResolutions,
				    i_bottomQuark, i_lightParton1, antiLepton, 
				    i_antiBottomQuark, i_lightParton2, lepton,
				    MET,
				    173, 2.0 ,
				    80.4, 2.09);
	  //if(system.unstretchedBJetPtExists()) 
	  //  {
	  //    system.minimizeLightJetChiSquareOnly();
	  //    system.plotEllipses(TString("plotCheck_noSmearing_")+jentry);
	  //  }
	  //system.minimizeLightJetChiSquareOnly();
	  //cout << "Checking the top and W masses at the intersections of the 'measured' ellipses" << endl;
	  //system.checkMasses();

	  chi2 = system.calcSolution();
	  //cout << "Checking the top and W masses at the intersections of the 'best' ellipses" << endl;
	  //system.checkMasses();
	  //if(chi2 > 0) system.plotEllipses(TString("plotCheck_")+jentry);
	  
	  minStatusCode = system.getMinimizerStatusCode();


	  //"best" 4-momentum from minimization
	  int iSol;
	  XYZTLorentzVector bestTopQuark, bestWPlus, bestBottomQuark, bestLightParton1;
	  XYZTLorentzVector bestAntiTopQuark, bestWMinus, bestAntiBottomQuark, bestLightParton2;
	  XYZTLorentzVector bestNeutrino, bestAntiNeutrino;
	  vector<XYZTLorentzVector> bestNeutrinos, bestAntiNeutrinos;
	  vector<int> bestNeutrinoPairingTruths, bestAntiNeutrinoPairingTruths;
	  XYZTLorentzVector bestMET;
	  vector<double> minDeltas;
	  
	  system.fillBestMomenta(bestBottomQuark,bestAntiBottomQuark,
				 bestLightParton1,bestLightParton2,
				 bestNeutrinos,bestNeutrinoPairingTruths,
				 bestAntiNeutrinos,bestAntiNeutrinoPairingTruths,
				 bestMET,bestPairing,minDeltas);


	  bestBottomQuarkPt  = bestBottomQuark.Pt();
          bestBottomQuarkEta = bestBottomQuark.Eta();
          bestBottomQuarkPhi = bestBottomQuark.Phi();
          bestBottomQuarkE   = bestBottomQuark.E();
          bestAntiBottomQuarkPt  = bestAntiBottomQuark.Pt();
          bestAntiBottomQuarkEta = bestAntiBottomQuark.Eta();
          bestAntiBottomQuarkPhi = bestAntiBottomQuark.Phi();
          bestAntiBottomQuarkE   = bestAntiBottomQuark.E();

          deltaPtGenBestBottomQuark = bestBottomQuark.Pt() - bottomQuark.Pt();
          deltaRGenBestBottomQuark  = deltaR(bestBottomQuark,bottomQuark);
          deltaPtGenBestAntiBottomQuark = bestAntiBottomQuark.Pt() - antiBottomQuark.Pt();
	  deltaRGenBestAntiBottomQuark  = deltaR(bestAntiBottomQuark,antiBottomQuark);

          deltaPtSmearedBestBottomQuark = bestBottomQuark.Pt() - smearedBottomQuark.Pt();
          deltaRSmearedBestBottomQuark  = deltaR(smearedBottomQuark,bestBottomQuark);
          deltaPtSmearedBestAntiBottomQuark = bestAntiBottomQuark.Pt() - smearedAntiBottomQuark.Pt();
          deltaRSmearedBestAntiBottomQuark  = deltaR(smearedAntiBottomQuark,bestAntiBottomQuark);

	  deltaRBestBottomQuarks = deltaR(bestBottomQuark,bestAntiBottomQuark);

	  bestLightParton1Pt  = bestLightParton1.Pt();
          bestLightParton1Eta = bestLightParton1.Eta();
          bestLightParton1Phi = bestLightParton1.Phi();
          bestLightParton1E   = bestLightParton1.E();
          bestLightParton2Pt  = bestLightParton2.Pt();
          bestLightParton2Eta = bestLightParton2.Eta();
          bestLightParton2Phi = bestLightParton2.Phi();
          bestLightParton2E   = bestLightParton2.E();

          deltaPtGenBestLightParton1 = bestLightParton1.Pt() - lightParton1.Pt();
	  deltaPxGenBestLightParton1 = bestLightParton1.Px() - lightParton1.Px();
          deltaPyGenBestLightParton1 = bestLightParton1.Py() - lightParton1.Py();
          deltaRGenBestLightParton1  = deltaR(bestLightParton1,lightParton1);
          deltaPtGenBestLightParton2 = bestLightParton2.Pt() - lightParton2.Pt();
	  deltaPxGenBestLightParton2 = bestLightParton2.Px() - lightParton2.Px();
	  deltaPyGenBestLightParton2 = bestLightParton2.Py() - lightParton2.Py();
          deltaRGenBestLightParton2  = deltaR(bestLightParton2,lightParton2);

          deltaPtSmearedBestLightParton1 = bestLightParton1.Pt() - smearedLightParton1.Pt();
	  deltaPxSmearedBestLightParton1 = bestLightParton1.Px() - smearedLightParton1.Px();
          deltaPySmearedBestLightParton1 = bestLightParton1.Py() - smearedLightParton1.Py();
          deltaRSmearedBestLightParton1  = deltaR(smearedLightParton1,bestLightParton1);
          deltaPtSmearedBestLightParton2 = bestLightParton2.Pt() - smearedLightParton2.Pt();
	  deltaPxSmearedBestLightParton2 = bestLightParton2.Px() - smearedLightParton2.Px();                                                                                
          deltaPySmearedBestLightParton2 = bestLightParton2.Py() - smearedLightParton2.Py();
          deltaRSmearedBestLightParton2  = deltaR(smearedLightParton2,bestLightParton2);

	  deltaRBestLightPartons = deltaR(bestLightParton1,bestLightParton2);


	  //solution 1
	  iSol = 0;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino1PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino1PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	  
	  bestTopQuark1Pt  = bestTopQuark.Pt();
	  bestTopQuark1Eta = bestTopQuark.Eta();
	  bestTopQuark1Phi = bestTopQuark.Phi();
	  bestTopQuark1E   = bestTopQuark.E();
	  bestTopQuark1M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark1Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark1Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark1Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark1E   = bestAntiTopQuark.E();
	  bestAntiTopQuark1M   = sqrt(bestAntiTopQuark.M2());
	  
	  bestWPlus1Pt  = bestWPlus.Pt();
	  bestWPlus1Eta = bestWPlus.Eta();
	  bestWPlus1Phi = bestWPlus.Phi();
	  bestWPlus1E   = bestWPlus.E();
	  bestWPlus1M   = sqrt(bestWPlus.M2());
	  bestWMinus1Pt  = bestWMinus.Pt();
	  bestWMinus1Eta = bestWMinus.Eta();
	  bestWMinus1Phi = bestWMinus.Phi();
	  bestWMinus1E   = bestWMinus.E();
	  bestWMinus1M   = sqrt(bestWMinus.M2());
	  
	  bestNeutrino1Pt  = bestNeutrino.Pt();
	  bestNeutrino1Eta = bestNeutrino.Eta();
	  bestNeutrino1Phi = bestNeutrino.Phi();
	  bestNeutrino1E   = bestNeutrino.E();
	  bestAntiNeutrino1Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino1Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino1Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino1E   = bestAntiNeutrino.E();

          deltaPtGenBestTopQuark1 = bestTopQuark.Pt() - topQuark.Pt();
	  deltaRGenBestTopQuark1 = deltaR(bestTopQuark,topQuark);
	  deltaMGenBestTopQuark1 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
	  deltaPtGenBestAntiTopQuark1 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
	  deltaRGenBestAntiTopQuark1 = deltaR(bestAntiTopQuark,antiTopQuark);
	  deltaMGenBestAntiTopQuark1 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus1 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus1  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus1  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus1 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus1  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus1  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino1 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino1  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino1 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino1  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark1 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark1  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark1  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark1 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark1  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark1  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus1  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus1   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus1   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus1 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus1  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus1  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks1 = deltaR(bestTopQuark,bestAntiTopQuark);
	  deltaRBestWBosons1   = deltaR(bestWPlus,bestWMinus);
	  deltaRBestNeutrinos1 = deltaR(bestNeutrino,bestAntiNeutrino);

	  //solution 2
	  iSol = 1;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino2PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino2PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark2Pt  = bestTopQuark.Pt();
	  bestTopQuark2Eta = bestTopQuark.Eta();
	  bestTopQuark2Phi = bestTopQuark.Phi();
	  bestTopQuark2E   = bestTopQuark.E();
	  bestTopQuark2M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark2Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark2Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark2Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark2E   = bestAntiTopQuark.E();
	  bestAntiTopQuark2M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus2Pt  = bestWPlus.Pt();
	  bestWPlus2Eta = bestWPlus.Eta();
	  bestWPlus2Phi = bestWPlus.Phi();
	  bestWPlus2E   = bestWPlus.E();
	  bestWPlus2M   = sqrt(bestWPlus.M2());
	  bestWMinus2Pt  = bestWMinus.Pt();
	  bestWMinus2Eta = bestWMinus.Eta();
	  bestWMinus2Phi = bestWMinus.Phi();
	  bestWMinus2E   = bestWMinus.E();
	  bestWMinus2M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino2Pt  = bestNeutrino.Pt();
	  bestNeutrino2Eta = bestNeutrino.Eta();
	  bestNeutrino2Phi = bestNeutrino.Phi();
	  bestNeutrino2E   = bestNeutrino.E();
	  bestAntiNeutrino2Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino2Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino2Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino2E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark2 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark2 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark2 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark2 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark2 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark2 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus2 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus2  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus2  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus2 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus2  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus2  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino2 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino2  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino2 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino2  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark2 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark2  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark2  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark2 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark2  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark2  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus2  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus2   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus2   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus2 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus2  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus2  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks2 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons2   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos2 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 3
	  iSol = 2;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino3PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino3PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark3Pt  = bestTopQuark.Pt();
	  bestTopQuark3Eta = bestTopQuark.Eta();
	  bestTopQuark3Phi = bestTopQuark.Phi();
	  bestTopQuark3E   = bestTopQuark.E();
	  bestTopQuark3M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark3Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark3Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark3Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark3E   = bestAntiTopQuark.E();
	  bestAntiTopQuark3M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus3Pt  = bestWPlus.Pt();
	  bestWPlus3Eta = bestWPlus.Eta();
	  bestWPlus3Phi = bestWPlus.Phi();
	  bestWPlus3E   = bestWPlus.E();
	  bestWPlus3M   = sqrt(bestWPlus.M2());
	  bestWMinus3Pt  = bestWMinus.Pt();
	  bestWMinus3Eta = bestWMinus.Eta();
	  bestWMinus3Phi = bestWMinus.Phi();
	  bestWMinus3E   = bestWMinus.E();
	  bestWMinus3M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino3Pt  = bestNeutrino.Pt();
	  bestNeutrino3Eta = bestNeutrino.Eta();
	  bestNeutrino3Phi = bestNeutrino.Phi();
	  bestNeutrino3E   = bestNeutrino.E();
	  bestAntiNeutrino3Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino3Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino3Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino3E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark3 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark3 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark3 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark3 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark3 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark3 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus3 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus3  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus3  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus3 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus3  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus3  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino3 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino3  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino3 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino3  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark3 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark3  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark3  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark3 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark3  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark3  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus3  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus3   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus3   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus3 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus3  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus3  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks3 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons3   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos3 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 4
	  iSol = 3;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino4PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino4PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark4Pt  = bestTopQuark.Pt();
	  bestTopQuark4Eta = bestTopQuark.Eta();
	  bestTopQuark4Phi = bestTopQuark.Phi();
	  bestTopQuark4E   = bestTopQuark.E();
	  bestTopQuark4M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark4Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark4Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark4Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark4E   = bestAntiTopQuark.E();
	  bestAntiTopQuark4M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus4Pt  = bestWPlus.Pt();
	  bestWPlus4Eta = bestWPlus.Eta();
	  bestWPlus4Phi = bestWPlus.Phi();
	  bestWPlus4E   = bestWPlus.E();
	  bestWPlus4M   = sqrt(bestWPlus.M2());
	  bestWMinus4Pt  = bestWMinus.Pt();
	  bestWMinus4Eta = bestWMinus.Eta();
	  bestWMinus4Phi = bestWMinus.Phi();
	  bestWMinus4E   = bestWMinus.E();
	  bestWMinus4M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino4Pt  = bestNeutrino.Pt();
	  bestNeutrino4Eta = bestNeutrino.Eta();
	  bestNeutrino4Phi = bestNeutrino.Phi();
	  bestNeutrino4E   = bestNeutrino.E();
	  bestAntiNeutrino4Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino4Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino4Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino4E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark4 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark4 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark4 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark4 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark4 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark4 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus4 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus4  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus4  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus4 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus4  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus4  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino4 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino4  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino4 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino4  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark4 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark4  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark4  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark4 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark4  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark4  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus4  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus4   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus4   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus4 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus4  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus4  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks4 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons4   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos4 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 5
	  iSol = 4;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino5PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino5PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark5Pt  = bestTopQuark.Pt();
	  bestTopQuark5Eta = bestTopQuark.Eta();
	  bestTopQuark5Phi = bestTopQuark.Phi();
	  bestTopQuark5E   = bestTopQuark.E();
	  bestTopQuark5M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark5Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark5Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark5Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark5E   = bestAntiTopQuark.E();
	  bestAntiTopQuark5M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus5Pt  = bestWPlus.Pt();
	  bestWPlus5Eta = bestWPlus.Eta();
	  bestWPlus5Phi = bestWPlus.Phi();
	  bestWPlus5E   = bestWPlus.E();
	  bestWPlus5M   = sqrt(bestWPlus.M2());
	  bestWMinus5Pt  = bestWMinus.Pt();
	  bestWMinus5Eta = bestWMinus.Eta();
	  bestWMinus5Phi = bestWMinus.Phi();
	  bestWMinus5E   = bestWMinus.E();
	  bestWMinus5M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino5Pt  = bestNeutrino.Pt();
	  bestNeutrino5Eta = bestNeutrino.Eta();
	  bestNeutrino5Phi = bestNeutrino.Phi();
	  bestNeutrino5E   = bestNeutrino.E();
	  bestAntiNeutrino5Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino5Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino5Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino5E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark5 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark5 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark5 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark5 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark5 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark5 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus5 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus5  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus5  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus5 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus5  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus5  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino5 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino5  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino5 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino5  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark5 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark5  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark5  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark5 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark5  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark5  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus5  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus5   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus5   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus5 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus5  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus5  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks5 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons5   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos5 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 6
	  iSol = 5;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino6PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino6PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark6Pt  = bestTopQuark.Pt();
	  bestTopQuark6Eta = bestTopQuark.Eta();
	  bestTopQuark6Phi = bestTopQuark.Phi();
	  bestTopQuark6E   = bestTopQuark.E();
	  bestTopQuark6M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark6Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark6Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark6Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark6E   = bestAntiTopQuark.E();
	  bestAntiTopQuark6M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus6Pt  = bestWPlus.Pt();
	  bestWPlus6Eta = bestWPlus.Eta();
	  bestWPlus6Phi = bestWPlus.Phi();
	  bestWPlus6E   = bestWPlus.E();
	  bestWPlus6M   = sqrt(bestWPlus.M2());
	  bestWMinus6Pt  = bestWMinus.Pt();
	  bestWMinus6Eta = bestWMinus.Eta();
	  bestWMinus6Phi = bestWMinus.Phi();
	  bestWMinus6E   = bestWMinus.E();
	  bestWMinus6M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino6Pt  = bestNeutrino.Pt();
	  bestNeutrino6Eta = bestNeutrino.Eta();
	  bestNeutrino6Phi = bestNeutrino.Phi();
	  bestNeutrino6E   = bestNeutrino.E();
	  bestAntiNeutrino6Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino6Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino6Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino6E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark6 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark6 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark6 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark6 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark6 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark6 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus6 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus6  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus6  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus6 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus6  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus6  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino6 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino6  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino6 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino6  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark6 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark6  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark6  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark6 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark6  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark6  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus6  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus6   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus6   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus6 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus6  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus6  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks6 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons6   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos6 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 7
	  iSol = 6;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino7PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino7PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark7Pt  = bestTopQuark.Pt();
	  bestTopQuark7Eta = bestTopQuark.Eta();
	  bestTopQuark7Phi = bestTopQuark.Phi();
	  bestTopQuark7E   = bestTopQuark.E();
	  bestTopQuark7M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark7Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark7Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark7Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark7E   = bestAntiTopQuark.E();
	  bestAntiTopQuark7M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus7Pt  = bestWPlus.Pt();
	  bestWPlus7Eta = bestWPlus.Eta();
	  bestWPlus7Phi = bestWPlus.Phi();
	  bestWPlus7E   = bestWPlus.E();
	  bestWPlus7M   = sqrt(bestWPlus.M2());
	  bestWMinus7Pt  = bestWMinus.Pt();
	  bestWMinus7Eta = bestWMinus.Eta();
	  bestWMinus7Phi = bestWMinus.Phi();
	  bestWMinus7E   = bestWMinus.E();
	  bestWMinus7M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino7Pt  = bestNeutrino.Pt();
	  bestNeutrino7Eta = bestNeutrino.Eta();
	  bestNeutrino7Phi = bestNeutrino.Phi();
	  bestNeutrino7E   = bestNeutrino.E();
	  bestAntiNeutrino7Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino7Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino7Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino7E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark7 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark7 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark7 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark7 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark7 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark7 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus7 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus7  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus7  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus7 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus7  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus7  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino7 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino7  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino7 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino7  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark7 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark7  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark7  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark7 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark7  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark7  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus7  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus7   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus7   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus7 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus7  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus7  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks7 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons7   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos7 = deltaR(bestNeutrino,bestAntiNeutrino);


	  //solution 8
	  iSol = 7;
	  bestNeutrino     = bestNeutrinos.at(iSol); 
	  bestAntiNeutrino = bestAntiNeutrinos.at(iSol);
	  bestNeutrino8PairingTruth     = bestNeutrinoPairingTruths.at(iSol);
	  bestAntiNeutrino8PairingTruth = bestAntiNeutrinoPairingTruths.at(iSol);

	  bestWPlus  = bestNeutrino + antiLepton;
	  bestWMinus = bestAntiNeutrino + lepton;
	  if(bestPairing)
	    {
	      bestTopQuark = bestWPlus + bestBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestAntiBottomQuark;
	    }
	  else
	    {
	      bestTopQuark = bestWPlus + bestAntiBottomQuark;
	      bestAntiTopQuark = bestWMinus + bestBottomQuark;
	    }
	    
	  bestTopQuark8Pt  = bestTopQuark.Pt();
	  bestTopQuark8Eta = bestTopQuark.Eta();
	  bestTopQuark8Phi = bestTopQuark.Phi();
	  bestTopQuark8E   = bestTopQuark.E();
	  bestTopQuark8M   = sqrt(bestTopQuark.M2());
	  bestAntiTopQuark8Pt  = bestAntiTopQuark.Pt();
	  bestAntiTopQuark8Eta = bestAntiTopQuark.Eta();
	  bestAntiTopQuark8Phi = bestAntiTopQuark.Phi();
	  bestAntiTopQuark8E   = bestAntiTopQuark.E();
	  bestAntiTopQuark8M   = sqrt(bestAntiTopQuark.M2());
	    
	  bestWPlus8Pt  = bestWPlus.Pt();
	  bestWPlus8Eta = bestWPlus.Eta();
	  bestWPlus8Phi = bestWPlus.Phi();
	  bestWPlus8E   = bestWPlus.E();
	  bestWPlus8M   = sqrt(bestWPlus.M2());
	  bestWMinus8Pt  = bestWMinus.Pt();
	  bestWMinus8Eta = bestWMinus.Eta();
	  bestWMinus8Phi = bestWMinus.Phi();
	  bestWMinus8E   = bestWMinus.E();
	  bestWMinus8M   = sqrt(bestWMinus.M2());
	    
	  bestNeutrino8Pt  = bestNeutrino.Pt();
	  bestNeutrino8Eta = bestNeutrino.Eta();
	  bestNeutrino8Phi = bestNeutrino.Phi();
	  bestNeutrino8E   = bestNeutrino.E();
	  bestAntiNeutrino8Pt  = bestAntiNeutrino.Pt();
	  bestAntiNeutrino8Eta = bestAntiNeutrino.Eta();
	  bestAntiNeutrino8Phi = bestAntiNeutrino.Phi();
	  bestAntiNeutrino8E   = bestAntiNeutrino.E();

	  deltaPtGenBestTopQuark8 = bestTopQuark.Pt() - topQuark.Pt();
          deltaRGenBestTopQuark8 = deltaR(bestTopQuark,topQuark);
          deltaMGenBestTopQuark8 = sqrt(bestTopQuark.M2()) - sqrt(topQuark.M2());
          deltaPtGenBestAntiTopQuark8 = bestAntiTopQuark.Pt() - antiTopQuark.Pt();
          deltaRGenBestAntiTopQuark8 = deltaR(bestAntiTopQuark,antiTopQuark);
          deltaMGenBestAntiTopQuark8 = sqrt(bestAntiTopQuark.M2()) - sqrt(antiTopQuark.M2());

          deltaPtGenBestWPlus8 = bestWPlus.Pt() - WPlus.Pt();
          deltaRGenBestWPlus8  = deltaR(bestWPlus,WPlus);
          deltaMGenBestWPlus8  = sqrt(bestWPlus.M2()) -sqrt(WPlus.M2());
          deltaPtGenBestWMinus8 = bestWMinus.Pt() - WMinus.Pt();
          deltaRGenBestWMinus8  = deltaR(bestWMinus,WMinus);
          deltaMGenBestWMinus8  = sqrt(bestWMinus.M2()) - sqrt(WMinus.M2());

	  deltaPtGenBestNeutrino8 = bestNeutrino.Pt() - neutrino.Pt();
          deltaRGenBestNeutrino8  = deltaR(bestNeutrino,neutrino);
          deltaPtGenBestAntiNeutrino8 = bestAntiNeutrino.Pt() - antiNeutrino.Pt();
          deltaRGenBestAntiNeutrino8  = deltaR(bestAntiNeutrino,antiNeutrino);

          deltaPtSmearedBestTopQuark8 = bestTopQuark.Pt() - recoTopQuark.Pt();
          deltaRSmearedBestTopQuark8  = deltaR(bestTopQuark,recoTopQuark);
          deltaMSmearedBestTopQuark8  = sqrt(bestTopQuark.M2()) - sqrt(recoTopQuark.M2());
          deltaPtSmearedBestAntiTopQuark8 = bestAntiTopQuark.Pt() - recoAntiTopQuark.Pt();
          deltaRSmearedBestAntiTopQuark8  = deltaR(bestAntiTopQuark,recoAntiTopQuark);
          deltaMSmearedBestAntiTopQuark8  = sqrt(bestAntiTopQuark.M2()) - sqrt(recoAntiTopQuark.M2());

          deltaPtSmearedBestWPlus8  = bestWPlus.Pt() - recoWPlus.Pt();
          deltaRSmearedBestWPlus8   = deltaR(bestWPlus,recoWPlus);
          deltaMSmearedBestWPlus8   = sqrt(bestWPlus.M2()) - sqrt(recoWPlus.M2());
          deltaPtSmearedBestWMinus8 = bestWMinus.Pt() - recoWMinus.Pt();
          deltaRSmearedBestWMinus8  = deltaR(bestWMinus,recoWMinus);
          deltaMSmearedBestWMinus8  = sqrt(bestWMinus.M2()) - sqrt(recoWMinus.M2());

	  deltaRBestTopQuarks8 = deltaR(bestTopQuark,bestAntiTopQuark);
          deltaRBestWBosons8   = deltaR(bestWPlus,bestWMinus);
          deltaRBestNeutrinos8 = deltaR(bestNeutrino,bestAntiNeutrino);

	  

	  if(minDeltas.size() == 14)
	    {
	      deltaPtBottomQuarkFromMin      = minDeltas.at(0);
	      deltaPhiBottomQuarkFromMin     = minDeltas.at(1);
	      deltaEtaBottomQuarkFromMin     = minDeltas.at(2);
	      deltaPtAntiBottomQuarkFromMin  = minDeltas.at(3);
	      deltaPhiAntiBottomQuarkFromMin = minDeltas.at(4);
	      deltaEtaAntiBottomQuarkFromMin = minDeltas.at(5);
	      deltaMTopQuarkFromMin          = minDeltas.at(6);
	      deltaMAntiTopQuarkFromMin      = minDeltas.at(7);
	      deltaMWPlusFromMin             = minDeltas.at(8);
	      deltaMWMinusFromMin            = minDeltas.at(9);
	      deltaPxLightParton1FromMin     = minDeltas.at(10);
	      deltaPyLightParton1FromMin     = minDeltas.at(11); 
	      deltaPxLightParton2FromMin     = minDeltas.at(12);
	      deltaPyLightParton2FromMin     = minDeltas.at(13);
	    }
	  
	  
	  outTree->Fill();

	} //end loop over smearings

   }

   outFile->cd();
   outFile->Write();
   outFile->Close();

}
