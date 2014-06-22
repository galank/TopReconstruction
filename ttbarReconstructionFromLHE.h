//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  8 15:44:02 2014 by ROOT version 5.34/04
// from TTree Physics/Physics
// found on file: lhe.root
//////////////////////////////////////////////////////////

#ifndef ttbarReconstructionFromLHE_h
#define ttbarReconstructionFromLHE_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "Math/GenVector/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class ttbarReconstructionFromLHE {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           n_particles;
   vector<int>     *PID;
   vector<double>  *P_X;
   vector<double>  *P_Y;
   vector<double>  *P_Z;
   vector<double>  *E;
   vector<double>  *M;
   vector<int>     *status;
   vector<int>     *particleID;
   vector<int>     *parent1ID;
   vector<int>     *parent2ID;

   // List of branches
   TBranch        *b_n_particles;   //!
   TBranch        *b_PID;   //!
   TBranch        *b_P_X;   //!
   TBranch        *b_P_Y;   //!
   TBranch        *b_P_Z;   //!
   TBranch        *b_E;   //!
   TBranch        *b_M;   //!
   TBranch        *b_status;   //!
   TBranch        *b_particleID;   //!
   TBranch        *b_parent1ID;   //!
   TBranch        *b_parent2ID;   //!

   ttbarReconstructionFromLHE(TString dirName, TTree *tree=0);
   virtual ~ttbarReconstructionFromLHE();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int nSmearings = 1, int whichLoop = 0, int maxLoops = 1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void printVector(XYZTLorentzVector&);

   double deltaR(XYZTLorentzVector& , XYZTLorentzVector&  );

   //store minimizer results in a tree;
   TString dirName_;
   TFile* outFile;
   TTree* outTree;
   void initOutput(int whichLoop);

   //from the minimization/quartic solver/smearing
   double chi2;
   bool bestPairing;
   int minStatusCode;
   int eventNumber;
   int smearingIndex;
   
   //compare true and smeared
   double deltaPtGenSmearedTopQuark     , deltaPtGenSmearedAntiTopQuark    ;
   double  deltaRGenSmearedTopQuark     ,  deltaRGenSmearedAntiTopQuark    ;
   double  deltaMGenSmearedTopQuark     ,  deltaMGenSmearedAntiTopQuark    ;
   double deltaPtGenSmearedBottomQuark  , deltaPtGenSmearedAntiBottomQuark ;
   double  deltaRGenSmearedBottomQuark  ,  deltaRGenSmearedAntiBottomQuark ;
   double deltaPtGenSmearedWPlus        , deltaPtGenSmearedWMinus          ;
   double  deltaRGenSmearedWPlus        ,  deltaRGenSmearedWMinus          ;
   double  deltaMGenSmearedWPlus        ,  deltaMGenSmearedWMinus          ;
   double deltaPtGenSmearedLightParton1 , deltaPtGenSmearedLightParton2    ;
   double deltaPxGenSmearedLightParton1 , deltaPxGenSmearedLightParton2    ;
   double deltaPyGenSmearedLightParton1 , deltaPyGenSmearedLightParton2    ;
   double  deltaRGenSmearedLightParton1 ,  deltaRGenSmearedLightParton2    ;

   //compare true and corrected
   double deltaPtGenBestTopQuark1    , deltaPtGenBestAntiTopQuark1   ;
   double  deltaRGenBestTopQuark1    ,  deltaRGenBestAntiTopQuark1   ;
   double  deltaMGenBestTopQuark1    ,  deltaMGenBestAntiTopQuark1   ;
   double deltaPtGenBestWPlus1       , deltaPtGenBestWMinus1         ;
   double  deltaRGenBestWPlus1       ,  deltaRGenBestWMinus1         ;
   double  deltaMGenBestWPlus1       ,  deltaMGenBestWMinus1         ;
   double deltaPtGenBestNeutrino1    , deltaPtGenBestAntiNeutrino1   ;
   double  deltaRGenBestNeutrino1    ,  deltaRGenBestAntiNeutrino1   ;
   double deltaPtGenBestTopQuark2    , deltaPtGenBestAntiTopQuark2   ;
   double  deltaRGenBestTopQuark2    ,  deltaRGenBestAntiTopQuark2   ;
   double  deltaMGenBestTopQuark2    ,  deltaMGenBestAntiTopQuark2   ;
   double deltaPtGenBestWPlus2       , deltaPtGenBestWMinus2         ;
   double  deltaRGenBestWPlus2       ,  deltaRGenBestWMinus2         ;
   double  deltaMGenBestWPlus2       ,  deltaMGenBestWMinus2         ;
   double deltaPtGenBestNeutrino2    , deltaPtGenBestAntiNeutrino2   ;
   double  deltaRGenBestNeutrino2    ,  deltaRGenBestAntiNeutrino2   ;
   double deltaPtGenBestTopQuark3    , deltaPtGenBestAntiTopQuark3   ;
   double  deltaRGenBestTopQuark3    ,  deltaRGenBestAntiTopQuark3   ;
   double  deltaMGenBestTopQuark3    ,  deltaMGenBestAntiTopQuark3   ;
   double deltaPtGenBestWPlus3       , deltaPtGenBestWMinus3         ;
   double  deltaRGenBestWPlus3       ,  deltaRGenBestWMinus3         ;
   double  deltaMGenBestWPlus3       ,  deltaMGenBestWMinus3         ;
   double deltaPtGenBestNeutrino3    , deltaPtGenBestAntiNeutrino3   ;
   double  deltaRGenBestNeutrino3    ,  deltaRGenBestAntiNeutrino3   ;
   double deltaPtGenBestTopQuark4    , deltaPtGenBestAntiTopQuark4   ;
   double  deltaRGenBestTopQuark4    ,  deltaRGenBestAntiTopQuark4   ;
   double  deltaMGenBestTopQuark4    ,  deltaMGenBestAntiTopQuark4   ;
   double deltaPtGenBestWPlus4       , deltaPtGenBestWMinus4         ;
   double  deltaRGenBestWPlus4       ,  deltaRGenBestWMinus4         ;
   double  deltaMGenBestWPlus4       ,  deltaMGenBestWMinus4         ;
   double deltaPtGenBestNeutrino4    , deltaPtGenBestAntiNeutrino4   ;
   double  deltaRGenBestNeutrino4    ,  deltaRGenBestAntiNeutrino4   ;
   double deltaPtGenBestTopQuark5    , deltaPtGenBestAntiTopQuark5   ;
   double  deltaRGenBestTopQuark5    ,  deltaRGenBestAntiTopQuark5   ;
   double  deltaMGenBestTopQuark5    ,  deltaMGenBestAntiTopQuark5   ;
   double deltaPtGenBestWPlus5       , deltaPtGenBestWMinus5         ;
   double  deltaRGenBestWPlus5       ,  deltaRGenBestWMinus5         ;
   double  deltaMGenBestWPlus5       ,  deltaMGenBestWMinus5         ;
   double deltaPtGenBestNeutrino5    , deltaPtGenBestAntiNeutrino5   ;
   double  deltaRGenBestNeutrino5    ,  deltaRGenBestAntiNeutrino5   ;
   double deltaPtGenBestTopQuark6    , deltaPtGenBestAntiTopQuark6   ;
   double  deltaRGenBestTopQuark6    ,  deltaRGenBestAntiTopQuark6   ;
   double  deltaMGenBestTopQuark6    ,  deltaMGenBestAntiTopQuark6   ;
   double deltaPtGenBestWPlus6       , deltaPtGenBestWMinus6         ;
   double  deltaRGenBestWPlus6       ,  deltaRGenBestWMinus6         ;
   double  deltaMGenBestWPlus6       ,  deltaMGenBestWMinus6         ;
   double deltaPtGenBestNeutrino6    , deltaPtGenBestAntiNeutrino6   ;
   double  deltaRGenBestNeutrino6    ,  deltaRGenBestAntiNeutrino6   ;
   double deltaPtGenBestTopQuark7    , deltaPtGenBestAntiTopQuark7   ;
   double  deltaRGenBestTopQuark7    ,  deltaRGenBestAntiTopQuark7   ;
   double  deltaMGenBestTopQuark7    ,  deltaMGenBestAntiTopQuark7   ;
   double deltaPtGenBestWPlus7       , deltaPtGenBestWMinus7         ;
   double  deltaRGenBestWPlus7       ,  deltaRGenBestWMinus7         ;
   double  deltaMGenBestWPlus7       ,  deltaMGenBestWMinus7         ;
   double deltaPtGenBestNeutrino7    , deltaPtGenBestAntiNeutrino7   ;
   double  deltaRGenBestNeutrino7    ,  deltaRGenBestAntiNeutrino7   ;
   double deltaPtGenBestTopQuark8    , deltaPtGenBestAntiTopQuark8   ;
   double  deltaRGenBestTopQuark8    ,  deltaRGenBestAntiTopQuark8   ;
   double  deltaMGenBestTopQuark8    ,  deltaMGenBestAntiTopQuark8   ;
   double deltaPtGenBestWPlus8       , deltaPtGenBestWMinus8         ;
   double  deltaRGenBestWPlus8       ,  deltaRGenBestWMinus8         ;
   double  deltaMGenBestWPlus8       ,  deltaMGenBestWMinus8         ;
   double deltaPtGenBestNeutrino8    , deltaPtGenBestAntiNeutrino8   ;
   double  deltaRGenBestNeutrino8    ,  deltaRGenBestAntiNeutrino8   ;
   double deltaPtGenBestBottomQuark  , deltaPtGenBestAntiBottomQuark ;
   double  deltaRGenBestBottomQuark  ,  deltaRGenBestAntiBottomQuark ;
   double deltaPtGenBestLightParton1 , deltaPtGenBestLightParton2    ;
   double deltaPxGenBestLightParton1 , deltaPxGenBestLightParton2    ;
   double deltaPyGenBestLightParton1 , deltaPyGenBestLightParton2    ;
   double  deltaRGenBestLightParton1 ,  deltaRGenBestLightParton2    ;

   //compare smeared and corrected
   double deltaPtSmearedBestTopQuark1    , deltaPtSmearedBestAntiTopQuark1   ;
   double  deltaRSmearedBestTopQuark1    ,  deltaRSmearedBestAntiTopQuark1   ;
   double  deltaMSmearedBestTopQuark1    ,  deltaMSmearedBestAntiTopQuark1   ;
   double deltaPtSmearedBestWPlus1       , deltaPtSmearedBestWMinus1         ;
   double  deltaRSmearedBestWPlus1       ,  deltaRSmearedBestWMinus1         ;
   double  deltaMSmearedBestWPlus1       ,  deltaMSmearedBestWMinus1         ;
   double deltaPtSmearedBestNeutrino1    , deltaPtSmearedBestAntiNeutrino1   ;
   double  deltaRSmearedBestNeutrino1    ,  deltaRSmearedBestAntiNeutrino1   ;
   double deltaPtSmearedBestTopQuark2    , deltaPtSmearedBestAntiTopQuark2   ;
   double  deltaRSmearedBestTopQuark2    ,  deltaRSmearedBestAntiTopQuark2   ;
   double  deltaMSmearedBestTopQuark2    ,  deltaMSmearedBestAntiTopQuark2   ;
   double deltaPtSmearedBestWPlus2       , deltaPtSmearedBestWMinus2         ;
   double  deltaRSmearedBestWPlus2       ,  deltaRSmearedBestWMinus2         ;
   double  deltaMSmearedBestWPlus2       ,  deltaMSmearedBestWMinus2         ;
   double deltaPtSmearedBestNeutrino2    , deltaPtSmearedBestAntiNeutrino2   ;
   double  deltaRSmearedBestNeutrino2    ,  deltaRSmearedBestAntiNeutrino2   ;
   double deltaPtSmearedBestTopQuark3    , deltaPtSmearedBestAntiTopQuark3   ;
   double  deltaRSmearedBestTopQuark3    ,  deltaRSmearedBestAntiTopQuark3   ;
   double  deltaMSmearedBestTopQuark3    ,  deltaMSmearedBestAntiTopQuark3   ;
   double deltaPtSmearedBestWPlus3       , deltaPtSmearedBestWMinus3         ;
   double  deltaRSmearedBestWPlus3       ,  deltaRSmearedBestWMinus3         ;
   double  deltaMSmearedBestWPlus3       ,  deltaMSmearedBestWMinus3         ;
   double deltaPtSmearedBestNeutrino3    , deltaPtSmearedBestAntiNeutrino3   ;
   double  deltaRSmearedBestNeutrino3    ,  deltaRSmearedBestAntiNeutrino3   ;
   double deltaPtSmearedBestTopQuark4    , deltaPtSmearedBestAntiTopQuark4   ;
   double  deltaRSmearedBestTopQuark4    ,  deltaRSmearedBestAntiTopQuark4   ;
   double  deltaMSmearedBestTopQuark4    ,  deltaMSmearedBestAntiTopQuark4   ;
   double deltaPtSmearedBestWPlus4       , deltaPtSmearedBestWMinus4         ;
   double  deltaRSmearedBestWPlus4       ,  deltaRSmearedBestWMinus4         ;
   double  deltaMSmearedBestWPlus4       ,  deltaMSmearedBestWMinus4         ;
   double deltaPtSmearedBestNeutrino4    , deltaPtSmearedBestAntiNeutrino4   ;
   double  deltaRSmearedBestNeutrino4    ,  deltaRSmearedBestAntiNeutrino4   ;
   double deltaPtSmearedBestTopQuark5    , deltaPtSmearedBestAntiTopQuark5   ;
   double  deltaRSmearedBestTopQuark5    ,  deltaRSmearedBestAntiTopQuark5   ;
   double  deltaMSmearedBestTopQuark5    ,  deltaMSmearedBestAntiTopQuark5   ;
   double deltaPtSmearedBestWPlus5       , deltaPtSmearedBestWMinus5         ;
   double  deltaRSmearedBestWPlus5       ,  deltaRSmearedBestWMinus5         ;
   double  deltaMSmearedBestWPlus5       ,  deltaMSmearedBestWMinus5         ;
   double deltaPtSmearedBestNeutrino5    , deltaPtSmearedBestAntiNeutrino5   ;
   double  deltaRSmearedBestNeutrino5    ,  deltaRSmearedBestAntiNeutrino5   ;
   double deltaPtSmearedBestTopQuark6    , deltaPtSmearedBestAntiTopQuark6   ;
   double  deltaRSmearedBestTopQuark6    ,  deltaRSmearedBestAntiTopQuark6   ;
   double  deltaMSmearedBestTopQuark6    ,  deltaMSmearedBestAntiTopQuark6   ;
   double deltaPtSmearedBestWPlus6       , deltaPtSmearedBestWMinus6         ;
   double  deltaRSmearedBestWPlus6       ,  deltaRSmearedBestWMinus6         ;
   double  deltaMSmearedBestWPlus6       ,  deltaMSmearedBestWMinus6         ;
   double deltaPtSmearedBestNeutrino6    , deltaPtSmearedBestAntiNeutrino6   ;
   double  deltaRSmearedBestNeutrino6    ,  deltaRSmearedBestAntiNeutrino6   ;
   double deltaPtSmearedBestTopQuark7    , deltaPtSmearedBestAntiTopQuark7   ;
   double  deltaRSmearedBestTopQuark7    ,  deltaRSmearedBestAntiTopQuark7   ;
   double  deltaMSmearedBestTopQuark7    ,  deltaMSmearedBestAntiTopQuark7   ;
   double deltaPtSmearedBestWPlus7       , deltaPtSmearedBestWMinus7         ;
   double  deltaRSmearedBestWPlus7       ,  deltaRSmearedBestWMinus7         ;
   double  deltaMSmearedBestWPlus7       ,  deltaMSmearedBestWMinus7         ;
   double deltaPtSmearedBestNeutrino7    , deltaPtSmearedBestAntiNeutrino7   ;
   double  deltaRSmearedBestNeutrino7    ,  deltaRSmearedBestAntiNeutrino7   ;
   double deltaPtSmearedBestTopQuark8    , deltaPtSmearedBestAntiTopQuark8   ;
   double  deltaRSmearedBestTopQuark8    ,  deltaRSmearedBestAntiTopQuark8   ;
   double  deltaMSmearedBestTopQuark8    ,  deltaMSmearedBestAntiTopQuark8   ;
   double deltaPtSmearedBestWPlus8       , deltaPtSmearedBestWMinus8         ;
   double  deltaRSmearedBestWPlus8       ,  deltaRSmearedBestWMinus8         ;
   double  deltaMSmearedBestWPlus8       ,  deltaMSmearedBestWMinus8         ;
   double deltaPtSmearedBestNeutrino8    , deltaPtSmearedBestAntiNeutrino8   ;
   double  deltaRSmearedBestNeutrino8    ,  deltaRSmearedBestAntiNeutrino8   ;
   double deltaPtSmearedBestBottomQuark  , deltaPtSmearedBestAntiBottomQuark ;
   double  deltaRSmearedBestBottomQuark  ,  deltaRSmearedBestAntiBottomQuark ;
   double deltaPtSmearedBestLightParton1 , deltaPtSmearedBestLightParton2    ;
   double deltaPxSmearedBestLightParton1 , deltaPxSmearedBestLightParton2    ;
   double deltaPySmearedBestLightParton1 , deltaPySmearedBestLightParton2    ;
   double  deltaRSmearedBestLightParton1 ,  deltaRSmearedBestLightParton2    ;

   //true momenta
   double        genTopQuarkPt ,        genTopQuarkEta ,        genTopQuarkPhi ,        genTopQuarkE ,     genTopQuarkM ;
   double    genAntiTopQuarkPt ,    genAntiTopQuarkEta ,    genAntiTopQuarkPhi ,    genAntiTopQuarkE , genAntiTopQuarkM ;
   double           genWPlusPt ,           genWPlusEta ,           genWPlusPhi ,           genWPlusE ,        genWPlusM ;
   double          genWMinusPt ,          genWMinusEta ,          genWMinusPhi ,          genWMinusE ,       genWMinusM ;
   double     genBottomQuarkPt ,     genBottomQuarkEta ,     genBottomQuarkPhi ,     genBottomQuarkE ;
   double genAntiBottomQuarkPt , genAntiBottomQuarkEta , genAntiBottomQuarkPhi , genAntiBottomQuarkE ;
   double        genNeutrinoPt ,        genNeutrinoEta ,        genNeutrinoPhi ,        genNeutrinoE ;
   double    genAntiNeutrinoPt ,    genAntiNeutrinoEta ,    genAntiNeutrinoPhi ,    genAntiNeutrinoE ;
   double    genLightParton1Pt ,    genLightParton1Eta ,    genLightParton1Phi ,    genLightParton1E ;
   double    genLightParton2Pt ,    genLightParton2Eta ,    genLightParton2Phi ,    genLightParton2E ;
   double          genLeptonPt ,          genLeptonEta ,          genLeptonPhi ,          genLeptonE ;
   double      genAntiLeptonPt ,      genAntiLeptonEta ,      genAntiLeptonPhi ,      genAntiLeptonE ;

   //smeared momenta
   double        smearedTopQuarkPt ,        smearedTopQuarkEta ,        smearedTopQuarkPhi ,        smearedTopQuarkE ,     smearedTopQuarkM ;
   double    smearedAntiTopQuarkPt ,    smearedAntiTopQuarkEta ,    smearedAntiTopQuarkPhi ,    smearedAntiTopQuarkE , smearedAntiTopQuarkM ;
   double           smearedWPlusPt ,           smearedWPlusEta ,           smearedWPlusPhi ,           smearedWPlusE ,        smearedWPlusM ;
   double          smearedWMinusPt ,          smearedWMinusEta ,          smearedWMinusPhi ,          smearedWMinusE ,       smearedWMinusM ;
   double     smearedBottomQuarkPt ,     smearedBottomQuarkEta ,     smearedBottomQuarkPhi ,     smearedBottomQuarkE ;
   double smearedAntiBottomQuarkPt , smearedAntiBottomQuarkEta , smearedAntiBottomQuarkPhi , smearedAntiBottomQuarkE ;
   double    smearedLightParton1Pt ,    smearedLightParton1Eta ,    smearedLightParton1Phi ,    smearedLightParton1E ;
   double    smearedLightParton2Pt ,    smearedLightParton2Eta ,    smearedLightParton2Phi ,    smearedLightParton2E ;

   //best momenta
   double        bestTopQuark1Pt ,       bestTopQuark1Eta ,       bestTopQuark1Phi ,       bestTopQuark1E ,     bestTopQuark1M ;
   double    bestAntiTopQuark1Pt ,   bestAntiTopQuark1Eta ,   bestAntiTopQuark1Phi ,   bestAntiTopQuark1E , bestAntiTopQuark1M ;
   double           bestWPlus1Pt ,          bestWPlus1Eta ,          bestWPlus1Phi ,          bestWPlus1E ,        bestWPlus1M ;
   double          bestWMinus1Pt ,         bestWMinus1Eta ,         bestWMinus1Phi ,         bestWMinus1E ,       bestWMinus1M ;
   double        bestNeutrino1Pt ,       bestNeutrino1Eta ,       bestNeutrino1Phi ,       bestNeutrino1E ,     bestNeutrino1PairingTruth ;
   double    bestAntiNeutrino1Pt ,   bestAntiNeutrino1Eta ,   bestAntiNeutrino1Phi ,   bestAntiNeutrino1E , bestAntiNeutrino1PairingTruth ;
   double        bestTopQuark2Pt ,       bestTopQuark2Eta ,       bestTopQuark2Phi ,       bestTopQuark2E ,     bestTopQuark2M ;
   double    bestAntiTopQuark2Pt ,   bestAntiTopQuark2Eta ,   bestAntiTopQuark2Phi ,   bestAntiTopQuark2E , bestAntiTopQuark2M ;
   double           bestWPlus2Pt ,          bestWPlus2Eta ,          bestWPlus2Phi ,          bestWPlus2E ,        bestWPlus2M ;
   double          bestWMinus2Pt ,         bestWMinus2Eta ,         bestWMinus2Phi ,         bestWMinus2E ,       bestWMinus2M ;
   double        bestNeutrino2Pt ,       bestNeutrino2Eta ,       bestNeutrino2Phi ,       bestNeutrino2E ,     bestNeutrino2PairingTruth ;
   double    bestAntiNeutrino2Pt ,   bestAntiNeutrino2Eta ,   bestAntiNeutrino2Phi ,   bestAntiNeutrino2E , bestAntiNeutrino2PairingTruth ;
   double        bestTopQuark3Pt ,       bestTopQuark3Eta ,       bestTopQuark3Phi ,       bestTopQuark3E ,     bestTopQuark3M ;
   double    bestAntiTopQuark3Pt ,   bestAntiTopQuark3Eta ,   bestAntiTopQuark3Phi ,   bestAntiTopQuark3E , bestAntiTopQuark3M ;
   double           bestWPlus3Pt ,          bestWPlus3Eta ,          bestWPlus3Phi ,          bestWPlus3E ,        bestWPlus3M ;
   double          bestWMinus3Pt ,         bestWMinus3Eta ,         bestWMinus3Phi ,         bestWMinus3E ,       bestWMinus3M ;
   double        bestNeutrino3Pt ,       bestNeutrino3Eta ,       bestNeutrino3Phi ,       bestNeutrino3E ,     bestNeutrino3PairingTruth ;
   double    bestAntiNeutrino3Pt ,   bestAntiNeutrino3Eta ,   bestAntiNeutrino3Phi ,   bestAntiNeutrino3E , bestAntiNeutrino3PairingTruth ;
   double        bestTopQuark4Pt ,       bestTopQuark4Eta ,       bestTopQuark4Phi ,       bestTopQuark4E ,     bestTopQuark4M ;
   double    bestAntiTopQuark4Pt ,   bestAntiTopQuark4Eta ,   bestAntiTopQuark4Phi ,   bestAntiTopQuark4E , bestAntiTopQuark4M ;
   double           bestWPlus4Pt ,          bestWPlus4Eta ,          bestWPlus4Phi ,          bestWPlus4E ,        bestWPlus4M ;
   double          bestWMinus4Pt ,         bestWMinus4Eta ,         bestWMinus4Phi ,         bestWMinus4E ,       bestWMinus4M ;
   double        bestNeutrino4Pt ,       bestNeutrino4Eta ,       bestNeutrino4Phi ,       bestNeutrino4E ,     bestNeutrino4PairingTruth ;
   double    bestAntiNeutrino4Pt ,   bestAntiNeutrino4Eta ,   bestAntiNeutrino4Phi ,   bestAntiNeutrino4E , bestAntiNeutrino4PairingTruth ;
   double        bestTopQuark5Pt ,       bestTopQuark5Eta ,       bestTopQuark5Phi ,       bestTopQuark5E ,     bestTopQuark5M ;
   double    bestAntiTopQuark5Pt ,   bestAntiTopQuark5Eta ,   bestAntiTopQuark5Phi ,   bestAntiTopQuark5E , bestAntiTopQuark5M ;
   double           bestWPlus5Pt ,          bestWPlus5Eta ,          bestWPlus5Phi ,          bestWPlus5E ,        bestWPlus5M ;
   double          bestWMinus5Pt ,         bestWMinus5Eta ,         bestWMinus5Phi ,         bestWMinus5E ,       bestWMinus5M ;
   double        bestNeutrino5Pt ,       bestNeutrino5Eta ,       bestNeutrino5Phi ,       bestNeutrino5E ,     bestNeutrino5PairingTruth ;
   double    bestAntiNeutrino5Pt ,   bestAntiNeutrino5Eta ,   bestAntiNeutrino5Phi ,   bestAntiNeutrino5E , bestAntiNeutrino5PairingTruth ;
   double        bestTopQuark6Pt ,       bestTopQuark6Eta ,       bestTopQuark6Phi ,       bestTopQuark6E ,     bestTopQuark6M ;
   double    bestAntiTopQuark6Pt ,   bestAntiTopQuark6Eta ,   bestAntiTopQuark6Phi ,   bestAntiTopQuark6E , bestAntiTopQuark6M ;
   double           bestWPlus6Pt ,          bestWPlus6Eta ,          bestWPlus6Phi ,          bestWPlus6E ,        bestWPlus6M ;
   double          bestWMinus6Pt ,         bestWMinus6Eta ,         bestWMinus6Phi ,         bestWMinus6E ,       bestWMinus6M ;
   double        bestNeutrino6Pt ,       bestNeutrino6Eta ,       bestNeutrino6Phi ,       bestNeutrino6E ,     bestNeutrino6PairingTruth ;
   double    bestAntiNeutrino6Pt ,   bestAntiNeutrino6Eta ,   bestAntiNeutrino6Phi ,   bestAntiNeutrino6E , bestAntiNeutrino6PairingTruth ;
   double        bestTopQuark7Pt ,       bestTopQuark7Eta ,       bestTopQuark7Phi ,       bestTopQuark7E ,     bestTopQuark7M ;
   double    bestAntiTopQuark7Pt ,   bestAntiTopQuark7Eta ,   bestAntiTopQuark7Phi ,   bestAntiTopQuark7E , bestAntiTopQuark7M ;
   double           bestWPlus7Pt ,          bestWPlus7Eta ,          bestWPlus7Phi ,          bestWPlus7E ,        bestWPlus7M ;
   double          bestWMinus7Pt ,         bestWMinus7Eta ,         bestWMinus7Phi ,         bestWMinus7E ,       bestWMinus7M ;
   double        bestNeutrino7Pt ,       bestNeutrino7Eta ,       bestNeutrino7Phi ,       bestNeutrino7E ,     bestNeutrino7PairingTruth ;
   double    bestAntiNeutrino7Pt ,   bestAntiNeutrino7Eta ,   bestAntiNeutrino7Phi ,   bestAntiNeutrino7E , bestAntiNeutrino7PairingTruth ;
   double        bestTopQuark8Pt ,       bestTopQuark8Eta ,       bestTopQuark8Phi ,       bestTopQuark8E ,     bestTopQuark8M ;
   double    bestAntiTopQuark8Pt ,   bestAntiTopQuark8Eta ,   bestAntiTopQuark8Phi ,   bestAntiTopQuark8E , bestAntiTopQuark8M ;
   double           bestWPlus8Pt ,          bestWPlus8Eta ,          bestWPlus8Phi ,          bestWPlus8E ,        bestWPlus8M ;
   double          bestWMinus8Pt ,         bestWMinus8Eta ,         bestWMinus8Phi ,         bestWMinus8E ,       bestWMinus8M ;
   double        bestNeutrino8Pt ,       bestNeutrino8Eta ,       bestNeutrino8Phi ,       bestNeutrino8E ,     bestNeutrino8PairingTruth ;
   double    bestAntiNeutrino8Pt ,   bestAntiNeutrino8Eta ,   bestAntiNeutrino8Phi ,   bestAntiNeutrino8E , bestAntiNeutrino8PairingTruth ;
   double     bestBottomQuarkPt  ,     bestBottomQuarkEta ,     bestBottomQuarkPhi ,     bestBottomQuarkE ;
   double bestAntiBottomQuarkPt  , bestAntiBottomQuarkEta , bestAntiBottomQuarkPhi , bestAntiBottomQuarkE ;
   double    bestLightParton1Pt  ,    bestLightParton1Eta ,    bestLightParton1Phi ,    bestLightParton1E ;
   double    bestLightParton2Pt  ,    bestLightParton2Eta ,    bestLightParton2Phi ,    bestLightParton2E ;

   //deltas from the minimization
   double  deltaPtBottomQuarkFromMin ,  deltaPtAntiBottomQuarkFromMin ;
   double deltaPhiBottomQuarkFromMin , deltaPhiAntiBottomQuarkFromMin ;
   double deltaEtaBottomQuarkFromMin , deltaEtaAntiBottomQuarkFromMin ;
   double deltaPxLightParton1FromMin ,     deltaPxLightParton2FromMin ;
   double deltaPyLightParton1FromMin ,     deltaPyLightParton2FromMin ;
   double      deltaMTopQuarkFromMin ,      deltaMAntiTopQuarkFromMin ;
   double         deltaMWPlusFromMin ,            deltaMWMinusFromMin ;

   //other
   double deltaRGenBottomQuarks , deltaRSmearedBottomQuarks , deltaRBestBottomQuarks ;
   double deltaRGenLightPartons , deltaRSmearedLightPartons , deltaRBestLightPartons ;
   double    deltaRGenTopQuarks ,    deltaRSmearedTopQuarks ;
   double      deltaRGenWBosons ,      deltaRSmearedWBosons ;
   double    deltaRGenNeutrinos ;
   double deltaRBestTopQuarks1, deltaRBestTopQuarks2, deltaRBestTopQuarks3, deltaRBestTopQuarks4, deltaRBestTopQuarks5, deltaRBestTopQuarks6, deltaRBestTopQuarks7, deltaRBestTopQuarks8;
   double deltaRBestWBosons1, deltaRBestWBosons2, deltaRBestWBosons3, deltaRBestWBosons4, deltaRBestWBosons5, deltaRBestWBosons6, deltaRBestWBosons7, deltaRBestWBosons8;
   double deltaRBestNeutrinos1, deltaRBestNeutrinos2, deltaRBestNeutrinos3, deltaRBestNeutrinos4, deltaRBestNeutrinos5, deltaRBestNeutrinos6, deltaRBestNeutrinos7, deltaRBestNeutrinos8;

};

#endif

#ifdef ttbarReconstructionFromLHE_cxx
ttbarReconstructionFromLHE::ttbarReconstructionFromLHE(TString dirName, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  dirName_ = dirName;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lhe.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("lhe.root");
      }
      f->GetObject("Physics",tree);

   }
   Init(tree);
}

ttbarReconstructionFromLHE::~ttbarReconstructionFromLHE()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

double ttbarReconstructionFromLHE::deltaR(XYZTLorentzVector& p2, XYZTLorentzVector& p1 )
{
  double deltaEta = p1.Eta()-p2.Eta();
  double deltaPhi = p1.Phi()-p2.Phi();
  while(deltaPhi> 3.14159265359)
    {
      deltaPhi -= 2*3.14159265359;
    }
  while(deltaPhi< -3.14159265359)
    {
      deltaPhi += 2*3.14159265359;
    }
  return sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
}

void ttbarReconstructionFromLHE::initOutput(int whichLoop)
{
  TString fileName = dirName_+"output_";
  fileName += whichLoop;
  fileName += ".root";
  outFile = new TFile(fileName,"RECREATE","tree");
  outTree = new TTree("tree","tree");

  outTree->Branch( "chi2"          ,  &chi2          ) ;
  outTree->Branch( "bestPairing"   ,  &bestPairing   ) ;
  outTree->Branch( "minStatusCode" ,  &minStatusCode ) ;
  outTree->Branch( "eventNumber"   ,  &eventNumber   ) ;
  outTree->Branch( "smearingIndex" ,  &smearingIndex ) ;

  outTree->Branch( "genTopQuarkPt"         ,  &genTopQuarkPt         ) ;
  outTree->Branch( "genTopQuarkEta"        ,  &genTopQuarkEta        ) ;    
  outTree->Branch( "genTopQuarkPhi"        ,  &genTopQuarkPhi        ) ;    
  outTree->Branch( "genTopQuarkE"          ,  &genTopQuarkE          ) ;      
  outTree->Branch( "genTopQuarkM"          ,  &genTopQuarkM          ) ;      
  outTree->Branch( "genAntiTopQuarkPt"     ,  &genAntiTopQuarkPt     ) ; 
  outTree->Branch( "genAntiTopQuarkEta"    ,  &genAntiTopQuarkEta    ) ; 
  outTree->Branch( "genAntiTopQuarkPhi"    ,  &genAntiTopQuarkPhi    ) ;
  outTree->Branch( "genAntiTopQuarkE"      ,  &genAntiTopQuarkE      ) ;  
  outTree->Branch( "genAntiTopQuarkM"      ,  &genAntiTopQuarkM      ) ;  
  outTree->Branch( "genWPlusPt"            ,  &genWPlusPt            ) ;
  outTree->Branch( "genWPlusEta"           ,  &genWPlusEta           ) ;
  outTree->Branch( "genWPlusPhi"           ,  &genWPlusPhi           ) ;
  outTree->Branch( "genWPlusE"             ,  &genWPlusE             ) ;
  outTree->Branch( "genWPlusM"             ,  &genWPlusM             ) ;
  outTree->Branch( "genWMinusPt"           ,  &genWMinusPt           ) ;
  outTree->Branch( "genWMinusEta"          ,  &genWMinusEta          ) ;
  outTree->Branch( "genWMinusPhi"          ,  &genWMinusPhi          ) ;
  outTree->Branch( "genWMinusE"            ,  &genWMinusE            ) ;
  outTree->Branch( "genWMinusM"            ,  &genWMinusM            ) ;
  outTree->Branch( "genBottomQuarkPt"      ,  &genBottomQuarkPt      ) ;
  outTree->Branch( "genBottomQuarkEta"     ,  &genBottomQuarkEta     ) ;
  outTree->Branch( "genBottomQuarkPhi"     ,  &genBottomQuarkPhi     ) ;
  outTree->Branch( "genBottomQuarkE"       ,  &genBottomQuarkE       ) ;
  outTree->Branch( "genAntiBottomQuarkPt"  ,  &genAntiBottomQuarkPt  ) ;
  outTree->Branch( "genAntiBottomQuarkEta" ,  &genAntiBottomQuarkEta ) ;
  outTree->Branch( "genAntiBottomQuarkPhi" ,  &genAntiBottomQuarkPhi ) ;
  outTree->Branch( "genAntiBottomQuarkE"   ,  &genAntiBottomQuarkE   ) ;
  outTree->Branch( "genLightParton1Pt"     ,  &genLightParton1Pt     ) ;
  outTree->Branch( "genLightParton1Eta"    ,  &genLightParton1Eta    ) ;
  outTree->Branch( "genLightParton1Phi"    ,  &genLightParton1Phi    ) ;
  outTree->Branch( "genLightParton1E"      ,  &genLightParton1E      ) ;
  outTree->Branch( "genLightParton2Pt"     ,  &genLightParton2Pt     ) ;
  outTree->Branch( "genLightParton2Eta"    ,  &genLightParton2Eta    ) ;
  outTree->Branch( "genLightParton2Phi"    ,  &genLightParton2Phi    ) ;
  outTree->Branch( "genLightParton2E"      ,  &genLightParton2E      ) ;
  outTree->Branch( "genNeutrinoPt"         ,  &genNeutrinoPt         ) ;
  outTree->Branch( "genNeutrinoEta"        ,  &genNeutrinoEta        ) ;
  outTree->Branch( "genNeutrinoPhi"        ,  &genNeutrinoPhi        ) ;
  outTree->Branch( "genNeutrinoE"          ,  &genNeutrinoE          ) ;
  outTree->Branch( "genAntiNeutrinoPt"     ,  &genAntiNeutrinoPt     ) ;
  outTree->Branch( "genAntiNeutrinoEta"    ,  &genAntiNeutrinoEta    ) ;
  outTree->Branch( "genAntiNeutrinoPhi"    ,  &genAntiNeutrinoPhi    ) ;
  outTree->Branch( "genAntiNeutrinoE"      ,  &genAntiNeutrinoE      ) ;
  outTree->Branch( "genLeptonPt"           ,  &genLeptonPt           ) ;
  outTree->Branch( "genLeptonEta"          ,  &genLeptonEta          ) ;
  outTree->Branch( "genLeptonPhi"          ,  &genLeptonPhi          ) ;
  outTree->Branch( "genLeptonE"            ,  &genLeptonE            ) ;
  outTree->Branch( "genAntiLeptonPt"       ,  &genAntiLeptonPt       ) ;
  outTree->Branch( "genAntiLeptonEta"      ,  &genAntiLeptonEta      ) ;
  outTree->Branch( "genAntiLeptonPhi"      ,  &genAntiLeptonPhi      ) ;
  outTree->Branch( "genAntiLeptonE"        ,  &genAntiLeptonE        ) ;

  outTree->Branch( "smearedBottomQuarkPt"         ,  &smearedBottomQuarkPt      ) ;
  outTree->Branch( "smearedBottomQuarkEta"        ,  &smearedBottomQuarkEta     ) ;
  outTree->Branch( "smearedBottomQuarkPhi"        ,  &smearedBottomQuarkPhi     ) ;
  outTree->Branch( "smearedBottomQuarkE"          ,  &smearedBottomQuarkE       ) ;
  outTree->Branch( "smearedAntiBottomQuarkPt"     ,  &smearedAntiBottomQuarkPt  ) ;
  outTree->Branch( "smearedAntiBottomQuarkEta"    ,  &smearedAntiBottomQuarkEta ) ;
  outTree->Branch( "smearedAntiBottomQuarkPhi"    ,  &smearedAntiBottomQuarkPhi ) ;
  outTree->Branch( "smearedAntiBottomQuarkE"      ,  &smearedAntiBottomQuarkE   ) ;
  outTree->Branch( "smearedLightParton1Pt"        ,  &smearedLightParton1Pt     ) ;
  outTree->Branch( "smearedLightParton1Eta"       ,  &smearedLightParton1Eta    ) ;
  outTree->Branch( "smearedLightParton1Phi"       ,  &smearedLightParton1Phi    ) ;
  outTree->Branch( "smearedLightParton1E"         ,  &smearedLightParton1E      ) ;
  outTree->Branch( "smearedLightParton2Pt"        ,  &smearedLightParton2Pt     ) ;
  outTree->Branch( "smearedLightParton2Eta"       ,  &smearedLightParton2Eta    ) ;
  outTree->Branch( "smearedLightParton2Phi"       ,  &smearedLightParton2Phi    ) ;
  outTree->Branch( "smearedLightParton2E"         ,  &smearedLightParton2E      ) ;
  outTree->Branch( "smearedTopQuarkPt"            ,  &smearedTopQuarkPt         ) ;
  outTree->Branch( "smearedTopQuarkEta"           ,  &smearedTopQuarkEta        ) ;
  outTree->Branch( "smearedTopQuarkPhi"           ,  &smearedTopQuarkPhi        ) ;
  outTree->Branch( "smearedTopQuarkE"             ,  &smearedTopQuarkE          ) ;
  outTree->Branch( "smearedTopQuarkM"             ,  &smearedTopQuarkM          ) ;
  outTree->Branch( "smearedAntiTopQuarkPt"        ,  &smearedAntiTopQuarkPt     ) ;
  outTree->Branch( "smearedAntiTopQuarkEta"       ,  &smearedAntiTopQuarkEta    ) ;
  outTree->Branch( "smearedAntiTopQuarkPhi"       ,  &smearedAntiTopQuarkPhi    ) ;
  outTree->Branch( "smearedAntiTopQuarkE"         ,  &smearedAntiTopQuarkE      ) ;
  outTree->Branch( "smearedAntiTopQuarkM"         ,  &smearedAntiTopQuarkM      ) ;
  outTree->Branch( "smearedWPlusPt"               ,  &smearedWPlusPt            ) ;
  outTree->Branch( "smearedWPlusEta"              ,  &smearedWPlusEta           ) ;
  outTree->Branch( "smearedWPlusPhi"              ,  &smearedWPlusPhi           ) ;
  outTree->Branch( "smearedWPlusE"                ,  &smearedWPlusE             ) ;
  outTree->Branch( "smearedWPlusM"                ,  &smearedWPlusM             ) ;
  outTree->Branch( "smearedWMinusPt"              ,  &smearedWMinusPt           ) ;
  outTree->Branch( "smearedWMinusEta"             ,  &smearedWMinusEta          ) ;
  outTree->Branch( "smearedWMinusPhi"             ,  &smearedWMinusPhi          ) ;
  outTree->Branch( "smearedWMinusE"               ,  &smearedWMinusE            ) ;
  outTree->Branch( "smearedWMinusM"               ,  &smearedWMinusM            ) ;

  outTree->Branch( "bestTopQuark1Pt"               ,  &bestTopQuark1Pt               ) ;
  outTree->Branch( "bestTopQuark1Eta"              ,  &bestTopQuark1Eta              ) ;
  outTree->Branch( "bestTopQuark1Phi"              ,  &bestTopQuark1Phi              ) ;
  outTree->Branch( "bestTopQuark1E"                ,  &bestTopQuark1E                ) ;
  outTree->Branch( "bestTopQuark1M"                ,  &bestTopQuark1M                ) ;
  outTree->Branch( "bestAntiTopQuark1Pt"           ,  &bestAntiTopQuark1Pt           ) ;
  outTree->Branch( "bestAntiTopQuark1Eta"          ,  &bestAntiTopQuark1Eta          ) ;
  outTree->Branch( "bestAntiTopQuark1Phi"          ,  &bestAntiTopQuark1Phi          ) ;
  outTree->Branch( "bestAntiTopQuark1E"            ,  &bestAntiTopQuark1E            ) ;
  outTree->Branch( "bestAntiTopQuark1M"            ,  &bestAntiTopQuark1M            ) ;
  outTree->Branch( "bestWPlus1Pt"                  ,  &bestWPlus1Pt                  ) ;
  outTree->Branch( "bestWPlus1Eta"                 ,  &bestWPlus1Eta                 ) ;
  outTree->Branch( "bestWPlus1Phi"                 ,  &bestWPlus1Phi                 ) ;
  outTree->Branch( "bestWPlus1E"                   ,  &bestWPlus1E                   ) ;
  outTree->Branch( "bestWPlus1M"                   ,  &bestWPlus1M                   ) ;
  outTree->Branch( "bestWMinus1Pt"                 ,  &bestWMinus1Pt                 ) ;
  outTree->Branch( "bestWMinus1Eta"                ,  &bestWMinus1Eta                ) ;
  outTree->Branch( "bestWMinus1Phi"                ,  &bestWMinus1Phi                ) ;
  outTree->Branch( "bestWMinus1E"                  ,  &bestWMinus1E                  ) ;
  outTree->Branch( "bestWMinus1M"                  ,  &bestWMinus1M                  ) ;
  outTree->Branch( "bestNeutrino1Pt"               ,  &bestNeutrino1Pt               ) ;
  outTree->Branch( "bestNeutrino1Eta"              ,  &bestNeutrino1Eta              ) ;
  outTree->Branch( "bestNeutrino1Phi"              ,  &bestNeutrino1Phi              ) ;
  outTree->Branch( "bestNeutrino1E"                ,  &bestNeutrino1E                ) ;
  outTree->Branch( "bestNeutrino1PairingTruth"     ,  &bestNeutrino1PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino1Pt"           ,  &bestAntiNeutrino1Pt           ) ;
  outTree->Branch( "bestAntiNeutrino1Eta"          ,  &bestAntiNeutrino1Eta          ) ;
  outTree->Branch( "bestAntiNeutrino1Phi"          ,  &bestAntiNeutrino1Phi          ) ;
  outTree->Branch( "bestAntiNeutrino1E"            ,  &bestAntiNeutrino1E            ) ;
  outTree->Branch( "bestAntiNeutrino1PairingTruth" ,  &bestAntiNeutrino1PairingTruth ) ;
  outTree->Branch( "bestTopQuark2Pt"               ,  &bestTopQuark2Pt               ) ;
  outTree->Branch( "bestTopQuark2Eta"              ,  &bestTopQuark2Eta              ) ;
  outTree->Branch( "bestTopQuark2Phi"              ,  &bestTopQuark2Phi              ) ;
  outTree->Branch( "bestTopQuark2E"                ,  &bestTopQuark2E                ) ;
  outTree->Branch( "bestTopQuark2M"                ,  &bestTopQuark2M                ) ;
  outTree->Branch( "bestAntiTopQuark2Pt"           ,  &bestAntiTopQuark2Pt           ) ;
  outTree->Branch( "bestAntiTopQuark2Eta"          ,  &bestAntiTopQuark2Eta          ) ;
  outTree->Branch( "bestAntiTopQuark2Phi"          ,  &bestAntiTopQuark2Phi          ) ;
  outTree->Branch( "bestAntiTopQuark2E"            ,  &bestAntiTopQuark2E            ) ;
  outTree->Branch( "bestAntiTopQuark2M"            ,  &bestAntiTopQuark2M            ) ;
  outTree->Branch( "bestWPlus2Pt"                  ,  &bestWPlus2Pt                  ) ;
  outTree->Branch( "bestWPlus2Eta"                 ,  &bestWPlus2Eta                 ) ;
  outTree->Branch( "bestWPlus2Phi"                 ,  &bestWPlus2Phi                 ) ;
  outTree->Branch( "bestWPlus2E"                   ,  &bestWPlus2E                   ) ;
  outTree->Branch( "bestWPlus2M"                   ,  &bestWPlus2M                   ) ;
  outTree->Branch( "bestWMinus2Pt"                 ,  &bestWMinus2Pt                 ) ;
  outTree->Branch( "bestWMinus2Eta"                ,  &bestWMinus2Eta                ) ;
  outTree->Branch( "bestWMinus2Phi"                ,  &bestWMinus2Phi                ) ;
  outTree->Branch( "bestWMinus2E"                  ,  &bestWMinus2E                  ) ;
  outTree->Branch( "bestWMinus2M"                  ,  &bestWMinus2M                  ) ;
  outTree->Branch( "bestNeutrino2Pt"               ,  &bestNeutrino2Pt               ) ;
  outTree->Branch( "bestNeutrino2Eta"              ,  &bestNeutrino2Eta              ) ;
  outTree->Branch( "bestNeutrino2Phi"              ,  &bestNeutrino2Phi              ) ;
  outTree->Branch( "bestNeutrino2E"                ,  &bestNeutrino2E                ) ;
  outTree->Branch( "bestNeutrino2PairingTruth"     ,  &bestNeutrino2PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino2Pt"           ,  &bestAntiNeutrino2Pt           ) ;
  outTree->Branch( "bestAntiNeutrino2Eta"          ,  &bestAntiNeutrino2Eta          ) ;
  outTree->Branch( "bestAntiNeutrino2Phi"          ,  &bestAntiNeutrino2Phi          ) ;
  outTree->Branch( "bestAntiNeutrino2E"            ,  &bestAntiNeutrino2E            ) ;
  outTree->Branch( "bestAntiNeutrino2PairingTruth" ,  &bestAntiNeutrino2PairingTruth ) ;
  outTree->Branch( "bestTopQuark3Pt"               ,  &bestTopQuark3Pt               ) ;
  outTree->Branch( "bestTopQuark3Eta"              ,  &bestTopQuark3Eta              ) ;
  outTree->Branch( "bestTopQuark3Phi"              ,  &bestTopQuark3Phi              ) ;
  outTree->Branch( "bestTopQuark3E"                ,  &bestTopQuark3E                ) ;
  outTree->Branch( "bestTopQuark3M"                ,  &bestTopQuark3M                ) ;
  outTree->Branch( "bestAntiTopQuark3Pt"           ,  &bestAntiTopQuark3Pt           ) ;
  outTree->Branch( "bestAntiTopQuark3Eta"          ,  &bestAntiTopQuark3Eta          ) ;
  outTree->Branch( "bestAntiTopQuark3Phi"          ,  &bestAntiTopQuark3Phi          ) ;
  outTree->Branch( "bestAntiTopQuark3E"            ,  &bestAntiTopQuark3E            ) ;
  outTree->Branch( "bestAntiTopQuark3M"            ,  &bestAntiTopQuark3M            ) ;
  outTree->Branch( "bestWPlus3Pt"                  ,  &bestWPlus3Pt                  ) ;
  outTree->Branch( "bestWPlus3Eta"                 ,  &bestWPlus3Eta                 ) ;
  outTree->Branch( "bestWPlus3Phi"                 ,  &bestWPlus3Phi                 ) ;
  outTree->Branch( "bestWPlus3E"                   ,  &bestWPlus3E                   ) ;
  outTree->Branch( "bestWPlus3M"                   ,  &bestWPlus3M                   ) ;
  outTree->Branch( "bestWMinus3Pt"                 ,  &bestWMinus3Pt                 ) ;
  outTree->Branch( "bestWMinus3Eta"                ,  &bestWMinus3Eta                ) ;
  outTree->Branch( "bestWMinus3Phi"                ,  &bestWMinus3Phi                ) ;
  outTree->Branch( "bestWMinus3E"                  ,  &bestWMinus3E                  ) ;
  outTree->Branch( "bestWMinus3M"                  ,  &bestWMinus3M                  ) ;
  outTree->Branch( "bestNeutrino3Pt"               ,  &bestNeutrino3Pt               ) ;
  outTree->Branch( "bestNeutrino3Eta"              ,  &bestNeutrino3Eta              ) ;
  outTree->Branch( "bestNeutrino3Phi"              ,  &bestNeutrino3Phi              ) ;
  outTree->Branch( "bestNeutrino3E"                ,  &bestNeutrino3E                ) ;
  outTree->Branch( "bestNeutrino3PairingTruth"     ,  &bestNeutrino3PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino3Pt"           ,  &bestAntiNeutrino3Pt           ) ;
  outTree->Branch( "bestAntiNeutrino3Eta"          ,  &bestAntiNeutrino3Eta          ) ;
  outTree->Branch( "bestAntiNeutrino3Phi"          ,  &bestAntiNeutrino3Phi          ) ;
  outTree->Branch( "bestAntiNeutrino3E"            ,  &bestAntiNeutrino3E            ) ;
  outTree->Branch( "bestAntiNeutrino3PairingTruth" ,  &bestAntiNeutrino3PairingTruth ) ;
  outTree->Branch( "bestTopQuark4Pt"               ,  &bestTopQuark4Pt               ) ;
  outTree->Branch( "bestTopQuark4Eta"              ,  &bestTopQuark4Eta              ) ;
  outTree->Branch( "bestTopQuark4Phi"              ,  &bestTopQuark4Phi              ) ;
  outTree->Branch( "bestTopQuark4E"                ,  &bestTopQuark4E                ) ;
  outTree->Branch( "bestTopQuark4M"                ,  &bestTopQuark4M                ) ;
  outTree->Branch( "bestAntiTopQuark4Pt"           ,  &bestAntiTopQuark4Pt           ) ;
  outTree->Branch( "bestAntiTopQuark4Eta"          ,  &bestAntiTopQuark4Eta          ) ;
  outTree->Branch( "bestAntiTopQuark4Phi"          ,  &bestAntiTopQuark4Phi          ) ;
  outTree->Branch( "bestAntiTopQuark4E"            ,  &bestAntiTopQuark4E            ) ;
  outTree->Branch( "bestAntiTopQuark4M"            ,  &bestAntiTopQuark4M            ) ;
  outTree->Branch( "bestWPlus4Pt"                  ,  &bestWPlus4Pt                  ) ;
  outTree->Branch( "bestWPlus4Eta"                 ,  &bestWPlus4Eta                 ) ;
  outTree->Branch( "bestWPlus4Phi"                 ,  &bestWPlus4Phi                 ) ;
  outTree->Branch( "bestWPlus4E"                   ,  &bestWPlus4E                   ) ;
  outTree->Branch( "bestWPlus4M"                   ,  &bestWPlus4M                   ) ;
  outTree->Branch( "bestWMinus4Pt"                 ,  &bestWMinus4Pt                 ) ;
  outTree->Branch( "bestWMinus4Eta"                ,  &bestWMinus4Eta                ) ;
  outTree->Branch( "bestWMinus4Phi"                ,  &bestWMinus4Phi                ) ;
  outTree->Branch( "bestWMinus4E"                  ,  &bestWMinus4E                  ) ;
  outTree->Branch( "bestWMinus4M"                  ,  &bestWMinus4M                  ) ;
  outTree->Branch( "bestNeutrino4Pt"               ,  &bestNeutrino4Pt               ) ;
  outTree->Branch( "bestNeutrino4Eta"              ,  &bestNeutrino4Eta              ) ;
  outTree->Branch( "bestNeutrino4Phi"              ,  &bestNeutrino4Phi              ) ;
  outTree->Branch( "bestNeutrino4E"                ,  &bestNeutrino4E                ) ;
  outTree->Branch( "bestNeutrino4PairingTruth"     ,  &bestNeutrino4PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino4Pt"           ,  &bestAntiNeutrino4Pt           ) ;
  outTree->Branch( "bestAntiNeutrino4Eta"          ,  &bestAntiNeutrino4Eta          ) ;
  outTree->Branch( "bestAntiNeutrino4Phi"          ,  &bestAntiNeutrino4Phi          ) ;
  outTree->Branch( "bestAntiNeutrino4E"            ,  &bestAntiNeutrino4E            ) ;
  outTree->Branch( "bestAntiNeutrino4PairingTruth" ,  &bestAntiNeutrino4PairingTruth ) ;
  outTree->Branch( "bestTopQuark5Pt"               ,  &bestTopQuark5Pt               ) ;
  outTree->Branch( "bestTopQuark5Eta"              ,  &bestTopQuark5Eta              ) ;
  outTree->Branch( "bestTopQuark5Phi"              ,  &bestTopQuark5Phi              ) ;
  outTree->Branch( "bestTopQuark5E"                ,  &bestTopQuark5E                ) ;
  outTree->Branch( "bestTopQuark5M"                ,  &bestTopQuark5M                ) ;
  outTree->Branch( "bestAntiTopQuark5Pt"           ,  &bestAntiTopQuark5Pt           ) ;
  outTree->Branch( "bestAntiTopQuark5Eta"          ,  &bestAntiTopQuark5Eta          ) ;
  outTree->Branch( "bestAntiTopQuark5Phi"          ,  &bestAntiTopQuark5Phi          ) ;
  outTree->Branch( "bestAntiTopQuark5E"            ,  &bestAntiTopQuark5E            ) ;
  outTree->Branch( "bestAntiTopQuark5M"            ,  &bestAntiTopQuark5M            ) ;
  outTree->Branch( "bestWPlus5Pt"                  ,  &bestWPlus5Pt                  ) ;
  outTree->Branch( "bestWPlus5Eta"                 ,  &bestWPlus5Eta                 ) ;
  outTree->Branch( "bestWPlus5Phi"                 ,  &bestWPlus5Phi                 ) ;
  outTree->Branch( "bestWPlus5E"                   ,  &bestWPlus5E                   ) ;
  outTree->Branch( "bestWPlus5M"                   ,  &bestWPlus5M                   ) ;
  outTree->Branch( "bestWMinus5Pt"                 ,  &bestWMinus5Pt                 ) ;
  outTree->Branch( "bestWMinus5Eta"                ,  &bestWMinus5Eta                ) ;
  outTree->Branch( "bestWMinus5Phi"                ,  &bestWMinus5Phi                ) ;
  outTree->Branch( "bestWMinus5E"                  ,  &bestWMinus5E                  ) ;
  outTree->Branch( "bestWMinus5M"                  ,  &bestWMinus5M                  ) ;
  outTree->Branch( "bestNeutrino5Pt"               ,  &bestNeutrino5Pt               ) ;
  outTree->Branch( "bestNeutrino5Eta"              ,  &bestNeutrino5Eta              ) ;
  outTree->Branch( "bestNeutrino5Phi"              ,  &bestNeutrino5Phi              ) ;
  outTree->Branch( "bestNeutrino5E"                ,  &bestNeutrino5E                ) ;
  outTree->Branch( "bestNeutrino5PairingTruth"     ,  &bestNeutrino5PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino5Pt"           ,  &bestAntiNeutrino5Pt           ) ;
  outTree->Branch( "bestAntiNeutrino5Eta"          ,  &bestAntiNeutrino5Eta          ) ;
  outTree->Branch( "bestAntiNeutrino5Phi"          ,  &bestAntiNeutrino5Phi          ) ;
  outTree->Branch( "bestAntiNeutrino5E"            ,  &bestAntiNeutrino5E            ) ;
  outTree->Branch( "bestAntiNeutrino5PairingTruth" ,  &bestAntiNeutrino5PairingTruth ) ;
  outTree->Branch( "bestTopQuark6Pt"               ,  &bestTopQuark6Pt               ) ;
  outTree->Branch( "bestTopQuark6Eta"              ,  &bestTopQuark6Eta              ) ;
  outTree->Branch( "bestTopQuark6Phi"              ,  &bestTopQuark6Phi              ) ;
  outTree->Branch( "bestTopQuark6E"                ,  &bestTopQuark6E                ) ;
  outTree->Branch( "bestTopQuark6M"                ,  &bestTopQuark6M                ) ;
  outTree->Branch( "bestAntiTopQuark6Pt"           ,  &bestAntiTopQuark6Pt           ) ;
  outTree->Branch( "bestAntiTopQuark6Eta"          ,  &bestAntiTopQuark6Eta          ) ;
  outTree->Branch( "bestAntiTopQuark6Phi"          ,  &bestAntiTopQuark6Phi          ) ;
  outTree->Branch( "bestAntiTopQuark6E"            ,  &bestAntiTopQuark6E            ) ;
  outTree->Branch( "bestAntiTopQuark6M"            ,  &bestAntiTopQuark6M            ) ;
  outTree->Branch( "bestWPlus6Pt"                  ,  &bestWPlus6Pt                  ) ;
  outTree->Branch( "bestWPlus6Eta"                 ,  &bestWPlus6Eta                 ) ;
  outTree->Branch( "bestWPlus6Phi"                 ,  &bestWPlus6Phi                 ) ;
  outTree->Branch( "bestWPlus6E"                   ,  &bestWPlus6E                   ) ;
  outTree->Branch( "bestWPlus6M"                   ,  &bestWPlus6M                   ) ;
  outTree->Branch( "bestWMinus6Pt"                 ,  &bestWMinus6Pt                 ) ;
  outTree->Branch( "bestWMinus6Eta"                ,  &bestWMinus6Eta                ) ;
  outTree->Branch( "bestWMinus6Phi"                ,  &bestWMinus6Phi                ) ;
  outTree->Branch( "bestWMinus6E"                  ,  &bestWMinus6E                  ) ;
  outTree->Branch( "bestWMinus6M"                  ,  &bestWMinus6M                  ) ;
  outTree->Branch( "bestNeutrino6Pt"               ,  &bestNeutrino6Pt               ) ;
  outTree->Branch( "bestNeutrino6Eta"              ,  &bestNeutrino6Eta              ) ;
  outTree->Branch( "bestNeutrino6Phi"              ,  &bestNeutrino6Phi              ) ;
  outTree->Branch( "bestNeutrino6E"                ,  &bestNeutrino6E                ) ;
  outTree->Branch( "bestNeutrino6PairingTruth"     ,  &bestNeutrino6PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino6Pt"           ,  &bestAntiNeutrino6Pt           ) ;
  outTree->Branch( "bestAntiNeutrino6Eta"          ,  &bestAntiNeutrino6Eta          ) ;
  outTree->Branch( "bestAntiNeutrino6Phi"          ,  &bestAntiNeutrino6Phi          ) ;
  outTree->Branch( "bestAntiNeutrino6E"            ,  &bestAntiNeutrino6E            ) ;
  outTree->Branch( "bestAntiNeutrino6PairingTruth" ,  &bestAntiNeutrino6PairingTruth ) ;
  outTree->Branch( "bestTopQuark7Pt"               ,  &bestTopQuark7Pt               ) ;
  outTree->Branch( "bestTopQuark7Eta"              ,  &bestTopQuark7Eta              ) ;
  outTree->Branch( "bestTopQuark7Phi"              ,  &bestTopQuark7Phi              ) ;
  outTree->Branch( "bestTopQuark7E"                ,  &bestTopQuark7E                ) ;
  outTree->Branch( "bestTopQuark7M"                ,  &bestTopQuark7M                ) ;
  outTree->Branch( "bestAntiTopQuark7Pt"           ,  &bestAntiTopQuark7Pt           ) ;
  outTree->Branch( "bestAntiTopQuark7Eta"          ,  &bestAntiTopQuark7Eta          ) ;
  outTree->Branch( "bestAntiTopQuark7Phi"          ,  &bestAntiTopQuark7Phi          ) ;
  outTree->Branch( "bestAntiTopQuark7E"            ,  &bestAntiTopQuark7E            ) ;
  outTree->Branch( "bestAntiTopQuark7M"            ,  &bestAntiTopQuark7M            ) ;
  outTree->Branch( "bestWPlus7Pt"                  ,  &bestWPlus7Pt                  ) ;
  outTree->Branch( "bestWPlus7Eta"                 ,  &bestWPlus7Eta                 ) ;
  outTree->Branch( "bestWPlus7Phi"                 ,  &bestWPlus7Phi                 ) ;
  outTree->Branch( "bestWPlus7E"                   ,  &bestWPlus7E                   ) ;
  outTree->Branch( "bestWPlus7M"                   ,  &bestWPlus7M                   ) ;
  outTree->Branch( "bestWMinus7Pt"                 ,  &bestWMinus7Pt                 ) ;
  outTree->Branch( "bestWMinus7Eta"                ,  &bestWMinus7Eta                ) ;
  outTree->Branch( "bestWMinus7Phi"                ,  &bestWMinus7Phi                ) ;
  outTree->Branch( "bestWMinus7E"                  ,  &bestWMinus7E                  ) ;
  outTree->Branch( "bestWMinus7M"                  ,  &bestWMinus7M                  ) ;
  outTree->Branch( "bestNeutrino7Pt"               ,  &bestNeutrino7Pt               ) ;
  outTree->Branch( "bestNeutrino7Eta"              ,  &bestNeutrino7Eta              ) ;
  outTree->Branch( "bestNeutrino7Phi"              ,  &bestNeutrino7Phi              ) ;
  outTree->Branch( "bestNeutrino7E"                ,  &bestNeutrino7E                ) ;
  outTree->Branch( "bestNeutrino7PairingTruth"     ,  &bestNeutrino7PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino7Pt"           ,  &bestAntiNeutrino7Pt           ) ;
  outTree->Branch( "bestAntiNeutrino7Eta"          ,  &bestAntiNeutrino7Eta          ) ;
  outTree->Branch( "bestAntiNeutrino7Phi"          ,  &bestAntiNeutrino7Phi          ) ;
  outTree->Branch( "bestAntiNeutrino7E"            ,  &bestAntiNeutrino7E            ) ;
  outTree->Branch( "bestAntiNeutrino7PairingTruth" ,  &bestAntiNeutrino7PairingTruth ) ;
  outTree->Branch( "bestTopQuark8Pt"               ,  &bestTopQuark8Pt               ) ;
  outTree->Branch( "bestTopQuark8Eta"              ,  &bestTopQuark8Eta              ) ;
  outTree->Branch( "bestTopQuark8Phi"              ,  &bestTopQuark8Phi              ) ;
  outTree->Branch( "bestTopQuark8E"                ,  &bestTopQuark8E                ) ;
  outTree->Branch( "bestTopQuark8M"                ,  &bestTopQuark8M                ) ;
  outTree->Branch( "bestAntiTopQuark8Pt"           ,  &bestAntiTopQuark8Pt           ) ;
  outTree->Branch( "bestAntiTopQuark8Eta"          ,  &bestAntiTopQuark8Eta          ) ;
  outTree->Branch( "bestAntiTopQuark8Phi"          ,  &bestAntiTopQuark8Phi          ) ;
  outTree->Branch( "bestAntiTopQuark8E"            ,  &bestAntiTopQuark8E            ) ;
  outTree->Branch( "bestAntiTopQuark8M"            ,  &bestAntiTopQuark8M            ) ;
  outTree->Branch( "bestWPlus8Pt"                  ,  &bestWPlus8Pt                  ) ;
  outTree->Branch( "bestWPlus8Eta"                 ,  &bestWPlus8Eta                 ) ;
  outTree->Branch( "bestWPlus8Phi"                 ,  &bestWPlus8Phi                 ) ;
  outTree->Branch( "bestWPlus8E"                   ,  &bestWPlus8E                   ) ;
  outTree->Branch( "bestWPlus8M"                   ,  &bestWPlus8M                   ) ;
  outTree->Branch( "bestWMinus8Pt"                 ,  &bestWMinus8Pt                 ) ;
  outTree->Branch( "bestWMinus8Eta"                ,  &bestWMinus8Eta                ) ;
  outTree->Branch( "bestWMinus8Phi"                ,  &bestWMinus8Phi                ) ;
  outTree->Branch( "bestWMinus8E"                  ,  &bestWMinus8E                  ) ;
  outTree->Branch( "bestWMinus8M"                  ,  &bestWMinus8M                  ) ;
  outTree->Branch( "bestNeutrino8Pt"               ,  &bestNeutrino8Pt               ) ;
  outTree->Branch( "bestNeutrino8Eta"              ,  &bestNeutrino8Eta              ) ;
  outTree->Branch( "bestNeutrino8Phi"              ,  &bestNeutrino8Phi              ) ;
  outTree->Branch( "bestNeutrino8E"                ,  &bestNeutrino8E                ) ;
  outTree->Branch( "bestNeutrino8PairingTruth"     ,  &bestNeutrino8PairingTruth     ) ;
  outTree->Branch( "bestAntiNeutrino8Pt"           ,  &bestAntiNeutrino8Pt           ) ;
  outTree->Branch( "bestAntiNeutrino8Eta"          ,  &bestAntiNeutrino8Eta          ) ;
  outTree->Branch( "bestAntiNeutrino8Phi"          ,  &bestAntiNeutrino8Phi          ) ;
  outTree->Branch( "bestAntiNeutrino8E"            ,  &bestAntiNeutrino8E            ) ;
  outTree->Branch( "bestAntiNeutrino8PairingTruth" ,  &bestAntiNeutrino8PairingTruth ) ;
  outTree->Branch( "bestBottomQuarkPt"             ,  &bestBottomQuarkPt             ) ;
  outTree->Branch( "bestBottomQuarkEta"            ,  &bestBottomQuarkEta            ) ;
  outTree->Branch( "bestBottomQuarkPhi"            ,  &bestBottomQuarkPhi            ) ;
  outTree->Branch( "bestBottomQuarkE"              ,  &bestBottomQuarkE              ) ;
  outTree->Branch( "bestAntiBottomQuarkPt"         ,  &bestAntiBottomQuarkPt         ) ;
  outTree->Branch( "bestAntiBottomQuarkEta"        ,  &bestAntiBottomQuarkEta        ) ;
  outTree->Branch( "bestAntiBottomQuarkPhi"        ,  &bestAntiBottomQuarkPhi        ) ;
  outTree->Branch( "bestAntiBottomQuarkE"          ,  &bestAntiBottomQuarkE          ) ;
  outTree->Branch( "bestLightParton1Pt"            ,  &bestLightParton1Pt            ) ;
  outTree->Branch( "bestLightParton1Eta"           ,  &bestLightParton1Eta           ) ;
  outTree->Branch( "bestLightParton1Phi"           ,  &bestLightParton1Phi           ) ;
  outTree->Branch( "bestLightParton1E"             ,  &bestLightParton1E             ) ;
  outTree->Branch( "bestLightParton2Pt "           ,  &bestLightParton2Pt            ) ;
  outTree->Branch( "bestLightParton2Eta"           ,  &bestLightParton2Eta           ) ;
  outTree->Branch( "bestLightParton2Phi"           ,  &bestLightParton2Phi           ) ;
  outTree->Branch( "bestLightParton2E  "           ,  &bestLightParton2E             ) ;


  outTree->Branch( "deltaPtGenSmearedTopQuark"        ,  &deltaPtGenSmearedTopQuark        ) ;
  outTree->Branch( "deltaPtGenSmearedAntiTopQuark"    ,  &deltaPtGenSmearedAntiTopQuark    ) ;
  outTree->Branch( "deltaRGenSmearedTopQuark"         ,  &deltaRGenSmearedTopQuark         ) ;
  outTree->Branch( "deltaRGenSmearedAntiTopQuark"     ,  &deltaRGenSmearedAntiTopQuark     ) ;
  outTree->Branch( "deltaMGenSmearedTopQuark"         ,  &deltaMGenSmearedTopQuark         ) ;
  outTree->Branch( "deltaMGenSmearedAntiTopQuark"     ,  &deltaMGenSmearedAntiTopQuark     ) ;
  outTree->Branch( "deltaPtGenSmearedWPlus"           ,  &deltaPtGenSmearedWPlus           ) ;
  outTree->Branch( "deltaPtGenSmearedWMinus"          ,  &deltaPtGenSmearedWMinus          ) ;
  outTree->Branch( "deltaRGenSmearedWPlus"            ,  &deltaRGenSmearedWPlus            ) ;
  outTree->Branch( "deltaRGenSmearedWMinus"           ,  &deltaRGenSmearedWMinus           ) ;
  outTree->Branch( "deltaMGenSmearedWPlus"            ,  &deltaMGenSmearedWPlus            ) ;
  outTree->Branch( "deltaMGenSmearedWMinus"           ,  &deltaMGenSmearedWMinus           ) ;
  outTree->Branch( "deltaPtGenSmearedBottomQuark"     ,  &deltaPtGenSmearedBottomQuark     ) ;
  outTree->Branch( "deltaPtGenSmearedAntiBottomQuark" ,  &deltaPtGenSmearedAntiBottomQuark ) ;
  outTree->Branch( "deltaRGenSmearedBottomQuark"      ,  &deltaRGenSmearedBottomQuark      ) ;
  outTree->Branch( "deltaRGenSmearedAntiBottomQuark"  ,  &deltaRGenSmearedAntiBottomQuark  ) ;
  outTree->Branch( "deltaPtGenSmearedLightParton1"    ,  &deltaPtGenSmearedLightParton1    ) ;
  outTree->Branch( "deltaPtGenSmearedLightParton2"    ,  &deltaPtGenSmearedLightParton2    ) ;
  outTree->Branch( "deltaPxGenSmearedLightParton1"    ,  &deltaPxGenSmearedLightParton1    ) ; 
  outTree->Branch( "deltaPxGenSmearedLightParton2"    ,  &deltaPxGenSmearedLightParton2    ) ;
  outTree->Branch( "deltaPyGenSmearedLightParton1"    ,  &deltaPyGenSmearedLightParton1    ) ; 
  outTree->Branch( "deltaPyGenSmearedLightParton2"    ,  &deltaPyGenSmearedLightParton2    ) ;
  outTree->Branch( "deltaRGenSmearedLightParton1"     ,  &deltaRGenSmearedLightParton1     ) ;
  outTree->Branch( "deltaRGenSmearedLightParton2"     ,  &deltaRGenSmearedLightParton2     ) ;

  outTree->Branch( "deltaPtGenBestBottomQuark"     ,  &deltaPtGenBestBottomQuark     ) ;
  outTree->Branch( "deltaPtGenBestAntiBottomQuark" ,  &deltaPtGenBestAntiBottomQuark ) ;
  outTree->Branch( "deltaRGenBestBottomQuark"      ,  &deltaRGenBestBottomQuark      ) ;
  outTree->Branch( "deltaRGenBestAntiBottomQuark"  ,  &deltaRGenBestAntiBottomQuark  ) ;
  outTree->Branch( "deltaPtGenBestLightParton1"    ,  &deltaPtGenBestLightParton1    ) ;
  outTree->Branch( "deltaPtGenBestLightParton2"    ,  &deltaPtGenBestLightParton2    ) ;
  outTree->Branch( "deltaPxGenBestLightParton1"    ,  &deltaPxGenBestLightParton1    ) ;
  outTree->Branch( "deltaPxGenBestLightParton2"    ,  &deltaPxGenBestLightParton2    ) ;
  outTree->Branch( "deltaPyGenBestLightParton1"    ,  &deltaPyGenBestLightParton1    ) ;
  outTree->Branch( "deltaPyGenBestLightParton2"    ,  &deltaPyGenBestLightParton2    ) ;
  outTree->Branch( "deltaRGenBestLightParton1"     ,  &deltaRGenBestLightParton1     ) ;
  outTree->Branch( "deltaRGenBestLightParton2"     ,  &deltaRGenBestLightParton2     ) ;
  outTree->Branch( "deltaPtGenBestTopQuark1"       ,  &deltaPtGenBestTopQuark1       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark1"   ,  &deltaPtGenBestAntiTopQuark1   ) ;
  outTree->Branch( "deltaRGenBestTopQuark1"        ,  &deltaRGenBestTopQuark1        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark1"    ,  &deltaRGenBestAntiTopQuark1    ) ;
  outTree->Branch( "deltaMGenBestTopQuark1"        ,  &deltaMGenBestTopQuark1        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark1"    ,  &deltaMGenBestAntiTopQuark1    ) ;
  outTree->Branch( "deltaPtGenBestWPlus1"          ,  &deltaPtGenBestWPlus1          ) ;
  outTree->Branch( "deltaPtGenBestWMinus1"         ,  &deltaPtGenBestWMinus1         ) ;
  outTree->Branch( "deltaRGenBestWPlus1"           ,  &deltaRGenBestWPlus1           ) ;
  outTree->Branch( "deltaRGenBestWMinus1"          ,  &deltaRGenBestWMinus1          ) ;
  outTree->Branch( "deltaMGenBestWPlus1"           ,  &deltaMGenBestWPlus1           ) ;
  outTree->Branch( "deltaMGenBestWMinus1"          ,  &deltaMGenBestWMinus1          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino1"       ,  &deltaPtGenBestNeutrino1       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino1"   ,  &deltaPtGenBestAntiNeutrino1   ) ;
  outTree->Branch( "deltaRGenBestNeutrino1"        ,  &deltaRGenBestNeutrino1        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino1"    ,  &deltaRGenBestAntiNeutrino1    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark2"       ,  &deltaPtGenBestTopQuark2       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark2"   ,  &deltaPtGenBestAntiTopQuark2   ) ;
  outTree->Branch( "deltaRGenBestTopQuark2"        ,  &deltaRGenBestTopQuark2        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark2"    ,  &deltaRGenBestAntiTopQuark2    ) ;
  outTree->Branch( "deltaMGenBestTopQuark2"        ,  &deltaMGenBestTopQuark2        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark2"    ,  &deltaMGenBestAntiTopQuark2    ) ;
  outTree->Branch( "deltaPtGenBestWPlus2"          ,  &deltaPtGenBestWPlus2          ) ;
  outTree->Branch( "deltaPtGenBestWMinus2"         ,  &deltaPtGenBestWMinus2         ) ;
  outTree->Branch( "deltaRGenBestWPlus2"           ,  &deltaRGenBestWPlus2           ) ;
  outTree->Branch( "deltaRGenBestWMinus2"          ,  &deltaRGenBestWMinus2          ) ;
  outTree->Branch( "deltaMGenBestWPlus2"           ,  &deltaMGenBestWPlus2           ) ;
  outTree->Branch( "deltaMGenBestWMinus2"          ,  &deltaMGenBestWMinus2          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino2"       ,  &deltaPtGenBestNeutrino2       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino2"   ,  &deltaPtGenBestAntiNeutrino2   ) ;
  outTree->Branch( "deltaRGenBestNeutrino2"        ,  &deltaRGenBestNeutrino2        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino2"    ,  &deltaRGenBestAntiNeutrino2    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark3"       ,  &deltaPtGenBestTopQuark3       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark3"   ,  &deltaPtGenBestAntiTopQuark3   ) ;
  outTree->Branch( "deltaRGenBestTopQuark3"        ,  &deltaRGenBestTopQuark3        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark3"    ,  &deltaRGenBestAntiTopQuark3    ) ;
  outTree->Branch( "deltaMGenBestTopQuark3"        ,  &deltaMGenBestTopQuark3        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark3"    ,  &deltaMGenBestAntiTopQuark3    ) ;
  outTree->Branch( "deltaPtGenBestWPlus3"          ,  &deltaPtGenBestWPlus3          ) ;
  outTree->Branch( "deltaPtGenBestWMinus3"         ,  &deltaPtGenBestWMinus3         ) ;
  outTree->Branch( "deltaRGenBestWPlus3"           ,  &deltaRGenBestWPlus3           ) ;
  outTree->Branch( "deltaRGenBestWMinus3"          ,  &deltaRGenBestWMinus3          ) ;
  outTree->Branch( "deltaMGenBestWPlus3"           ,  &deltaMGenBestWPlus3           ) ;
  outTree->Branch( "deltaMGenBestWMinus3"          ,  &deltaMGenBestWMinus3          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino3"       ,  &deltaPtGenBestNeutrino3       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino3"   ,  &deltaPtGenBestAntiNeutrino3   ) ;
  outTree->Branch( "deltaRGenBestNeutrino3"        ,  &deltaRGenBestNeutrino3        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino3"    ,  &deltaRGenBestAntiNeutrino3    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark4"       ,  &deltaPtGenBestTopQuark4       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark4"   ,  &deltaPtGenBestAntiTopQuark4   ) ;
  outTree->Branch( "deltaRGenBestTopQuark4"        ,  &deltaRGenBestTopQuark4        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark4"    ,  &deltaRGenBestAntiTopQuark4    ) ;
  outTree->Branch( "deltaMGenBestTopQuark4"        ,  &deltaMGenBestTopQuark4        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark4"    ,  &deltaMGenBestAntiTopQuark4    ) ;
  outTree->Branch( "deltaPtGenBestWPlus4"          ,  &deltaPtGenBestWPlus4          ) ;
  outTree->Branch( "deltaPtGenBestWMinus4"         ,  &deltaPtGenBestWMinus4         ) ;
  outTree->Branch( "deltaRGenBestWPlus4"           ,  &deltaRGenBestWPlus4           ) ;
  outTree->Branch( "deltaRGenBestWMinus4"          ,  &deltaRGenBestWMinus4          ) ;
  outTree->Branch( "deltaMGenBestWPlus4"           ,  &deltaMGenBestWPlus4           ) ;
  outTree->Branch( "deltaMGenBestWMinus4"          ,  &deltaMGenBestWMinus4          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino4"       ,  &deltaPtGenBestNeutrino4       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino4"   ,  &deltaPtGenBestAntiNeutrino4   ) ;
  outTree->Branch( "deltaRGenBestNeutrino4"        ,  &deltaRGenBestNeutrino4        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino4"    ,  &deltaRGenBestAntiNeutrino4    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark5"       ,  &deltaPtGenBestTopQuark5       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark5"   ,  &deltaPtGenBestAntiTopQuark5   ) ;
  outTree->Branch( "deltaRGenBestTopQuark5"        ,  &deltaRGenBestTopQuark5        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark5"    ,  &deltaRGenBestAntiTopQuark5    ) ;
  outTree->Branch( "deltaMGenBestTopQuark5"        ,  &deltaMGenBestTopQuark5        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark5"    ,  &deltaMGenBestAntiTopQuark5    ) ;
  outTree->Branch( "deltaPtGenBestWPlus5"          ,  &deltaPtGenBestWPlus5          ) ;
  outTree->Branch( "deltaPtGenBestWMinus5"         ,  &deltaPtGenBestWMinus5         ) ;
  outTree->Branch( "deltaRGenBestWPlus5"           ,  &deltaRGenBestWPlus5           ) ;
  outTree->Branch( "deltaRGenBestWMinus5"          ,  &deltaRGenBestWMinus5          ) ;
  outTree->Branch( "deltaMGenBestWPlus5"           ,  &deltaMGenBestWPlus5           ) ;
  outTree->Branch( "deltaMGenBestWMinus5"          ,  &deltaMGenBestWMinus5          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino5"       ,  &deltaPtGenBestNeutrino5       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino5"   ,  &deltaPtGenBestAntiNeutrino5   ) ;
  outTree->Branch( "deltaRGenBestNeutrino5"        ,  &deltaRGenBestNeutrino5        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino5"    ,  &deltaRGenBestAntiNeutrino5    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark6"       ,  &deltaPtGenBestTopQuark6       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark6"   ,  &deltaPtGenBestAntiTopQuark6   ) ;
  outTree->Branch( "deltaRGenBestTopQuark6"        ,  &deltaRGenBestTopQuark6        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark6"    ,  &deltaRGenBestAntiTopQuark6    ) ;
  outTree->Branch( "deltaMGenBestTopQuark6"        ,  &deltaMGenBestTopQuark6        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark6"    ,  &deltaMGenBestAntiTopQuark6    ) ;
  outTree->Branch( "deltaPtGenBestWPlus6"          ,  &deltaPtGenBestWPlus6          ) ;
  outTree->Branch( "deltaPtGenBestWMinus6"         ,  &deltaPtGenBestWMinus6         ) ;
  outTree->Branch( "deltaRGenBestWPlus6"           ,  &deltaRGenBestWPlus6           ) ;
  outTree->Branch( "deltaRGenBestWMinus6"          ,  &deltaRGenBestWMinus6          ) ;
  outTree->Branch( "deltaMGenBestWPlus6"           ,  &deltaMGenBestWPlus6           ) ;
  outTree->Branch( "deltaMGenBestWMinus6"          ,  &deltaMGenBestWMinus6          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino6"       ,  &deltaPtGenBestNeutrino6       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino6"   ,  &deltaPtGenBestAntiNeutrino6   ) ;
  outTree->Branch( "deltaRGenBestNeutrino6"        ,  &deltaRGenBestNeutrino6        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino6"    ,  &deltaRGenBestAntiNeutrino6    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark7"       ,  &deltaPtGenBestTopQuark7       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark7"   ,  &deltaPtGenBestAntiTopQuark7   ) ;
  outTree->Branch( "deltaRGenBestTopQuark7"        ,  &deltaRGenBestTopQuark7        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark7"    ,  &deltaRGenBestAntiTopQuark7    ) ;
  outTree->Branch( "deltaMGenBestTopQuark7"        ,  &deltaMGenBestTopQuark7        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark7"    ,  &deltaMGenBestAntiTopQuark7    ) ;
  outTree->Branch( "deltaPtGenBestWPlus7"          ,  &deltaPtGenBestWPlus7          ) ;
  outTree->Branch( "deltaPtGenBestWMinus7"         ,  &deltaPtGenBestWMinus7         ) ;
  outTree->Branch( "deltaRGenBestWPlus7"           ,  &deltaRGenBestWPlus7           ) ;
  outTree->Branch( "deltaRGenBestWMinus7"          ,  &deltaRGenBestWMinus7          ) ;
  outTree->Branch( "deltaMGenBestWPlus7"           ,  &deltaMGenBestWPlus7           ) ;
  outTree->Branch( "deltaMGenBestWMinus7"          ,  &deltaMGenBestWMinus7          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino7"       ,  &deltaPtGenBestNeutrino7       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino7"   ,  &deltaPtGenBestAntiNeutrino7   ) ;
  outTree->Branch( "deltaRGenBestNeutrino7"        ,  &deltaRGenBestNeutrino7        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino7"    ,  &deltaRGenBestAntiNeutrino7    ) ;
  outTree->Branch( "deltaPtGenBestTopQuark8"       ,  &deltaPtGenBestTopQuark8       ) ;
  outTree->Branch( "deltaPtGenBestAntiTopQuark8"   ,  &deltaPtGenBestAntiTopQuark8   ) ;
  outTree->Branch( "deltaRGenBestTopQuark8"        ,  &deltaRGenBestTopQuark8        ) ;
  outTree->Branch( "deltaRGenBestAntiTopQuark8"    ,  &deltaRGenBestAntiTopQuark8    ) ;
  outTree->Branch( "deltaMGenBestTopQuark8"        ,  &deltaMGenBestTopQuark8        ) ;
  outTree->Branch( "deltaMGenBestAntiTopQuark8"    ,  &deltaMGenBestAntiTopQuark8    ) ;
  outTree->Branch( "deltaPtGenBestWPlus8"          ,  &deltaPtGenBestWPlus8          ) ;
  outTree->Branch( "deltaPtGenBestWMinus8"         ,  &deltaPtGenBestWMinus8         ) ;
  outTree->Branch( "deltaRGenBestWPlus8"           ,  &deltaRGenBestWPlus8           ) ;
  outTree->Branch( "deltaRGenBestWMinus8"          ,  &deltaRGenBestWMinus8          ) ;
  outTree->Branch( "deltaMGenBestWPlus8"           ,  &deltaMGenBestWPlus8           ) ;
  outTree->Branch( "deltaMGenBestWMinus8"          ,  &deltaMGenBestWMinus8          ) ;
  outTree->Branch( "deltaPtGenBestNeutrino8"       ,  &deltaPtGenBestNeutrino8       ) ;
  outTree->Branch( "deltaPtGenBestAntiNeutrino8"   ,  &deltaPtGenBestAntiNeutrino8   ) ;
  outTree->Branch( "deltaRGenBestNeutrino8"        ,  &deltaRGenBestNeutrino8        ) ;
  outTree->Branch( "deltaRGenBestAntiNeutrino8"    ,  &deltaRGenBestAntiNeutrino8    ) ;


  outTree->Branch( "deltaPtSmearedBestTopQuark1"       ,  &deltaPtSmearedBestTopQuark1       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark1"   ,  &deltaPtSmearedBestAntiTopQuark1   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark1"        ,  &deltaRSmearedBestTopQuark1        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark1"    ,  &deltaRSmearedBestAntiTopQuark1    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark1"        ,  &deltaMSmearedBestTopQuark1        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark1"    ,  &deltaMSmearedBestAntiTopQuark1    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus1"          ,  &deltaPtSmearedBestWPlus1          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus1"         ,  &deltaPtSmearedBestWMinus1         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus1"           ,  &deltaRSmearedBestWPlus1           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus1"          ,  &deltaRSmearedBestWMinus1          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus1"           ,  &deltaMSmearedBestWPlus1           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus1"          ,  &deltaMSmearedBestWMinus1          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino1"       ,  &deltaPtSmearedBestNeutrino1       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino1"   ,  &deltaPtSmearedBestAntiNeutrino1   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino1"        ,  &deltaRSmearedBestNeutrino1        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino1"    ,  &deltaRSmearedBestAntiNeutrino1    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark2"       ,  &deltaPtSmearedBestTopQuark2       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark2"   ,  &deltaPtSmearedBestAntiTopQuark2   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark2"        ,  &deltaRSmearedBestTopQuark2        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark2"    ,  &deltaRSmearedBestAntiTopQuark2    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark2"        ,  &deltaMSmearedBestTopQuark2        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark2"    ,  &deltaMSmearedBestAntiTopQuark2    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus2"          ,  &deltaPtSmearedBestWPlus2          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus2"         ,  &deltaPtSmearedBestWMinus2         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus2"           ,  &deltaRSmearedBestWPlus2           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus2"          ,  &deltaRSmearedBestWMinus2          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus2"           ,  &deltaMSmearedBestWPlus2           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus2"          ,  &deltaMSmearedBestWMinus2          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino2"       ,  &deltaPtSmearedBestNeutrino2       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino2"   ,  &deltaPtSmearedBestAntiNeutrino2   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino2"        ,  &deltaRSmearedBestNeutrino2        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino2"    ,  &deltaRSmearedBestAntiNeutrino2    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark3"       ,  &deltaPtSmearedBestTopQuark3       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark3"   ,  &deltaPtSmearedBestAntiTopQuark3   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark3"        ,  &deltaRSmearedBestTopQuark3        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark3"    ,  &deltaRSmearedBestAntiTopQuark3    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark3"        ,  &deltaMSmearedBestTopQuark3        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark3"    ,  &deltaMSmearedBestAntiTopQuark3    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus3"          ,  &deltaPtSmearedBestWPlus3          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus3"         ,  &deltaPtSmearedBestWMinus3         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus3"           ,  &deltaRSmearedBestWPlus3           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus3"          ,  &deltaRSmearedBestWMinus3          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus3"           ,  &deltaMSmearedBestWPlus3           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus3"          ,  &deltaMSmearedBestWMinus3          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino3"       ,  &deltaPtSmearedBestNeutrino3       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino3"   ,  &deltaPtSmearedBestAntiNeutrino3   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino3"        ,  &deltaRSmearedBestNeutrino3        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino3"    ,  &deltaRSmearedBestAntiNeutrino3    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark4"       ,  &deltaPtSmearedBestTopQuark4       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark4"   ,  &deltaPtSmearedBestAntiTopQuark4   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark4"        ,  &deltaRSmearedBestTopQuark4        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark4"    ,  &deltaRSmearedBestAntiTopQuark4    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark4"        ,  &deltaMSmearedBestTopQuark4        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark4"    ,  &deltaMSmearedBestAntiTopQuark4    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus4"          ,  &deltaPtSmearedBestWPlus4          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus4"         ,  &deltaPtSmearedBestWMinus4         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus4"           ,  &deltaRSmearedBestWPlus4           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus4"          ,  &deltaRSmearedBestWMinus4          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus4"           ,  &deltaMSmearedBestWPlus4           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus4"          ,  &deltaMSmearedBestWMinus4          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino4"       ,  &deltaPtSmearedBestNeutrino4       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino4"   ,  &deltaPtSmearedBestAntiNeutrino4   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino4"        ,  &deltaRSmearedBestNeutrino4        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino4"    ,  &deltaRSmearedBestAntiNeutrino4    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark5"       ,  &deltaPtSmearedBestTopQuark5       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark5"   ,  &deltaPtSmearedBestAntiTopQuark5   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark5"        ,  &deltaRSmearedBestTopQuark5        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark5"    ,  &deltaRSmearedBestAntiTopQuark5    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark5"        ,  &deltaMSmearedBestTopQuark5        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark5"    ,  &deltaMSmearedBestAntiTopQuark5    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus5"          ,  &deltaPtSmearedBestWPlus5          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus5"         ,  &deltaPtSmearedBestWMinus5         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus5"           ,  &deltaRSmearedBestWPlus5           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus5"          ,  &deltaRSmearedBestWMinus5          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus5"           ,  &deltaMSmearedBestWPlus5           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus5"          ,  &deltaMSmearedBestWMinus5          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino5"       ,  &deltaPtSmearedBestNeutrino5       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino5"   ,  &deltaPtSmearedBestAntiNeutrino5   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino5"        ,  &deltaRSmearedBestNeutrino5        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino5"    ,  &deltaRSmearedBestAntiNeutrino5    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark6"       ,  &deltaPtSmearedBestTopQuark6       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark6"   ,  &deltaPtSmearedBestAntiTopQuark6   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark6"        ,  &deltaRSmearedBestTopQuark6        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark6"    ,  &deltaRSmearedBestAntiTopQuark6    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark6"        ,  &deltaMSmearedBestTopQuark6        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark6"    ,  &deltaMSmearedBestAntiTopQuark6    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus6"          ,  &deltaPtSmearedBestWPlus6          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus6"         ,  &deltaPtSmearedBestWMinus6         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus6"           ,  &deltaRSmearedBestWPlus6           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus6"          ,  &deltaRSmearedBestWMinus6          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus6"           ,  &deltaMSmearedBestWPlus6           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus6"          ,  &deltaMSmearedBestWMinus6          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino6"       ,  &deltaPtSmearedBestNeutrino6       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino6"   ,  &deltaPtSmearedBestAntiNeutrino6   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino6"        ,  &deltaRSmearedBestNeutrino6        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino6"    ,  &deltaRSmearedBestAntiNeutrino6    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark7"       ,  &deltaPtSmearedBestTopQuark7       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark7"   ,  &deltaPtSmearedBestAntiTopQuark7   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark7"        ,  &deltaRSmearedBestTopQuark7        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark7"    ,  &deltaRSmearedBestAntiTopQuark7    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark7"        ,  &deltaMSmearedBestTopQuark7        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark7"    ,  &deltaMSmearedBestAntiTopQuark7    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus7"          ,  &deltaPtSmearedBestWPlus7          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus7"         ,  &deltaPtSmearedBestWMinus7         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus7"           ,  &deltaRSmearedBestWPlus7           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus7"          ,  &deltaRSmearedBestWMinus7          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus7"           ,  &deltaMSmearedBestWPlus7           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus7"          ,  &deltaMSmearedBestWMinus7          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino7"       ,  &deltaPtSmearedBestNeutrino7       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino7"   ,  &deltaPtSmearedBestAntiNeutrino7   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino7"        ,  &deltaRSmearedBestNeutrino7        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino7"    ,  &deltaRSmearedBestAntiNeutrino7    ) ;
  outTree->Branch( "deltaPtSmearedBestTopQuark8"       ,  &deltaPtSmearedBestTopQuark8       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiTopQuark8"   ,  &deltaPtSmearedBestAntiTopQuark8   ) ;
  outTree->Branch( "deltaRSmearedBestTopQuark8"        ,  &deltaRSmearedBestTopQuark8        ) ;
  outTree->Branch( "deltaRSmearedBestAntiTopQuark8"    ,  &deltaRSmearedBestAntiTopQuark8    ) ;
  outTree->Branch( "deltaMSmearedBestTopQuark8"        ,  &deltaMSmearedBestTopQuark8        ) ;
  outTree->Branch( "deltaMSmearedBestAntiTopQuark8"    ,  &deltaMSmearedBestAntiTopQuark8    ) ;
  outTree->Branch( "deltaPtSmearedBestWPlus8"          ,  &deltaPtSmearedBestWPlus8          ) ;
  outTree->Branch( "deltaPtSmearedBestWMinus8"         ,  &deltaPtSmearedBestWMinus8         ) ;
  outTree->Branch( "deltaRSmearedBestWPlus8"           ,  &deltaRSmearedBestWPlus8           ) ;
  outTree->Branch( "deltaRSmearedBestWMinus8"          ,  &deltaRSmearedBestWMinus8          ) ;
  outTree->Branch( "deltaMSmearedBestWPlus8"           ,  &deltaMSmearedBestWPlus8           ) ;
  outTree->Branch( "deltaMSmearedBestWMinus8"          ,  &deltaMSmearedBestWMinus8          ) ;
  outTree->Branch( "deltaPtSmearedBestNeutrino8"       ,  &deltaPtSmearedBestNeutrino8       ) ;
  outTree->Branch( "deltaPtSmearedBestAntiNeutrino8"   ,  &deltaPtSmearedBestAntiNeutrino8   ) ;
  outTree->Branch( "deltaRSmearedBestNeutrino8"        ,  &deltaRSmearedBestNeutrino8        ) ;
  outTree->Branch( "deltaRSmearedBestAntiNeutrino8"    ,  &deltaRSmearedBestAntiNeutrino8    ) ;
  outTree->Branch( "deltaPtSmearedBestBottomQuark"     ,  &deltaPtSmearedBestBottomQuark     ) ;
  outTree->Branch( "deltaPtSmearedBestAntiBottomQuark" ,  &deltaPtSmearedBestAntiBottomQuark ) ;
  outTree->Branch( "deltaRSmearedBestBottomQuark"      ,  &deltaRSmearedBestBottomQuark      ) ;
  outTree->Branch( "deltaRSmearedBestAntiBottomQuark"  ,  &deltaRSmearedBestAntiBottomQuark  ) ;
  outTree->Branch( "deltaPtSmearedBestLightParton1"    ,  &deltaPtSmearedBestLightParton1    ) ;
  outTree->Branch( "deltaPtSmearedBestLightParton2"    ,  &deltaPtSmearedBestLightParton2    ) ;
  outTree->Branch( "deltaPxSmearedBestLightParton1"    ,  &deltaPxSmearedBestLightParton1    ) ;
  outTree->Branch( "deltaPxSmearedBestLightParton2"    ,  &deltaPxSmearedBestLightParton2    ) ;
  outTree->Branch( "deltaPySmearedBestLightParton1"    ,  &deltaPySmearedBestLightParton1    ) ;
  outTree->Branch( "deltaPySmearedBestLightParton2"    ,  &deltaPySmearedBestLightParton2    ) ;
  outTree->Branch( "deltaRSmearedBestLightParton1"     ,  &deltaRSmearedBestLightParton1     ) ;
  outTree->Branch( "deltaRSmearedBestLightParton2"     ,  &deltaRGenSmearedLightParton2      ) ;

  outTree->Branch( "deltaPtBottomQuarkFromMin"      ,  &deltaPtBottomQuarkFromMin      ) ;  
  outTree->Branch( "deltaPtAntiBottomQuarkFromMin"  ,  &deltaPtAntiBottomQuarkFromMin  ) ;
  outTree->Branch( "deltaPhiBottomQuarkFromMin"     ,  &deltaPhiBottomQuarkFromMin     ) ;
  outTree->Branch( "deltaPhiAntiBottomQuarkFromMin" ,  &deltaPhiAntiBottomQuarkFromMin ) ;
  outTree->Branch( "deltaEtaBottomQuarkFromMin"     ,  &deltaEtaBottomQuarkFromMin     ) ;
  outTree->Branch( "deltaPhiAntiBottomQuarkFromMin" ,  &deltaPhiAntiBottomQuarkFromMin ) ;
  outTree->Branch( "deltaPxLightParton1FromMin"     ,  &deltaPxLightParton1FromMin     ) ;
  outTree->Branch( "deltaPxLightParton2FromMin"     ,  &deltaPxLightParton2FromMin     ) ;
  outTree->Branch( "deltaPyLightParton1FromMin"     ,  &deltaPyLightParton1FromMin     ) ;
  outTree->Branch( "deltaPyLightParton2FromMin"     ,  &deltaPyLightParton2FromMin     ) ;
  outTree->Branch( "deltaMTopQuarkFromMin"          ,  &deltaMTopQuarkFromMin          ) ;
  outTree->Branch( "deltaMAntiTopQuarkFromMin"      ,  &deltaMAntiTopQuarkFromMin      ) ;
  outTree->Branch( "deltaMWPlusFromMin"             ,  &deltaMWPlusFromMin             ) ;
  outTree->Branch( "deltaMWMinusFromMin"            ,  &deltaMWMinusFromMin            ) ;

  outTree->Branch( "deltaRGenBottomQuarks"     ,  &deltaRGenBottomQuarks     ) ;
  outTree->Branch( "deltaRSmearedBottomQuarks" ,  &deltaRSmearedBottomQuarks ) ;
  outTree->Branch( "deltaRBestBottomQuarks"    ,  &deltaRBestBottomQuarks    ) ;
  outTree->Branch( "deltaRGenLightPartons"     ,  &deltaRGenLightPartons     ) ;
  outTree->Branch( "deltaRSmearedLightPartons" ,  &deltaRSmearedLightPartons ) ;
  outTree->Branch( "deltaRBestLightPartons"    ,  &deltaRBestLightPartons    ) ;
  outTree->Branch( "deltaRGenTopQuarks"        ,  &deltaRGenTopQuarks        ) ;
  outTree->Branch( "deltaRSmearedTopQuarks"    ,  &deltaRSmearedTopQuarks    ) ;
  outTree->Branch( "deltaRBestTopQuarks1"      ,  &deltaRBestTopQuarks1      ) ;
  outTree->Branch( "deltaRBestTopQuarks2"      ,  &deltaRBestTopQuarks2      ) ;
  outTree->Branch( "deltaRBestTopQuarks3"      ,  &deltaRBestTopQuarks3      ) ;
  outTree->Branch( "deltaRBestTopQuarks4"      ,  &deltaRBestTopQuarks4      ) ;
  outTree->Branch( "deltaRBestTopQuarks5"      ,  &deltaRBestTopQuarks5      ) ;
  outTree->Branch( "deltaRBestTopQuarks6"      ,  &deltaRBestTopQuarks6      ) ;
  outTree->Branch( "deltaRBestTopQuarks7"      ,  &deltaRBestTopQuarks7      ) ;
  outTree->Branch( "deltaRBestTopQuarks8"      ,  &deltaRBestTopQuarks8      ) ;
  outTree->Branch( "deltaRGenWBosons"          ,  &deltaRGenWBosons          ) ;
  outTree->Branch( "deltaRSmearedWBosons"      ,  &deltaRSmearedWBosons      ) ;
  outTree->Branch( "deltaRBestWBosons1"        ,  &deltaRBestWBosons1        ) ;
  outTree->Branch( "deltaRBestWBosons2"        ,  &deltaRBestWBosons2        ) ;
  outTree->Branch( "deltaRBestWBosons3"        ,  &deltaRBestWBosons3        ) ;
  outTree->Branch( "deltaRBestWBosons4"        ,  &deltaRBestWBosons4        ) ;
  outTree->Branch( "deltaRBestWBosons5"        ,  &deltaRBestWBosons5        ) ;
  outTree->Branch( "deltaRBestWBosons6"        ,  &deltaRBestWBosons6        ) ;
  outTree->Branch( "deltaRBestWBosons7"        ,  &deltaRBestWBosons7        ) ;
  outTree->Branch( "deltaRBestWBosons8"        ,  &deltaRBestWBosons8        ) ;
  outTree->Branch( "deltaRGenNeutrinos"        ,  &deltaRGenNeutrinos        ) ;
  outTree->Branch( "deltaRBestNeutrinos1"      ,  &deltaRBestNeutrinos1      ) ;
  outTree->Branch( "deltaRBestNeutrinos2"      ,  &deltaRBestNeutrinos2      ) ;
  outTree->Branch( "deltaRBestNeutrinos3"      ,  &deltaRBestNeutrinos3      ) ;
  outTree->Branch( "deltaRBestNeutrinos4"      ,  &deltaRBestNeutrinos4      ) ;
  outTree->Branch( "deltaRBestNeutrinos5"      ,  &deltaRBestNeutrinos5      ) ;
  outTree->Branch( "deltaRBestNeutrinos6"      ,  &deltaRBestNeutrinos6      ) ;
  outTree->Branch( "deltaRBestNeutrinos7"      ,  &deltaRBestNeutrinos7      ) ;
  outTree->Branch( "deltaRBestNeutrinos8"      ,  &deltaRBestNeutrinos8      ) ;

}

Int_t ttbarReconstructionFromLHE::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ttbarReconstructionFromLHE::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ttbarReconstructionFromLHE::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PID = 0;
   P_X = 0;
   P_Y = 0;
   P_Z = 0;
   E = 0;
   M = 0;
   status = 0;
   particleID = 0;
   parent1ID = 0;
   parent2ID = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n_particles", &n_particles, &b_n_particles);
   fChain->SetBranchAddress("PID", &PID, &b_PID);
   fChain->SetBranchAddress("P_X", &P_X, &b_P_X);
   fChain->SetBranchAddress("P_Y", &P_Y, &b_P_Y);
   fChain->SetBranchAddress("P_Z", &P_Z, &b_P_Z);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("M", &M, &b_M);
   fChain->SetBranchAddress("status", &status, &b_status);
   fChain->SetBranchAddress("particleID", &particleID, &b_particleID);
   fChain->SetBranchAddress("parent1ID", &parent1ID, &b_parent1ID);
   fChain->SetBranchAddress("parent2ID", &parent2ID, &b_parent2ID);
   Notify();
}

Bool_t ttbarReconstructionFromLHE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ttbarReconstructionFromLHE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ttbarReconstructionFromLHE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ttbarReconstructionFromLHE_cxx
