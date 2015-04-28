//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  8 15:44:02 2014 by ROOT version 5.34/04
// from TTree Physics/Physics
// found on file: lhe.root
//////////////////////////////////////////////////////////

#ifndef topReconstructionFromLHE_h
#define topReconstructionFromLHE_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "Math/GenVector/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class topReconstructionFromLHE {
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

   topReconstructionFromLHE(TTree *tree=0);
   virtual ~topReconstructionFromLHE();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString, int, int);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void printVector(XYZTLorentzVector&);

   double deltaR(XYZTLorentzVector& , XYZTLorentzVector&  );
   double deltaPhi(XYZTLorentzVector& , XYZTLorentzVector&  );

   //store minimizer results in a tree;
   TFile* outFile;
   TTree* outTree;
   void initOutput(TString, int);


   //Branches

   int eventNumber;

   int innerMinStatus;
   int outerMinStatus;

   double totalChi2;
   double topSystemChi2;
   double topMassChi2;
   double hadronicChi2;
   double nonTopChi2;

   double outerMinEdm;


   //deltas from the minimization
   double bJet1PtDelta, bJet1PhiDelta, bJet1EtaDelta;
   double bJet2PtDelta, bJet2PhiDelta, bJet2EtaDelta;

   double W1Daughter1PtDelta, W1Daughter1PhiDelta, W1Daughter1EtaDelta;
   double W2Daughter1PtDelta, W2Daughter1PhiDelta, W2Daughter1EtaDelta;

   double W1Daughter2PtDelta, W1Daughter2PhiDelta, W1Daughter2EtaDelta;
   double W2Daughter2PtDelta, W2Daughter2PhiDelta, W2Daughter2EtaDelta;

   double top1MassDelta, W1MassDelta;
   double top2MassDelta, W2MassDelta;

   double lightJet1PxDelta, lightJet1PyDelta;
   double lightJet2PxDelta, lightJet2PyDelta;


   //double METx_, METy_, recoMETx_, recoMETy_;


   //compare true and smeared
   double top1DeltaPtTrueSmeared        , top1DeltaRTrueSmeared         , top1DeltaMTrueSmeared        ;
   double W1DeltaPtTrueSmeared          , W1DeltaRTrueSmeared           , W1DeltaMTrueSmeared          ;
   double bJet1DeltaPtTrueSmeared       , bJet1DeltaRTrueSmeared        , bJet1DeltaMTrueSmeared       ;
   double W1Daughter1DeltaPtTrueSmeared , W1Daughter1DeltaRTrueSmeared  , W1Daughter1DeltaMTrueSmeared ;
   double W1Daughter2DeltaPtTrueSmeared , W1Daughter2DeltaRTrueSmeared  , W1Daughter2DeltaMTrueSmeared ;
   double lightJet1DeltaPxTrueSmeared   , lightJet1DeltaPyTrueSmeared   , lightJet1DeltaRTrueSmeared   ;


   double top2DeltaPtTrueSmeared        , top2DeltaRTrueSmeared         , top2DeltaMTrueSmeared        ;
   double W2DeltaPtTrueSmeared          , W2DeltaRTrueSmeared           , W2DeltaMTrueSmeared          ;
   double bJet2DeltaPtTrueSmeared       , bJet2DeltaRTrueSmeared        , bJet2DeltaMTrueSmeared       ;
   double W2Daughter1DeltaPtTrueSmeared , W2Daughter1DeltaRTrueSmeared  , W2Daughter1DeltaMTrueSmeared ;
   double W2Daughter2DeltaPtTrueSmeared , W2Daughter2DeltaRTrueSmeared  , W2Daughter2DeltaMTrueSmeared ;
   double lightJet2DeltaPxTrueSmeared   , lightJet2DeltaPyTrueSmeared   , lightJet2DeltaRTrueSmeared   ;


   //compare true and corrected
   double top1DeltaPtTrueBest        , top1DeltaRTrueBest         , top1DeltaMTrueBest        ;
   double W1DeltaPtTrueBest          , W1DeltaRTrueBest           , W1DeltaMTrueBest          ;
   double bJet1DeltaPtTrueBest       , bJet1DeltaRTrueBest        , bJet1DeltaMTrueBest       ;
   double W1Daughter1DeltaPtTrueBest , W1Daughter1DeltaRTrueBest  , W1Daughter1DeltaMTrueBest ;
   double W1Daughter2DeltaPtTrueBest , W1Daughter2DeltaRTrueBest  , W1Daughter2DeltaMTrueBest ;
   double lightJet1DeltaPxTrueBest   , lightJet1DeltaPyTrueBest   , lightJet1DeltaRTrueBest   ;


   double top2DeltaPtTrueBest        , top2DeltaRTrueBest         , top2DeltaMTrueBest        ;
   double W2DeltaPtTrueBest          , W2DeltaRTrueBest           , W2DeltaMTrueBest          ;
   double bJet2DeltaPtTrueBest       , bJet2DeltaRTrueBest        , bJet2DeltaMTrueBest       ;
   double W2Daughter1DeltaPtTrueBest , W2Daughter1DeltaRTrueBest  , W2Daughter1DeltaMTrueBest ;
   double W2Daughter2DeltaPtTrueBest , W2Daughter2DeltaRTrueBest  , W2Daughter2DeltaMTrueBest ;
   double lightJet2DeltaPxTrueBest   , lightJet2DeltaPyTrueBest   , lightJet2DeltaRTrueBest   ;


   //compare smeared and corrected
   double top1DeltaPtSmearedBest        , top1DeltaRSmearedBest         , top1DeltaMSmearedBest        ;
   double W1DeltaPtSmearedBest          , W1DeltaRSmearedBest           , W1DeltaMSmearedBest          ;
   double bJet1DeltaPtSmearedBest       , bJet1DeltaRSmearedBest        , bJet1DeltaMSmearedBest       ;
   double W1Daughter1DeltaPtSmearedBest , W1Daughter1DeltaRSmearedBest  , W1Daughter1DeltaMSmearedBest ;
   double W1Daughter2DeltaPtSmearedBest , W1Daughter2DeltaRSmearedBest  , W1Daughter2DeltaMSmearedBest ;
   double lightJet1DeltaPxSmearedBest   , lightJet1DeltaPySmearedBest   , lightJet1DeltaRSmearedBest   ;


   double top2DeltaPtSmearedBest        , top2DeltaRSmearedBest         , top2DeltaMSmearedBest        ;
   double W2DeltaPtSmearedBest          , W2DeltaRSmearedBest           , W2DeltaMSmearedBest          ;
   double bJet2DeltaPtSmearedBest       , bJet2DeltaRSmearedBest        , bJet2DeltaMSmearedBest       ;
   double W2Daughter1DeltaPtSmearedBest , W2Daughter1DeltaRSmearedBest  , W2Daughter1DeltaMSmearedBest ;
   double W2Daughter2DeltaPtSmearedBest , W2Daughter2DeltaRSmearedBest  , W2Daughter2DeltaMSmearedBest ;
   double lightJet2DeltaPxSmearedBest   , lightJet2DeltaPySmearedBest   , lightJet2DeltaRSmearedBest   ;



   //true momenta
   double top1TruePt        , top1TrueEta       , top1TruePhi       , top1TrueE        ;
   double W1TruePt          , W1TrueEta         , W1TruePhi         , W1TrueE          ;
   double bJet1TruePt       , bJet1TrueEta      , bJet1TruePhi      , bJet1TrueE       ;
   double W1Daughter1TruePt , W1Daughter1TrueEta, W1Daughter1TruePhi, W1Daughter1TrueE ;
   double W1Daughter2TruePt , W1Daughter2TrueEta, W1Daughter2TruePhi, W1Daughter2TrueE ;
   double lightJet1TruePx   , lightJet1TruePy   , lightJet1TruePz   , lightJet1TrueE   ;
   double top1TrueMass      , W1TrueMass        ;
   double lightJet1TruePt   ;

   double top2TruePt        , top2TrueEta       , top2TruePhi       , top2TrueE        ;
   double W2TruePt          , W2TrueEta         , W2TruePhi         , W2TrueE          ;
   double bJet2TruePt       , bJet2TrueEta      , bJet2TruePhi      , bJet2TrueE       ;
   double W2Daughter1TruePt , W2Daughter1TrueEta, W2Daughter1TruePhi, W2Daughter1TrueE ;
   double W2Daughter2TruePt , W2Daughter2TrueEta, W2Daughter2TruePhi, W2Daughter2TrueE ;
   double lightJet2TruePx   , lightJet2TruePy   , lightJet2TruePz   , lightJet2TrueE   ;
   double top2TrueMass      , W2TrueMass        ;
   double lightJet2TruePt   ;


   //smeared momenta
   double top1SmearedPt        , top1SmearedEta       , top1SmearedPhi       , top1SmearedE        ;
   double W1SmearedPt          , W1SmearedEta         , W1SmearedPhi         , W1SmearedE          ;
   double bJet1SmearedPt       , bJet1SmearedEta      , bJet1SmearedPhi      , bJet1SmearedE       ;
   double W1Daughter1SmearedPt , W1Daughter1SmearedEta, W1Daughter1SmearedPhi, W1Daughter1SmearedE ;
   double W1Daughter2SmearedPt , W1Daughter2SmearedEta, W1Daughter2SmearedPhi, W1Daughter2SmearedE ;
   double lightJet1SmearedPx   , lightJet1SmearedPy   , lightJet1SmearedPz   , lightJet1SmearedE   ;
   double top1SmearedMass      , W1SmearedMass        ;
   double lightJet1SmearedPt   ;

   double top2SmearedPt        , top2SmearedEta       , top2SmearedPhi       , top2SmearedE        ;
   double W2SmearedPt          , W2SmearedEta         , W2SmearedPhi         , W2SmearedE          ;
   double bJet2SmearedPt       , bJet2SmearedEta      , bJet2SmearedPhi      , bJet2SmearedE       ;
   double W2Daughter1SmearedPt , W2Daughter1SmearedEta, W2Daughter1SmearedPhi, W2Daughter1SmearedE ;
   double W2Daughter2SmearedPt , W2Daughter2SmearedEta, W2Daughter2SmearedPhi, W2Daughter2SmearedE ;
   double lightJet2SmearedPx   , lightJet2SmearedPy   , lightJet2SmearedPz   , lightJet2SmearedE   ;
   double top2SmearedMass      , W2SmearedMass        ;
   double lightJet2SmearedPt   ;



   //best momenta
   double top1BestPx        , top1BestPy       , top1BestPz       ;
   double W1BestPx          , W1BestPy         , W1BestPz         ;
   double bJet1BestPx       , bJet1BestPy      , bJet1BestPz      ;
   double W1Daughter1BestPx , W1Daughter1BestPy, W1Daughter1BestPz;
   double W1Daughter2BestPx , W1Daughter2BestPy, W1Daughter2BestPz;

   double top2BestPx        , top2BestPy       , top2BestPz       ;
   double W2BestPx          , W2BestPy         , W2BestPz         ;
   double bJet2BestPx       , bJet2BestPy      , bJet2BestPz      ;
   double W2Daughter1BestPx , W2Daughter1BestPy, W2Daughter1BestPz;
   double W2Daughter2BestPx , W2Daughter2BestPy, W2Daughter2BestPz;

   double top1BestPt        , top1BestEta       , top1BestPhi       , top1BestE        ;
   double W1BestPt          , W1BestEta         , W1BestPhi         , W1BestE          ;
   double bJet1BestPt       , bJet1BestEta      , bJet1BestPhi      , bJet1BestE       ;
   double W1Daughter1BestPt , W1Daughter1BestEta, W1Daughter1BestPhi, W1Daughter1BestE ;
   double W1Daughter2BestPt , W1Daughter2BestEta, W1Daughter2BestPhi, W1Daughter2BestE ;
   double lightJet1BestPx   , lightJet1BestPy   , lightJet1BestPz   , lightJet1BestE   ;
   double top1BestMass      , W1BestMass        ;
   double lightJet1BestPt   ;

   double top2BestPt        , top2BestEta       , top2BestPhi       , top2BestE        ;
   double W2BestPt          , W2BestEta         , W2BestPhi         , W2BestE          ;
   double bJet2BestPt       , bJet2BestEta      , bJet2BestPhi      , bJet2BestE       ;
   double W2Daughter1BestPt , W2Daughter1BestEta, W2Daughter1BestPhi, W2Daughter1BestE ;
   double W2Daughter2BestPt , W2Daughter2BestEta, W2Daughter2BestPhi, W2Daughter2BestE ;
   double lightJet2BestPx   , lightJet2BestPy   , lightJet2BestPz   , lightJet2BestE   ;
   double top2BestMass      , W2BestMass        ;
   double lightJet2BestPt   ;


   bool smearingSwitchedLightJetOrdering;

};

#endif

#ifdef topReconstructionFromLHE_cxx
topReconstructionFromLHE::topReconstructionFromLHE(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lhe.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("lhe.root");
      }
      f->GetObject("Physics",tree);

   }
   Init(tree);
}

topReconstructionFromLHE::~topReconstructionFromLHE()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

double topReconstructionFromLHE::deltaR(XYZTLorentzVector& p2, XYZTLorentzVector& p1 )
{
  double thisDeltaEta = p1.Eta()-p2.Eta();
  double thisDeltaPhi = p1.Phi()-p2.Phi();
  while(thisDeltaPhi> 3.14159265359)
    {
      thisDeltaPhi -= 2*3.14159265359;
    }
  while(thisDeltaPhi< -3.14159265359)
    {
      thisDeltaPhi += 2*3.14159265359;
    }
  return sqrt(thisDeltaPhi*thisDeltaPhi+thisDeltaEta*thisDeltaEta);
}

double topReconstructionFromLHE::deltaPhi(XYZTLorentzVector& p2, XYZTLorentzVector& p1 )
{
  double thisDeltaPhi = p1.Phi()-p2.Phi();
  while(thisDeltaPhi> 3.14159265359)
    {
      thisDeltaPhi -= 2*3.14159265359;
    }
  while(thisDeltaPhi< -3.14159265359)
    {
      thisDeltaPhi += 2*3.14159265359;
    }
  return thisDeltaPhi;
}

void topReconstructionFromLHE::initOutput(TString dir, int whichLoop)
{
  TString fileName=dir; 
  if(dir!="") fileName+="/";
  fileName+="output_";
  fileName+=whichLoop;
  fileName+=".root";
  outFile = new TFile(fileName,"RECREATE","tree");
  outTree = new TTree("tree","tree");

  outTree->Branch( "eventNumber"  ,  &eventNumber  ) ;

  outTree->Branch( "innerMinStatus"  ,  &innerMinStatus  ) ;
  outTree->Branch( "outerMinStatus"  ,  &outerMinStatus  ) ;
  outTree->Branch( "outerMinEdm"     ,  &outerMinEdm     );

  outTree->Branch( "totalChi2"     ,  &totalChi2     ) ;
  outTree->Branch( "topSystemChi2" ,  &topSystemChi2 ) ;
  outTree->Branch( "topMassChi2"   ,  &topMassChi2   ) ;
  outTree->Branch( "hadronicChi2"  ,  &hadronicChi2  );
  outTree->Branch( "nonTopChi2"    ,  &nonTopChi2    ) ;

  outTree->Branch( "bJet1PtDelta"   ,  &bJet1PtDelta  ) ;
  outTree->Branch( "bJet1PhiDelta"  ,  &bJet1PhiDelta ) ;
  outTree->Branch( "bJet1EtaDelta"  ,  &bJet1EtaDelta ) ;
  outTree->Branch( "bJet2PtDelta"   ,  &bJet2PtDelta  ) ;
  outTree->Branch( "bJet2PhiDelta"  ,  &bJet2PhiDelta ) ;
  outTree->Branch( "bJet2EtaDelta"  ,  &bJet2EtaDelta ) ;

  outTree->Branch( "W1Daughter1PtDelta"   ,  &W1Daughter1PtDelta  ) ;
  outTree->Branch( "W1Daughter1PhiDelta"  ,  &W1Daughter1PhiDelta ) ;
  outTree->Branch( "W1Daughter1EtaDelta"  ,  &W1Daughter1EtaDelta ) ;
  outTree->Branch( "W2Daughter1PtDelta"   ,  &W2Daughter1PtDelta  ) ;
  outTree->Branch( "W2Daughter1PhiDelta"  ,  &W2Daughter1PhiDelta ) ;
  outTree->Branch( "W2Daughter1EtaDelta"  ,  &W2Daughter1EtaDelta ) ;

  outTree->Branch( "W1Daughter2PtDelta"   ,  &W1Daughter2PtDelta  ) ;
  outTree->Branch( "W1Daughter2PhiDelta"  ,  &W1Daughter2PhiDelta ) ;
  outTree->Branch( "W1Daughter2EtaDelta"  ,  &W1Daughter2EtaDelta ) ;
  outTree->Branch( "W2Daughter2PtDelta"   ,  &W2Daughter2PtDelta  ) ;
  outTree->Branch( "W2Daughter2PhiDelta"  ,  &W2Daughter2PhiDelta ) ;
  outTree->Branch( "W2Daughter2EtaDelta"  ,  &W2Daughter2EtaDelta ) ;

  outTree->Branch( "top1MassDelta"  ,  &top1MassDelta  ) ;
  outTree->Branch( "W1MassDelta"    ,  &W1MassDelta    ) ;
  outTree->Branch( "top2MassDelta"  ,  &top2MassDelta  ) ;
  outTree->Branch( "W2MassDelta"    ,  &W2MassDelta    ) ;

  outTree->Branch( "lightJet1PxDelta"  ,  &lightJet1PxDelta  ) ;
  outTree->Branch( "lightJet1PyDelta"  ,  &lightJet1PyDelta  ) ;
  outTree->Branch( "lightJet2PxDelta"  ,  &lightJet2PxDelta  ) ;
  outTree->Branch( "lightJet2PyDelta"  ,  &lightJet2PyDelta  ) ;

  outTree->Branch( "top1TruePt"          ,  &top1TruePt         ) ;  
  outTree->Branch( "top1TrueEta"         ,  &top1TrueEta        ) ; 
  outTree->Branch( "top1TruePhi"         ,  &top1TruePhi        ) ; 
  outTree->Branch( "top1TrueE"           ,  &top1TrueE          ) ;
  outTree->Branch( "top1TrueMass"        ,  &top1TrueMass       ) ; 
  outTree->Branch( "W1TruePt"            ,  &W1TruePt           ) ;
  outTree->Branch( "W1TrueEta"           ,  &W1TrueEta          ) ;
  outTree->Branch( "W1TruePhi"           ,  &W1TruePhi          ) ;
  outTree->Branch( "W1TrueE"             ,  &W1TrueE            ) ;
  outTree->Branch( "W1TrueMass"          ,  &W1TrueMass         ) ;
  outTree->Branch( "bJet1TruePt"         ,  &bJet1TruePt        ) ;
  outTree->Branch( "bJet1TrueEta"        ,  &bJet1TrueEta       ) ;
  outTree->Branch( "bJet1TruePhi"        ,  &bJet1TruePhi       ) ;
  outTree->Branch( "bJet1TrueE"          ,  &bJet1TrueE         ) ;
  outTree->Branch( "W1Daughter1TruePt"   ,  &W1Daughter1TruePt  ) ;
  outTree->Branch( "W1Daughter1TrueEta"  ,  &W1Daughter1TrueEta ) ;
  outTree->Branch( "W1Daughter1TruePhi"  ,  &W1Daughter1TruePhi ) ;
  outTree->Branch( "W1Daughter1TrueE"    ,  &W1Daughter1TrueE   ) ;
  outTree->Branch( "W1Daughter2TruePt"   ,  &W1Daughter2TruePt  ) ; 
  outTree->Branch( "W1Daughter2TrueEta"  ,  &W1Daughter2TrueEta ) ;
  outTree->Branch( "W1Daughter2TruePhi"  ,  &W1Daughter2TruePhi ) ;
  outTree->Branch( "W1Daughter2TrueE"    ,  &W1Daughter2TrueE   ) ;
  outTree->Branch( "lightJet1TruePx"     ,  &lightJet1TruePx    ) ;
  outTree->Branch( "lightJet1TruePy"     ,  &lightJet1TruePy    ) ;
  outTree->Branch( "lightJet1TruePz"     ,  &lightJet1TruePz    ) ;
  outTree->Branch( "lightJet1TrueE"      ,  &lightJet1TrueE     ) ;
  outTree->Branch( "lightJet1TruePt"     ,  &lightJet1TruePt    ) ;

  outTree->Branch( "top2TruePt"          ,  &top2TruePt         ) ;
  outTree->Branch( "top2TrueEta"         ,  &top2TrueEta        ) ;
  outTree->Branch( "top2TruePhi"         ,  &top2TruePhi        ) ;
  outTree->Branch( "top2TrueE"           ,  &top2TrueE          ) ;
  outTree->Branch( "top2TrueMass"        ,  &top2TrueMass       ) ;
  outTree->Branch( "W2TruePt"            ,  &W2TruePt           ) ;
  outTree->Branch( "W2TrueEta"           ,  &W2TrueEta          ) ;
  outTree->Branch( "W2TruePhi"           ,  &W2TruePhi          ) ;
  outTree->Branch( "W2TrueE"             ,  &W2TrueE            ) ;
  outTree->Branch( "W2TrueMass"          ,  &W2TrueMass         ) ;
  outTree->Branch( "bJet2TruePt"         ,  &bJet2TruePt        ) ;
  outTree->Branch( "bJet2TrueEta"        ,  &bJet2TrueEta       ) ;
  outTree->Branch( "bJet2TruePhi"        ,  &bJet2TruePhi       ) ;
  outTree->Branch( "bJet2TrueE"          ,  &bJet2TrueE         ) ;
  outTree->Branch( "W2Daughter1TruePt"   ,  &W2Daughter1TruePt  ) ;
  outTree->Branch( "W2Daughter1TrueEta"  ,  &W2Daughter1TrueEta ) ;
  outTree->Branch( "W2Daughter1TruePhi"  ,  &W2Daughter1TruePhi ) ;
  outTree->Branch( "W2Daughter1TrueE"    ,  &W2Daughter1TrueE   ) ;
  outTree->Branch( "W2Daughter2TruePt"   ,  &W2Daughter2TruePt  ) ;
  outTree->Branch( "W2Daughter2TrueEta"  ,  &W2Daughter2TrueEta ) ;
  outTree->Branch( "W2Daughter2TruePhi"  ,  &W2Daughter2TruePhi ) ;
  outTree->Branch( "W2Daughter2TrueE"    ,  &W2Daughter2TrueE   ) ;
  outTree->Branch( "lightJet2TruePx"     ,  &lightJet2TruePx    ) ;
  outTree->Branch( "lightJet2TruePy"     ,  &lightJet2TruePy    ) ;
  outTree->Branch( "lightJet2TruePz"     ,  &lightJet2TruePz    ) ;
  outTree->Branch( "lightJet2TrueE"      ,  &lightJet2TrueE     ) ;
  outTree->Branch( "lightJet2TruePt"     ,  &lightJet2TruePt    ) ;


  outTree->Branch( "top1SmearedPt"          ,  &top1SmearedPt         ) ;
  outTree->Branch( "top1SmearedEta"         ,  &top1SmearedEta        ) ;
  outTree->Branch( "top1SmearedPhi"         ,  &top1SmearedPhi        ) ;
  outTree->Branch( "top1SmearedE"           ,  &top1SmearedE          ) ;
  outTree->Branch( "top1SmearedMass"        ,  &top1SmearedMass       ) ;
  outTree->Branch( "W1SmearedPt"            ,  &W1SmearedPt           ) ;
  outTree->Branch( "W1SmearedEta"           ,  &W1SmearedEta          ) ;
  outTree->Branch( "W1SmearedPhi"           ,  &W1SmearedPhi          ) ;
  outTree->Branch( "W1SmearedE"             ,  &W1SmearedE            ) ;
  outTree->Branch( "W1SmearedMass"          ,  &W1SmearedMass         ) ;
  outTree->Branch( "bJet1SmearedPt"         ,  &bJet1SmearedPt        ) ;
  outTree->Branch( "bJet1SmearedEta"        ,  &bJet1SmearedEta       ) ;
  outTree->Branch( "bJet1SmearedPhi"        ,  &bJet1SmearedPhi       ) ;
  outTree->Branch( "bJet1SmearedE"          ,  &bJet1SmearedE         ) ;
  outTree->Branch( "W1Daughter1SmearedPt"   ,  &W1Daughter1SmearedPt  ) ;
  outTree->Branch( "W1Daughter1SmearedEta"  ,  &W1Daughter1SmearedEta ) ;
  outTree->Branch( "W1Daughter1SmearedPhi"  ,  &W1Daughter1SmearedPhi ) ;
  outTree->Branch( "W1Daughter1SmearedE"    ,  &W1Daughter1SmearedE   ) ;
  outTree->Branch( "W1Daughter2SmearedPt"   ,  &W1Daughter2SmearedPt  ) ;
  outTree->Branch( "W1Daughter2SmearedEta"  ,  &W1Daughter2SmearedEta ) ;
  outTree->Branch( "W1Daughter2SmearedPhi"  ,  &W1Daughter2SmearedPhi ) ;
  outTree->Branch( "W1Daughter2SmearedE"    ,  &W1Daughter2SmearedE   ) ;
  outTree->Branch( "lightJet1SmearedPx"     ,  &lightJet1SmearedPx    ) ;
  outTree->Branch( "lightJet1SmearedPy"     ,  &lightJet1SmearedPy    ) ;
  outTree->Branch( "lightJet1SmearedPz"     ,  &lightJet1SmearedPz    ) ;
  outTree->Branch( "lightJet1SmearedE"      ,  &lightJet1SmearedE     ) ;
  outTree->Branch( "lightJet1SmearedPt"     ,  &lightJet1SmearedPt    ) ;

  outTree->Branch( "top2SmearedPt"          ,  &top2SmearedPt         ) ;
  outTree->Branch( "top2SmearedEta"         ,  &top2SmearedEta        ) ;
  outTree->Branch( "top2SmearedPhi"         ,  &top2SmearedPhi        ) ;
  outTree->Branch( "top2SmearedE"           ,  &top2SmearedE          ) ;
  outTree->Branch( "top2SmearedMass"        ,  &top2SmearedMass       ) ;
  outTree->Branch( "W2SmearedPt"            ,  &W2SmearedPt           ) ;
  outTree->Branch( "W2SmearedEta"           ,  &W2SmearedEta          ) ;
  outTree->Branch( "W2SmearedPhi"           ,  &W2SmearedPhi          ) ;
  outTree->Branch( "W2SmearedE"             ,  &W2SmearedE            ) ;
  outTree->Branch( "W2SmearedMass"          ,  &W2SmearedMass         ) ;
  outTree->Branch( "bJet2SmearedPt"         ,  &bJet2SmearedPt        ) ;
  outTree->Branch( "bJet2SmearedEta"        ,  &bJet2SmearedEta       ) ;
  outTree->Branch( "bJet2SmearedPhi"        ,  &bJet2SmearedPhi       ) ;
  outTree->Branch( "bJet2SmearedE"          ,  &bJet2SmearedE         ) ;
  outTree->Branch( "W2Daughter1SmearedPt"   ,  &W2Daughter1SmearedPt  ) ;
  outTree->Branch( "W2Daughter1SmearedEta"  ,  &W2Daughter1SmearedEta ) ;
  outTree->Branch( "W2Daughter1SmearedPhi"  ,  &W2Daughter1SmearedPhi ) ;
  outTree->Branch( "W2Daughter1SmearedE"    ,  &W2Daughter1SmearedE   ) ;
  outTree->Branch( "W2Daughter2SmearedPt"   ,  &W2Daughter2SmearedPt  ) ;
  outTree->Branch( "W2Daughter2SmearedEta"  ,  &W2Daughter2SmearedEta ) ;
  outTree->Branch( "W2Daughter2SmearedPhi"  ,  &W2Daughter2SmearedPhi ) ;
  outTree->Branch( "W2Daughter2SmearedE"    ,  &W2Daughter2SmearedE   ) ;
  outTree->Branch( "lightJet2SmearedPx"     ,  &lightJet2SmearedPx    ) ;
  outTree->Branch( "lightJet2SmearedPy"     ,  &lightJet2SmearedPy    ) ;
  outTree->Branch( "lightJet2SmearedPz"     ,  &lightJet2SmearedPz    ) ;
  outTree->Branch( "lightJet2SmearedE"      ,  &lightJet2SmearedE     ) ;
  outTree->Branch( "lightJet2SmearedPt"     ,  &lightJet2SmearedPt    ) ;


  outTree->Branch( "top1BestPt"          ,  &top1BestPt         ) ;
  outTree->Branch( "top1BestEta"         ,  &top1BestEta        ) ;
  outTree->Branch( "top1BestPhi"         ,  &top1BestPhi        ) ;
  outTree->Branch( "top1BestE"           ,  &top1BestE          ) ;
  outTree->Branch( "top1BestMass"        ,  &top1BestMass       ) ;
  outTree->Branch( "W1BestPt"            ,  &W1BestPt           ) ;
  outTree->Branch( "W1BestEta"           ,  &W1BestEta          ) ;
  outTree->Branch( "W1BestPhi"           ,  &W1BestPhi          ) ;
  outTree->Branch( "W1BestE"             ,  &W1BestE            ) ;
  outTree->Branch( "W1BestMass"          ,  &W1BestMass         ) ;
  outTree->Branch( "bJet1BestPt"         ,  &bJet1BestPt        ) ;
  outTree->Branch( "bJet1BestEta"        ,  &bJet1BestEta       ) ;
  outTree->Branch( "bJet1BestPhi"        ,  &bJet1BestPhi       ) ;
  outTree->Branch( "bJet1BestE"          ,  &bJet1BestE         ) ;
  outTree->Branch( "W1Daughter1BestPt"   ,  &W1Daughter1BestPt  ) ;
  outTree->Branch( "W1Daughter1BestEta"  ,  &W1Daughter1BestEta ) ;
  outTree->Branch( "W1Daughter1BestPhi"  ,  &W1Daughter1BestPhi ) ;
  outTree->Branch( "W1Daughter1BestE"    ,  &W1Daughter1BestE   ) ;
  outTree->Branch( "W1Daughter2BestPt"   ,  &W1Daughter2BestPt  ) ;
  outTree->Branch( "W1Daughter2BestEta"  ,  &W1Daughter2BestEta ) ;
  outTree->Branch( "W1Daughter2BestPhi"  ,  &W1Daughter2BestPhi ) ;
  outTree->Branch( "W1Daughter2BestE"    ,  &W1Daughter2BestE   ) ;
  outTree->Branch( "lightJet1BestPx"     ,  &lightJet1BestPx    ) ;
  outTree->Branch( "lightJet1BestPy"     ,  &lightJet1BestPy    ) ;
  outTree->Branch( "lightJet1BestPz"     ,  &lightJet1BestPz    ) ;
  outTree->Branch( "lightJet1BestE"      ,  &lightJet1BestE     ) ;
  outTree->Branch( "lightJet1BestPt"     ,  &lightJet1BestPt    ) ;

  outTree->Branch( "top2BestPt"          ,  &top2BestPt         ) ;
  outTree->Branch( "top2BestEta"         ,  &top2BestEta        ) ;
  outTree->Branch( "top2BestPhi"         ,  &top2BestPhi        ) ;
  outTree->Branch( "top2BestE"           ,  &top2BestE          ) ;
  outTree->Branch( "top2BestMass"        ,  &top2BestMass       ) ;
  outTree->Branch( "W2BestPt"            ,  &W2BestPt           ) ;
  outTree->Branch( "W2BestEta"           ,  &W2BestEta          ) ;
  outTree->Branch( "W2BestPhi"           ,  &W2BestPhi          ) ;
  outTree->Branch( "W2BestE"             ,  &W2BestE            ) ;
  outTree->Branch( "W2BestMass"          ,  &W2BestMass         ) ;
  outTree->Branch( "bJet2BestPt"         ,  &bJet2BestPt        ) ;
  outTree->Branch( "bJet2BestEta"        ,  &bJet2BestEta       ) ;
  outTree->Branch( "bJet2BestPhi"        ,  &bJet2BestPhi       ) ;
  outTree->Branch( "bJet2BestE"          ,  &bJet2BestE         ) ;
  outTree->Branch( "W2Daughter1BestPt"   ,  &W2Daughter1BestPt  ) ;
  outTree->Branch( "W2Daughter1BestEta"  ,  &W2Daughter1BestEta ) ;
  outTree->Branch( "W2Daughter1BestPhi"  ,  &W2Daughter1BestPhi ) ;
  outTree->Branch( "W2Daughter1BestE"    ,  &W2Daughter1BestE   ) ;
  outTree->Branch( "W2Daughter2BestPt"   ,  &W2Daughter2BestPt  ) ;
  outTree->Branch( "W2Daughter2BestEta"  ,  &W2Daughter2BestEta ) ;
  outTree->Branch( "W2Daughter2BestPhi"  ,  &W2Daughter2BestPhi ) ;
  outTree->Branch( "W2Daughter2BestE"    ,  &W2Daughter2BestE   ) ;
  outTree->Branch( "lightJet2BestPx"     ,  &lightJet2BestPx    ) ;
  outTree->Branch( "lightJet2BestPy"     ,  &lightJet2BestPy    ) ;
  outTree->Branch( "lightJet2BestPz"     ,  &lightJet2BestPz    ) ;
  outTree->Branch( "lightJet2BestE"      ,  &lightJet2BestE     ) ;
  outTree->Branch( "lightJet2BestPt"     ,  &lightJet2BestPt    ) ;



  outTree->Branch( "top1DeltaPtTrueSmeared"         ,  &top1DeltaPtTrueSmeared        ) ;
  outTree->Branch( "top1DeltaRTrueSmeared"          ,  &top1DeltaRTrueSmeared         ) ;
  outTree->Branch( "top1DeltaMTrueSmeared"          ,  &top1DeltaMTrueSmeared         ) ;
  outTree->Branch( "W1DeltaPtTrueSmeared"           ,  &W1DeltaPtTrueSmeared          ) ;
  outTree->Branch( "W1DeltaRTrueSmeared"            ,  &W1DeltaRTrueSmeared           ) ;
  outTree->Branch( "W1DeltaMTrueSmeared"            ,  &W1DeltaMTrueSmeared           ) ;
  outTree->Branch( "bJet1DeltaPtTrueSmeared"        ,  &bJet1DeltaPtTrueSmeared       ) ;
  outTree->Branch( "bJet1DeltaRTrueSmeared"         ,  &bJet1DeltaRTrueSmeared        ) ;
  //outTree->Branch( "bJet1DeltaMTrueSmeared"         ,  &bJet1DeltaMTrueSmeared        ) ;
  outTree->Branch( "W1Daughter1DeltaPtTrueSmeared"  ,  &W1Daughter1DeltaPtTrueSmeared ) ;
  outTree->Branch( "W1Daughter1DeltaRTrueSmeared"   ,  &W1Daughter1DeltaRTrueSmeared  ) ;
  //outTree->Branch( "W1Daughter1DeltaMTrueSmeared"   ,  &W1Daughter1DeltaMTrueSmeared  ) ;
  outTree->Branch( "W1Daughter2DeltaPtTrueSmeared"  ,  &W1Daughter2DeltaPtTrueSmeared ) ;
  outTree->Branch( "W1Daughter2DeltaRTrueSmeared"   ,  &W1Daughter2DeltaRTrueSmeared  ) ;
  //outTree->Branch( "W1Daughter2DeltaMTrueSmeared"   ,  &W1Daughter2DeltaMTrueSmeared  ) ;
  outTree->Branch( "lightJet1DeltaPxTrueSmeared"    ,  &lightJet1DeltaPxTrueSmeared   ) ;
  outTree->Branch( "lightJet1DeltaPyTrueSmeared"    ,  &lightJet1DeltaPyTrueSmeared   ) ;
  outTree->Branch( "lightJet1DeltaRTrueSmeared"     ,  &lightJet1DeltaRTrueSmeared    ) ;

  outTree->Branch( "top2DeltaPtTrueSmeared"         ,  &top2DeltaPtTrueSmeared        ) ;
  outTree->Branch( "top2DeltaRTrueSmeared"          ,  &top2DeltaRTrueSmeared         ) ;
  outTree->Branch( "top2DeltaMTrueSmeared"          ,  &top2DeltaMTrueSmeared         ) ;
  outTree->Branch( "W2DeltaPtTrueSmeared"           ,  &W2DeltaPtTrueSmeared          ) ;
  outTree->Branch( "W2DeltaRTrueSmeared"            ,  &W2DeltaRTrueSmeared           ) ;
  outTree->Branch( "W2DeltaMTrueSmeared"            ,  &W2DeltaMTrueSmeared           ) ;
  outTree->Branch( "bJet2DeltaPtTrueSmeared"        ,  &bJet2DeltaPtTrueSmeared       ) ;
  outTree->Branch( "bJet2DeltaRTrueSmeared"         ,  &bJet2DeltaRTrueSmeared        ) ;
  //outTree->Branch( "bJet2DeltaMTrueSmeared"         ,  &bJet2DeltaMTrueSmeared        ) ;
  outTree->Branch( "W2Daughter1DeltaPtTrueSmeared"  ,  &W2Daughter1DeltaPtTrueSmeared ) ;
  outTree->Branch( "W2Daughter1DeltaRTrueSmeared"   ,  &W2Daughter1DeltaRTrueSmeared  ) ;
  //outTree->Branch( "W2Daughter1DeltaMTrueSmeared"   ,  &W2Daughter1DeltaMTrueSmeared  ) ;
  outTree->Branch( "W2Daughter2DeltaPtTrueSmeared"  ,  &W2Daughter2DeltaPtTrueSmeared ) ;
  outTree->Branch( "W2Daughter2DeltaRTrueSmeared"   ,  &W2Daughter2DeltaRTrueSmeared  ) ;
  //outTree->Branch( "W2Daughter2DeltaMTrueSmeared"   ,  &W2Daughter2DeltaMTrueSmeared  ) ;
  outTree->Branch( "lightJet2DeltaPxTrueSmeared"    ,  &lightJet2DeltaPxTrueSmeared   ) ;
  outTree->Branch( "lightJet2DeltaPyTrueSmeared"    ,  &lightJet2DeltaPyTrueSmeared   ) ;
  outTree->Branch( "lightJet2DeltaRTrueSmeared"     ,  &lightJet2DeltaRTrueSmeared    ) ;


  outTree->Branch( "top1DeltaPtTrueBest"         ,  &top1DeltaPtTrueBest        ) ;
  outTree->Branch( "top1DeltaRTrueBest"          ,  &top1DeltaRTrueBest         ) ;
  outTree->Branch( "top1DeltaMTrueBest"          ,  &top1DeltaMTrueBest         ) ;
  outTree->Branch( "W1DeltaPtTrueBest"           ,  &W1DeltaPtTrueBest          ) ;
  outTree->Branch( "W1DeltaRTrueBest"            ,  &W1DeltaRTrueBest           ) ;
  outTree->Branch( "W1DeltaMTrueBest"            ,  &W1DeltaMTrueBest           ) ;
  outTree->Branch( "bJet1DeltaPtTrueBest"        ,  &bJet1DeltaPtTrueBest       ) ;
  outTree->Branch( "bJet1DeltaRTrueBest"         ,  &bJet1DeltaRTrueBest        ) ;
  //outTree->Branch( "bJet1DeltaMTrueBest"         ,  &bJet1DeltaMTrueBest        ) ;
  outTree->Branch( "W1Daughter1DeltaPtTrueBest"  ,  &W1Daughter1DeltaPtTrueBest ) ;
  outTree->Branch( "W1Daughter1DeltaRTrueBest"   ,  &W1Daughter1DeltaRTrueBest  ) ;
  //outTree->Branch( "W1Daughter1DeltaMTrueBest"   ,  &W1Daughter1DeltaMTrueBest  ) ;
  outTree->Branch( "W1Daughter2DeltaPtTrueBest"  ,  &W1Daughter2DeltaPtTrueBest ) ;
  outTree->Branch( "W1Daughter2DeltaRTrueBest"   ,  &W1Daughter2DeltaRTrueBest  ) ;
  //outTree->Branch( "W1Daughter2DeltaMTrueBest"   ,  &W1Daughter2DeltaMTrueBest  ) ;
  outTree->Branch( "lightJet1DeltaPxTrueBest"    ,  &lightJet1DeltaPxTrueBest   ) ;
  outTree->Branch( "lightJet1DeltaPyTrueBest"    ,  &lightJet1DeltaPyTrueBest   ) ;
  outTree->Branch( "lightJet1DeltaRTrueBest"     ,  &lightJet1DeltaRTrueBest    ) ;

  outTree->Branch( "top2DeltaPtTrueBest"         ,  &top2DeltaPtTrueBest        ) ;
  outTree->Branch( "top2DeltaRTrueBest"          ,  &top2DeltaRTrueBest         ) ;
  outTree->Branch( "top2DeltaMTrueBest"          ,  &top2DeltaMTrueBest         ) ;
  outTree->Branch( "W2DeltaPtTrueBest"           ,  &W2DeltaPtTrueBest          ) ;
  outTree->Branch( "W2DeltaRTrueBest"            ,  &W2DeltaRTrueBest           ) ;
  outTree->Branch( "W2DeltaMTrueBest"            ,  &W2DeltaMTrueBest           ) ;
  outTree->Branch( "bJet2DeltaPtTrueBest"        ,  &bJet2DeltaPtTrueBest       ) ;
  outTree->Branch( "bJet2DeltaRTrueBest"         ,  &bJet2DeltaRTrueBest        ) ;
  //outTree->Branch( "bJet2DeltaMTrueBest"         ,  &bJet2DeltaMTrueBest        ) ;
  outTree->Branch( "W2Daughter1DeltaPtTrueBest"  ,  &W2Daughter1DeltaPtTrueBest ) ;
  outTree->Branch( "W2Daughter1DeltaRTrueBest"   ,  &W2Daughter1DeltaRTrueBest  ) ;
  //outTree->Branch( "W2Daughter1DeltaMTrueBest"   ,  &W2Daughter1DeltaMTrueBest  ) ;
  outTree->Branch( "W2Daughter2DeltaPtTrueBest"  ,  &W2Daughter2DeltaPtTrueBest ) ;
  outTree->Branch( "W2Daughter2DeltaRTrueBest"   ,  &W2Daughter2DeltaRTrueBest  ) ;
  //outTree->Branch( "W2Daughter2DeltaMTrueBest"   ,  &W2Daughter2DeltaMTrueBest  ) ;
  outTree->Branch( "lightJet2DeltaPxTrueBest"    ,  &lightJet2DeltaPxTrueBest   ) ;
  outTree->Branch( "lightJet2DeltaPyTrueBest"    ,  &lightJet2DeltaPyTrueBest   ) ;
  outTree->Branch( "lightJet2DeltaRTrueBest"     ,  &lightJet2DeltaRTrueBest    ) ;


  outTree->Branch( "top1DeltaPtSmearedBest"         ,  &top1DeltaPtSmearedBest        ) ;
  outTree->Branch( "top1DeltaRSmearedBest"          ,  &top1DeltaRSmearedBest         ) ;
  outTree->Branch( "top1DeltaMSmearedBest"          ,  &top1DeltaMSmearedBest         ) ;
  outTree->Branch( "W1DeltaPtSmearedBest"           ,  &W1DeltaPtSmearedBest          ) ;
  outTree->Branch( "W1DeltaRSmearedBest"            ,  &W1DeltaRSmearedBest           ) ;
  outTree->Branch( "W1DeltaMSmearedBest"            ,  &W1DeltaMSmearedBest           ) ;
  outTree->Branch( "bJet1DeltaPtSmearedBest"        ,  &bJet1DeltaPtSmearedBest       ) ;
  outTree->Branch( "bJet1DeltaRSmearedBest"         ,  &bJet1DeltaRSmearedBest        ) ;
  outTree->Branch( "bJet1DeltaMSmearedBest"         ,  &bJet1DeltaMSmearedBest        ) ;
  outTree->Branch( "W1Daughter1DeltaPtSmearedBest"  ,  &W1Daughter1DeltaPtSmearedBest ) ;
  outTree->Branch( "W1Daughter1DeltaRSmearedBest"   ,  &W1Daughter1DeltaRSmearedBest  ) ;
  outTree->Branch( "W1Daughter1DeltaMSmearedBest"   ,  &W1Daughter1DeltaMSmearedBest  ) ;
  outTree->Branch( "W1Daughter2DeltaPtSmearedBest"  ,  &W1Daughter2DeltaPtSmearedBest ) ;
  outTree->Branch( "W1Daughter2DeltaRSmearedBest"   ,  &W1Daughter2DeltaRSmearedBest  ) ;
  outTree->Branch( "W1Daughter2DeltaMSmearedBest"   ,  &W1Daughter2DeltaMSmearedBest  ) ;
  outTree->Branch( "lightJet1DeltaPxSmearedBest"    ,  &lightJet1DeltaPxSmearedBest   ) ;
  outTree->Branch( "lightJet1DeltaPySmearedBest"    ,  &lightJet1DeltaPySmearedBest   ) ;
  outTree->Branch( "lightJet1DeltaRSmearedBest"     ,  &lightJet1DeltaRSmearedBest    ) ;

  outTree->Branch( "top2DeltaPtSmearedBest"         ,  &top2DeltaPtSmearedBest        ) ;
  outTree->Branch( "top2DeltaRSmearedBest"          ,  &top2DeltaRSmearedBest         ) ;
  outTree->Branch( "top2DeltaMSmearedBest"          ,  &top2DeltaMSmearedBest         ) ;
  outTree->Branch( "W2DeltaPtSmearedBest"           ,  &W2DeltaPtSmearedBest          ) ;
  outTree->Branch( "W2DeltaRSmearedBest"            ,  &W2DeltaRSmearedBest           ) ;
  outTree->Branch( "W2DeltaMSmearedBest"            ,  &W2DeltaMSmearedBest           ) ;
  outTree->Branch( "bJet2DeltaPtSmearedBest"        ,  &bJet2DeltaPtSmearedBest       ) ;
  outTree->Branch( "bJet2DeltaRSmearedBest"         ,  &bJet2DeltaRSmearedBest        ) ;
  outTree->Branch( "bJet2DeltaMSmearedBest"         ,  &bJet2DeltaMSmearedBest        ) ;
  outTree->Branch( "W2Daughter1DeltaPtSmearedBest"  ,  &W2Daughter1DeltaPtSmearedBest ) ;
  outTree->Branch( "W2Daughter1DeltaRSmearedBest"   ,  &W2Daughter1DeltaRSmearedBest  ) ;
  outTree->Branch( "W2Daughter1DeltaMSmearedBest"   ,  &W2Daughter1DeltaMSmearedBest  ) ;
  outTree->Branch( "W2Daughter2DeltaPtSmearedBest"  ,  &W2Daughter2DeltaPtSmearedBest ) ;
  outTree->Branch( "W2Daughter2DeltaRSmearedBest"   ,  &W2Daughter2DeltaRSmearedBest  ) ;
  outTree->Branch( "W2Daughter2DeltaMSmearedBest"   ,  &W2Daughter2DeltaMSmearedBest  ) ;
  outTree->Branch( "lightJet2DeltaPxSmearedBest"    ,  &lightJet2DeltaPxSmearedBest   ) ;
  outTree->Branch( "lightJet2DeltaPySmearedBest"    ,  &lightJet2DeltaPySmearedBest   ) ;
  outTree->Branch( "lightJet2DeltaRSmearedBest"     ,  &lightJet2DeltaRSmearedBest    ) ;

  //outTree->Branch( "smearingSwitchedLightJetOrdering" ,  &smearingSwitchedLightJetOrdering ) ;
}

Int_t topReconstructionFromLHE::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t topReconstructionFromLHE::LoadTree(Long64_t entry)
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

void topReconstructionFromLHE::Init(TTree *tree)
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

Bool_t topReconstructionFromLHE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void topReconstructionFromLHE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t topReconstructionFromLHE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef topReconstructionFromLHE_cxx
