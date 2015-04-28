#include "TROOT.h"

void loadFiles(TString dir, int whichLoop, int maxLoops)
{
  gROOT->LoadMacro("neutrinoSolutions.cxx+g");
  gROOT->LoadMacro("WDaughterEllipseCalculator.cxx+g");
  gROOT->LoadMacro("topSystemChiSquare.cxx+g");
  gROOT->LoadMacro("leptonicTopSystemChiSquare.cxx+g");
  gROOT->LoadMacro("hadronicTopSystemChiSquare.cxx+g");
  gROOT->LoadMacro("lightJetChiSquareMinimumSolver.cxx+g");
  gROOT->LoadMacro("topEventMinimizer.cxx+g");
  gROOT->LoadMacro("topReconstructionFromLHE.C+g");
  topReconstructionFromLHE t;
  t.Loop(dir,whichLoop,maxLoops);
}
