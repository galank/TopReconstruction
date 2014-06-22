#include "TROOT.h"

void loadFiles(TString dirName, int nSmearings, int whichLoop, int maxLoops)
{
  gROOT->LoadMacro("neutrinoSolutions.cxx+g");
  gROOT->LoadMacro("NeutrinoEllipseCalculator.cxx+g");
  gROOT->LoadMacro("lightJetChiSquareMinimumSolver.cxx+g");
  gROOT->LoadMacro("topSystemChiSquare.cxx+g");
  gROOT->LoadMacro("ttbarReconstructionFromLHE.C+g");
  ttbarReconstructionFromLHE m(dirName);
  m.Loop(nSmearings,whichLoop,maxLoops);
}
