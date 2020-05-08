#ifndef timefitSIM_h
#define timefitSIM_h

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRandom.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Math/GSLMinimizer.h"

#include <vector>
#include <iostream>

class timefitSIM {
 public:

  std::vector<float> cost1;
  std::vector<float> mean1;
  std::vector<float> sigma1;
  std::vector<float> cost2;
  std::vector<float> cost3;

  float shift_23;
  float sigma_23;

  std::vector<float> cost1_err;
  std::vector<float> mean1_err;
  std::vector<float> sigma1_err;
  std::vector<float> cost2_err;
  std::vector<float> cost3_err;

  float shift_23_err;
  float sigma_23_err;

  timefitSIM();
  virtual ~timefitSIM();
  virtual timefitSIM fitSimultaneous (std::vector<TH1F*> histos, double xMin, double xMax);
  virtual void  FitSIM();
  virtual void  Fit3Peaks();
  
};

#endif
