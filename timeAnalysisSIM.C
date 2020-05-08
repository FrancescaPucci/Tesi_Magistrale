#define timeAnalysisSIM_cxx
#include "timeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include "timeAnalysisSIM.h"

#include <iostream>

//DA QUI

#include "TH1F.h"
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

#include <vector>
#include <iostream>

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


using std::begin;
using std::end;
//A QUI

void timeAnalysisSIM::Loop()
{
//   In a ROOT session, you can do:
//      root> .L timeAnalysis.C
//      root> timeAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int thresholds[] = {10, 20, 40};
  float overvoltages[] = {4.0, 6.0, 8.0};
  float mip_213[] = {4.84, 10.30, 19.69};
  float mip_219[] = {3.39, 8.14, 16.14};
  
  
  for (auto& th: thresholds){
    for (auto& ov: overvoltages){
      
      profs[Form("enratioVStDiff_ov%.1f_th%d",ov,th)] = new TProfile(Form("enratioVStDifff_ov%.1f_th%d",ov,th), Form("enratioVStDiff_ov%.1f_th%d",ov,th), 60, 0.85, 1.15);
      
      histos[Form("tDiff_ov%.1f_th%d_corrected",ov,th)] = new TH1F(Form("tDiff_ov%.1f_th%d_corrected",ov,th), Form("tDiff_ov%.1f_th%d_corrected",ov,th), 120, -2000, 2000);

      fitfunc[Form("f_ov%.1f_th%d",ov,th)] =new TF1(Form("f_ov%.1f_th%d",ov,th), "pol2", -2000, 2000);

    }
  }

  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain -> GetEntry(jentry);
    nbytes += nb;
    int ich = -1;
    int nch = 0;
    int th1 = int (step2/10000)%100 - 1;

    if (step1!=4 && step1!=6 && step1!=8)
      continue;
    if (th1!=10 && th1!=20 && th1!=40)
      continue;
    if (ntracks!=1)
      continue;
    if (x_dut<xMin)
      continue;
    if (x_dut>xMax)
      continue;

    for (auto& ch: channels)
      {
	++ich;
	if (qfine[ch]<qfineMin[ich])
	  continue;
	if (qfine[ch]>qfineMax[ich])
	  continue;
	if (energy[ch]<energyMin[ich][Form("%.1f", step1)])
	  continue;
	if (energy[ch]>energyMax[ich][Form("%.1f", step1)])
	  continue;
	++nch;
      }

    if(nch<2)
      continue;

    //histos[Form("tDiff_ov%.1f_th%d", step1, th1)]->Fill((time[channels[0]] - time[channels[1]]));

     if(step1 == 4.0){
       int i=0;
       profs[Form("enratioVStDiff_ov%.1f_th%d", step1, th1)]->Fill(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])), (time[channels[0]] - time[channels[1]]));
     }else if (step1==6.0){
       int i=1;
       profs[Form("enratioVStDiff_ov%.1f_th%d", step1, th1)]->Fill(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])), (time[channels[0]] - time[channels[1]]));
     }else if (step1==8.0){
       int i=2;
       profs[Form("enratioVStDiff_ov%.1f_th%d", step1, th1)]->Fill(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])), (time[channels[0]] - time[channels[1]]));
     }
 
  }
   
  for (auto& th: thresholds){
    for (auto& ov: overvoltages){
      double p0, p1, p2;
      p0 = profs[Form("enratioVStDiff_ov%.1f_th%d", ov, th)]->Fit("pol2", "S", "", 0.85, 1.15)->Parameter(0);
      p1 = profs[Form("enratioVStDiff_ov%.1f_th%d", ov, th)]->Fit("pol2", "S", "", 0.85, 1.15)->Parameter(1);
      p2 = profs[Form("enratioVStDiff_ov%.1f_th%d", ov, th)]->Fit("pol2", "S", "", 0.85, 1.15)->Parameter(2);
      fitfunc[Form("f_ov%.1f_th%d",ov,th)]->SetParameter(0, p0);
      fitfunc[Form("f_ov%.1f_th%d",ov,th)]->SetParameter(1, p1);
      fitfunc[Form("f_ov%.1f_th%d",ov,th)]->SetParameter(2, p2);
    }
  }
 
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain -> GetEntry(jentry);
    nbytes += nb;
    int ich = -1;
    int nch = 0;
    int th1 = int (step2/10000)%100 - 1;

    if (step1!=4 && step1!=6 && step1!=8)
      continue;
    if (th1!=10 && th1!=20 && th1!=40)
      continue;
    if (ntracks!=1)
      continue;
    if (x_dut<xMin)
      continue;
    if (x_dut>xMax)
      continue;

    for (auto& ch: channels)
      {
	++ich;
	if (qfine[ch]<qfineMin[ich])
	  continue;
	if (qfine[ch]>qfineMax[ich])
	  continue;
	if (energy[ch]<energyMin[ich][Form("%.1f", step1)])
	  continue;
	if (energy[ch]>energyMax[ich][Form("%.1f", step1)])
	  continue;
	++nch;
      }

    if(nch<2)
      continue;

     if(step1 == 4.0){
       int i=0;
       double v = fitfunc[Form("f_ov%.1f_th%d",4.0,10)]->Eval(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])));
       histos[Form("tDiff_ov%.1f_th%d_corrected", step1, th1)]->Fill((time[channels[0]] - time[channels[1]]) - v);
     }else if (step1==6.0){
       int i=1;
      double v = fitfunc[Form("f_ov%.1f_th%d",4.0,10)]->Eval(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])));
       histos[Form("tDiff_ov%.1f_th%d_corrected", step1, th1)]->Fill((time[channels[0]] - time[channels[1]]) - v);
     }else if (step1==8.0){
       int i=2;
       double v = fitfunc[Form("f_ov%.1f_th%d",4.0,10)]->Eval(((energy[channels[0]]/mip_213[i])/(energy[channels[1]]/mip_219[i])));
       histos[Form("tDiff_ov%.1f_th%d_corrected", step1, th1)]->Fill((time[channels[0]] - time[channels[1]]) - v);
     }
  }

  
  fOut = new TFile(outFileName.c_str(), "RECREATE");
  for (auto& h: histos)
   h.second->Write();
  
  fOut->Write();
  fOut->Save();
  fOut->Close();
      
}



