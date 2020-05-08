#define timefitSIM_cxx
#include "timefitSIM.h"

timefitSIM::timefitSIM() {
}
timefitSIM::~timefitSIM() {
}
struct GlobalChi2 {

  GlobalChi2 (std::vector<ROOT::Math::IMultiGenFunction*> f):
    f_(f) {};
  double operator() (const double *par) const {
    double val = 0;
    int Nfun = f_.size();
    for (int ifun=0; ifun < Nfun; ifun++) {
      double p1[7];
      p1[0] = par[5*ifun];
      p1[1] = par[5*ifun+1];
      p1[2] = par[5*ifun+2];
      p1[3] = par[5*ifun+3];
      p1[4] = par[5*ifun+4];
      p1[5] = par[5*Nfun];
      p1[6] = par[5*Nfun+1];
      val += (*(f_[ifun]))(p1);
    }
    return val;
  }

  std::vector<ROOT::Math::IMultiGenFunction*> f_;

};

Double_t gausSum (Double_t *x, Double_t *par) {
  Double_t G = par[0]*(TMath::Gaus(x[0], par[1], par[2]));
  G += par[3]*(TMath::Gaus(x[0], par[1]+par[5], par[6]*par[2]));
  G += par[4]*(TMath::Gaus(x[0], par[1]-par[5], par[6]*par[2]));
  return G;
}

Double_t G (Double_t *x, Double_t *par) {
  return par[0]*(TMath::Gaus(x[0], par[1], par[2]));
}


timefitSIM timefitSIM::fitSimultaneous (std::vector<TH1F*> histos, double xMin, double xMax) {
  timefitSIM fr;
  int Nhist = histos.size();

  std::vector<TF1*> func;
  std::vector<ROOT::Math::WrappedMultiTF1*> wf;
  std::vector<ROOT::Fit::BinData*> data;
  std::vector<ROOT::Math::IMultiGenFunction*> chi2;

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range;
  range.SetRange(xMin, xMax);

  const int Npar = 5*Nhist + 2;
  double par0[Npar];
  ROOT::Fit::Fitter fitter;

  int nBins = 0;
  for (int ifun=0; ifun<Nhist; ifun++) {
    TF1* f = new TF1 (Form("ftDiff_%d", ifun), gausSum, xMin, xMax, 7);
    func.push_back(f);
    ROOT::Math::WrappedMultiTF1* wf1= new ROOT::Math::WrappedMultiTF1(*f,1);
    wf.push_back(wf1);
    ROOT::Fit::BinData* data1=new ROOT::Fit::BinData(opt,range);
    data.push_back(data1);
    ROOT::Fit::FillData(*(data[ifun]), histos[ifun]);
    ROOT::Fit::PoissonLLFunction* chi2_1= new ROOT::Fit::PoissonLLFunction(*(data[ifun]), *(wf[ifun]));
    chi2.push_back((ROOT::Math::IMultiGenFunction*) chi2_1);

    nBins += (*(data[ifun])).Size();
  }

  GlobalChi2 globalChi2(chi2);

  
  for (int ifun = 0; ifun<Nhist; ifun++){
    par0[5*ifun] = 100;
    par0[5*ifun+1] = 100;
    par0[5*ifun+2] = 100;
    par0[5*ifun+3] = 50;
    par0[5*ifun+4] = 50;
  }
  par0[5*Nhist] = 350;
  par0[5*Nhist+1] = 2.;

  fitter.Config().SetParamsSettings(Npar, par0);
  
  for (int ifun=0; ifun<Nhist; ifun++){
    double hmax = histos[ifun]->GetXaxis()->GetBinCenter(histos[ifun]->GetMaximumBin());
    double hmin = hmax-200;
    hmax+=200;
    fitter.Config().ParSettings(5*ifun).SetLimits( 0, 999999);
    fitter.Config().ParSettings(5*ifun+1).SetLimits( hmin, hmax);
    fitter.Config().ParSettings(5*ifun+2).SetLimits( 50, 150);
    fitter.Config().ParSettings(5*ifun+3).SetLimits( 0, 999999);
    fitter.Config().ParSettings(5*ifun+4).SetLimits( 0, 999999);
  }
  fitter.Config().ParSettings(5*Nhist).SetLimits( 200, 500);
  fitter.Config().ParSettings(5*Nhist+1).SetLimits( 0., 5.);
  
    ROOT::Math::MinimizerOptions fitopt;
    fitopt.SetTolerance(0.001);
    fitopt.SetErrorDef(0.5);
    fitter.Config().SetMinimizerOptions(fitopt);

    fitter.FitFCN(Npar, globalChi2, par0, nBins, kTRUE);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

    for (int ifun = 0; ifun<Nhist; ifun++) {
      fr.cost1.push_back(result.Value(5*ifun));
      fr.mean1.push_back(result.Value(5*ifun+1));
      fr.sigma1.push_back(result.Value(5*ifun+2));
      fr.cost2.push_back(result.Value(5*ifun+3));
      fr.cost3.push_back(result.Value(5*ifun+4));
      fr.cost1_err.push_back(result.Error(5*ifun));
      fr.mean1_err.push_back(result.Error(5*ifun+1));
      fr.sigma1_err.push_back(result.Error(5*ifun+2));
      fr.cost2_err.push_back(result.Error(5*ifun+3));
      fr.cost3_err.push_back(result.Error(5*ifun+4));
    }
    fr.shift_23 = result.Value(Npar-2);
    fr.sigma_23 = result.Value(Npar-1);
    fr.shift_23_err = result.Error(Npar-2);
    fr.sigma_23_err = result.Error(Npar-1);
    return fr;
}

void timefitSIM::FitSIM () {

  TCanvas *c = new TCanvas("c", "c", 800,600);
  TFile* out = TFile::Open("tDiff_sim.root", "RECREATE");

  int i = 0;
  std::vector<TH1F*> histograms;
  TFile *in = TFile::Open("hist.root");
  TIter next(in->GetListOfKeys());
  TKey *key;
  while((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1F")) continue;
    histograms.push_back((TH1F*)key->ReadObj());
    }
  
  std::cout << "# histos = " << histograms.size() << std::endl;
  int Nhist = histograms.size();
  
  timefitSIM fr = fitSimultaneous(histograms,-1000,1000);
  for (int j=0; j<Nhist; j++){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    histograms[j]->SetMarkerStyle(20);
    histograms[j]->SetMarkerSize(0.6);
    histograms[j]->SetMarkerColor(kBlack);
    histograms[j]->Draw();
    TF1* f = new TF1(Form("Histo_%d", j), gausSum, -1000, 1000, 100);
    f->SetParameter(0, fr.cost1[j]);
    f->SetParameter(1, fr.mean1[j]);
    f->SetParameter(2, fr.sigma1[j]);
    f->SetParameter(3, fr.cost2[j]);
    f->SetParameter(4, fr.cost3[j]);
    f->SetParameter(5, fr.shift_23);
    f->SetParameter(6, fr.sigma_23);
    f->SetLineColor(kYellow);
    f->SetLineWidth(3);
 
    f->Draw("SAME");
    c->SaveAs(Form("FitSIM_%d.png", j+1));

  }
  /*TF1* f1 = new TF1("f1", G, -1000, 1000);
    f1->SetParameter(0, 585.977 );
    f1->SetParameter(1, 51.6622);
    f1->SetParameter(2, 92.4277);
    f1->Draw();
    c1->SaveAs("FitSIM_A.png");
  std::cout << "AAAAAAAAAAAAAAAAAAA" << std::endl;
  
    
  */
    /*
       auto f2 = new TF1("f2", gaus , -1000, 1000);
    f2->SetParameter(0, fr.cost2[j]);
    f2->SetParameter(1, ((fr.mean1[j])+(fr.shift_23)));
    f2->SetParameter(2, ((fr.sigma_23)*(fr.sigma1[j])));
    f2->SetLineColor(kGreen);
    f2->Draw("SAME");
auto f3 = new TF1("f3","[cost] *TMath::Gaus(x, [sigma], [mean]) ", -1000, 1000);
    f3->SetParameter("cost", fr.cost3[j]);
    f3->SetParameter("mean", ((fr.mean1[j])-(shift_23)));
    f3->SetParameter("sigma", ((fr.sigma_23)*(fr.sigma1[j])));
    f3->SetLineColor(kBlue);
    f3->Draw("SAME");
    */
    

  out->Write();
  out->Save();
}



void timefitSIM::Fit3Peaks () {
  TCanvas *c = new TCanvas("c", "c", 800,600);
  TFile* out = TFile::Open("tDiff_sim.root", "RECREATE");

  int i = 0;
  std::vector<TH1F*> histograms1;
  std::vector<TH1F*> histograms2;
  std::vector<TH1F*> histograms3;
  TFile *in = TFile::Open("timeAnalysisSIM.root");
  if(!in->IsOpen() || !out->IsOpen()){
    std::cout << "Problems opening one of the root files. Exiting." << std::endl;
    exit(-1);
  }
  auto hi1 = in->Get("tDiff_ov4.0_th10_corrected");
  auto h1 = dynamic_cast<TH1F*>(hi1);
  auto hi2 = in->Get("tDiff_ov4.0_th20_corrected");
  auto h2 = dynamic_cast<TH1F*>(hi2);
  auto hi3 = in->Get("tDiff_ov4.0_th40_corrected");
  auto h3 = dynamic_cast<TH1F*>(hi3);
  auto hi4 = in->Get("tDiff_ov6.0_th10_corrected");
  auto h4 = dynamic_cast<TH1F*>(hi4);
  auto hi5 = in->Get("tDiff_ov6.0_th20_corrected");
  auto h5 = dynamic_cast<TH1F*>(hi5);
  auto hi6 = in->Get("tDiff_ov6.0_th40_corrected");
  auto h6 = dynamic_cast<TH1F*>(hi6);
  auto hi7 = in->Get("tDiff_ov8.0_th10_corrected");
  auto h7 = dynamic_cast<TH1F*>(hi7);
  auto hi8 = in->Get("tDiff_ov8.0_th20_corrected");
  auto h8 = dynamic_cast<TH1F*>(hi8);
  auto hi9 = in->Get("tDiff_ov8.0_th40_corrected");
  auto h9 = dynamic_cast<TH1F*>(hi9);

  histograms1.push_back(h1);
  histograms1.push_back(h2);
  histograms1.push_back(h3);
  histograms2.push_back(h4);
  histograms2.push_back(h5);
  histograms2.push_back(h6);
  histograms3.push_back(h7);
  histograms3.push_back(h8);
  histograms3.push_back(h9);

  int Nhist = histograms1.size();
  
  timefitSIM fr1 = fitSimultaneous(histograms1,-1000,1000);
  for (int j=0; j<Nhist; j++){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    histograms1[j]->SetMarkerStyle(20);
    histograms1[j]->SetMarkerSize(0.6);
    histograms1[j]->SetMarkerColor(kBlack);
    histograms1[j]->GetXaxis()->SetTitle("tDiff [ps]");
    histograms1[j]->GetYaxis()->SetTitle("Entries");
    histograms1[j]->Draw();

    TF1* f = new TF1(Form("Histo_%d", j), gausSum, -1000, 1000, 100);
    TF1* f1 = new TF1(Form("Histo1_%d", j), G, -1000, 1000, 100);
    TF1* f2 = new TF1(Form("Histo2_%d", j), G, -1000, 1000, 100);
    TF1* f3 = new TF1(Form("Histo3_%d", j), G, -1000, 1000, 100);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.027);
    auto legend = new TLegend(0.53, 0.75, 0.89, 0.89);
    f->SetParameter(0, fr1.cost1[j]);
    f->SetParameter(1, fr1.mean1[j]);
    f->SetParameter(2, fr1.sigma1[j]);
    f->SetParameter(3, fr1.cost2[j]);
    f->SetParameter(4, fr1.cost3[j]);
    f->SetParameter(5, fr1.shift_23);
    f->SetParameter(6, fr1.sigma_23);
    f1->SetParameter(0, fr1.cost1[j]);
    f1->SetParameter(1, fr1.mean1[j]);
    f1->SetParameter(2, fr1.sigma1[j]);
    f2->SetParameter(0, fr1.cost2[j]);
    f2->SetParameter(1, fr1.mean1[j]+fr1.shift_23);
    f2->SetParameter(2, fr1.sigma1[j]*fr1.sigma_23);
    f3->SetParameter(0, fr1.cost3[j]);
    f3->SetParameter(1, fr1.mean1[j]-fr1.shift_23);
    f3->SetParameter(2, fr1.sigma1[j]*fr1.sigma_23);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->SetLineStyle(2);
    f2->SetLineColor(kGreen+2);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f3->SetLineColor(kGreen+2);
    f3->SetLineWidth(2);
    f3->SetLineStyle(2);
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
    
    legend->AddEntry(f, Form("Central peak width: #sigma = %.2f", fr1.sigma1[j]), "");
    legend->AddEntry(f, Form("#splitline{Fraction of events in }{the central peak = %.2f}", ((fr1.cost1[j]*sqrt(2.*M_PI)*fr1.sigma1[j])/((fr1.cost1[j]*sqrt(2.*M_PI)*fr1.sigma1[j])+(fr1.cost2[j]*sqrt(2.*M_PI)*fr1.sigma_23*fr1.sigma1[j])+(fr1.cost3[j]*sqrt(2.*M_PI)*fr1.sigma_23*fr1.sigma1[j])))), "");  //(fr1.cost1[j]/histograms1[j]->GetEntries()))
 
    f2->Draw("SAME");
    f3->Draw("SAME");
    f1->Draw("SAME");
    f->Draw("SAME");
    legend->Draw();
    c->SaveAs(Form("FitSIM_ov4.0_%d.png", j+1));
  }

  timefitSIM fr2 = fitSimultaneous(histograms2,-1000,1000);
  for (int j=0; j<Nhist; j++){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    histograms2[j]->SetMarkerStyle(20);
    histograms2[j]->SetMarkerSize(0.6);
    histograms2[j]->SetMarkerColor(kBlack);
    histograms2[j]->GetXaxis()->SetTitle("tDiff [ps]");
    histograms2[j]->GetYaxis()->SetTitle("Entries");
    histograms2[j]->Draw();

    TF1* f = new TF1(Form("Histo_%d", j), gausSum, -1000, 1000, 100);
    TF1* f1 = new TF1(Form("Histo1_%d", j), G, -1000, 1000, 100);
    TF1* f2 = new TF1(Form("Histo2_%d", j), G, -1000, 1000, 100);
    TF1* f3 = new TF1(Form("Histo3_%d", j), G, -1000, 1000, 100);
    auto legend = new TLegend(0.53, 0.75, 0.89, 0.89);
    f->SetParameter(0, fr2.cost1[j]);
    f->SetParameter(1, fr2.mean1[j]);
    f->SetParameter(2, fr2.sigma1[j]);
    f->SetParameter(3, fr2.cost2[j]);
    f->SetParameter(4, fr2.cost3[j]);
    f->SetParameter(5, fr2.shift_23);
    f->SetParameter(6, fr2.sigma_23);
    f1->SetParameter(0, fr2.cost1[j]);
    f1->SetParameter(1, fr2.mean1[j]);
    f1->SetParameter(2, fr2.sigma1[j]);
    f2->SetParameter(0, fr2.cost2[j]);
    f2->SetParameter(1, fr2.mean1[j]+fr2.shift_23);
    f2->SetParameter(2, fr2.sigma1[j]*fr2.sigma_23);
    f3->SetParameter(0, fr2.cost3[j]);
    f3->SetParameter(1, fr2.mean1[j]-fr2.shift_23);
    f3->SetParameter(2, fr2.sigma1[j]*fr2.sigma_23);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->SetLineStyle(2);
    f2->SetLineColor(kGreen+2);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f3->SetLineColor(kGreen+2);
    f3->SetLineWidth(2);
    f3->SetLineStyle(2);
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
    
    legend->AddEntry(f, Form("Central peak width: #sigma = %.2f", fr2.sigma1[j]), "");
    legend->AddEntry(f, Form("#splitline{Fraction of events in }{the central peak = %.2f}", ((fr2.cost1[j]*sqrt(2.*M_PI)*fr2.sigma1[j])/((fr2.cost1[j]*sqrt(2.*M_PI)*fr2.sigma1[j])+(fr2.cost2[j]*sqrt(2.*M_PI)*fr2.sigma_23*fr2.sigma1[j])+(fr2.cost3[j]*sqrt(2.*M_PI)*fr2.sigma_23*fr2.sigma1[j])))), "");       //(fr2.cost1[j]*sqrt(2.*M_PI)/histograms2[j]->GetEntries()))
 
    f2->Draw("SAME");
    f3->Draw("SAME");
    f1->Draw("SAME");
    f->Draw("SAME");
    legend->Draw();
    c->SaveAs(Form("FitSIM_ov6.0_%d.png", j+1));
  }
    timefitSIM fr3 = fitSimultaneous(histograms3,-1000,1000);
  for (int j=0; j<Nhist; j++){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    histograms3[j]->SetMarkerStyle(20);
    histograms3[j]->SetMarkerSize(0.6);
    histograms3[j]->SetMarkerColor(kBlack);
    histograms3[j]->GetXaxis()->SetTitle("tDiff [ps]");
    histograms3[j]->GetYaxis()->SetTitle("Entries");
    histograms3[j]->Draw();
      
    TF1* f = new TF1(Form("Histo_%d", j), gausSum, -1000, 1000, 100);
    TF1* f1 = new TF1(Form("Histo1_%d", j), G, -1000, 1000, 100);
    TF1* f2 = new TF1(Form("Histo2_%d", j), G, -1000, 1000, 100);
    TF1* f3 = new TF1(Form("Histo3_%d", j), G, -1000, 1000, 100);
    auto legend = new TLegend(0.53, 0.75, 0.89, 0.89);
    f->SetParameter(0, fr3.cost1[j]);
    f->SetParameter(1, fr3.mean1[j]);
    f->SetParameter(2, fr3.sigma1[j]);
    f->SetParameter(3, fr3.cost2[j]);
    f->SetParameter(4, fr3.cost3[j]);
    f->SetParameter(5, fr3.shift_23);
    f->SetParameter(6, fr3.sigma_23);
    f1->SetParameter(0, fr3.cost1[j]);
    f1->SetParameter(1, fr3.mean1[j]);
    f1->SetParameter(2, fr3.sigma1[j]);
    f2->SetParameter(0, fr3.cost2[j]);
    f2->SetParameter(1, fr3.mean1[j]+fr3.shift_23);
    f2->SetParameter(2, fr3.sigma1[j]*fr3.sigma_23);
    f3->SetParameter(0, fr3.cost3[j]);
    f3->SetParameter(1, fr3.mean1[j]-fr3.shift_23);
    f3->SetParameter(2, fr3.sigma1[j]*fr3.sigma_23);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->SetLineStyle(2);
    f2->SetLineColor(kGreen+2);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f3->SetLineColor(kGreen+2);
    f3->SetLineWidth(2);
    f3->SetLineStyle(2);
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
std::cout << "ENTRIES: " << histograms3[j]->GetEntries() << "  NORM: " << (fr3.cost1[j]*sqrt(2.*M_PI)*fr3.sigma1[j]) << std::endl;
    legend->AddEntry(f, Form("Central peak width: #sigma = %.2f", fr3.sigma1[j]), "");
    legend->AddEntry(f, Form("#splitline{Fraction of events in }{the central peak = %.2f}", (fr3.cost1[j]*sqrt(2.*M_PI)*fr3.sigma1[j])/((fr3.cost1[j]*sqrt(2.*M_PI)*fr3.sigma1[j])+(fr3.cost2[j]*sqrt(2.*M_PI)*fr3.sigma_23*fr3.sigma1[j])+(fr3.cost3[j]*sqrt(2.*M_PI)*fr3.sigma_23*fr3.sigma1[j]))), "");                    // (fr3.cost1[j]/histograms3[j]->GetEntries()))

    
    f2->Draw("SAME");
    f3->Draw("SAME");
    f1->Draw("SAME");
    f->Draw("SAME");
    legend->Draw();
    c->SaveAs(Form("FitSIM_ov8.0_%d.png", j+1));
  }
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  auto gr1 = new TGraphErrors();
  auto gr2 = new TGraphErrors();
  auto gr3 = new TGraphErrors();
  auto leg = new TLegend(0.68, 0.69, 0.88, 0.89);
  auto l = new TLatex();
  TH2F* a = new TH2F("a", "a", 10, 0, 50, 10, 0, 100);
  a->GetXaxis()->SetTitle("Thresholds [DAC counts]");
  a->GetYaxis()->SetTitle("#sigma_{bar} [ps]");
  a->SetStats(0);
  a->Draw();
  l->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.);
  gr1->SetPoint(0, 10, fr1.sigma1[0]/2);
  gr1->SetPoint(1, 20, fr1.sigma1[1]/2);
  gr1->SetPoint(2, 40, fr1.sigma1[2]/2);
  gr1->SetPointError(0, 0, fr1.sigma1_err[0]/2);
  gr1->SetPointError(1, 0, fr1.sigma1_err[1]/2);
  gr1->SetPointError(2, 0, fr1.sigma1_err[2]/2);
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.7);
  gr1->SetLineColor(kRed);
  gr1->SetLineStyle(1);
  gr1->SetLineWidth(1);
  gr1->Draw("LPSAME");
  leg->AddEntry(gr1, "OV = 4.0V", "LP");
  gr2->SetPoint(0, 10, fr2.sigma1[0]/2);
  gr2->SetPoint(1, 20, fr2.sigma1[1]/2);
  gr2->SetPoint(2, 40, fr2.sigma1[2]/2);
  gr2->SetPointError(0, 0, fr2.sigma1_err[0]/2);
  gr2->SetPointError(1, 0, fr2.sigma1_err[1]/2);
  gr2->SetPointError(2, 0, fr2.sigma1_err[2]/2);
  gr2->SetMarkerColor(kBlue+2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.7);
  gr2->SetLineColor(kBlue+2);
  gr2->SetLineStyle(1);
  gr2->SetLineWidth(1);
  gr2->Draw("LPSAME");
  leg->AddEntry(gr2, "OV = 6.0V", "LP");
  gr3->SetPoint(0, 10, fr3.sigma1[0]/2);
  gr3->SetPoint(1, 20, fr3.sigma1[1]/2);
  gr3->SetPoint(2, 40, fr3.sigma1[2]/2);
  gr3->SetPointError(0, 0, fr3.sigma1_err[0]/2);
  gr3->SetPointError(1, 0, fr3.sigma1_err[1]/2);
  gr3->SetPointError(2, 0, fr3.sigma1_err[2]/2);
  gr3->SetMarkerColor(kGreen+2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.7);
  gr3->SetLineColor(kGreen+2);
  gr3->SetLineStyle(1);
  gr3->SetLineWidth(1);
  leg->AddEntry(gr3, "OV = 8.0V", "LP");
  gr3->Draw("LPSAME");

  leg->Draw();
  l->DrawLatexNDC(0.12, 0.93, "S12572-015C - protons 120 GeV");
  c1->SaveAs("SIM_sigmaVSth.png");
  
  
}


#undef timefitSIM_cxx   
   
  
			  
	 
			 

    
