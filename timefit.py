import ROOT as R


def Map(tf):
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m


def gausSum (x, par):
    return (par[0]*R.TMath.Gaus(x[0], par[1], par[2]) + par[3]*R.TMath.Gaus(x[0], par[1]+par[4], par[5]*par[2]) + par[6]*R.TMath.Gaus(x[0], par[1]-par[4], par[5]*par[2]))

def G (x,par):
    return (par[0]*R.TMath.Gaus(x[0], par[1], par[2]))



def fitGausSum(h):
    f = R.TF1(h.GetName()+"_gausSumFit", gausSum, -1000, 1000, 7)
    hmax = h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    hrms = h.GetRMS()
    f.SetParameter(0, 100)
    f.SetParLimits(0, 0, 999999)
    f.SetParameter(1, 100)
    f.SetParLimits(1, hmax-200, hmax+200)
    f.SetParameter(2, 100)
    f.SetParLimits(2, 50, 180)
    f.SetParameter(3, 50)
    f.SetParLimits(3, 0, 999999)
    f.SetParameter(4, 350)
    f.SetParLimits(4, 50, 500)
    f.SetParameter(5, 2.)
    f.SetParLimits(5, 0., 5.)
    f.SetParameter(6, 50)
    f.SetParLimits(6, 0, 999999)

    f.Print()
    h.Fit(f, "RB", "", -1000, 1000)
    return f

ov = ['4.0', '6.0', '8.0']
th = [10, 20, 40]

R.gROOT.SetBatch(1)

c1 = R.TCanvas("c1", "c1", 800, 600)
R.gStyle.SetOptFit(111111)
R.gStyle.SetOptTitle(0)

l = R.TLatex()
l.SetTextSize(0.05)

a = R.TH2F("a", "a", 10, 0, 50, 10, 0, 100)
a.GetXaxis().SetTitle("Thresholds [DAC counts]")
a.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a.SetStats(0)

histosOut = {}
histosOut['ctr_FNAL']=R.TGraphErrors()
histosOut['ctr_FNAL'].SetName('ctr_FNAL')
histosOut['ctrRel_FNAL']=R.TGraphErrors()
histosOut['ctrRel_FNAL'].SetName('ctrRel_FNAL')

f = R.TFile("timeAnalysis.root")
histos = Map(f)

for v in ov:
    histosOut['tDiff_ov'+v] = R.TGraphErrors()
    histosOut['tDiff_ov'+v].SetName('tDiff_ov'+v)
    histosOut['tDiff_ov'+v+'_corrected'] = R.TGraphErrors()
    histosOut['tDiff_ov'+v+'_corrected'].SetName('tDiff_ov'+v+'_corrected')

    for it, t in enumerate(th):
        key = "ov{}_th{}".format(v, str(t))
        ff = fitGausSum(histos["tDiff_"+key])
        histosOut['tDiff_ov'+v].SetPoint(it, t, ff.GetParameter(2)/2.)
        histosOut['tDiff_ov'+v].SetPointError(it, 0, ff.GetParError(2)/2.)
        c1.SaveAs('tDiff_'+key+'.png')

a.Draw()
leg = R.TLegend(0.15, 0.65, 0.4, 0.85)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(0, 0)
leg.SetTextSize(0.05)

for v in ov:
    for it, t in enumerate(th):
        key = "ov{}_th{}_corrected".format(v, str(t))
        ff = fitGausSum(histos["tDiff_"+key])
        histosOut['tDiff_ov'+v+'_corrected'].SetPoint(it, t, ff.GetParameter(2)/2.)
        histosOut['tDiff_ov'+v+'_corrected'].SetPointError(it, 0, ff.GetParError(2)/2.)
        c1.SaveAs('tDiff_'+key+'.png')

p1 = [0, 0, 0]
p2 = [0, 0, 0]
p3 = [0, 0, 0]
for v in ov:
    for it, t in enumerate(th):
         key = "ov{}_th{}_corrected".format(v, str(t))
         ff = fitGausSum(histos["tDiff_"+key])
         p1[0]=ff.GetParameter(0)
         p1[1]=ff.GetParameter(1)
         p1[2]=ff.GetParameter(2)
         p2[0]=ff.GetParameter(3)
         p2[1]=ff.GetParameter(1)-ff.GetParameter(1)
         p2[2]=ff.GetParameter(2)*ff.GetParameter(5)
         p3[0]=ff.GetParameter(6)
         p3[1]=ff.GetParameter(4)+ff.GetParameter(1)
         p3[2]= ff.GetParameter(2)*ff.GetParameter(5)
   
         f1 = R.TF1("_G1", G, -1000, 1000, 7)
         f1.SetParameter(0, p1[0])
         f1.SetParameter(1, p1[1])
         f1.SetParameter(2, p1[2])
         f1.SetLineColor(R.kOrange)
         f2 = R.TF1("_G2", G, -1000, 1000, 7)
         f2.SetParameter(0, p2[0])
         f2.SetParameter(1, p2[1])
         f2.SetParameter(2, p2[2])
         f2.SetLineColor(R.kBlue)
         f3 = R.TF1("_G2", G, -1000, 1000, 7)
         f3.SetParameter(0, p3[0])
         f3.SetParameter(1, p3[1])
         f3.SetParameter(2, p3[2])
         f3.SetLineColor(R.kGreen)
         f1.Draw()
         f2.Draw("SAME")
         f3.Draw("SAME")
         c1.SaveAs("3Peaks_"+key+".png")
   
         

        

a.Draw()
leg = R.TLegend(0.15, 0.65, 0.4, 0.85)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(0, 0)
leg.SetTextSize(0.05)

for iv, v in enumerate(ov):
    histosOut['tDiff_ov'+v].SetMarkerColor(R.kBlack+iv)
    histosOut['tDiff_ov'+v].SetMarkerStyle(24)
    histosOut['tDiff_ov'+v].SetMarkerSize(1.2)
    histosOut['tDiff_ov'+v].SetLineColor(R.kBlack+iv)
    histosOut['tDiff_ov'+v].SetLineStyle(3)
    histosOut['tDiff_ov'+v].Draw("LPSAME")
    histosOut['tDiff_ov'+v+'_corrected'].SetMarkerColor(R.kBlack+iv)
    histosOut['tDiff_ov'+v+'_corrected'].SetMarkerStyle(20)
    histosOut['tDiff_ov'+v+'_corrected'].SetMarkerSize(1.2)
    histosOut['tDiff_ov'+v+'_corrected'].SetLineColor(R.kBlack+iv)
    histosOut['tDiff_ov'+v+'_corrected'].SetLineStyle(1)
    histosOut['tDiff_ov'+v+'_corrected'].Draw("LPSAME")
    leg.AddEntry(histosOut['tDiff_ov'+v+'_corrected'], "OV=%s V"%v, "LP")
leg.Draw()
l.DrawLatexNDC(0.12, 0.93, 'S12572-015C - protons 120 GeV')
c1.SaveAs("tDiff.png")

a = R.TH2F("a", "a", 10, 0, 10, 10, 0, 100)
a.GetXaxis().SetTitle("Overvoltage [V]")
a.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a.SetStats(0)

histosOut['tDiff_TH20'] = R.TGraphErrors()
histosOut['tDiff_TH20'].SetName('tDiff_TH20')
histosOut['tDiff_TH20_NC'] = R.TGraphErrors()
histosOut['tDiff_TH20_NC'].SetName('tDiff_TH20_NC')

for iv,v in enumerate(ov):
    for it, t in enumerate(th):
        if (t==20):
            key = "ov{}_th{}_corrected".format(v, str(t))
            ff = fitGausSum(histos["tDiff_"+key])
            print t, v, iv
            histosOut['tDiff_TH20'].SetPoint(iv, float(v), ff.GetParameter(2)/2.)

for iv,v in enumerate(ov):
    for it, t in enumerate(th):
        if (t==20):
            key = "ov{}_th{}".format(v,str(t))
            ff = fitGausSum(histos["tDiff_"+key])
            print t, v, iv
            histosOut['tDiff_TH20_NC'].SetPoint(iv, float(v), ff.GetParameter(2)/2.)


a.Draw()
leg = R.TLegend(0.15, 0.65, 0.4, 0.85)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(0, 0)
leg.SetTextSize(0.05)

histosOut['tDiff_TH20'].SetMarkerColor(R.kOrange)
histosOut['tDiff_TH20'].SetMarkerStyle(20)
histosOut['tDiff_TH20'].SetMarkerSize(1.2)
histosOut['tDiff_TH20'].SetLineColor(R.kOrange)
histosOut['tDiff_TH20'].SetLineStyle(1)
histosOut['tDiff_TH20'].Draw("LPSAME")
histosOut['tDiff_TH20_NC'].SetMarkerColor(R.kOrange)
histosOut['tDiff_TH20_NC'].SetMarkerStyle(24)
histosOut['tDiff_TH20_NC'].SetMarkerSize(1.2)
histosOut['tDiff_TH20_NC'].SetLineColor(R.kOrange)
histosOut['tDiff_TH20_NC'].SetLineStyle(3)
histosOut['tDiff_TH20_NC'].Draw("LPSAME")
leg.AddEntry(histosOut['tDiff_TH20'], "TH=20", "LP")
leg.Draw()
l.DrawLatexNDC(0.12, 0.93, 'S12572-015C - protons 120 GeV')
c1.SaveAs("tDiff_TH20.png")


        
            
