import ROOT as R


def Map(tf):
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m


def gausSum (x, par):
    return (par[0]*R.TMath.Gaus(x[0], par[1], par[2]) + par[3]*R.TMath.Gaus(x[0], par[1]+par[4], par[5]*par[2]) + par[6]*R.TMath.Gaus(x[0], par[1]-par[4], par[5]*par[2]))


ov = ['4.0', '6.0', '8.0']
th = [10, 20, 40]

R.gROOT.SetBatch(1)

c1 = R.TCanvas("c1", "c1", 800, 600)

fil1 = R.TFile("timeAnalysisSIM.root")
histos_2D = Map(fil1)
fil2 = R.TFile("timeAnalysis.root")
histos = Map(fil2)
chi = 0


for v in ov:
    for t in th:
        key = "ov{}_th{}".format(v, t)
        h = histos_2D["tDiffVSenRatio_"+key]
        histos['tDiffProjection_'+key] = histos_2D["tDiffVSenRatio_"+key].ProjectionY('tDiffProjection_'+key, 0, -1, "e")
        f = R.TF1(h.GetName()+"_gausSumFit", gausSum, -1000, 1000, 7)
        hmax = h.GetXaxis().GetBinCenter(h.GetMaximumBin())
        hrms = h.GetRMS()
        f.SetParameter(0, 100)
        f.SetParLimits(0, 0, 999999)
        f.SetParameter(1, 100)
        f.SetParLimits(1, hmax-50, hmax+50)
        f.SetParameter(2, 100)
        f.SetParLimits(2, 50, 180)
        f.SetParameter(3, 50)
        f.SetParLimits(3, 0, 999999)
        f.SetParameter(4, 350)
        f.SetParLimits(4, 50, 500)
        f.SetParameter(5, 2.)
        f.SetParLimits(5, 1., 5.)
        f.SetParameter(6, 50)
        f.SetParLimits(6, 0, 999999)

        histos_2D["tDiffVSenRatio_"+key].FitSlicesY(f, 0, -1, 0, "R", 0)
        if (histos_2D["tDiffVSenRatio_"+key].GetFunction('f') != NULL):
            chi= chi + histos_2D["tDiffVSenRatio_"+key].GetFunction('f').GetChisquare()



fit_res = {}
f = []
for i in range(0, 100):
    f = R.Math.IMultiGenFunction()

def Chi2(par):
    retval = 0
    Nfun = f.size()
    p1 = []
    for ifun in range(0, Nfun):
        p1[0] = par[2*ifun]
        p1[1] = par[2*ifun+1]
        p1[2] = par[2*ifun+2]
        p1[3] = par[2*Nfun+3]
        p1[4] = par[2*Nfun+4]
        p1[5] = par[2*Nfun]
        p1[6] = par[2*Nfun+1]
        retval = retval + p1[0]
        





     




