import ROOT as R
import os

chain = R.TChain("data")
files = os.listdir(".")

for file in files:
    chain.Add(file)

nent=chain.GetEntries()
channels=[213, 219]
histos={}

histos['STEP1']=R.TH1F('STEP1', 'STEP1', 10, 0, 10)
histos['TH1']=R.TH1F('TH1','TH1', 100, 0, 100.5)
for ch in channels:
    histos['qfineVStot_ch'+str(ch)]=R.TProfile('qfineVStot_ch'+str(ch), 'qfineVStot_ch'+str(ch), 1000, 0, 1000.5)

histos['E_'+str(channels[0])+'VSE_'+str(channels[1])]=R.TProfile('E_'+str(channels[0])+'VSE_'+str(channels[1]), 'E_'+str(channels[0])+'VSE_'+str(channels[1]), 500, 5, 25.5)

ranges={}

for ev in range(0,nent):
    chain.GetEntry(ev)
    th1=(chain.step2/10000)%100-1
    histos['STEP1'].Fill(chain.step1)
    histos['TH1'].Fill(th1)
    for ch in channels:
        histos['qfineVStot_ch'+str(ch)].Fill(chain.tot[ch]/1000, chain.qfine[ch])

c1=R.TCanvas("c1", "c1", 600, 400)
   
for ch in channels:
    ranges['ch'+str(ch)]={}
    histos['qfineVStot_ch'+str(ch)].Fit("pol0","R","", 50.5, 196.5)
    ranges['ch'+str(ch)]['qfineMin']=int(histos['qfineVStot_ch'+str(ch)].GetFunction("pol0").GetParameter(0))+2
    histos['qfineVStot_ch'+str(ch)].Fit("pol0","R","", 650.5, 800.5)
    ranges['ch'+str(ch)]['qfineMax']=int(histos['qfineVStot_ch'+str(ch)].GetFunction("pol0").GetParameter(0))-20
    histos['qfineVStot_ch'+str(ch)].Fit("pol1","R","", 250.5, 575.5)
    ranges['ch'+str(ch)]['qfineSlope']=histos['qfineVStot_ch'+str(ch)].GetFunction("pol1").GetParameter(1)

ov=[]
for i in range(1,histos['STEP1'].GetNbinsX()+1):
    if (histos['STEP1'].GetBinContent(i)>0):
        ov.append(histos['STEP1'].GetXaxis().GetBinLowEdge(i))
print ov

for v in ov:           
    histos['EnergyTot'+'_ov%.1f'%v]=R.TH1F('EnergyTot'+'_ov%.1f'%v, 'EnergyTot'+'_ov%.1f'%v, 400, 0, 100.5)
    for ch in channels:
        histos['Energy_ch'+str(ch)+'_ov%.1f'%v]=R.TH1F('Energy_ch'+str(ch)+'_ov%.1f'%v, 'Energy_ch'+str(ch)+'_ov%.1f'%v, 400, 0, 100.5)
        histos['Tot_ch'+str(ch)+'_ov%.1f'%v]=R.TH1F('Tot_ch'+str(ch)+'_ov%.1f'%v, 'Tot_ch'+str(ch)+'_ov%.1f'%v, 200, 0, 1000.5)



for i in range(0, nent/4):
    chain.GetEntry(i)
    o='_ov%.1f'%(chain.step1)
    etot = 0
    ech = 0

    for ch in channels:
        if (chain.qfine[ch] < ranges['ch'+str(ch)]['qfineMin']):
            continue
        if (chain.qfine[ch] > ranges['ch'+str(ch)]['qfineMax']):
            continue
        if (chain.energy[ch]<0):
            continue
            
        ech=ech+1
        etot=etot+chain.energy[ch]
    
        histos['Energy_ch'+str(ch)+o].Fill(chain.energy[ch])
        histos['Tot_ch'+str(ch)+o].Fill(chain.tot[ch])
    if(ech>1):
        histos['EnergyTot'+o].Fill(etot)
        histos['E_'+str(channels[0])+'VSE_'+str(channels[1])].Fill(chain.energy[channels[0]], chain.energy[channels[1]])


histos['E_'+str(channels[0])+'VSE_'+str(channels[1])].Draw()
c1.SaveAs('E_'+str(channels[0])+'VSE_'+str(channels[1])+'.png')

#Definisco una funzione che sia convoluzione di una distribuzione di Landau e di una Gaussiana. Utilizzo i seguenti parametri per il fit:
#par[0] = larghezza (scala) della densita di Landau
#par[1] = most probable (MP, location) parametro della distribuzione di landau
#par[2] = area totale (integrale da inf a +inf), costante di normalizzazione
#par[3] = larghezza (SIGMA) della funzione Gaussiana convoluta

#Nella distribuzione di Landau il massimo si trova a x=-0.22278298 con il par[1] pari a 0. In questa funzione lo shift viene corretto, cosi che il massimo reale sia identico al parametro piu probabile.

def LanGaufun(x,par):

    invsq2pi = 0.3989422804014  #(2pi)^(-1/2)
    mpshift = -0.22278298     
    np = 100.0                 #numero di passi
    sc = 5.0                   #la convoluzione si estende tra + e - sc sigma gaussiane
    ssum = 0.0

    mpc = par[1] - mpshift * par[0] #Correzone di MP shift
    #il range dela convoluzione ha come estremi xlow e xupp
    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]
    step = (xupp - xlow) / np

    #integrale di convoluzione tra una Landau e una Gaussiana tramite la somma
    for i in range(1, int(np/2)+1):
        xx = xlow + (float(i)-0.5) * step
        fland = R.TMath.Landau(xx, mpc, par[0]) / par[0]
        ssum = ssum + fland * R.TMath.Gaus(x[0], xx, par[3])

        xx = xupp - (float(i)-0.5) * step
        fland = R.TMath.Landau(xx, mpc, par[0]) / par[0]
        ssum = ssum + fland * R.TMath.Gaus(x[0], xx, par[3])

    return (par[2] * step * ssum * invsq2pi / par[3])


#Voglio utilizzarela mia funzione per fare dei fit, devo implementare un metodo per farlo. I parametri par[i] sono gli stessi definiti nella precedente funzione
def fitLanGaus(h):
    histos[h.GetName()+"_lanGausFit"]=R.TF1(h.GetName()+"_lanGausFit", LanGaufun, 0, 100, 4)
    h.GetXaxis().SetRangeUser(0, 100)
    hmax = h.GetXaxis().GetBinLowEdge(h.GetMaximumBin())
    hrms = h.GetRMS()
    histos[h.GetName()+"_lanGausFit"].SetParameter(0, 1)
    histos[h.GetName()+"_lanGausFit"].SetParLimits(0, 0., 2.)
    histos[h.GetName()+"_lanGausFit"].SetParameter(1,hmax)
    histos[h.GetName()+"_lanGausFit"].SetParameter(2, histos['Energy_ch'+str(ch)+ovS].GetMaximum())
    histos[h.GetName()+"_lanGausFit"].SetParameter(3, 0.001)
    histos[h.GetName()+"_lanGausFit"].SetParLimits(3, 0.001, 2)

    print(hmax, hrms)

    h.Fit(histos[h.GetName()+"_lanGausFit"], "RB", "", hmax-hmax*0.15, hmax+0.8*hrms)
    h.GetXaxis().SetRangeUser(hmax-hmax*0.15, hmax+0.8*hrms)
    c1.SaveAs(h.GetName()+"_lanGausFit.png")


R.gStyle.SetOptFit(111111)
for v in ov:
    ovS='_ov%.1f'%(v)
    for ch in channels:
        fitLanGaus(histos['Energy_ch'+str(ch)+ovS])
    fitLanGaus(histos['EnergyTot'+ovS])


histos['EnergyVSOV']=R.TGraphErrors()
histos['EnergyVSOV'].SetName('EnergyVSOV')
for ch in channels:
    histos['Energy_ch'+str(ch)+'VSOV']=R.TGraphErrors()
    histos['Energy_ch'+str(ch)+'VSOV'].SetName('Energy_ch'+str(ch)+'VSOV')

mip = {}
for v in ov:
    ovS='_ov%.1f'%(v)
    i=histos['EnergyVSOV'].GetN()
    histos['EnergyVSOV'].SetPoint(i, v , histos['EnergyTot'+ovS+'_lanGausFit'].GetParameter(1))
    for ch in channels:
        i=histos['Energy_ch'+str(ch)+'VSOV'].GetN()
        histos['Energy_ch'+str(ch)+'VSOV'].SetPoint(i, v , histos['Energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1))
        mip[str(ch)+ovS] = histos['Energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1)

print "MIPPPPPPPPPP", mip

R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)

a=R.TH2F("a", "a", 10, 2, 10, 10, 0, 50)
a.GetXaxis().SetTitle("OverVoltage [V]")
a.GetYaxis().SetTitle("MIP Peak [ADC]")
a.SetStats(0)

colors={0:R.kRed, 1:R.kBlue}
labels={0:'Left', 1:'Right'}

l=R.TLatex()
l.SetTextSize(0.05)

a.Draw()
leg=R.TLegend(0.15, 0.65, 0.4, 0.85)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(0, 0)
leg.SetTextSize(0.05)

histos['EnergyVSOV'].SetLineColor(R.kBlack)
histos['EnergyVSOV'].SetMarkerColor(R.kBlack)
histos['EnergyVSOV'].SetMarkerSize(1.2)
histos['EnergyVSOV'].SetMarkerStyle(20)
histos['EnergyVSOV'].Draw("LPSAME")
leg.AddEntry(histos['EnergyVSOV'], "Bar Energy SUM", "LP")
for ich, ch in enumerate(channels):
    histos['Energy_ch'+str(ch)+'VSOV'].SetLineColor(colors[ich])
    histos['Energy_ch'+str(ch)+'VSOV'].SetMarkerColor(colors[ich])
    histos['Energy_ch'+str(ch)+'VSOV'].SetMarkerSize(1.2)
    histos['Energy_ch'+str(ch)+'VSOV'].SetMarkerStyle(20)
    histos['Energy_ch'+str(ch)+'VSOV'].Draw("LPSAME")
    histos['Energy_ch'+str(ch)+'VSOV'].Draw("LPSAME")
    leg.AddEntry(histos['Energy_ch'+str(ch)+'VSOV'], "Bar Energy "+labels[ich], "LP")

leg.Draw()
l.DrawLatexNDC(0.12, 0.93, 'S12572-015C - protons 120 GeV')
c1.SaveAs("EnergyVSOV.png")


for ch in channels:
    for v in ov:
        ovS='_ov%.1f'%(v)
        ranges['ch'+str(ch)]['eMin'+ovS]=histos['Energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1)*0.9
        ranges['ch'+str(ch)]['eMax'+ovS]=histos['Energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1)*2
print "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR", ranges

R.gROOT.ProcessLine(".L timeAnalysis.C+")

def Map(tf):
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m

t=R.timeAnalysis(chain)
t.outFileName="timeAnalysis.root"
t.xMin=19.5-3;
t.xMax=19.5+3;
t.energyMin.resize(len(channels))
t.energyMax.resize(len(channels))
for ich,ch in enumerate(channels):
    t.channels.push_back(ch)
    t.qfineMin.push_back(ranges["ch"+str(ch)]['qfineMin'])
    t.qfineMax.push_back(ranges["ch"+str(ch)]['qfineMax'])
    for v in ov:
        ovS='%.1f'%(v)
        t.energyMin[ich][ovS]=ranges["ch"+str(ch)]['eMin_ov'+ovS]
        t.energyMax[ich][ovS]=ranges["ch"+str(ch)]['eMax_ov'+ovS]
t.Loop()

fTime=R.TFile("timeAnalysis.root")
                        
fOut=R.TFile("histos.root", "RECREATE")
for hN,h in histos.iteritems():
    h.Write()
fOut.Write()
fOut.Close()
                         

               

        


        
