import ROOT as R


R.gROOT.ProcessLine(".L timefitSIM.C")
s=R.timefitSIM()
s.Fit3Peaks()
