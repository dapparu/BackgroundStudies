import ROOT

import subprocess
import os,subprocess,sys


listFiles = ["outfile_ias25_pt60_ih0_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_pThresholdVarBins80_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_massThresholdVarBins500_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso1_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET1_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"]

c1 = ROOT.TCanvas()
c1.SetLogy()
it = 0
'''
for f in listFiles:
    command = ["rootprint -d dir_"+f+" -D colz "+f+".root:ih_p_all"+" &"]
    #process = subprocess.Popen(command,shell=True)
    fi = ROOT.TFile(f+".root")
    h = fi.Get("massFrom1DTemplatesEtaBinning_all")
    c1.cd()
    h.SetLineColor(2*it)
    h.Draw("same")
    it+=1

c1.SaveAs("c1.pdf")
'''
inputfile = ROOT.TFile("outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root")
h3 = inputfile.Get("ih_p_eta_all")
#h3.RebinX(2)
c2=ROOT.TCanvas()
c2.cd()
'''
for x in range(0,h3.GetNbinsX()+1):
    h3.GetXaxis().SetRange(x,x+1)
    h2 = h3.Project3D("zy")
    h2.GetXaxis().SetRangeUser(0,2000)
    h2.GetYaxis().SetRangeUser(0,6)
    h2.Draw("colz")
    prof2=h2.ProfileX()
    prof2.SetLineColor(2)
    prof2.Draw("same")
#    c2.SaveAs("ih_p_eta/eta"+str(x)+".pdf")
'''

h_p_c=inputfile.Get("ih_p_eta_regionC")
h_p_d=inputfile.Get("ih_p_eta_regionD")
c3=ROOT.TCanvas()
c3.cd()
h_p_c.RebinX(2)
h_p_d.RebinX(2)
for x in range(0,h_p_c.GetNbinsX()+1):
    h_p_c.GetXaxis().SetRange(x,x+1)
    h_p_d.GetXaxis().SetRange(x,x+1)
    p_c=h_p_c.Project3D("y")
    p_d=h_p_d.Project3D("y")
    p_c.GetXaxis().SetRangeUser(0,2000)
    p_d.GetXaxis().SetRangeUser(0,2000)
    if p_c.Integral()>0:
        p_c.Scale(1./p_c.Integral())
    if p_d.Integral()>0:
        p_d.Scale(1./p_d.Integral())
    p_c.Draw()
    p_d.SetLineColor(2)
    p_d.Draw("same")
    c3.SetLogy()
    c3.SaveAs("ih_p_eta/p_c_d_"+str(x)+".pdf")
