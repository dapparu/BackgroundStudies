import subprocess
import os,subprocess,sys

#sys.argv.append('-b-')
import ROOT

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
#inputfile = ROOT.TFile("outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root")
#inputfile = ROOT.TFile("outfile_ias25_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root")
#inputfile = ROOT.TFile("outfile_ias25_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root")
inputfile = ROOT.TFile("outfile-50eta50_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso1_invMET0_TOF.root")
h3 = inputfile.Get("ih_p_eta_regionB")
regB = inputfile.Get("ih_p_eta_regionB")
regC = inputfile.Get("ih_p_eta_regionC")
regD = inputfile.Get("ih_p_eta_regionD")
regCias = inputfile.Get("ias_p_eta_regionC")
regDias = inputfile.Get("ias_p_eta_regionD")
regB.Rebin3D(2,50,10)
regC.Rebin3D(2,50,10)
regD.Rebin3D(2,50,10)
regCias.Rebin3D(2,50,10)
regDias.Rebin3D(2,50,10)
#h3.RebinX(2)
c2=ROOT.TCanvas()
c2.SetLogy()
c3=ROOT.TCanvas()
c2.Divide(2,2)
ROOT.gStyle.SetOptStat(0)

def chi2(h,b):
    chi2=0
    df=0
    for i in range(1,b):
        if h.GetBinContent(i)<=0:
            continue
        df+=1
        chi2+=((h.GetBinContent(i)-1)*(h.GetBinContent(i)-1))/(h.GetBinError(i)*h.GetBinError(i))
    return ROOT.TMath.Prob(chi2,df)
    #return format(ROOT.TMath.Prob(chi2,df),".3f")
    #return format(chi2/df,".3f")
    #return chi2,df,chi2/df

def uptobin(h,val):
    for i in range(1,h.GetNbinsX()+1):
        if h.GetBinContent(i)<val and h.GetBinContent(i)>0:
            return i

def fit(h,b):
    fr=h.Fit("pol1","SR","same",0,b)
    print ("fr:",fr)
    pv=ROOT.TPaveText(0.6,0.6,0.88,0.88,"ARC,NDC")
    if fr==None:
        return pv
    pv.AddText("par0: "+format(fr.Parameter(0),".3f")+" +/- "+format(fr.ParError(0),".3f"))
    pv.AddText("par1: "+format(fr.Parameter(1),".3f")+" +/- "+format(fr.ParError(1),".3f"))
    pv.AddText("chi2/dof: "+format(fr.Chi2(),".3f")+" / "+str(fr.Ndf()))
    return pv

def functIasP(x):
    eta=int(100*regCias.GetXaxis().GetBinLowEdge(x))
    regCias.GetXaxis().SetRange(x,x+1)
    regDias.GetXaxis().SetRange(x,x+1)
    if regCias.Integral()<=0:
        return
    c3.cd()
    C_ias_p = regCias.Project3D("zy").Clone()
    D_ias_p = regDias.Project3D("zy").Clone()
    C_ias_p_pfy=C_ias_p.ProfileY().Clone()
    D_ias_p_pfy=D_ias_p.ProfileY().Clone()
    C_ias_p_pfy.Draw()
    D_ias_p_pfy.Draw("same")
    D_ias_p_pfy.SetLineColor(2)
    C_ias_p_pfy.GetYaxis().SetRangeUser(100,200)
    c3.SaveAs("templatesEta/ias_p_pfy"+str(x)+"_eta"+str(eta)+".pdf")

def funct(x):
    eta=int(100*regB.GetXaxis().GetBinLowEdge(x))
    leg1=ROOT.TLegend(0.7,0.7,0.9,0.9)
    leg2=ROOT.TLegend(0.7,0.7,0.9,0.9)
    regB.GetXaxis().SetRange(x,x+1)
    if regB.Integral()<=0:
        return 0,0,0,0
    regC.GetXaxis().SetRange(x,x+1)
    regD.GetXaxis().SetRange(x,x+1)
    c2.cd(1)
    c2.SetLogy()
    ROOT.gPad.SetLogy()
    C_p=regC.Project3D("y").Clone()
    D_p=regD.Project3D("y").Clone()
    a=C_p.KolmogorovTest(D_p)
    text1="Kolmogorov test between C and D: "+format(C_p.KolmogorovTest(D_p),".3f")
    C_p.SetTitle(text1)
    bintotestP=uptobin(D_p,30)
#    if regC.Integral()<=0:
#        continue
    C_p.Scale(1./C_p.Integral())
    D_p.Scale(1./D_p.Integral())
    C_p.Sumw2()
    D_p.Sumw2()
    leg1.AddEntry(C_p,"Region C","lep")
    leg1.AddEntry(D_p,"Region D","lep")
    C_p.Draw("E0")
    C_p.GetXaxis().SetRangeUser(0,4000)
    D_p.Draw("E0,SAME")
    leg1.Draw("SAME")
    C_p.SetLineColor(46)
    D_p.SetLineColor(1)
    c2.cd(2)
    c2.SetLogy()
    ROOT.gPad.SetLogy()
    B_ih=regB.Project3D("z").Clone()
    D_ih=regD.Project3D("z").Clone()
    b=B_ih.KolmogorovTest(D_ih)
    text2="Kolmogorov test between B and D: "+format(B_ih.KolmogorovTest(D_ih),".3f")
    B_ih.SetTitle(text2)
    bintotestIh=uptobin(D_ih,30)
    B_ih.Scale(1./B_ih.Integral())
    D_ih.Scale(1./D_ih.Integral())
    B_ih.Sumw2()
    D_ih.Sumw2()
    leg2.AddEntry(B_ih,"Region B","lep")
    leg2.AddEntry(D_ih,"Region D","lep")
    B_ih.Draw("E0")
    B_ih.GetXaxis().SetRangeUser(0,10)
    D_ih.Draw("E0,SAME")
    leg2.Draw("SAME")
    B_ih.SetLineColor(8)
    D_ih.SetLineColor(1)
    c2.cd(3)
    C_pR = C_p.Clone()
    C_pR.Divide(D_p)
    C_pR.Draw()
    C_pR.GetYaxis().SetTitle("C/D")
    C_pR.GetYaxis().SetRangeUser(0,2)
    #p1=fit(C_pR,bintotestP)
    #p1.Draw("same")
    p1=ROOT.TPaveText(0.6,0.6,0.8,0.8,"NDCARC")
    c=chi2(C_pR,bintotestP)
    text3="Chi2 probability compatibility with 1: "+format(chi2(C_pR,bintotestP),".3f")
    p1.AddText(text3)
    p1.SetAllWith("Test","size",12)
    C_pR.SetTitle(text3)
    #p1.Draw("same")
    c2.Update()
    c2.cd(4)
    B_ihR = B_ih.Clone()
    B_ihR.Divide(D_ih)
    B_ihR.Draw()
    B_ihR.GetYaxis().SetTitle("B/D")
    B_ihR.GetYaxis().SetRangeUser(0,2)
    d=chi2(B_ihR,bintotestIh)
    text4="Chi2 probability compatibility with 1: "+format(chi2(B_ihR,bintotestIh),".3f")
    B_ihR.SetTitle(text4)
    p2=fit(B_ihR,bintotestIh)
    p2.Draw("same")
    c2.Update()
    c2.SaveAs("templatesEta/largerbinning_bin"+str(x)+"_eta"+str(eta)+".C")
    c2.SaveAs("templatesEta/largerbinning_bin"+str(x)+"_eta"+str(eta)+".pdf")
    del C_p
    del D_p
    del B_ih
    del D_ih
    print "func:",a,b,c,d
    return a,b,c,d

ha=ROOT.TH1F("",";#eta;Kolmogorov test C & D",60,-3,3)
hb=ROOT.TH1F("",";#eta;Kolmogorov test B & D",60,-3,3)
hc=ROOT.TH1F("",";#eta;Chi2 compatibility C & D at 95%",60,-3,3)
hd=ROOT.TH1F("",";#eta;Chi2 compatibility B & D at 95%",60,-3,3)

for x in range(0,regB.GetNbinsX()+1):

    functIasP(x)


'''    a,b,c,d=funct(x)
    print a,b,c,d
    ha.SetBinContent(x,a)
    hb.SetBinContent(x,b)
    if c>0.05:
        hc.SetBinContent(x,1)
    else:
        hc.SetBinContent(x,0)
    if d>0.05:
        hd.SetBinContent(x,1)
    else:
        hd.SetBinContent(x,0)

hc.SetFillColor(1)
hc.SetLineColor(1)
hd.SetFillColor(1)
hd.SetLineColor(1)
c5=ROOT.TCanvas("","",600,600)
ha.Draw()
c5.SaveAs("ha.pdf")
hb.Draw()
c5.SaveAs("hb.pdf")
hc.Draw()
c5.SaveAs("hc.pdf")
hd.Draw()
c5.SaveAs("hd.pdf")'''



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
    c2.SaveAs("ih_p_eta_regionB/regB_eta"+str(x)+".pdf")
'''
'''
c_1D_all = inputfile.Get("plotting_mass1D_all")
c_1D_all->cd(1)
h_obs_all = c_1D_all.GetPrimitive("massFromTree_all")
h_pred_all = c_1D_all.GetPrimitive("massFrom1DTemplatesEtaBinning_all")
c_1D_all->cd(2)
h_ratio = c_1D_all.GetPrimitive("massFrom1DTemplatesEtaBinning_all")
'''

'''
ih_p_c=inputfile.Get("ih_p_eta_regionC")
ih_p_d=inputfile.Get("ih_p_eta_regionD")
c3=ROOT.TCanvas()
c4=ROOT.TCanvas()
ih_p_c.RebinX(2)
ih_p_d.RebinX(2)
for x in range(0,ih_p_c.GetNbinsX()+1):
    ih_p_c.GetXaxis().SetRange(x,x+1)
    ih_p_d.GetXaxis().SetRange(x,x+1)
    p_c=ih_p_c.Project3D("y")
    p_d=ih_p_d.Project3D("y")
    ih_c=ih_p_c.Project3D("z")
    ih_d=ih_p_d.Project3D("z")
    p_c.GetXaxis().SetRangeUser(0,2000)
    p_d.GetXaxis().SetRangeUser(0,2000)
    ih_c.GetXaxis().SetRangeUser(0,10)
    ih_d.GetXaxis().SetRangeUser(0,10)
    if p_c.Integral()>0:
        p_c.Scale(1./p_c.Integral())
    if p_d.Integral()>0:
        p_d.Scale(1./p_d.Integral())
    if ih_c.Integral()>0:
        ih_c.Scale(1./ih_c.Integral())
    if ih_d.Integral()>0:
        ih_d.Scale(1./ih_d.Integral())
    c3.cd()
    p_c.Draw()
    p_d.SetLineColor(2)
    p_d.Draw("same")
    c3.SetLogy()
    #c3.SaveAs("ih_p_eta/p_c_d_"+str(x)+".pdf")
    c4.cd()
    ih_c.Draw()
    ih_d.SetLineColor(2)
    ih_d.Draw("same")
    c4.SetLogy()
    c4.SaveAs("ih_p_eta/ih_c_d_"+str(x)+".pdf")
    '''
