#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include "HscpCandidates.h"

void makePlotMass()
{

    bool reb = false;
    std::string quan = "90";
    
    bool blind=false;
    if(quan=="90") blind=true;
    if(quan=="50_100") blind=true;
    blind = false;

    bool signal=false;
    bool gluino=true;
    bool stau=false;

    //std::string inputfilename = "outfile_SingleMu_TkOnly_2017C_Ih3p24_-50eta50_ias100_pt60_ih32_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_2017C_Ih3p24_slicesIasOuter_-50eta50_ias100_pt60_ih32_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    
    
   // std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p17_slicesIasAll_-50eta50_ias50_pt60_ih31_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p24_slicesIasAll_-50eta50_ias50_pt60_ih32_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_slicesIasAll_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_slicesIasAll_-50eta50_ias100_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_slicesIasAll_-50eta50_ias200_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    
    
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_slicesIasAll_-50eta50_ias25_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_PtErrCut2_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_slicesIasAll_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_WJetsToLNu_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_QCD_Pt-20_MuEnrichedPt15_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    
    
    
    //std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    //std::string inputfilename = "outfile_WJetsToLNu_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    std::string inputfilename = "outfile_TTTo2L2Nu_TkOnly_2017_Ih3p47_CutPtErrQ95_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins200_pbins200_massbins50_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";

    TFile* ifile1 = new TFile((inputfilename+".root").c_str());
    //TFile* ifile2 = new TFile("outfile_Gluino_M-2000_TkOnly_2017C_Ih3p24_slicesIasAll_-50eta50_ias25_pt60_ih32_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    
    TFile* ifile2 = new TFile("outfile_Stau_M-1599_TkOnly_2017C_slicesIasAll_-50eta50_ias25_pt60_ih32_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    
    //17 mars 2022
    /*TFile* ifile_Gluino1000 = new TFile("outfile_Gluino_M-1000_TkOnly_2017C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1400 = new TFile("outfile_Gluino_M-1400_TkOnly_2017C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino2000 = new TFile("outfile_Gluino_M-2000_TkOnly_2017C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_ppStau871 = new TFile("outfile_ppStau_M-871_TkOnly_2017C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_ppStau1029 = new TFile("outfile_ppStau_M-1029_TkOnly_2017C_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    */
   
    //24 mars 2022
    /*TFile* ifile_Gluino1000 = new TFile("outfile_Gluino_M-1000_TkOnly_2017_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1400 = new TFile("outfile_Gluino_M-1400_TkOnly_2017_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1600 = new TFile("outfile_Gluino_M-1600_TkOnly_2017_Ih3p47_pT500_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino2000 = new TFile("outfile_Gluino_M-2000_TkOnly_2017_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino2600 = new TFile("outfile_Gluino_M-2600_TkOnly_2017_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    */
    TFile* ifile_ppStau871 = new TFile("outfile_ppStau_M-871_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_ppStau1029 = new TFile("outfile_ppStau_M-1029_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");

    //26 avril 2022
    TFile* ifile_Gluino1000 = new TFile("outfile_Gluino_M-1000_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1400 = new TFile("outfile_Gluino_M-1400_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1600 = new TFile("outfile_Gluino_M-1600_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino2000 = new TFile("outfile_Gluino_M-2000_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino2600 = new TFile("outfile_Gluino_M-2600_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    TFile* ifile_Gluino1600on = new TFile("outfile_GluinoOnlyNeutral_M-1600_TkOnly_2017_Ih3p47_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root");
    
    //TFile* ifile2 = new TFile("ooutfile_SingleMu_TkOnly_2017C_Ih2p9_-50eta50_ias100_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed.root");
//    TFile* ifile1 = new TFile("outfile_ias50_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
//    TFile* ifile1 = new TFile("outfile_ias50_pt60_ih29_p13000_etabins120_ihbins1000_pbins2000_massbins6500_invIso0_invMET0_eta50_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
/* 
    TH1F* h1 = (TH1F*)ifile1->Get("mass_obs_boundedIas");
    TH1F* h2 = (TH1F*)ifile1->Get("mass_pred1D_boundedIas");
    TH1F* h3 = (TH1F*)ifile1->Get("mass_predBC_boundedIas");
    TH1F* h2R = (TH1F*)ifile1->Get("mass_pred1DR_boundedIas");
    TH1F* h3R = (TH1F*)ifile1->Get("mass_predBCR_boundedIas");
*/

  //  ifile1=ifile2;

    //std::string quan = "50_90";
    //quan = "";

    //if(reb)quan+="_rebinned";

    /*std::string st1 = "mass_obs_"+quan+"ias_50";
    std::string st2 = "mass_predBC_"+quan+"ias_50";
    std::string st3 = "mass_predBCR_"+quan+"ias_50";
    std::string stSignal = "massObs_q"+quan;*/

    quan = "50_100";

    std::string st1 = "mass_obs_"+quan;
    std::string st2 = "mass_predBC_"+quan;
    std::string st3 = "mass_predBCR_"+quan;
    std::string stSignal = "massObs_q"+quan;
   
if(quan==""){
    st1 = "mass_obs";
    st2 = "mass_predBC";
    st3 = "mass_predBCR";
    stSignal = "massObs";
}


    //TH1F* h1 = (TH1F*)ifile1->Get("mass_obs");
    TH1F* h1 = (TH1F*)ifile1->Get(st1.c_str());
    //TH1F* h2 = (TH1F*)ifile1->Get("mass_pred1D");
    //TH1F* h3 = (TH1F*)ifile1->Get("mass_predBC");
    TH1F* h3 = (TH1F*)ifile1->Get(st2.c_str());
    //TH1F* h4 = (TH1F*)ifile1->Get("mass_predDB");
    //TH1F* h5 = (TH1F*)ifile1->Get("mass_predDC");
    //TH1F* h2R = (TH1F*)ifile1->Get("mass_pred1DR");
    //TH1F* h3R = (TH1F*)ifile1->Get("mass_predBCR");
    TH1F* h3R = (TH1F*)ifile1->Get(st3.c_str());
    //TH1F* h4R = (TH1F*)ifile1->Get("mass_predDBR");
    //TH1F* h5R = (TH1F*)ifile1->Get("mass_predDCR");

    /*overflowLastBin(h1);
    overflowLastBin(h3);
    h1=rebinHisto(h1);
    h3=rebinHisto(h3);
    h3R = ratioIntegral(h3,h1);*/

    TH1F* hSignal = (TH1F*)ifile2->Get("massObs");

    TH1F* h_Gluino1000 = (TH1F*) ifile_Gluino1000->Get(stSignal.c_str());
    TH1F* h_Gluino1400 = (TH1F*) ifile_Gluino1400->Get(stSignal.c_str());
    TH1F* h_Gluino1600 = (TH1F*) ifile_Gluino1600->Get(stSignal.c_str());
    TH1F* h_Gluino2000 = (TH1F*) ifile_Gluino2000->Get(stSignal.c_str());
    TH1F* h_Gluino2600 = (TH1F*) ifile_Gluino2600->Get(stSignal.c_str());
    TH1F* h_ppStau871 = (TH1F*) ifile_ppStau871->Get(stSignal.c_str());
    TH1F* h_ppStau1029 = (TH1F*) ifile_ppStau1029->Get(stSignal.c_str());
    TH1F* h_Gluino1600on = (TH1F*) ifile_Gluino1600on->Get(stSignal.c_str());

    /*overflowLastBin(h1);
    overflowLastBin(h3);
    //overflowLastBin(h3R);
    overflowLastBin(h_Gluino1000);
    overflowLastBin(h_Gluino1400);
    overflowLastBin(h_Gluino1600);
    overflowLastBin(h_Gluino2000);
    overflowLastBin(h_Gluino2600);
    overflowLastBin(h_Gluino1600on);*/

    int max_mass=2000;
    if(reb){
        h1 = rebinHisto(h1);
        h3 = rebinHisto(h3);
        //h3R = rebinHisto(h3R);
        if(signal){
        h_Gluino1000 = rebinHisto(h_Gluino1000);
        h_Gluino1400 = rebinHisto(h_Gluino1400);
        h_Gluino1600 = rebinHisto(h_Gluino1600);
        h_Gluino2000 = rebinHisto(h_Gluino2000);
        h_Gluino2600 = rebinHisto(h_Gluino2600);
        h_Gluino1600on = rebinHisto(h_Gluino1600on);}
        quan+="_rebinned";
    }
    else{
        overflowLastBin(h1,max_mass);
        overflowLastBin(h3,max_mass);
        if(signal){
        overflowLastBin(h_Gluino1000,max_mass);
        overflowLastBin(h_Gluino1400,max_mass);
        overflowLastBin(h_Gluino1600,max_mass);
        overflowLastBin(h_Gluino2000,max_mass);
        overflowLastBin(h_Gluino2600,max_mass);
        overflowLastBin(h_Gluino1600on,max_mass);}
    }
    h3R = ratioIntegral(h3,h1);

    double SystErr = 0.5;
    float Min_bin = 0.005;
    h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    TH1F* h1Err = (TH1F*)h1->Clone();
    TH1F* h3Err = (TH1F*)h3->Clone();
    for(int j=0;j<h3->GetNbinsX()+1;j++){
        double err= sqrt(pow(h3->GetBinError(j),2) + pow(h3->GetBinContent(j)*SystErr,2));
        h3Err->SetBinError(j,err);
        if(h3Err->GetBinContent(j)<Min_bin && j>5){for(int i=j+1;i<h3Err->GetNbinsX();i++){h3Err->SetBinError(i,0);}}
    }

 /*   TH1F* h1_2 = (TH1F*)ifile2->Get("mass_obs");
    TH1F* h3_2 = (TH1F*)ifile2->Get("mass_predBC");
    TH1F* h3R_2 = (TH1F*)ifile2->Get("mass_predBCR");
*/
    //h1->Scale(h1->GetEntries());
    //h1_2->Scale(h1_2->GetEntries());


    float min_entries = 5e-6;
    float max_entries = 1e6;
    

    TH1F* h1b = (TH1F*)ifile1->Get("mass_obs");

    gStyle->SetOptStat(0);
    //gStyle->SetErrorX(0);

    TLegend* leg;
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    TPad* t1 = new TPad("t1","t1", 0.0, 0.30, 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    t1->SetTopMargin(0.06);
    t1->SetBottomMargin(0.001);
    c1->cd();
    TPad* t2 = new TPad("t2","t2", 0.0, 0.2, 1.0, 0.3); 
    t2->Draw();
    //t2->cd();
    t2->SetGridy(true);
    //t2->SetPad(0,0.15,1.0,0.3);
    t2->SetTopMargin(0.1);
    t2->SetBottomMargin(0.05);
    TPad* t3 = new TPad("t3","t3", 0.0, 0.0, 1.0, 0.2); 
    t3->Draw();
    //t3->cd();
    t3->SetGridy(true);
    //t3->SetPad(0,0.0,1.0,0.15);
    t3->SetTopMargin(0.1);
    t3->SetBottomMargin(0.6);
    t1->cd();
    c1->SetLogy(true);


    h1->SetBinContent(h1->GetNbinsX(),h1->GetBinContent(h1->GetNbinsX())+h1->GetBinContent(h1->GetNbinsX()+1));
    //h1->GetYaxis()->SetRangeUser(1e-8,1);
    h3Err->GetYaxis()->SetTitle("Tracks / bin");
    h3Err->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
    h3Err->GetYaxis()->SetLabelSize(20); //font size
    h3Err->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
    h3Err->GetYaxis()->SetTitleSize(20); //font size
    //h3Err->GetYaxis()->SetNdivisions(505);


    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(1);
    h1->SetMarkerSize(1.0);
    h1->SetLineColor(1);
    h1->SetFillColor(0);
    h1->GetXaxis()->SetRangeUser(0,max_mass);
    //h1->GetYaxis()->SetRangeUser(1,1e7);

/*    h1_2->SetMarkerStyle(33);
    h1_2->SetMarkerColor(4);
    h1_2->SetMarkerSize(0.8);
    h1_2->SetLineColor(4);
    h1_2->SetFillColor(0);
    //h1_2->Draw("same HIST P");
*/
/*
    h1b->SetMarkerStyle(20);
    h1b->SetMarkerColor(30);
    h1b->SetMarkerSize(0.8);
    h1b->SetLineColor(30);
    h1b->SetFillColor(0);
    //h1b->Draw("same HIST P");
*/

    /*h2->SetMarkerStyle(21);
    h2->SetMarkerColor(9);
    h2->SetMarkerSize(0.8);
    h2->SetLineColor(9);
    h2->SetFillColor(0);
    //h2->Draw("same P");*/

    if(signal){
        hSignal->SetMarkerStyle(21);
        hSignal->SetMarkerColor(8);
        hSignal->GetXaxis()->SetRangeUser(0,max_mass);
        hSignal->GetYaxis()->SetRangeUser(1e-3,1e6);
        hSignal->SetMarkerSize(1.5);
        hSignal->SetLineColor(1);
        hSignal->SetFillColor(8);
        //hSignal->Draw("same HIST");
       
        if(gluino){
        h_Gluino1000->GetYaxis()->SetTitle("Tracks / bin");
        h_Gluino1000->SetFillStyle(0);
        h_Gluino1000->SetFillColor(8);
        h_Gluino1000->SetLineColor(8);
        h_Gluino1000->SetLineWidth(2);
        h_Gluino1000->SetLineStyle(2);
        h_Gluino1000->GetXaxis()->SetRangeUser(0,max_mass);
        h_Gluino1000->GetYaxis()->SetRangeUser(min_entries,max_entries);
        h_Gluino1000->Draw("same HIST");

        h_Gluino1400->SetFillStyle(0);
        h_Gluino1400->SetLineStyle(1);
        h_Gluino1400->SetLineWidth(2);
        h_Gluino1400->SetFillColor(30);
        h_Gluino1400->SetLineColor(30);
        h_Gluino1400->Draw("same HIST");

        h_Gluino1600->SetFillStyle(0);
        h_Gluino1600->SetLineStyle(3);
        h_Gluino1600->SetLineWidth(2);
        h_Gluino1600->SetFillColor(38);
        h_Gluino1600->SetLineColor(38);
        h_Gluino1600->Draw("same HIST");

        h_Gluino2000->SetFillStyle(0);
        h_Gluino2000->SetLineStyle(7);
        h_Gluino2000->SetLineWidth(2);
        h_Gluino2000->SetFillColor(32);
        h_Gluino2000->SetLineColor(32);
        h_Gluino2000->Draw("same HIST");

        h_Gluino2600->SetFillStyle(0);
        h_Gluino2600->SetLineStyle(9);
        h_Gluino2600->SetLineWidth(2);
        h_Gluino2600->SetFillColor(34);
        h_Gluino2600->SetLineColor(34);
        h_Gluino2600->Draw("same HIST");

        h_Gluino1600on->SetFillStyle(0);
        h_Gluino1600on->SetLineStyle(4);
        h_Gluino1600on->SetLineWidth(2);
        h_Gluino1600on->SetFillColor(46);
        h_Gluino1600on->SetLineColor(46);
        h_Gluino1600on->Draw("same HIST");



        }
        if(stau){
        h_ppStau871->GetXaxis()->SetRangeUser(0,max_mass);
        h_ppStau871->GetYaxis()->SetRangeUser(min_entries,max_entries);
        h_ppStau871->SetFillStyle(0);
        h_ppStau871->SetFillColor(9);
        h_ppStau871->SetLineColor(9);
        h_ppStau871->SetLineWidth(2);
        h_ppStau871->Draw("same HIST");


        h_ppStau1029->SetFillStyle(0);
        h_ppStau1029->SetLineStyle(2);
        h_ppStau1029->SetLineWidth(2);
        h_ppStau1029->SetFillColor(38);
        h_ppStau1029->SetLineColor(38);
        h_ppStau1029->Draw("same HIST");
        }



    }

    h3Err->SetMarkerStyle(22);
    h3Err->SetMarkerColor(5);
    h3Err->SetMarkerSize(1.0);
    h3Err->SetLineColor(5);
    h3Err->SetFillColor(5);
    h3Err->SetFillStyle(1001);
    h3Err->GetXaxis()->SetRangeUser(0,max_mass);
    h3Err->GetYaxis()->SetRangeUser(min_entries,max_entries);
    h3Err->GetXaxis()->SetTitle("");
    h3Err->Draw("same E5");


    h3->SetMarkerStyle(21);
    h3->SetMarkerColor(2);
    h3->SetMarkerSize(1);
    h3->SetLineColor(2);
    h3->SetFillColor(0);
    h3->Draw("same HIST P");


    if(!blind)h1->Draw("same E1");


    


 /*   h3_2->SetMarkerStyle(34);
    h3_2->SetMarkerColor(42);
    h3_2->SetMarkerSize(0.8);
    h3_2->SetLineColor(42);
    h3_2->SetFillColor(0);
    //h3_2->Draw("same P");
*/
    
    /*h4->SetMarkerStyle(21);
    h4->SetMarkerColor(42);
    h4->SetMarkerSize(0.8);
    h4->SetLineColor(42);
    h4->SetFillColor(0);
    //h4->Draw("same P");

    h5->SetMarkerStyle(21);
    h5->SetMarkerColor(12);
    h5->SetMarkerSize(0.8);
    h5->SetLineColor(12);
    h5->SetFillColor(0);
    //h5->Draw("same P");*/




    //leg = new TLegend(0.82,0.85,0.47,0.69);
    leg = new TLegend(0.82,0.9,0.4,0.6);
    //leg = new TLegend(0.82,0.85,0.3,0.5);

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);

    TH1F* h3leg = (TH1F*)h3->Clone();
    h3leg->SetFillColor(h3Err->GetFillColor());
    h3leg->SetFillStyle(h3Err->GetFillStyle());

    if(!blind)leg->AddEntry(h1,"Observed","PE1");
    //leg->AddEntry(h1_2,"Observed, new Ias calculation","PE1");
    //leg->AddEntry(h1b,"Observed, 0.05<Ias","PE1");
    //leg->AddEntry(h2,"Prediction templates from D","PE1");
    //leg->AddEntry(h3,"Prediction templates from B & C","PE1");
    leg->AddEntry(h3leg,"Data-based SM prediction","PF");
    //leg->AddEntry(hSignal,"Gluino (M = 2000 GeV)","F");
    if(signal){
        if(gluino)leg->AddEntry(h_Gluino1000,"Gluino (M = 1000 GeV)","L");
        if(gluino)leg->AddEntry(h_Gluino1400,"Gluino (M = 1400 GeV)","L");
        if(gluino)leg->AddEntry(h_Gluino1600,"Gluino (M = 1600 GeV)","L");
        if(gluino)leg->AddEntry(h_Gluino2000,"Gluino (M = 2000 GeV)","L");
        if(gluino)leg->AddEntry(h_Gluino2600,"Gluino (M = 2600 GeV)","L");
        if(gluino)leg->AddEntry(h_Gluino1600on,"Gluino CS (M = 1600 GeV)","L");
        if(stau)leg->AddEntry(h_ppStau871,"ppStau (M = 871 GeV)","L");
        if(stau)leg->AddEntry(h_ppStau1029,"ppStau (M = 1029 GeV)","L");
        
        //leg->AddEntry(hSignal,"Stau (M = 1599 GeV)","F");
    }
    //leg->AddEntry(h3_2,"Prediction templates from B & C, new Ias calculation","PE1");
    //leg->AddEntry(h4,"Prediction templates from B & D","PE1");
    //leg->AddEntry(h5,"Prediction templates from C & D","PE1");

    leg->Draw("same");
    //if(reb) {TLine* LineLastBin = new TLine(h1->GetBinLowEdge(h1->FindBin(max_mass)),0,h1->GetBinLowEdge(h1->FindBin(max_mass)),max_entries); LineLastBin->SetLineStyle(3); LineLastBin->SetLineColor(1); LineLastBin->Draw("same");}
    //else 
{TLine* LineLastBin = new TLine(h1->GetBinLowEdge(h1->FindBin(max_mass)-1),0,h1->GetBinLowEdge(h1->FindBin(max_mass)-1),max_entries); LineLastBin->SetLineStyle(3); LineLastBin->SetLineColor(1); LineLastBin->Draw("same");}
    c1->cd();
    t2->cd();

   TH1D* frameR = new TH1D("frameR", "frameR", 1,0, max_mass);
   frameR->GetXaxis()->SetNdivisions(505);
   frameR->SetTitle("");
   frameR->SetStats(kFALSE);
   frameR->GetXaxis()->SetTitle("");
   //frameR->GetXaxis()->SetTitle("Mass (GeV)");
   frameR->GetYaxis()->SetTitle("Ratio - R");
   frameR->SetMaximum(2.5);
   frameR->SetMinimum(0.0);
   frameR->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetYaxis()->SetLabelSize(12); //font size
   frameR->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetYaxis()->SetTitleSize(14); //font size
   frameR->GetYaxis()->SetNdivisions(503);
   frameR->GetXaxis()->SetNdivisions(505);
   frameR->GetXaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetXaxis()->SetLabelSize(20); //font size
   frameR->GetXaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetXaxis()->SetTitleSize(20); //font size
   frameR->GetXaxis()->SetTitleOffset(3.75);
   frameR->Draw("AXIS");


    /*h2R->SetMarkerStyle(21);
    h2R->SetMarkerColor(9);
    h2R->SetMarkerSize(0.8);
    h2R->SetLineColor(9);
    h2R->SetFillColor(0);
    //h2R->Draw("same E1");*/

    h3R->SetMarkerStyle(21);
    h3R->SetMarkerColor(2);
    h3R->SetMarkerSize(0.7);
    h3R->SetLineColor(2);
    h3R->SetFillColor(0);
    if(!blind)h3R->Draw("same E0");
    h3R->GetXaxis()->SetRangeUser(0,max_mass);

    
/*    h3R_2->SetMarkerStyle(34);
    h3R_2->SetMarkerColor(42);
    h3R_2->SetMarkerSize(0.8);
    h3R_2->SetLineColor(42);
    h3R_2->SetFillColor(0);
    //h3R_2->Draw("same E1");



    h4R->SetMarkerStyle(21);
    h4R->SetMarkerColor(42);
    h4R->SetMarkerSize(0.8);
    h4R->SetLineColor(42);
    h4R->SetFillColor(0);
    //h4R->Draw("same E1");

    h5R->SetMarkerStyle(21);
    h5R->SetMarkerColor(12);
    h5R->SetMarkerSize(0.8);
    h5R->SetLineColor(12);
    h5R->SetFillColor(0);
    //h5R->Draw("same E1");*/

    TH1F* h1Cl = (TH1F*) h1->Clone();
    //h1Cl->Divide(h1_2);
    //h1Cl->Draw("same E1");

    TLine* LineAtOne = new TLine(0,1,max_mass,1); LineAtOne->SetLineStyle(3); LineAtOne->SetLineColor(1); LineAtOne->Draw("same");
    TLine* LineAt1p2 = new TLine(0,1.2,max_mass,1.2); LineAt1p2->SetLineStyle(4); LineAt1p2->SetLineColor(9); LineAt1p2->Draw("same");
    TLine* LineAt0p8 = new TLine(0,0.8,max_mass,0.8); LineAt0p8->SetLineStyle(4); LineAt0p8->SetLineColor(9); LineAt0p8->Draw("same");

    TLatex T;
    T.SetTextAlign(33);
    T.SetTextSize (0.15);
    //T.DrawLatex (1580, 2.4, "R = #int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");

    c1->cd();
    t3->cd();

    TH1D* frameR2 = new TH1D("frameR2", "frameR2", 1,0, max_mass);
   frameR2->GetXaxis()->SetNdivisions(505);
   frameR2->SetTitle("");
   frameR2->SetStats(kFALSE);
   frameR2->GetXaxis()->SetTitle("");
   frameR2->GetXaxis()->SetTitle("Mass (GeV)");
   //frameR2->GetYaxis()->SetTitle("Ratio - R");
   frameR2->SetMaximum(1.5);
   frameR2->SetMinimum(0.5);
   frameR2->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR2->GetYaxis()->SetLabelSize(12); //font size
   frameR2->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR2->GetYaxis()->SetTitleSize(20); //font size
   frameR2->GetYaxis()->SetNdivisions(503);
   frameR2->GetXaxis()->SetNdivisions(510);
   frameR2->GetXaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR2->GetXaxis()->SetLabelSize(20); //font size
   frameR2->GetXaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR2->GetXaxis()->SetTitleSize(20); //font size
   frameR2->GetXaxis()->SetTitleOffset(3.75);
   frameR2->Draw("AXIS");

   if(!blind)h3R->Draw("same E0");

   LineAtOne->Draw("same");
   LineAt1p2->Draw("same");
   LineAt0p8->Draw("same");

   c1->cd();

   std::string directory = "plotsMassApril28/";

    //c1->SaveAs(("plotsMassPdf/MassPlot_"+inputfilename+"_quantile"+quan+".pdf").c_str());
    if(quan==""){
        c1->SaveAs((directory+"MassPlot_"+inputfilename+".pdf").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+".png").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+".C").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+".root").c_str());
    }
    else {
        c1->SaveAs((directory+"MassPlot_"+inputfilename+"_quantile_"+quan+".pdf").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+"_quantile_"+quan+".png").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+"_quantile_"+quan+".C").c_str());
        c1->SaveAs((directory+"MassPlot_"+inputfilename+"_quantile_"+quan+".root").c_str());
    }
   
}
