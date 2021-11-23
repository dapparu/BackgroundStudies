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

    //std::string inputfilename = "outfile-9eta9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-21eta-9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-15eta-9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile9eta15_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-50eta50_ias25_pt60_ih29_p2000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-50eta50_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso1_invMET0_TOF_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    std::string inputfilename = "outfile-50eta50_ias200_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso1_invMET0_TOF_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-9eta9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_2_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-9eta9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_2_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";
    //std::string inputfilename = "outfile-9eta9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed";

    //TFile* ifile1 = new TFile("outfile_ias25_pt60_ih0_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    //TFile* ifile1 = new TFile("outfile_ias100_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_eta50_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    //TFile* ifile1 = new TFile("outfile_ias25_pt60_ih0_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    //TFile* ifile2 = new TFile("outfile_ias25_pt60_ih0_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");

    //TFile* ifile1 = new TFile("outfile_ias50_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
//    TFile* ifile1 = new TFile("outfile-21eta-9_ias50_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    TFile* ifile1 = new TFile((inputfilename+".root").c_str());
//    TFile* ifile1 = new TFile("outfile_ias50_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
//    TFile* ifile1 = new TFile("outfile_ias50_pt60_ih29_p13000_etabins120_ihbins1000_pbins2000_massbins6500_invIso0_invMET0_eta50_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
/* 
    TH1F* h1 = (TH1F*)ifile1->Get("mass_obs_boundedIas");
    TH1F* h2 = (TH1F*)ifile1->Get("mass_pred1D_boundedIas");
    TH1F* h3 = (TH1F*)ifile1->Get("mass_predBC_boundedIas");
    TH1F* h2R = (TH1F*)ifile1->Get("mass_pred1DR_boundedIas");
    TH1F* h3R = (TH1F*)ifile1->Get("mass_predBCR_boundedIas");
*/

    TH1F* h1 = (TH1F*)ifile1->Get("mass_obs");
    TH1F* h2 = (TH1F*)ifile1->Get("mass_pred1D");
    TH1F* h3 = (TH1F*)ifile1->Get("mass_predBC");
    TH1F* h4 = (TH1F*)ifile1->Get("mass_predDB");
    TH1F* h5 = (TH1F*)ifile1->Get("mass_predDC");
    TH1F* h2R = (TH1F*)ifile1->Get("mass_pred1DR");
    TH1F* h3R = (TH1F*)ifile1->Get("mass_predBCR");
    TH1F* h4R = (TH1F*)ifile1->Get("mass_predDBR");
    TH1F* h5R = (TH1F*)ifile1->Get("mass_predDCR");


    TH1F* h1b = (TH1F*)ifile1->Get("mass_obs");

    gStyle->SetOptStat(0);

    TLegend* leg;
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    t1->SetTopMargin(0.06);
    t1->SetBottomMargin(0.001);
    c1->cd();
    TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2); 
    t2->Draw();
    t2->cd();
    t2->SetGridy(true);
    t2->SetPad(0,0.0,1.0,0.2);
    t2->SetTopMargin(0.001);
    t2->SetBottomMargin(0.5);
    t1->cd();
    c1->SetLogy(true);

    h1->GetYaxis()->SetRangeUser(1e-5,1);
    h1->GetYaxis()->SetTitle("Fraction of tracks / bin");
    h1->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
    h1->GetYaxis()->SetLabelSize(20); //font size
    h1->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
    h1->GetYaxis()->SetTitleSize(20); //font size
    h1->GetYaxis()->SetNdivisions(503);


    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(1);
    h1->SetMarkerSize(0.8);
    h1->SetLineColor(1);
    h1->SetFillColor(0);
    h1->Draw("same HIST P");

    h1b->SetMarkerStyle(20);
    h1b->SetMarkerColor(30);
    h1b->SetMarkerSize(0.8);
    h1b->SetLineColor(30);
    h1b->SetFillColor(0);
    //h1b->Draw("same HIST P");


    h2->SetMarkerStyle(21);
    h2->SetMarkerColor(9);
    h2->SetMarkerSize(0.8);
    h2->SetLineColor(9);
    h2->SetFillColor(0);
    h2->Draw("same P");

    h3->SetMarkerStyle(21);
    h3->SetMarkerColor(46);
    h3->SetMarkerSize(0.8);
    h3->SetLineColor(46);
    h3->SetFillColor(0);
    h3->Draw("same P");

    
    h4->SetMarkerStyle(21);
    h4->SetMarkerColor(42);
    h4->SetMarkerSize(0.8);
    h4->SetLineColor(42);
    h4->SetFillColor(0);
    h4->Draw("same P");

    h5->SetMarkerStyle(21);
    h5->SetMarkerColor(12);
    h5->SetMarkerSize(0.8);
    h5->SetLineColor(12);
    h5->SetFillColor(0);
    h5->Draw("same P");




    //leg = new TLegend(0.82,0.85,0.47,0.69);
    leg = new TLegend(0.82,0.85,0.3,0.5);

    leg->AddEntry(h1,"Observed","PE1");
    //leg->AddEntry(h1b,"Observed, 0.05<Ias","PE1");
    leg->AddEntry(h2,"Prediction templates from D","PE1");
    leg->AddEntry(h3,"Prediction templates from B & C","PE1");
    leg->AddEntry(h4,"Prediction templates from B & D","PE1");
    leg->AddEntry(h5,"Prediction templates from C & D","PE1");

    leg->Draw("same");

    c1->cd();
    t2->cd();

   TH1D* frameR = new TH1D("frameR", "frameR", 1,0, 4000);
   frameR->GetXaxis()->SetNdivisions(505);
   frameR->SetTitle("");
   frameR->SetStats(kFALSE);
   frameR->GetXaxis()->SetTitle("");
   frameR->GetXaxis()->SetTitle("Mass (GeV)");
   frameR->GetYaxis()->SetTitle("Ratio - R");
   frameR->SetMaximum(1.5);
   frameR->SetMinimum(0.5);
   frameR->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetYaxis()->SetLabelSize(20); //font size
   frameR->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetYaxis()->SetTitleSize(20); //font size
   frameR->GetYaxis()->SetNdivisions(503);
   frameR->GetXaxis()->SetNdivisions(505);
   frameR->GetXaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetXaxis()->SetLabelSize(20); //font size
   frameR->GetXaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
   frameR->GetXaxis()->SetTitleSize(20); //font size
   frameR->GetXaxis()->SetTitleOffset(3.75);
   frameR->Draw("AXIS");


    h2R->SetMarkerStyle(21);
    h2R->SetMarkerColor(9);
    h2R->SetMarkerSize(0.8);
    h2R->SetLineColor(9);
    h2R->SetFillColor(0);
    h2R->Draw("same E1");

    h3R->SetMarkerStyle(21);
    h3R->SetMarkerColor(46);
    h3R->SetMarkerSize(0.8);
    h3R->SetLineColor(46);
    h3R->SetFillColor(0);
    h3R->Draw("same E1");

    h4R->SetMarkerStyle(21);
    h4R->SetMarkerColor(42);
    h4R->SetMarkerSize(0.8);
    h4R->SetLineColor(42);
    h4R->SetFillColor(0);
    h4R->Draw("same E1");

    h5R->SetMarkerStyle(21);
    h5R->SetMarkerColor(12);
    h5R->SetMarkerSize(0.8);
    h5R->SetLineColor(12);
    h5R->SetFillColor(0);
    h5R->Draw("same E1");


    TLine* LineAtOne = new TLine(0,1,4000,1); LineAtOne->SetLineStyle(3); LineAtOne->Draw("same");

    TLatex T;
    T.SetTextAlign(33);
    T.SetTextSize (0.15);
    //T.DrawLatex (1580, 2.4, "R = #int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");

    c1->cd();

    c1->SaveAs(("plotsMassPdf/MassPlot_200p_"+inputfilename+".pdf").c_str());
   
}
