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

    //TFile* ifile1 = new TFile("outfile_ias25_pt60_ih0_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    //TFile* ifile1 = new TFile("outfile_ias100_pt60_ih29_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_eta50_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    TFile* ifile1 = new TFile("outfile_ias25_pt60_ih0_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");
    TFile* ifile2 = new TFile("outfile_ias25_pt60_ih0_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed.root");

    TH1F* h1 = (TH1F*)ifile1->Get("mass_obs");
    TH1F* h2 = (TH1F*)ifile2->Get("mass_obs");
    TH1F* h3 = (TH1F*)ifile1->Get("mass_pred2D");
    TH1F* h2R = (TH1F*)ifile1->Get("mass_pred2DR");
    TH1F* h3R = (TH1F*)ifile2->Get("mass_pred2DR");

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

    h1->GetYaxis()->SetTitle("Fraction of tracks/bin");
    h1->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
    h1->GetYaxis()->SetLabelSize(20); //font size
    h1->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
    h1->GetYaxis()->SetTitleSize(20); //font size
    h1->GetYaxis()->SetNdivisions(503);


    h1->SetMarkerStyle(22);
    h1->SetMarkerColor(2);
    h1->SetMarkerSize(1.0);
    h1->SetLineColor(2);
    h1->SetFillColor(0);
    h1->Draw("same HIST P");

    h2->SetMarkerStyle(23);
    h2->SetMarkerColor(9);
    h2->SetMarkerSize(1.0);
    h2->SetLineColor(9);
    h2->SetFillColor(0);
    h2->Draw("same HIST P");

    h3->SetMarkerStyle(20);
    h3->SetMarkerColor(1);
    h3->SetMarkerSize(1.0);
    h3->SetLineColor(1);
    h3->SetFillColor(0);
    h3->Draw("same P");


    //leg = new TLegend(0.82,0.85,0.47,0.69);
    leg = new TLegend(0.82,0.85,0.3,0.5);

    leg->AddEntry(h1,"Observed","PE1");
    leg->AddEntry(h2,"Observed with p cut","PE1");
    leg->AddEntry(h3,"Prediction (ih,p) template","PE1");

    leg->Draw("same");

    c1->cd();
    t2->cd();

   TH1D* frameR = new TH1D("frameR", "frameR", 1,0, 1600);
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


    h2R->SetMarkerStyle(22);
    h2R->SetMarkerColor(2);
    h2R->SetMarkerSize(1.);
    h2R->SetLineColor(2);
    h2R->SetFillColor(0);
    h2R->Draw("same E1");

    h3R->SetMarkerStyle(23);
    h3R->SetMarkerColor(9);
    h3R->SetMarkerSize(1.);
    h3R->SetLineColor(9);
    h3R->SetFillColor(0);
    h3R->Draw("same E1");

    TLine* LineAtOne = new TLine(0,1,1600,1); LineAtOne->SetLineStyle(3); LineAtOne->Draw("same");

    TLatex T;
    T.SetTextAlign(33);
    T.SetTextSize (0.15);
    //T.DrawLatex (1580, 2.4, "R = #int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");

    c1->cd();

    c1->SaveAs("plotMass.pdf");
   
}
