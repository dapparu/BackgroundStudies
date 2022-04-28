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

using namespace std;

void massSlices(){
    
    std::string inputfilename = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_pTerrCutEtaBin_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh2_rebinP2_rebinMass25_analysed";
    TFile* ifile1 = new TFile((inputfilename+".root").c_str());

    string prefix = "mass_predBC_";

    string suffixSave = "pt40.pdf";

    string listNameIas50[5] = {
    prefix+"50ias_50",
    prefix+"60ias_50",
    prefix+"70ias_50",
    prefix+"80ias_50",
    prefix+"90ias_50"
    };


    string listNameIas40[6] = {
    prefix+"40ias_40",
    prefix+"50ias_40",
    prefix+"60ias_40",
    prefix+"70ias_40",
    prefix+"80ias_40",
    prefix+"90ias_40"
    };


    string listNamePt50[5] = {
    prefix+"50pt_50",
    prefix+"60pt_50",
    prefix+"70pt_50",
    prefix+"80pt_50",
    prefix+"90pt_50"
    };


    string listNamePt40[6] = {
    prefix+"40pt_40",
    prefix+"50pt_40",
    prefix+"60pt_40",
    prefix+"70pt_40",
    prefix+"80pt_40",
    prefix+"90pt_40"
    };



    string listLegend40[6] = {"40","50","60","70","80","90"};
    string listLegend50[5] = {"50","60","70","80","90"};
    int listColor[6] = {30,38,43,46,12,41};

    TLegend leg(0.7,0.7,0.9,0.9);

    TH1F* h[6]={NULL,NULL,NULL,NULL,NULL,NULL};

    int j=6;

    for(int i=0;i<j;i++){
        h[i]=(TH1F*)ifile1->Get(listNamePt40[i].c_str());
    }

    TCanvas c1;
    gStyle->SetOptStat(0);

    h[0]->GetXaxis()->SetRangeUser(0,2000);
    h[0]->GetYaxis()->SetRangeUser(1e-3,1e5);

    for(int i=0;i<j;i++){
        //scale(h[i]);
        h[i]->SetMarkerColor(listColor[i]);
        h[i]->SetLineColor(listColor[i]);
        h[i]->Draw("same");
        leg.AddEntry(h[i],listLegend40[i].c_str(),"L");
    }

    leg.Draw("same");
    c1.SetLogy();
    
    c1.SaveAs(("massPlotSlices_"+suffixSave).c_str());

}
