#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include "HscpCandidates.h"

using namespace std;

void fitQuantiles(ofstream& ofile, TH2F* h, TH2F* hbis, const float& xq){
    //h->Rebin2D(2,1);
    hbis->Rebin2D(10,10);
    TH1F* h2 = (TH1F*)h->ProfileX()->Clone(); h2->Reset();
    TH1F* h2bis = (TH1F*)hbis->ProfileX()->Clone(); h2bis->Reset();
    double q[1]={xq};
    double p[1];
    double qbis[1]={xq};
    double pbis[1];
    vector<double> x, y;
    vector<double> xbis, ybis;
    for(int i=0;i<h->GetNbinsX()+1;i++){
        h->ProjectionY("",i,i)->GetQuantiles(1,p,q);
        //cout << "x: " << h->GetXaxis()->GetBinCenter(i) << " quantile: " << p[0] << endl;
        h2->SetBinContent(i,h->GetYaxis()->FindBin(p[0]));
        x.push_back(h->GetXaxis()->GetBinCenter(i));
        y.push_back(p[0]);
    }
    for(int i=0;i<hbis->GetNbinsX()+1;i++){
        hbis->ProjectionY("",i,i)->GetQuantiles(1,p,q);
        //cout << "x: " << h->GetXaxis()->GetBinCenter(i) << " quantile: " << p[0] << endl;
        h2bis->SetBinContent(i,hbis->GetYaxis()->FindBin(p[0]));
        xbis.push_back(hbis->GetXaxis()->GetBinCenter(i));
        ybis.push_back(p[0]);
    }

    TCanvas c1;
    TGraph gr(h->GetNbinsX()+1,x.data(),y.data());
    TGraph grbis(hbis->GetNbinsX()+1,xbis.data(),ybis.data());
    hbis->Draw("colz");
    hbis->GetXaxis()->SetRangeUser(0,3000);
    hbis->GetYaxis()->SetRangeUser(0,0.5);
    h->Rebin2D(2,2);
    h->Draw("colz,same");
    gr.Draw("P*");
    gr.GetXaxis()->SetRangeUser(0,3000);
    h->GetXaxis()->SetRangeUser(0,3000);
    gr.GetYaxis()->SetRangeUser(0,0.2);
    h->GetYaxis()->SetRangeUser(0,0.2);
    grbis.Draw("P*");
    grbis.SetMarkerColor(2);
    TFitResultPtr res = gr.Fit("pol1","S","same",50,150);
    string formula = to_string(res->Parameter(0))+"+"+to_string(res->Parameter(1))+"*x";
    ofile << "quantile: " << xq << " formula: " << formula << std::endl;
    TF1 f("f",formula.c_str(),0,6000);
    f.Draw("same");
    string tit = (h->GetName());
    tit+="_"+to_string((float)xq)+"quantile.pdf";
    //c1.SetLogz();
    c1.SaveAs(("fitQuantiles/"+tit).c_str());
}

void extractLineAtXsigma(ofstream& ofile, TH2F* h, const float& x){
    TCanvas c1;
    h->Draw("colz");
    //h->GetXaxis()->SetRangeUser(0,600);
    //h->GetYaxis()->SetRangeUser(0,0.2);
    TFitResultPtr res = h->ProfileX()->Fit("pol1","S","same",50,500);
    float sigma = h->ProjectionY()->GetStdDev();
    string formula = to_string(res->Parameter(0)+x*sigma)+"+"+to_string(res->Parameter(1))+"*x";
    ofile << "name: " << h->GetName() << " p0: " << res->Parameter(0) << " p1: " << res->Parameter(1) << " p0+" << x << "sig: " << res->Parameter(0)+x*sigma << " formula: " << formula << std::endl;
    TF1 f("f",formula.c_str(),50,500);
    string tit = (h->GetName());
    tit+="_"+to_string((int)x)+"sig.pdf";
    f.Draw("same");
    c1.SetLogz();
    c1.SaveAs(("cut_pTerr_eta_pdf/"+tit).c_str()); 
}

TF1 f_eta_0p8("f_eta_0p8","0.006260+0.000258*x",0,10000);
TF1 f_0p8_eta_1p7("f_0p8_eta_1p7","0.018252+0.000372*x",0,10000);
TF1 f_1p7_eta_2p1("f_1p7_eta_2p1","-0.008672+0.001016*x",0,10000);

void draw2H(TH2F* h1, TH2F* h2, TF1 f){
    TCanvas c1;
    //h1->Rebin2D(10,10);
    //h2->Rebin2D(10,10);
    h1->Draw("colz");
    //h1->Draw("colz,same");
    f.Draw("same");
    string tit = f.GetName();
    tit+=".pdf";
    c1.SaveAs(("fitQuantiles/"+tit).c_str());
}

void drawH(TH2F* h1, TF1 f){
    TCanvas c1;
    h1->Rebin2D(10,10);
    h1->Draw("colz");
    f.Draw("same");
    string tit = f.GetName();
    tit+=".pdf";
    c1.SaveAs(("fitQuantiles/"+tit).c_str());
}


void cutPtErr(){
    ofstream ofile("fitQuantiles/cut_pTerr_eta.txt");
    string stfile = "outfile_SingleMu_TkOnly_20UL17C_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root";
    string stfile2 = "outfile_Gluino_M-1600_TkOnly_2017_Ih3p47_pTerrCutVar_IasPtSlices_-50eta50_ias50_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0.root";
    TFile* f = new TFile(stfile.c_str());
    TFile* f2 = new TFile(stfile2.c_str());
    
    TH2F* h_eta_0p8 = (TH2F*) f->Get("pT_pTerrOverpT_eta_0p8");
    TH2F* h_0p8_eta_1p7 = (TH2F*) f->Get("pT_pTerrOverpT_0p8_eta_1p7");
    TH2F* h_1p7_eta_2p1 = (TH2F*) f->Get("pT_pTerrOverpT_1p7_eta_2p1");

    TH2F* h_eta_0p8_2 = (TH2F*) f2->Get("pT_pTerrOverpT_eta_0p8");
    TH2F* h_0p8_eta_1p7_2 = (TH2F*) f2->Get("pT_pTerrOverpT_0p8_eta_1p7");
    TH2F* h_1p7_eta_2p1_2 = (TH2F*) f2->Get("pT_pTerrOverpT_1p7_eta_2p1");


    fitQuantiles(ofile, h_eta_0p8, h_eta_0p8_2, 0.5);
    fitQuantiles(ofile, h_0p8_eta_1p7, h_0p8_eta_1p7_2, 0.5);
    fitQuantiles(ofile, h_1p7_eta_2p1, h_1p7_eta_2p1_2, 0.5);
    /*extractLineAtXsigma(ofile, h_eta_0p8, 2);
    extractLineAtXsigma(ofile, h_eta_0p8, 3);
    extractLineAtXsigma(ofile, h_0p8_eta_1p7, 2);
    extractLineAtXsigma(ofile, h_0p8_eta_1p7, 3);
    extractLineAtXsigma(ofile, h_1p7_eta_2p1, 2);
    extractLineAtXsigma(ofile, h_1p7_eta_2p1, 3);*/

    /*draw2H(h_eta_0p8, h_eta_0p8_2, f_eta_0p8);
    draw2H(h_0p8_eta_1p7, h_0p8_eta_1p7_2, f_0p8_eta_1p7);
    draw2H(h_1p7_eta_2p1, h_1p7_eta_2p1_2, f_1p7_eta_2p1);*/

    ofile.close();
     
}
