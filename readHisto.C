#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "HscpCandidates.h"

// neta=120, np=2000, nih=500, nmass=2000


void loadHistograms(region& r, TFile* f, const std::string& regionName, bool bool_rebin=true, int rebineta=1, int rebinp=1, int rebinih=1, int rebinmass=1)
{
    r.eta_p     = (TH2F*)f->Get(("eta_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.eta_p->Rebin2D(rebinp,rebineta);
    r.ih_eta    = (TH2F*)f->Get(("ih_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_eta->Rebin2D(rebineta,rebinih);
    r.ih_p      = (TH2F*)f->Get(("ih_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p->Rebin2D(rebinp,rebinih);
    r.mass      = (TH1F*)f->Get(("massFromTree_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    r.massFrom1DTemplatesEtaBinning      = (TH1F*)f->Get(("massFrom1DTemplatesEtaBinning_"+regionName).c_str())->Clone(); r.massFrom1DTemplatesEtaBinning->Reset(); if(bool_rebin) r.massFrom1DTemplatesEtaBinning->Rebin(rebinmass);
}

TH1F* massFrom2D(const region& r,const std::string& name,int binx=1,int biny=1)
{
    TH2F* h2 = (TH2F*)r.ih_p->Clone();
    h2->Rebin2D(binx,biny);
    int nbin=r.mass->GetNbinsX(); 
    float lowbin=r.mass->GetBinLowEdge(1), upbin=r.mass->GetBinLowEdge(nbin+1);
    TH1F* h1 = (TH1F*) r.mass->Clone("_massFrom2D_"); h1->Reset();
    //TH1F* h1 = new TH1F(("massFrom2D_"+name).c_str(),";Mass [GeV]",nbin,lowbin,upbin);
    for(int i=1;i<h2->GetNbinsX()+1;i++)
    {
        for(int j=1;j<h2->GetNbinsY()+1;j++)
        {
            float mom=h2->GetXaxis()->GetBinCenter(i), dedx=h2->GetYaxis()->GetBinCenter(j);
            h1->Fill(GetMass(mom,dedx,K,C),h2->GetBinContent(i,j));
        }
    }
    return h1;
}

void readHisto()
{
    TFile* ifile = new TFile("outfile_FullStat.root");

    
    region rall;
    region rd;
    bool bool_rebin=true;
    loadHistograms(rall,ifile,"all",bool_rebin,2,2,5,10); //rebin eta,p,ih,mass

//    ifile->Close(); delete ifile;

    invScale(rall.mass);
    std::vector<double> vectOfBins;
    std::vector<double> vectOfBins_P;

    //rebinning(rall.mass,80,vectOfBins);
    rebinning((TH1F*)rall.eta_p->ProjectionX(),500,vectOfBins_P);

    //rall.rebinEtaP(vectOfBins_P);

    rall.VectOfBins_P_ = vectOfBins_P;

    for(int i=0;i<21;i++) vectOfBins.push_back(i*200);

    //for(int i=0;i<vectOfBins.size();i++) std::cout << vectOfBins[i] << std::endl;

    TH1F* mhhh = (TH1F*) rall.mass->Rebin(vectOfBins.size()-1,"variableBins",vectOfBins.data());

    rall.mass = mhhh;


    rall.massFrom1DTemplatesEtaBinning = (TH1F*)rall.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());
   
    rall.fillMassFrom1DTemplatesEtaBinning();

    TH1F* h_massFrom2D = (TH1F*) massFrom2D(rall,"all");

    TFile* ofile = new TFile("analysed.root","RECREATE");

    rall.mass->Write();
    rall.eta_p->Write();

    mhhh->Write();
    rall.massFrom1DTemplatesEtaBinning->Write();
    rall.errMass->Write();

    plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,false,"mass1D_all","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,true,"mass1D_all_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,h_massFrom2D,false,"mass2D_all","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D,true,"mass2D_all_simpleRatio","Observed","Prediction from 2D template")->Write();
    /*plotting(rall.mass,massFrom2D(rall,"all"),false,"mass2D_all","Observed","Prediction from 2D template - p: 10 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,2),false,"mass2D_all_IhP_rebin2_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.2 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,10),false,"mass2D_all_IhP_rebin10_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,1),false,"mass2D_all_IhP_rebin1_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",5,1),false,"mass2D_all_IhP_rebin1_5","Observed","Prediction from 2D template - p: 50 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,1),false,"mass2D_all_IhP_rebin1_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 0.1 MeV/cm")->Write();*/

    rall.ih_p->Write();

    ofile->Close(); delete ofile;

}
