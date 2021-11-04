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
    r.ias_p     = (TH2F*)f->Get(("ias_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_p->Rebin2D(rebinp,rebinih);
    r.ias_pt    = (TH2F*)f->Get(("ias_pt_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_pt->Rebin2D(rebinp,rebinih);
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

    ifstream infile;
    infile.open("./configFile_readHist.txt");
    std::string line;
    std::string filename;
    int rebineta,rebinih,rebinp,rebinmass,thresholdP,thresholdMass;
    bool rebin,varBinsP,varBinsMass;
    while(std::getline(infile,line))
    {
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass >> varBinsP >> varBinsMass >> thresholdP >> thresholdMass;
    }

    std::string outfilename_;
    if(!varBinsMass && !varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_analysed";
    if(varBinsMass && varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_massThresholdVarBins"+to_string(thresholdMass)+"_pThresholdVarBins"+to_string(thresholdP)+"_analysed";
    if(varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_pThresholdVarBins"+to_string(thresholdP)+"_analysed";
    if(varBinsMass) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_massThresholdVarBins"+to_string(thresholdMass)+"_analysed";


    TFile* ifile = new TFile((filename+".root").c_str());
    //TFile* ifile = new TFile("outfile_FullStat.root");
    //TFile* ifile = new TFile("outfile_invIso.root");
    //TFile* ifile = new TFile("outfile_invMET.root");

    
    region rall;
    region rb;
    region rc;
    region rd;
    region rbc;
    
    bool bool_rebin=rebin;
    
    /*loadHistograms(rall,ifile,"all",bool_rebin,2,2,5,10); //rebin eta,p,ih,mass
    loadHistograms(rb,ifile,"regionB",bool_rebin,2,2,5,10); //rebin eta,p,ih,mass
    loadHistograms(rc,ifile,"regionC",bool_rebin,2,2,5,10); //rebin eta,p,ih,mass
    loadHistograms(rd,ifile,"regionD",bool_rebin,2,2,5,10); //rebin eta,p,ih,mass
    */
    
    loadHistograms(rall,ifile,"all",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rb,ifile,"regionB",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rc,ifile,"regionC",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rd,ifile,"regionD",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass



    invScale(rall.mass);
    std::vector<double> vectOfBins;
    std::vector<double> vectOfBins_P;

    if(varBinsMass) rebinning(rall.mass,thresholdMass,vectOfBins);
    if(varBinsP) rebinning((TH1F*)rall.eta_p->ProjectionX(),thresholdP,vectOfBins_P);

    /*rall.eta_p->RebinX(5);
    rb.eta_p->RebinX(5);
    rc.eta_p->RebinX(5);
    rd.eta_p->RebinX(5);*/
   
    rall.fillStdDev();
    rb.fillStdDev();
    rc.fillStdDev();
    rd.fillStdDev();
    
    rall.fillQuantile();
    rb.fillQuantile();
    rc.fillQuantile();
    rd.fillQuantile();

    rall.rebinQuantiles(5);
    rb.rebinQuantiles(5);
    rc.rebinQuantiles(5);
    rd.rebinQuantiles(5);


    //rall.rebinEtaP(vectOfBins_P);

    if(varBinsP)
    {
        rall.VectOfBins_P_ = vectOfBins_P;
        rb.VectOfBins_P_ = vectOfBins_P;
        rc.VectOfBins_P_ = vectOfBins_P;
        rd.VectOfBins_P_ = vectOfBins_P;
    }

    if(varBinsMass)
    {
        rall.mass = (TH1F*) rall.mass->Rebin(vectOfBins.size()-1,"variableBins",vectOfBins.data());
        rall.massFrom1DTemplatesEtaBinning = (TH1F*) rall.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());
    }
   
    rall.fillMassFrom1DTemplatesEtaBinning();

    etaReweighingP(rc.eta_p,rb.eta_p); 

    //rd.eta_p = rc.eta_p;
    //rd.ih_eta = rb.ih_eta;
    

    //rd.mass = (TH1F*) rd.mass->Rebin(vectOfBins.size()-1,"variableBins",vectOfBins.data());

    //rd.massFrom1DTemplatesEtaBinning = (TH1F*) rd.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());

    rd.fillMassFrom1DTemplatesEtaBinning();


    rbc = rd;
    rbc.eta_p = rc.eta_p;
    rbc.ih_eta = rb.ih_eta;

    rbc.massFrom1DTemplatesEtaBinning->Reset();

    rbc.fillMassFrom1DTemplatesEtaBinning();

    TProfile* profXAll = (TProfile*)rall.ias_p->ProfileX();
    TProfile* profXB = (TProfile*)rb.ias_p->ProfileX();
    TProfile* profXC = (TProfile*)rc.ias_p->ProfileX();
    TProfile* profXD = (TProfile*)rd.ias_p->ProfileX();

    TProfile* profYAll = (TProfile*)rall.ias_p->ProfileY();
    TProfile* profYB = (TProfile*)rb.ias_p->ProfileY();
    TProfile* profYC = (TProfile*)rc.ias_p->ProfileY();
    TProfile* profYD = (TProfile*)rd.ias_p->ProfileY();

    TH1F* h_massFrom2D = (TH1F*) massFrom2D(rall,"all");
               
    rall.Mass_errMass = (TH2F*)rall.Mass_errMass->Rebin2D(10,10);

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");
    //TFile* ofile = new TFile("analysed_FullStat.root","RECREATE");
    //TFile* ofile = new TFile("analysed_invIso.root","RECREATE");
    //TFile* ofile = new TFile("analysed_invMET.root","RECREATE");

    rall.mass->Write();
    rall.eta_p->Write();

    rall.massFrom1DTemplatesEtaBinning->Write();
    rall.errMass->Write();
    rall.Mass_errMass->Write();


    rall.stdDevIas_p_y->Write();
    rall.stdDevIas_p->Write();
    profXAll->Write();
    profYAll->Write();
    rall.ias_p->Write();
    rall.ias_pt->Write();


    profXB->Write();
    profYB->Write();

    rc.stdDevIas_p_y->Write();
    rc.stdDevIas_p->Write();
    profXC->Write();
    profYC->Write();
    rc.ias_p->Write();
    rc.ias_pt->Write();
    
    rd.stdDevIas_p_y->Write();
    rd.stdDevIas_p->Write();
    profXD->Write();
    profYD->Write();
    rd.ias_p->Write();
    rd.ias_pt->Write();

    plotting((TH1F*)rc.quantile01Ias_p,(TH1F*)rd.quantile01Ias_p,true,"quantile01_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile10Ias_p,(TH1F*)rd.quantile10Ias_p,true,"quantile10_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile30Ias_p,(TH1F*)rd.quantile30Ias_p,true,"quantile30_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile50Ias_p,(TH1F*)rd.quantile50Ias_p,true,"quantile50_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile70Ias_p,(TH1F*)rd.quantile70Ias_p,true,"quantile70_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile90Ias_p,(TH1F*)rd.quantile90Ias_p,true,"quantile90_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.quantile99Ias_p,(TH1F*)rd.quantile99Ias_p,true,"quantile99_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.stdDevIas_p,(TH1F*)rd.stdDevIas_p,true,"stddev_ias_p_c_d","region C","region D")->Write();
    //plotting((TH1F*)profC,(TH1F*)profD,true,"profile_ias_p_c_d","region C","region D")->Write();
    plotting((TH1F*)rc.eta_p->ProjectionX(),(TH1F*)rd.eta_p->ProjectionX(),true,"p_c_d","region C","region D")->Write();

    plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,false,"mass1D_all","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,true,"mass1D_all_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,h_massFrom2D,false,"mass2D_all","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D,true,"mass2D_all_simpleRatio","Observed","Prediction from 2D template")->Write();
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,false,"mass1D_regionD","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,true,"mass1D_regionD_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC","Observed","Prediction from 1D templates in B and C")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,true,"mass1D_regionBC_simpleRatio","Observed","Prediction from 1D templates in B and C")->Write();


    /*plotting(rall.mass,massFrom2D(rall,"all"),false,"mass2D_all","Observed","Prediction from 2D template - p: 10 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,2),false,"mass2D_all_IhP_rebin2_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.2 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,10),false,"mass2D_all_IhP_rebin10_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,1),false,"mass2D_all_IhP_rebin1_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",5,1),false,"mass2D_all_IhP_rebin1_5","Observed","Prediction from 2D template - p: 50 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,1),false,"mass2D_all_IhP_rebin1_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 0.1 MeV/cm")->Write();*/

    rall.ih_p->Write();

    std::cout << 600*GetMassErr(1000,20,3.8,0.1,600,K,C) << std::endl;
    std::cout << 600*GetMassErr(1000,4,3.8,0.1,600,K,C) << std::endl;
    std::cout << 600*GetMassErr(1000,100,3.8,0.5,600,K,C) << std::endl;
    ofile->Close(); delete ofile;

}
