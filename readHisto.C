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
    r.ih_p_eta  = (TH3F*)f->Get(("ih_p_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p_eta->Rebin3D(rebineta,rebinp,rebinih);
    r.eta_p     = (TH2F*)f->Get(("eta_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.eta_p->Rebin2D(rebinp,rebineta);
    r.ih_eta    = (TH2F*)f->Get(("ih_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_eta->Rebin2D(rebineta,rebinih);
    r.ih_p      = (TH2F*)f->Get(("ih_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p->Rebin2D(rebinp,rebinih);
    r.ias_p     = (TH2F*)f->Get(("ias_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_p->Rebin2D(rebinp,rebinih);
    r.ias_pt    = (TH2F*)f->Get(("ias_pt_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_pt->Rebin2D(rebinp,rebinih);
    r.mass      = (TH1F*)f->Get(("massFromTree_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    //r.mass      = TH1F*((TH2F*)f->Get(("Mass_errMass_"+regionName).c_str()))->ProjectionX(); if(bool_rebin) r.mass->Rebin(rebinmass);
    //r.massFrom1DTemplatesEtaBinning      = (TH1F*)(((TH2F*)f->Get(("Mass_errMass_"+regionName).c_str()))->ProjectionX())->Clone(); r.massFrom1DTemplatesEtaBinning->Reset(); if(bool_rebin) r.massFrom1DTemplatesEtaBinning->Rebin(rebinmass);
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

void massNormalisation(TH1F* h, const float& normalisation){

    for(int k=0;k<h->GetNbinsX()+1;k++){
        h->SetBinContent(k,h->GetBinContent(k)*normalisation);
        h->SetBinError(k,h->GetBinError(k)*normalisation);
    }

}

void saveHistoRatio(TH1F* h1,TH1F* h2,std::string st1,std::string st2,std::string st3,bool rebin=false){
    h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    h2->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    h1->SetName(st1.c_str());
    h2->SetName(st2.c_str());
    if(rebin){
        h1=rebinHisto(h1);
        h2=rebinHisto(h2);
    }
    h1->Write();
    h2->Write();
    TH1F* R = (TH1F*) ratioIntegral(h2,h1)->Clone();
    if(rebin) st3+="_rebinned";
    R->SetName(st3.c_str());
    R->Write();
}

TH1F meanHistoPE(std::vector<TH1F> vPE, float systErr=0.2){
    TH1F h(vPE[0]);
    //h = (TH1F)vPE[0].Clone(); 
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    for(int i=0;i<h.GetNbinsX()+1;i++){
        float mean=0, err=0;
        for(int pe=0;pe<vPE.size();pe++){
            mean += vPE[pe].GetBinContent(i);
        }
        mean /= vPE.size();
        for(int pe=0;pe<vPE.size();pe++){
            err += pow(mean - vPE[pe].GetBinContent(i),2);
        }
        err = sqrt(err/(vPE.size()-1));
        h.SetBinContent(i,mean);
        h.SetBinError(i,err);
    }
    for(int i=0;i<h.GetNbinsX()+1;i++){
        float err = sqrt( pow(h.GetBinError(i),2) + pow(h.GetBinContent(i)*systErr,2) );
        h.SetBinError(i,err);
    }
    
    return h;
}

void poissonHisto(TH2F &h){
    for(int i=0;i<h.GetNbinsX()+1;i++){
        for(int j=0;j<h.GetNbinsY()+1;j++){
            std::cout << "content: " << h.GetBinContent(i,j) << std::endl;
            h.SetBinContent(i,j,RNG->Poisson(h.GetBinContent(i,j)));
        }
    }
}

std::ofstream ofile("normalisation_vs_regionD.txt");

void doAll(region& b, region& c, region& bc, region& a, region& d, std::string st, int nPE=100, float systErr=0.2){
    std::vector<TH1F> vPE;
    std::cout << st << std::endl;
    d.mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    for(int pe=0;pe<nPE;pe++){
        std::cout << "pe: " << pe << std::endl;
        TH2F a_ih_eta(*a.ih_eta);
        TH2F b_ih_eta(*b.ih_eta);
        TH2F c_ih_eta(*c.ih_eta);
        TH2F b_eta_p(*b.eta_p);
        TH2F c_eta_p(*c.eta_p);
        poissonHisto(a_ih_eta);
        poissonHisto(b_ih_eta);
        poissonHisto(c_ih_eta);
        poissonHisto(b_eta_p);
        poissonHisto(c_eta_p);
        etaReweighingP(&c_eta_p,&b_eta_p);
        bc.eta_p = &c_eta_p;    bc.ih_eta = &b_ih_eta;
        /*float A = RNG->Poisson(a.ih_eta->GetEntries());
        float B = RNG->Poisson(b.ih_eta->GetEntries());
        float C = RNG->Poisson(c.ih_eta->GetEntries());*/
        float A = a_ih_eta.Integral();
        float B = b_ih_eta.Integral();
        float C = c_ih_eta.Integral();
        float normalisationABC = B*C/A;
        ofile << st << " PE: " << pe << " A: " << A << " B: " << B << " C: " << C << " --> "<< " normalisation: " << normalisationABC << " regionD: " << d.mass->GetEntries() << std::endl;
        bc.fillMassFrom1DTemplatesEtaBinning();
        scale(bc.massFrom1DTemplatesEtaBinning);
        massNormalisation(bc.massFrom1DTemplatesEtaBinning,normalisationABC);
        vPE.push_back(*bc.massFrom1DTemplatesEtaBinning);
    }
    TH1F h_temp = meanHistoPE(vPE, systErr);
    if(nPE>1) bc.massFrom1DTemplatesEtaBinning = &h_temp;
   
    bool rebin=false;

    saveHistoRatio(d.mass,bc.massFrom1DTemplatesEtaBinning,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),rebin=false);
    saveHistoRatio(d.mass,bc.massFrom1DTemplatesEtaBinning,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),rebin=true);
    
    overflowLastBin(d.mass);
    overflowLastBin(bc.massFrom1DTemplatesEtaBinning);
    
    plotting(d.mass,bc.massFrom1DTemplatesEtaBinning,false,("mass1D_regionBC_"+st).c_str(),"Observed","Prediction",rebin=false)->Write();
    plotting(d.mass,bc.massFrom1DTemplatesEtaBinning,false,("mass1D_regionBC_"+st).c_str(),"Observed","Prediction",rebin=true)->Write();
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
        ss >> filename >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass >> varBinsMass >> varBinsP >> thresholdMass >> thresholdP;
    }

    std::string outfilename_;
    if(!varBinsMass && !varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_analysed";
    else if(varBinsMass && varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_massThresholdVarBins"+to_string(thresholdMass)+"_pThresholdVarBins"+to_string(thresholdP)+"_analysed";
    else if(varBinsP) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_pThresholdVarBins"+to_string(thresholdP)+"_analysed";
    else if(varBinsMass) outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_massThresholdVarBins"+to_string(thresholdMass)+"_analysed";


    std::cout << outfilename_ << std::endl;

    TFile* ifile = new TFile((filename+".root").c_str());

    
    region rall;
    region ra;
    region rb;
    region rc;
    region rd;
    region rbc;
    region rdb;
    region rdc;
/*    region rb_boundedIas;
    region rc_boundedIas;
    region rd_boundedIas;
    region rbc_boundedIas;
    region rb_boundedPt;
    region rc_boundedPt;
    region rd_boundedPt;
    region rbc_boundedPt;*/

    region rb_40;
    region rb_50;
    region rb_60;
    region rb_70;
    region rb_80;
    region rb_90;
    region rb_50_90;
    region rb_50_100;

/*    region rc_40_pt;
    region rc_50_pt;
    region rc_60_pt;
    region rc_70_pt;
    region rc_80_pt;
    region rc_90_pt;*/

    region rbc_40;
    region rbc_50;
    region rbc_60;
    region rbc_70;
    region rbc_80;
    region rbc_90;
    region rbc_50_90;
    region rbc_50_100;
    
/*    region rbc_40_pt;
    region rbc_50_pt;
    region rbc_60_pt;
    region rbc_70_pt;
    region rbc_80_pt;
    region rbc_90_pt;*/

    region rd_40;
    region rd_50;
    region rd_60;
    region rd_70;
    region rd_80;
    region rd_90;
    region rd_50_90;
    region rd_50_100;
    
/*    region rd_40_pt;
    region rd_50_pt;
    region rd_60_pt;
    region rd_70_pt;
    region rd_80_pt;
    region rd_90_pt;*/

//    region ra_40;
//    region ra_40_pt;
    region ra_med;
//    region ra_med_pt;
//    region rc_40;
    region rc_med;
//    region rb_40_pt;
//    region rb_med_pt;
    
    bool bool_rebin=rebin;

    region rA_sc1;
    region rA_sc2;
    region rA_sc3;
    region rA_sc4;
    region rA_sc5;
    region rA_sc6;
    region rA_sc7;
    region rA_sc8;

    region rB_sc1;
    region rB_sc2;
    region rB_sc3;
    region rB_sc4;
    region rB_sc5;
    region rB_sc6;
    region rB_sc7;
    region rB_sc8;

    region rC_sc1;
    region rC_sc2;
    region rC_sc3;
    region rC_sc4;
    region rC_sc5;
    region rC_sc6;
    region rC_sc7;
    region rC_sc8;
 
    region rD_sc1;
    region rD_sc2;
    region rD_sc3;
    region rD_sc4;  
    region rD_sc5;  
    region rD_sc6;
    region rD_sc7;
    region rD_sc8;

    region rBC_sc1;
    region rBC_sc2;
    region rBC_sc3;
    region rBC_sc4;
    region rBC_sc5;
    region rBC_sc6;
    region rBC_sc7;
    region rBC_sc8;
    
    loadHistograms(rA_sc1,ifile,"regA_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc2,ifile,"regA_sc2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc3,ifile,"regA_sc3",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc4,ifile,"regA_sc4",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc5,ifile,"regA_sc5",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc6,ifile,"regA_sc6",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc7,ifile,"regA_sc7",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rA_sc8,ifile,"regA_sc8",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rB_sc1,ifile,"regB_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_sc2,ifile,"regB_sc2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_sc3,ifile,"regB_sc3",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_sc4,ifile,"regB_sc4",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rB_sc5,ifile,"regB_sc5",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rB_sc6,ifile,"regB_sc6",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rB_sc7,ifile,"regB_sc7",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rB_sc8,ifile,"regB_sc8",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    loadHistograms(rC_sc1,ifile,"regC_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_sc2,ifile,"regC_sc2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_sc3,ifile,"regC_sc3",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_sc4,ifile,"regC_sc4",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rC_sc5,ifile,"regC_sc5",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rC_sc6,ifile,"regC_sc6",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rC_sc7,ifile,"regC_sc7",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rC_sc8,ifile,"regC_sc8",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    loadHistograms(rD_sc1,ifile,"regD_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_sc2,ifile,"regD_sc2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_sc3,ifile,"regD_sc3",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_sc4,ifile,"regD_sc4",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rD_sc5,ifile,"regD_sc5",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rD_sc6,ifile,"regD_sc6",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rD_sc7,ifile,"regD_sc7",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rD_sc8,ifile,"regD_sc8",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    loadHistograms(rBC_sc1,ifile,"regD_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rBC_sc2,ifile,"regD_sc2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rBC_sc3,ifile,"regD_sc3",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rBC_sc4,ifile,"regD_sc4",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_sc5,ifile,"regD_sc5",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_sc6,ifile,"regD_sc6",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_sc7,ifile,"regD_sc7",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_sc8,ifile,"regD_sc8",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    
    loadHistograms(rall,ifile,"all",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(ra,ifile,"regionA",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rb,ifile,"regionB",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rc,ifile,"regionC",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rd,ifile,"regionD",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rbc,ifile,"regionD",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rdb,ifile,"regionD",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rdc,ifile,"regionD",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
/*    loadHistograms(rb_boundedIas,ifile,"regionB_boundedIas",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rc_boundedIas,ifile,"regionC",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rd_boundedIas,ifile,"regionD_boundedIas",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rbc_boundedIas,ifile,"regionD_boundedIas",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rb_boundedPt,ifile,"regionB",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rc_boundedPt,ifile,"regionC_boundedPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rd_boundedPt,ifile,"regionD_boundedPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass
    loadHistograms(rbc_boundedPt,ifile,"regionD_boundedPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); //rebin eta,p,ih,mass*/

//    loadHistograms(rb_40,ifile,"regionB_40",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50,ifile,"regionB_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60,ifile,"regionB_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70,ifile,"regionB_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80,ifile,"regionB_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90,ifile,"regionB_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50_90,ifile,"regionB_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50_100,ifile,"regionB_50_100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
/*    loadHistograms(rc_40_pt,ifile,"regionC_40pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_50_pt,ifile,"regionC_50pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_60_pt,ifile,"regionC_60pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_70_pt,ifile,"regionC_70pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_80_pt,ifile,"regionC_80pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_90_pt,ifile,"regionC_90pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); */

//    loadHistograms(rbc_40,ifile,"regionD_40",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50_90,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50_100,ifile,"regionD_50_100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
/*    loadHistograms(rbc_40_pt,ifile,"regionD_40pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50_pt,ifile,"regionD_50pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60_pt,ifile,"regionD_60pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70_pt,ifile,"regionD_70pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80_pt,ifile,"regionD_80pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90_pt,ifile,"regionD_90pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); */

//    loadHistograms(rd_40,ifile,"regionD_40",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50_90,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50_100,ifile,"regionD_50_100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
/*    loadHistograms(rd_40_pt,ifile,"regionD_40pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50_pt,ifile,"regionD_50pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60_pt,ifile,"regionD_60pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70_pt,ifile,"regionD_70pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80_pt,ifile,"regionD_80pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90_pt,ifile,"regionD_90pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); */
    
//    loadHistograms(ra_40,ifile,"regionA_40",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
//    loadHistograms(ra_40_pt,ifile,"regionA_40pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(ra_med,ifile,"regionA_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
//    loadHistograms(ra_med_pt,ifile,"regionA_med_pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
//    loadHistograms(rc_40,ifile,"regionC_40",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_med,ifile,"regionC_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
//    loadHistograms(rb_40_pt,ifile,"regionB_40pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
//    loadHistograms(rb_med_pt,ifile,"regionB_med_pt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    region rA_ias50_eta08;
    region rC_ias50_eta08;
    region rB_50ias90_eta08;
    region rD_50ias90_eta08;
    region rBC_50ias90_eta08;
    loadHistograms(rA_ias50_eta08,ifile,"regA_ias50_eta08",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_ias50_eta08,ifile,"regC_ias50_eta08",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_50ias90_eta08,ifile,"regB_50ias90_eta08",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_50ias90_eta08,ifile,"regD_50ias90_eta08",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_50ias90_eta08,ifile,"regD_50ias90_eta08",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    region rA_ias50_08eta17;
    region rC_ias50_08eta17;
    region rB_50ias90_08eta17;
    region rD_50ias90_08eta17;
    region rBC_50ias90_08eta17;
    loadHistograms(rA_ias50_08eta17,ifile,"regA_ias50_08eta17",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_ias50_08eta17,ifile,"regC_ias50_08eta17",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_50ias90_08eta17,ifile,"regB_50ias90_08eta17",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_50ias90_08eta17,ifile,"regD_50ias90_08eta17",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_50ias90_08eta17,ifile,"regD_50ias90_08eta17",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    region rA_ias50_17eta21;
    region rC_ias50_17eta21;
    region rB_50ias90_17eta21;
    region rD_50ias90_17eta21;
    region rBC_50ias90_17eta21;
    loadHistograms(rA_ias50_17eta21,ifile,"regA_ias50_17eta21",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_ias50_17eta21,ifile,"regC_ias50_17eta21",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_50ias90_17eta21,ifile,"regB_50ias90_17eta21",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_50ias90_17eta21,ifile,"regD_50ias90_17eta21",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_50ias90_17eta21,ifile,"regD_50ias90_17eta21",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    
    region rA_005ias01;
    region rB_005ias01;
    region rC_005ias01;
    region rD_005ias01;
    region rBC_005ias01;
    loadHistograms(rA_005ias01,ifile,"regA_005ias01",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_005ias01,ifile,"regB_005ias01",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_005ias01,ifile,"regC_005ias01",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_005ias01,ifile,"regD_005ias01",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_005ias01,ifile,"regD_005ias01",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    region rA_005ias015;
    region rB_005ias015;
    region rC_005ias015;
    region rD_005ias015;
    region rBC_005ias015;
    loadHistograms(rA_005ias015,ifile,"regA_005ias015",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rB_005ias015,ifile,"regB_005ias015",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rC_005ias015,ifile,"regC_005ias015",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rD_005ias015,ifile,"regD_005ias015",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rBC_005ias015,ifile,"regD_005ias015",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    std::cout << "Regions loaded" << std::endl;

    TH2F* h_cross1D_all = (TH2F*) rall.ih_p->Clone("cross1D"); h_cross1D_all->Reset();
    TH2F* mapDiff = (TH2F*) rall.ih_p->Clone("differencesMap"); mapDiff->Reset();

    invScale(rall.mass);
    std::vector<double> vectOfBins;
    std::vector<double> vectOfBins_P;

    if(varBinsMass) rebinning(rall.mass,thresholdMass,vectOfBins);
    //if(varBinsP) rebinning((TH1F*)rall.eta_p->ProjectionX(),thresholdP,vectOfBins_P);
    if(varBinsP) rebinning((TH1F*)rd.eta_p->ProjectionX(),thresholdP,vectOfBins_P);



    /*rall.fillStdDev();
    rb.fillStdDev();
    rc.fillStdDev();
    rd.fillStdDev();
    
    rall.fillQuantile();
    rb.fillQuantile();
    rc.fillQuantile();
    rd.fillQuantile();*/

    /*rall.rebinQuantiles(5);
    rb.rebinQuantiles(5);
    rc.rebinQuantiles(5);
    rd.rebinQuantiles(5);*/
   
    //crossHistos(h_cross1D_all,(TH1F*)rall.ih_p->ProjectionX(),(TH1F*)rall.ih_p->ProjectionY("",0,rall.ih_p->GetXaxis()->FindBin(150)));
    //crossHistos(h_cross1D_all,(TH1F*)rall.ih_p->ProjectionX(),(TH1F*)rall.ih_p->ProjectionY("",0,10));
    //crossHistos(h_cross1D_all,(TH1F*)rall.eta_p->ProjectionX(),(TH1F*)rall.ih_eta->ProjectionY());
    //crossHistosEtaBinning(h_cross1D_all,rall.eta_p,rall.ih_eta);
    //mapOfDifferences(mapDiff,rall.ih_p,h_cross1D_all);
    
    /*rall.cross1D();
    rb.cross1D();
    rc.cross1D();
    rd.cross1D();*/

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
        rd.mass = (TH1F*) rd.mass->Rebin(vectOfBins.size()-1,"variableBins",vectOfBins.data());
        rall.massFrom1DTemplatesEtaBinning = (TH1F*) rall.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());
        rd.massFrom1DTemplatesEtaBinning = (TH1F*) rd.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());
    }
/*    
    etaReweighingP(rc.eta_p,rb.eta_p); 
    //etaReweighingP(rc_boundedIas.eta_p,rb_boundedIas.eta_p); 
    //etaReweighingP(rc_boundedPt.eta_p,rb_boundedPt.eta_p); 
 
    etaReweighingP(rc_40.eta_p,rb_40.eta_p);
    etaReweighingP(rc_40.eta_p,rb_50.eta_p);
    etaReweighingP(rc_40.eta_p,rb_60.eta_p);
    etaReweighingP(rc_40.eta_p,rb_70.eta_p);
    etaReweighingP(rc_40.eta_p,rb_80.eta_p);
    etaReweighingP(rc_40.eta_p,rb_90.eta_p);

    //rbc = rd;
    rbc.eta_p = rc.eta_p;
    rbc.ih_eta = rb.ih_eta;

    rbc_boundedIas.eta_p = rc_boundedIas.eta_p;
    rbc_boundedIas.ih_eta = rb_boundedIas.ih_eta;

    rbc_boundedPt.eta_p = rc_boundedPt.eta_p;
    rbc_boundedPt.ih_eta = rb_boundedPt.ih_eta;
    
    rdc.eta_p = rc.eta_p;
    rdb.ih_eta = rb.ih_eta;


    rbc_40.eta_p = rc_40.eta_p;    rbc_40.ih_eta = rb_40.ih_eta;
    rbc_50.eta_p = rc_40.eta_p;    rbc_50.ih_eta = rb_50.ih_eta;
    rbc_60.eta_p = rc_40.eta_p;    rbc_60.ih_eta = rb_60.ih_eta;
    rbc_70.eta_p = rc_40.eta_p;    rbc_70.ih_eta = rb_70.ih_eta;
    rbc_80.eta_p = rc_40.eta_p;    rbc_80.ih_eta = rb_80.ih_eta;
    rbc_90.eta_p = rc_40.eta_p;    rbc_90.ih_eta = rb_90.ih_eta;

    rall.fillMassFrom1DTemplatesEtaBinning();

    float normalisationABC = rb.ih_eta->GetEntries()*rc.ih_eta->GetEntries()/ra.ih_eta->GetEntries();
    
    float normalisationABC_40 = rb_40.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    float normalisationABC_50 = rb_50.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    float normalisationABC_60 = rb_60.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    float normalisationABC_70 = rb_70.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    float normalisationABC_80 = rb_80.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    float normalisationABC_90 = rb_90.ih_eta->GetEntries()*rc_40.ih_eta->GetEntries()/ra_40.ih_eta->GetEntries();
    std::cout << " normaABD40: " << normalisationABC_40 << std::endl;
    std::cout << " normaABD50: " << normalisationABC_50 << std::endl;
    std::cout << " normaABD60: " << normalisationABC_60 << std::endl;
    std::cout << " normaABD70: " << normalisationABC_70 << std::endl;
    std::cout << " normaABD80: " << normalisationABC_80 << std::endl;
    std::cout << " normaABD90: " << normalisationABC_90 << std::endl;
//    std::cout << "normalisationABC: " << normalisationABC << std::endl;

    //rd.mass = (TH1F*) rd.mass->Rebin(vectOfBins.size()-1,"variableBins",vectOfBins.data());

    //rd.massFrom1DTemplatesEtaBinning = (TH1F*) rd.massFrom1DTemplatesEtaBinning->Rebin(vectOfBins.size()-1,"",vectOfBins.data());

    rd.fillMassFrom1DTemplatesEtaBinning();
    rbc.fillMassFrom1DTemplatesEtaBinning();

    rbc_40.fillMassFrom1DTemplatesEtaBinning();
    rbc_50.fillMassFrom1DTemplatesEtaBinning();
    rbc_60.fillMassFrom1DTemplatesEtaBinning();
    rbc_70.fillMassFrom1DTemplatesEtaBinning();
    rbc_80.fillMassFrom1DTemplatesEtaBinning();
    rbc_90.fillMassFrom1DTemplatesEtaBinning();
*/
    


    //rbc.fillMassFrom1DTemplatesEtaBinning(normalisationABC);

    /*rd_boundedIas.fillMassFrom1DTemplatesEtaBinning();
    rbc_boundedIas.fillMassFrom1DTemplatesEtaBinning();
    rd_boundedPt.fillMassFrom1DTemplatesEtaBinning();
    rbc_boundedPt.fillMassFrom1DTemplatesEtaBinning();
*/
    rdc.fillMassFrom1DTemplatesEtaBinning();
    rdb.fillMassFrom1DTemplatesEtaBinning();

    TProfile* ihprofXAll = (TProfile*)rall.ih_p->ProfileX();
    TProfile* ihprofXB = (TProfile*)rb.ih_p->ProfileX();
    TProfile* ihprofXC = (TProfile*)rc.ih_p->ProfileX();
    TProfile* ihprofXD = (TProfile*)rd.ih_p->ProfileX();


    TProfile* profXAll = (TProfile*)rall.ias_p->ProfileX();
    TProfile* profXB = (TProfile*)rb.ias_p->ProfileX();
    TProfile* profXC = (TProfile*)rc.ias_p->ProfileX();
    TProfile* profXD = (TProfile*)rd.ias_p->ProfileX();

    TProfile* profYAll = (TProfile*)rall.ias_p->ProfileY();
    TProfile* profYB = (TProfile*)rb.ias_p->ProfileY();
    TProfile* profYC = (TProfile*)rc.ias_p->ProfileY();
    TProfile* profYD = (TProfile*)rd.ias_p->ProfileY();

    TH1F* h_massFrom2D = (TH1F*) massFrom2D(rall,"all");
    TH1F* h_massFrom2D_D = (TH1F*) massFrom2D(rd,"regD");
               
    //rall.Mass_errMass = (TH2F*)rall.Mass_errMass->Rebin2D(10,10);

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving... " << std::endl;

    float systematicErr = 0.2;
    int pseudoExperiments = 100;

    doAll(rB_sc1,rC_sc1,rBC_sc1,rA_sc1,rD_sc1,"sc1",pseudoExperiments,systematicErr);
    //doAll(rB_sc2,rC_sc2,rBC_sc2,rA_sc2,rD_sc2,"sc2",100,systematicErr);
    doAll(rB_sc3,rC_sc3,rBC_sc3,rA_sc3,rD_sc3,"sc3",pseudoExperiments,systematicErr);
    //doAll(rB_sc4,rC_sc4,rBC_sc4,rA_sc4,rD_sc4,"sc4",100,systematicErr);
    doAll(rB_sc5,rC_sc5,rBC_sc5,rA_sc5,rD_sc5,"sc5",pseudoExperiments,systematicErr);
    doAll(rB_sc6,rC_sc6,rBC_sc6,rA_sc6,rD_sc6,"sc6",pseudoExperiments,systematicErr);
    doAll(rB_sc7,rC_sc7,rBC_sc7,rA_sc7,rD_sc7,"sc7",pseudoExperiments,systematicErr);
    doAll(rB_sc8,rC_sc8,rBC_sc8,rA_sc8,rD_sc8,"sc8",pseudoExperiments,systematicErr);

    /*doAll(rb,rc,rbc,ra,rd,"pt60_ias025");

    doAll(rb_50_90,rc_med,rbc_50_90,ra_med,rd_50_90,"50_90");
    doAll(rb_50_100,rc_med,rbc_50_100,ra_med,rd_50_100,"50_100");
    //doAll(rb_50_90,rc_med,rbc_50_90,ra_med,rd_50_90,"50_90",1);

    //doAll(rb_40,rc_40,rbc_40,ra_40,rd_40,"40",);
   
    doAll(rb_50,rc_med,rbc_50,ra_med,rd_50,"50ias_50");
    //doAll(rb_50,rc_med,rbc_50,ra_med,rd_50,"50ias_50",1);
    doAll(rb_60,rc_med,rbc_60,ra_med,rd_60,"60ias_50");
    //doAll(rb_60,rc_med,rbc_60,ra_med,rd_60,"60ias_50",1);
    doAll(rb_70,rc_med,rbc_70,ra_med,rd_70,"70ias_50");
    //doAll(rb_70,rc_med,rbc_70,ra_med,rd_70,"70ias_50",1);
    doAll(rb_80,rc_med,rbc_80,ra_med,rd_80,"80ias_50");
    //doAll(rb_80,rc_med,rbc_80,ra_med,rd_80,"80ias_50",1);
    doAll(rb_90,rc_med,rbc_90,ra_med,rd_90,"90ias_50");

    doAll(rB_50ias90_eta08, rC_ias50_eta08, rBC_50ias90_eta08, rA_ias50_eta08, rD_50ias90_eta08, "50ias90_eta08"); 
    doAll(rB_50ias90_08eta17, rC_ias50_08eta17, rBC_50ias90_08eta17, rA_ias50_08eta17, rD_50ias90_08eta17, "50ias90_08eta17"); 
    doAll(rB_50ias90_17eta21, rC_ias50_17eta21, rBC_50ias90_17eta21, rA_ias50_17eta21, rD_50ias90_17eta21, "50ias90_17eta21"); 

    doAll(rB_005ias01, rC_005ias01, rBC_005ias01, rA_005ias01, rD_005ias01, "005ias01");
    doAll(rB_005ias015, rC_005ias015, rBC_005ias015, rA_005ias015, rD_005ias015, "005ias015");*/
    
    /*doAll(rb_40,rc_40,rbc_40,ra_40,rd_40,"40ias_40");
    doAll(rb_50,rc_40,rbc_50,ra_40,rd_50,"50ias_40");
    doAll(rb_60,rc_40,rbc_60,ra_40,rd_60,"60ias_40");
    doAll(rb_70,rc_40,rbc_70,ra_40,rd_70,"70ias_40");
    doAll(rb_80,rc_40,rbc_80,ra_40,rd_80,"80ias_40");
    doAll(rb_90,rc_40,rbc_90,ra_40,rd_90,"90ias_40");*/
    
/*    doAll(rb_med_pt,rc_50_pt,rbc_50_pt,ra_med_pt,rd_50_pt,"50pt_50");
    doAll(rb_med_pt,rc_60_pt,rbc_60_pt,ra_med_pt,rd_60_pt,"60pt_50");
    doAll(rb_med_pt,rc_70_pt,rbc_70_pt,ra_med_pt,rd_70_pt,"70pt_50");
    doAll(rb_med_pt,rc_80_pt,rbc_80_pt,ra_med_pt,rd_80_pt,"80pt_50");
    doAll(rb_med_pt,rc_90_pt,rbc_90_pt,ra_med_pt,rd_90_pt,"90pt_50");
    
    doAll(rb_40_pt,rc_40_pt,rbc_40_pt,ra_40_pt,rd_40_pt,"40pt_40");
    doAll(rb_40_pt,rc_50_pt,rbc_50_pt,ra_40_pt,rd_50_pt,"50pt_40");
    doAll(rb_40_pt,rc_60_pt,rbc_60_pt,ra_40_pt,rd_60_pt,"60pt_40");
    doAll(rb_40_pt,rc_70_pt,rbc_70_pt,ra_40_pt,rd_70_pt,"70pt_40");
    doAll(rb_40_pt,rc_80_pt,rbc_80_pt,ra_40_pt,rd_80_pt,"80pt_40");
    doAll(rb_40_pt,rc_90_pt,rbc_90_pt,ra_40_pt,rd_90_pt,"90pt_40");*/

/*    rc.eta_p->Write();
    rd.eta_p->Write();
    //rc_boundedIas.eta_p->SetName("eta_p_regionC_boundedIas");
    //rc_boundedIas.eta_p->Write();
    //rd_boundedIas.eta_p->Write();

    rall.mass->Write();
    rall.eta_p->Write();

    rall.massFrom1DTemplatesEtaBinning->Write();
    //rall.errMass->Write();
    //rall.Mass_errMass->Write();


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

    rall.stdDevIh_p->Write();
    rall.quantile01Ih_p->Write();
    rall.quantile10Ih_p->Write();
    rall.quantile30Ih_p->Write();
    rall.quantile50Ih_p->Write();
    rall.quantile70Ih_p->Write();
    rall.quantile90Ih_p->Write();
    rall.quantile99Ih_p->Write();

    ihprofXAll->Write();
    ihprofXB->Write();
    ihprofXC->Write();
    ihprofXD->Write();
    
    rall.ih_p->Write();
    h_cross1D_all->Write();
    mapDiff->Write();
    rall.cross1Dtemplates->Write();
    rb.cross1Dtemplates->Write();
    rc.cross1Dtemplates->Write();
    rd.cross1Dtemplates->Write();

    rall.ih_used->Write();
    rd.ih_used->Write();
    //rall.mapM800->Write();
    rd.mapM800->Write();

    rd.momentumDistribM1000->Write();
    rd.dedxDistribM1000->Write();
*/
/*
    std::cout << "here1" << std::endl;

    //scale(rd.mass);
    scale(rd.massFrom1DTemplatesEtaBinning);
    //rbc.massFrom1DTemplatesEtaBinning->Scale(normalisationABC);
    scale(rbc.massFrom1DTemplatesEtaBinning);
    scale(h_massFrom2D_D);

    scale(rbc_40.massFrom1DTemplatesEtaBinning);
    scale(rbc_50.massFrom1DTemplatesEtaBinning);
    scale(rbc_60.massFrom1DTemplatesEtaBinning);
    scale(rbc_70.massFrom1DTemplatesEtaBinning);
    scale(rbc_80.massFrom1DTemplatesEtaBinning);
    scale(rbc_90.massFrom1DTemplatesEtaBinning);

    std::cout << "here2" << std::endl;
    for(int k=0;k<rbc.massFrom1DTemplatesEtaBinning->GetNbinsX()+1;k++){
        rbc.massFrom1DTemplatesEtaBinning->SetBinContent(k,rbc.massFrom1DTemplatesEtaBinning->GetBinContent(k)*normalisationABC);
    }

    massNormalisation(rbc_40.massFrom1DTemplatesEtaBinning,normalisationABC_40);
    massNormalisation(rbc_50.massFrom1DTemplatesEtaBinning,normalisationABC_50);
    massNormalisation(rbc_60.massFrom1DTemplatesEtaBinning,normalisationABC_60);
    massNormalisation(rbc_70.massFrom1DTemplatesEtaBinning,normalisationABC_70);
    massNormalisation(rbc_80.massFrom1DTemplatesEtaBinning,normalisationABC_80);
    massNormalisation(rbc_90.massFrom1DTemplatesEtaBinning,normalisationABC_90);
    std::cout << "here2" << std::endl;

    scale(rd_boundedIas.mass);
    scale(rd_boundedIas.massFrom1DTemplatesEtaBinning);
    scale(rbc_boundedIas.massFrom1DTemplatesEtaBinning);
    scale(rd_boundedPt.mass);
    scale(rd_boundedPt.massFrom1DTemplatesEtaBinning);
    scale(rbc_boundedPt.massFrom1DTemplatesEtaBinning);
    scale(rdb.massFrom1DTemplatesEtaBinning);
    scale(rdc.massFrom1DTemplatesEtaBinning);

//    std::cout << "D integral: " << rd.mass->Integral() << " entries: " << rd.mass->GetEntries() << std::endl;
//    std::cout << "BC integral: " << rbc.massFrom1DTemplatesEtaBinning->Integral() << " entries: " << rbc.massFrom1DTemplatesEtaBinning->GetEntries() << std::endl;
    
    rd.mass->SetName("mass_obs");
    rd.mass->Write();
    
    h_massFrom2D_D->SetName("mass_pred2D");
    h_massFrom2D_D->Write();

    TH1F* h2DR = (TH1F*) ratioIntegral(h_massFrom2D_D,rd.mass)->Clone();
    h2DR->SetName("mass_pred2DR");
    h2DR->Write();

    rd.massFrom1DTemplatesEtaBinning->SetName("mass_pred1D");
    rd.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DR = (TH1F*) ratioIntegral(rd.massFrom1DTemplatesEtaBinning,rd.mass)->Clone();
    h1DR->SetName("mass_pred1DR");
    h1DR->Write();

    rbc.massFrom1DTemplatesEtaBinning->SetName("mass_predBC");
    rbc.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DBCR = (TH1F*) ratioIntegral(rbc.massFrom1DTemplatesEtaBinning,rd.mass)->Clone();
    h1DBCR->SetName("mass_predBCR");
    h1DBCR->Write();

    std::cout << "here3" << std::endl;

    saveHistoRatio(rd_40.mass,rbc_40.massFrom1DTemplatesEtaBinning,"mass_obs_40","mass_predBC_40","mass_predBCR_40");
    saveHistoRatio(rd_50.mass,rbc_50.massFrom1DTemplatesEtaBinning,"mass_obs_50","mass_predBC_50","mass_predBCR_50");
    saveHistoRatio(rd_60.mass,rbc_60.massFrom1DTemplatesEtaBinning,"mass_obs_60","mass_predBC_60","mass_predBCR_60");
    saveHistoRatio(rd_70.mass,rbc_70.massFrom1DTemplatesEtaBinning,"mass_obs_70","mass_predBC_70","mass_predBCR_70");
    saveHistoRatio(rd_80.mass,rbc_80.massFrom1DTemplatesEtaBinning,"mass_obs_80","mass_predBC_80","mass_predBCR_80");
    saveHistoRatio(rd_90.mass,rbc_90.massFrom1DTemplatesEtaBinning,"mass_obs_90","mass_predBC_90","mass_predBCR_90");


    

    std::cout << "here4" << std::endl;




    rd_boundedIas.mass->SetName("mass_obs_boundedIas");
    rd_boundedIas.mass->Write();

    rd_boundedIas.massFrom1DTemplatesEtaBinning->SetName("mass_pred1D_boundedIas");
    rd_boundedIas.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DR_boundedIas = (TH1F*) ratioIntegral(rd_boundedIas.massFrom1DTemplatesEtaBinning,rd_boundedIas.mass)->Clone();
    h1DR_boundedIas->SetName("mass_pred1DR_boundedIas");
    h1DR_boundedIas->Write();

    rbc_boundedIas.massFrom1DTemplatesEtaBinning->SetName("mass_predBC_boundedIas");
    rbc_boundedIas.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DBCR_boundedIas = (TH1F*) ratioIntegral(rbc_boundedIas.massFrom1DTemplatesEtaBinning,rd_boundedIas.mass)->Clone();
    h1DBCR_boundedIas->SetName("mass_predBCR_boundedIas");
    h1DBCR_boundedIas->Write();


    rd_boundedPt.mass->SetName("mass_obs_boundedPt");
    rd_boundedPt.mass->Write();

    rd_boundedPt.massFrom1DTemplatesEtaBinning->SetName("mass_pred1D_boundedPt");
    rd_boundedPt.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DR_boundedPt = (TH1F*) ratioIntegral(rd_boundedPt.massFrom1DTemplatesEtaBinning,rd_boundedPt.mass)->Clone();
    h1DR_boundedPt->SetName("mass_pred1DR_boundedPt");
    h1DR_boundedPt->Write();

    rbc_boundedPt.massFrom1DTemplatesEtaBinning->SetName("mass_predBC_boundedPt");
    rbc_boundedPt.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DBCR_boundedPt = (TH1F*) ratioIntegral(rbc_boundedPt.massFrom1DTemplatesEtaBinning,rd_boundedPt.mass)->Clone();
    h1DBCR_boundedPt->SetName("mass_predBCR_boundedPt");
    h1DBCR_boundedPt->Write();



    rdc.massFrom1DTemplatesEtaBinning->SetName("mass_predDC");
    rdc.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DCR = (TH1F*) ratioIntegral(rdc.massFrom1DTemplatesEtaBinning,rd.mass)->Clone();
    h1DCR->SetName("mass_predDCR");
    h1DCR->Write();

    rdb.massFrom1DTemplatesEtaBinning->SetName("mass_predDB");
    rdb.massFrom1DTemplatesEtaBinning->Write();

    TH1F* h1DBR = (TH1F*) ratioIntegral(rdb.massFrom1DTemplatesEtaBinning,rd.mass)->Clone();
    h1DBR->SetName("mass_predDBR");
    h1DBR->Write();



    plotting((TH1F*)rb.quantile01Ih_p,(TH1F*)rd.quantile01Ih_p,true,"quantile01_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile10Ih_p,(TH1F*)rd.quantile10Ih_p,true,"quantile10_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile30Ih_p,(TH1F*)rd.quantile30Ih_p,true,"quantile30_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile50Ih_p,(TH1F*)rd.quantile50Ih_p,true,"quantile50_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile70Ih_p,(TH1F*)rd.quantile70Ih_p,true,"quantile70_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile90Ih_p,(TH1F*)rd.quantile90Ih_p,true,"quantile90_ih_p_b_d","region B","region D")->Write();
    plotting((TH1F*)rb.quantile99Ih_p,(TH1F*)rd.quantile99Ih_p,true,"quantile99_ih_p_b_d","region B","region D")->Write();

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
    plotting(rall.mass,h_massFrom2D_D,false,"mass2D_regionD","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D_D,true,"mass2D_regionD_simpleRatio","Observed","Prediction from 2D template")->Write();
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,false,"mass1D_regionD","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,true,"mass1D_regionD_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC","Observed","Prediction from 1D templates in B and C")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,true,"mass1D_regionBC_simpleRatio","Observed","Prediction from 1D templates in B and C")->Write();

    plotting(rd_40.mass,rbc_40.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_40","Observed","Prediction")->Write();
    plotting(rd_50.mass,rbc_50.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_50","Observed","Prediction")->Write();
    plotting(rd_60.mass,rbc_60.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_60","Observed","Prediction")->Write();
    plotting(rd_70.mass,rbc_70.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_70","Observed","Prediction")->Write();
    plotting(rd_80.mass,rbc_80.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_80","Observed","Prediction")->Write();
    plotting(rd_90.mass,rbc_90.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_90","Observed","Prediction")->Write();


    std::cout << "here5" << std::endl;
    plotting(rall.mass,massFrom2D(rall,"all"),false,"mass2D_all","Observed","Prediction from 2D template - p: 10 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,2),false,"mass2D_all_IhP_rebin2_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.2 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,10),false,"mass2D_all_IhP_rebin10_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,1),false,"mass2D_all_IhP_rebin1_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",5,1),false,"mass2D_all_IhP_rebin1_5","Observed","Prediction from 2D template - p: 50 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,1),false,"mass2D_all_IhP_rebin1_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 0.1 MeV/cm")->Write();




    TCanvas* c4 = new TCanvas();
    c4->cd();
    for(int i=0;i<rall.eta_p->GetNbinsX();i++)
    {
        //TH1F* hih1=(TH1F*)((TH2F*)rall.ih_p_eta->Project3D("zx"))->ProjectionY("",i,i);
        TH1F* hih1=(TH1F*)rall.ih_p->ProjectionY("",i,i+1);
        TH1F* hih2=(TH1F*)h_cross1D_all->ProjectionY("",i,i+1);
        if(hih1->Integral()<=0) continue;
        scale(hih1);
        scale(hih2);
        hih1->Draw();
        hih1->GetXaxis()->SetRangeUser(0,10);
        hih1->GetYaxis()->SetRangeUser(1e-8,1);

        hih2->SetLineColor(2);
        hih2->Draw("same");
        c4->SetLogy();
        //c4->SaveAs(("ih_cross/ih_p"+to_string(i)+".pdf").c_str());
        
    }
    */
}
