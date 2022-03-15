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
    }

}

void saveHistoRatio(TH1F* h1,TH1F* h2,std::string st1,std::string st2,std::string st3){
        h1->SetName(st1.c_str());
        h1->Write();
        h2->SetName(st2.c_str());
        h2->Write();
        TH1F* R = (TH1F*) ratioIntegral(h2,h1)->Clone();
        R->SetName(st3.c_str());
        R->Write();

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

    region rb_50;
    region rb_60;
    region rb_70;
    region rb_80;
    region rb_90;

    region rbc_50;
    region rbc_60;
    region rbc_70;
    region rbc_80;
    region rbc_90;

    region rd_50;
    region rd_60;
    region rd_70;
    region rd_80;
    region rd_90;

    region ra_med;
    region rc_med;
    
    bool bool_rebin=rebin;
    
    
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

    loadHistograms(rb_50,ifile,"regionB_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60,ifile,"regionB_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70,ifile,"regionB_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80,ifile,"regionB_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90,ifile,"regionB_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rbc_50,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(ra_med,ifile,"regionA_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_med,ifile,"regionC_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 


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
    
    etaReweighingP(rc.eta_p,rb.eta_p); 
    //etaReweighingP(rc_boundedIas.eta_p,rb_boundedIas.eta_p); 
    //etaReweighingP(rc_boundedPt.eta_p,rb_boundedPt.eta_p); 
  
    etaReweighingP(rc_med.eta_p,rb_50.eta_p);
    etaReweighingP(rc_med.eta_p,rb_60.eta_p);
    etaReweighingP(rc_med.eta_p,rb_70.eta_p);
    etaReweighingP(rc_med.eta_p,rb_80.eta_p);
    etaReweighingP(rc_med.eta_p,rb_90.eta_p);

    //rbc = rd;
    rbc.eta_p = rc.eta_p;
    rbc.ih_eta = rb.ih_eta;

    /*rbc_boundedIas.eta_p = rc_boundedIas.eta_p;
    rbc_boundedIas.ih_eta = rb_boundedIas.ih_eta;

    rbc_boundedPt.eta_p = rc_boundedPt.eta_p;
    rbc_boundedPt.ih_eta = rb_boundedPt.ih_eta;
    */
    rdc.eta_p = rc.eta_p;
    rdb.ih_eta = rb.ih_eta;

    rbc_50.eta_p = rc_med.eta_p;    rbc_50.ih_eta = rb_50.ih_eta;
    rbc_60.eta_p = rc_med.eta_p;    rbc_60.ih_eta = rb_60.ih_eta;
    rbc_70.eta_p = rc_med.eta_p;    rbc_70.ih_eta = rb_70.ih_eta;
    rbc_80.eta_p = rc_med.eta_p;    rbc_80.ih_eta = rb_80.ih_eta;
    rbc_90.eta_p = rc_med.eta_p;    rbc_90.ih_eta = rb_90.ih_eta;

    rall.fillMassFrom1DTemplatesEtaBinning();

    float normalisationABC = rb.ih_eta->GetEntries()*rc.ih_eta->GetEntries()/ra.ih_eta->GetEntries();
    
    float normalisationABC_50 = rb_50.ih_eta->GetEntries()*rc_med.ih_eta->GetEntries()/ra_med.ih_eta->GetEntries();
    float normalisationABC_60 = rb_60.ih_eta->GetEntries()*rc_med.ih_eta->GetEntries()/ra_med.ih_eta->GetEntries();
    float normalisationABC_70 = rb_70.ih_eta->GetEntries()*rc_med.ih_eta->GetEntries()/ra_med.ih_eta->GetEntries();
    float normalisationABC_80 = rb_80.ih_eta->GetEntries()*rc_med.ih_eta->GetEntries()/ra_med.ih_eta->GetEntries();
    float normalisationABC_90 = rb_90.ih_eta->GetEntries()*rc_med.ih_eta->GetEntries()/ra_med.ih_eta->GetEntries();
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

    rbc_50.fillMassFrom1DTemplatesEtaBinning();
    rbc_60.fillMassFrom1DTemplatesEtaBinning();
    rbc_70.fillMassFrom1DTemplatesEtaBinning();
    rbc_80.fillMassFrom1DTemplatesEtaBinning();
    rbc_90.fillMassFrom1DTemplatesEtaBinning();

    


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


    //scale(rd.mass);
    scale(rd.massFrom1DTemplatesEtaBinning);
    //rbc.massFrom1DTemplatesEtaBinning->Scale(normalisationABC);
    scale(rbc.massFrom1DTemplatesEtaBinning);
    scale(h_massFrom2D_D);

    scale(rbc_50.massFrom1DTemplatesEtaBinning);
    scale(rbc_60.massFrom1DTemplatesEtaBinning);
    scale(rbc_70.massFrom1DTemplatesEtaBinning);
    scale(rbc_80.massFrom1DTemplatesEtaBinning);
    scale(rbc_90.massFrom1DTemplatesEtaBinning);

    for(int k=0;k<rbc.massFrom1DTemplatesEtaBinning->GetNbinsX()+1;k++){
        rbc.massFrom1DTemplatesEtaBinning->SetBinContent(k,rbc.massFrom1DTemplatesEtaBinning->GetBinContent(k)*normalisationABC);
    }

    massNormalisation(rbc_50.massFrom1DTemplatesEtaBinning,normalisationABC_50);
    massNormalisation(rbc_60.massFrom1DTemplatesEtaBinning,normalisationABC_60);
    massNormalisation(rbc_70.massFrom1DTemplatesEtaBinning,normalisationABC_70);
    massNormalisation(rbc_80.massFrom1DTemplatesEtaBinning,normalisationABC_80);
    massNormalisation(rbc_90.massFrom1DTemplatesEtaBinning,normalisationABC_90);

    /*scale(rd_boundedIas.mass);
    scale(rd_boundedIas.massFrom1DTemplatesEtaBinning);
    scale(rbc_boundedIas.massFrom1DTemplatesEtaBinning);
    scale(rd_boundedPt.mass);
    scale(rd_boundedPt.massFrom1DTemplatesEtaBinning);
    scale(rbc_boundedPt.massFrom1DTemplatesEtaBinning);*/
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


    saveHistoRatio(rd_50.mass,rbc_50.massFrom1DTemplatesEtaBinning,"mass_obs_50","mass_predBC_50","mass_predBCR_50");
    saveHistoRatio(rd_60.mass,rbc_60.massFrom1DTemplatesEtaBinning,"mass_obs_60","mass_predBC_60","mass_predBCR_60");
    saveHistoRatio(rd_70.mass,rbc_70.massFrom1DTemplatesEtaBinning,"mass_obs_70","mass_predBC_70","mass_predBCR_70");
    saveHistoRatio(rd_80.mass,rbc_80.massFrom1DTemplatesEtaBinning,"mass_obs_80","mass_predBC_80","mass_predBCR_80");
    saveHistoRatio(rd_90.mass,rbc_90.massFrom1DTemplatesEtaBinning,"mass_obs_90","mass_predBC_90","mass_predBCR_90");


    




/*
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
    h1DBCR_boundedPt->Write();*/



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


/*
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
*/
    /*plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,false,"mass1D_all","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,rall.massFrom1DTemplatesEtaBinning,true,"mass1D_all_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rall.mass,h_massFrom2D,false,"mass2D_all","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D,true,"mass2D_all_simpleRatio","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D_D,false,"mass2D_regionD","Observed","Prediction from 2D template")->Write();
    plotting(rall.mass,h_massFrom2D_D,true,"mass2D_regionD_simpleRatio","Observed","Prediction from 2D template")->Write();*/
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,false,"mass1D_regionD","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rd.massFrom1DTemplatesEtaBinning,true,"mass1D_regionD_simpleRatio","Observed","Prediction from 1D templates")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC","Observed","Prediction from 1D templates in B and C")->Write();
    plotting(rd.mass,rbc.massFrom1DTemplatesEtaBinning,true,"mass1D_regionBC_simpleRatio","Observed","Prediction from 1D templates in B and C")->Write();

    plotting(rd_50.mass,rbc_50.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_50","Observed","Prediction")->Write();
    plotting(rd_60.mass,rbc_60.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_60","Observed","Prediction")->Write();
    plotting(rd_70.mass,rbc_70.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_70","Observed","Prediction")->Write();
    plotting(rd_80.mass,rbc_80.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_80","Observed","Prediction")->Write();
    plotting(rd_90.mass,rbc_90.massFrom1DTemplatesEtaBinning,false,"mass1D_regionBC_90","Observed","Prediction")->Write();


    /*plotting(rall.mass,massFrom2D(rall,"all"),false,"mass2D_all","Observed","Prediction from 2D template - p: 10 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,2),false,"mass2D_all_IhP_rebin2_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.2 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,10),false,"mass2D_all_IhP_rebin10_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",2,1),false,"mass2D_all_IhP_rebin1_2","Observed","Prediction from 2D template - p: 20 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",5,1),false,"mass2D_all_IhP_rebin1_5","Observed","Prediction from 2D template - p: 50 GeV - dEdx: 0.1 MeV/cm")->Write();
    plotting(rall.mass,massFrom2D(rall,"all",10,1),false,"mass2D_all_IhP_rebin1_10","Observed","Prediction from 2D template - p: 100 GeV - dEdx: 0.1 MeV/cm")->Write();*/




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

}
