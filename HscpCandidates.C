#define HscpCandidates_cxx
#include "HscpCandidates.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>


TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
    size_t pos = Path.find("/");
    if(pos < 256)
    {
        std::string firstPart   = Path.substr(0,pos);
        std::string endPart     = Path.substr(pos+1,Path.length());
        TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
        if(TMP!=nullptr)return GetObjectFromPath(TMP,endPart,GetACopy);
        printf("ObjectNotFound: %s::%s\n",File->GetName(), Path.c_str());
        return nullptr;
    }
    else
    {
        if(GetACopy)
        {
            return (File->Get(Path.c_str()))->Clone();
        }
        else
        {
            File->Get(Path.c_str());
        }
    }
}

void predMass(TH1F* res, TH1F* h_p, TH1F* h_ih, float norm=0)
{
    //double norm = h_p->Integral();
    scale(h_p); scale(h_ih);
    for(int i=0;i<h_p->GetNbinsX();i++)
    {
        for(int j=0;j<h_ih->GetNbinsX();j++)
        {
            float p = h_p->GetBinCenter(i);
            float ih = h_ih->GetBinCenter(j);
            float prob = h_p->GetBinContent(i) * h_ih->GetBinContent(j);
            float mass = GetMass(p,ih,K,C);
            res->Fill(mass,prob*norm);
        }
    }
}

float HscpCandidates::PFIsolationMuon(int i){
    float iso = MuonPFIsolationR03_sumChargedHadronPt->at(i) + std::max(0.0,MuonPFIsolationR03_sumPhotonPt->at(i)+MuonPFIsolationR03_sumNeutralHadronPt->at(i)-0.5*MuonPFIsolationR03_sumPUPt->at(i));
    return iso/Pt->at(i);
}

float HscpCandidates::PFIsolationTrack(int i){
    float iso = TrackPFIsolationR03_sumChargedHadronPt->at(i) + std::max(0.0,TrackPFIsolationR03_sumPhotonPt->at(i)+TrackPFIsolationR03_sumNeutralHadronPt->at(i)-0.5*TrackPFIsolationR03_sumPUPt->at(i));
    return iso/Pt->at(i);
}

float pTerrOverpT_forPt(const float& pt){
    return 2.5e-4*pt+3e-2;
}

void HscpCandidates::Loop()
{
//   In a ROOT session, you can do:
//      root> .L HscpCandidates.C
//      root> HscpCandidates t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
  
   Long64_t nentries = fChain->GetEntriesFast();


   std::cout << "nentries: " << nentries << std::endl;

   int ntot=0, nA=0, nB=0, nC=0, nD=0, nB_boundedIas=0, nD_boundedIas=0, nC_boundedPt=0, nD_boundedPt=0;

   region rAll("_all",etabins_,ihbins_,pbins_,massbins_);
   region rA("_regionA",etabins_,ihbins_,pbins_,massbins_);
   region rB("_regionB",etabins_,ihbins_,pbins_,massbins_);
   region rC("_regionC",etabins_,ihbins_,pbins_,massbins_);
   region rD("_regionD",etabins_,ihbins_,pbins_,massbins_);

   region rA_med("_regionA_med",etabins_,ihbins_,pbins_,massbins_);
   region rA_med_pt("_regionA_med_pt",etabins_,ihbins_,pbins_,massbins_);
   region rB_med_pt("_regionB_med_pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_med("_regionC_med",etabins_,ihbins_,pbins_,massbins_);


   region rA_40("_regionA_40",etabins_,ihbins_,pbins_,massbins_);
   region rA_40pt("_regionA_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rB_40pt("_regionB_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_40("_regionC_40",etabins_,ihbins_,pbins_,massbins_);

   
   region rB_40("_regionB_40",etabins_,ihbins_,pbins_,massbins_);
   region rB_50("_regionB_50",etabins_,ihbins_,pbins_,massbins_);
   region rB_60("_regionB_60",etabins_,ihbins_,pbins_,massbins_);
   region rB_70("_regionB_70",etabins_,ihbins_,pbins_,massbins_);
   region rB_80("_regionB_80",etabins_,ihbins_,pbins_,massbins_);
   region rB_90("_regionB_90",etabins_,ihbins_,pbins_,massbins_);
   region rB_50_90("_regionB_50_90",etabins_,ihbins_,pbins_,massbins_);

   region rD_40("_regionD_40",etabins_,ihbins_,pbins_,massbins_);
   region rD_50("_regionD_50",etabins_,ihbins_,pbins_,massbins_);
   region rD_60("_regionD_60",etabins_,ihbins_,pbins_,massbins_);
   region rD_70("_regionD_70",etabins_,ihbins_,pbins_,massbins_);
   region rD_80("_regionD_80",etabins_,ihbins_,pbins_,massbins_);
   region rD_90("_regionD_90",etabins_,ihbins_,pbins_,massbins_);
   region rD_50_90("_regionD_50_90",etabins_,ihbins_,pbins_,massbins_);

   region rC_40pt("_regionC_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_50pt("_regionC_50pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_60pt("_regionC_60pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_70pt("_regionC_70pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_80pt("_regionC_80pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_90pt("_regionC_90pt",etabins_,ihbins_,pbins_,massbins_);

   region rD_40pt("_regionD_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rD_50pt("_regionD_50pt",etabins_,ihbins_,pbins_,massbins_);
   region rD_60pt("_regionD_60pt",etabins_,ihbins_,pbins_,massbins_);
   region rD_70pt("_regionD_70pt",etabins_,ihbins_,pbins_,massbins_);
   region rD_80pt("_regionD_80pt",etabins_,ihbins_,pbins_,massbins_);
   region rD_90pt("_regionD_90pt",etabins_,ihbins_,pbins_,massbins_);



/*   //Ias_outer quantiles
   float quan50 = 0.052;
   float quan60 = 0.062;
   float quan70 = 0.075;
   float quan80 = 0.093;
   float quan90 = 0.117;
*/

/*   //Ias_all quantiles //data
   float quan50 = 0.039;
   float quan60 = 0.045;
   float quan70 = 0.053;
   float quan80 = 0.064;
   float quan90 = 0.082;
*/


   //Ias_all quantiles //MC WJets
   //double q[6]={ 0.029754446, 0.034736341, 0.040516857, 0.047670289, 0.057816411, 0.075263733 }; //MC WJets
   //double q_pt[6]={ 57.294373, 60.057283, 63.757415, 68.985768, 77.181910, 93.585251 }; //MC WJets
   //double q[5]={ 0.039, 0.045, 0.053, 0.064, 0.082 }; //data or signal
   double q[6]={ 0.033863945, 0.039237098, 0.045422508, 0.053164483, 0.063867635, 0.082214175 }; //data or signal
   double q_pt[6]={ 57.013944, 59.648194, 63.123116, 67.930817, 75.501106, 90.672101 }; //data or signal
   
   float quan40 = q[0];
   float quan50 = q[1];
   float quan60 = q[2];
   float quan70 = q[3];
   float quan80 = q[4];
   float quan90 = q[5];


   float quan40_pt = q_pt[0];
   float quan50_pt = q_pt[1];
   float quan60_pt = q_pt[2];
   float quan70_pt = q_pt[3];
   float quan80_pt = q_pt[4];
   float quan90_pt = q_pt[5];


  


   TH1F* h_massObs_q40 = new TH1F("massObs_q40",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50 = new TH1F("massObs_q50",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q60 = new TH1F("massObs_q60",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q70 = new TH1F("massObs_q70",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q80 = new TH1F("massObs_q80",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q90 = new TH1F("massObs_q90",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_90 = new TH1F("massObs_q50_90",";Mass (GeV)",80,0,4000);


   TH1F* h_massObs_q40_pt = new TH1F("massObs_q40_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_pt = new TH1F("massObs_q50_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q60_pt = new TH1F("massObs_q60_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q70_pt = new TH1F("massObs_q70_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q80_pt = new TH1F("massObs_q80_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q90_pt = new TH1F("massObs_q90_pt",";Mass (GeV)",80,0,4000);



   h_massObs_q40->Sumw2();
   h_massObs_q50->Sumw2();
   h_massObs_q60->Sumw2();
   h_massObs_q70->Sumw2();
   h_massObs_q80->Sumw2();
   h_massObs_q90->Sumw2();
   h_massObs_q50_90->Sumw2();
   
   h_massObs_q40_pt->Sumw2();
   h_massObs_q50_pt->Sumw2();
   h_massObs_q60_pt->Sumw2();
   h_massObs_q70_pt->Sumw2();
   h_massObs_q80_pt->Sumw2();
   h_massObs_q90_pt->Sumw2();


   TH1F* h_mT = new TH1F("mT",";m_{T} [GeV];tracks/bin",200,0,200);
   TH1F* h_mT_regD = new TH1F("mT_regD",";m_{T} [GeV];tracks/bin",200,0,200);
   TH1F* h_probQ = new TH1F("probQ",";probQ;tracks/bin",100,0,1);
   TH1F* h_isoPFMuon = new TH1F("isoPFMuon",";Iso/p_{T};",100,0,2);
   TH1F* h_njets = new TH1F("njets",";njets;",20,0,20);
   TH1F* h_nhscp = new TH1F("nhscp",";nhscp",20,0,20);
   TH1F* h_massObs = new TH1F("massObs",";Mass (GeV);tracks/bin",80,0,4000);
   TH1F* h_isoPFTk = new TH1F("isoPFTk",";p_{T} relative tracker PF-based isolation (GeV);tracks/bin",100,0,1);
   TH1F* h_ias_outer = new TH1F("ias_outer",";I_{as} outer;tracks/bin",100,0,1);
   TH1F* h_ias = new TH1F("ias",";I_{as};tracks/bin",100,0,1);
   TH1F* h_ECAL = new TH1F("ECAL",";E_{ecal} (GeV);tracks/bin",100,0,20); 
   TH1F* h_HCAL = new TH1F("HCAL",";E_{hcal} (GeV);tracks/bin",100,0,20); 

   TH1F* h_pT_isHighPurity = new TH1F("pT_isHighPurity",";p_{T} (GeV);tracks/bin",200,0,4000);
   TH1F* h_pT_QualityMask = new TH1F("pT_QualityMask",";p_{T} (GeV);tracks/bin",200,0,4000);
   TH2F* h_pT_isHighPurity_QualityMask = new TH2F("pT_isHighPurity_QualityMask",";p_{T} high purity (GeV);p_{T} quality mask (GeV);tracks/bin",200,0,4000,200,0,4000);

   TH1F* h_mass_isHighPurity = new TH1F("mass_isHighPurity",";Mass (GeV);tracks/bin",80,0,4000);
   TH1F* h_mass_QualityMask = new TH1F("mass_QualityMask",";Mass (GeV);tracks/bin",80,0,4000);
   TH2F* h_mass_isHighPurity_QualityMask = new TH2F("mass_isHighPurity_QualityMask",";Mass high purity (GeV);Mass quality mask (GeV);tracks/bin",80,0,4000,80,0,4000);
   
   TH1F* h_pT_pTerrCut = new TH1F("pT_pTerrCut","#frac{#sigma_{pT}}{p_{T}}<0.25;p_{T} (GeV);tracks/bin",200,0,4000);
   TH1F* h_pT_pTerrCut2 = new TH1F("pT_pTerrCut2","#frac{#sigma_{pT}}{p_{T}^{2}}<0.001;p_{T} (GeV);tracks/bin",200,0,4000);
   TH2F* h_pT_pTerrCut_pTerrCut2 = new TH2F("pT_pTerrCut_pTerrCut2",";p_{T} (GeV) #frac{#sigma_{pT}}{p_{T}}<0.25;p_{T} (GeV) #frac{#sigma_{pT}}{p_{T}^{2}}<0.001;tracks/bin",200,0,4000,200,0,4000);
   TH2F* h_pT_noCut_pTerrCut = new TH2F("pT_noCut_pTerrCut",";p_{T} (GeV);p_{T} (GeV) #frac{#sigma_{pT}}{p_{T}}<0.25;tracks/bin",200,0,4000,200,0,4000);

   TH1F* h_mass_pTerrCut = new TH1F("mass_pTerrCut","#frac{#sigma_{pT}}{p_{T}}<0.25;Mass (GeV);tracks/bin",80,0,4000);
   TH1F* h_mass_pTerrCut2 = new TH1F("mass_pTerrCut2","#frac{#sigma_{pT}}{p_{T}^{2}}<0.001;Mass (GeV);tracks/bin",80,0,4000);
   TH2F* h_mass_pTerrCut_pTerrCut2 = new TH2F("mass_pTerrCut_pTerrCut2",";Mass (GeV) #frac{#sigma_{pT}}{p_{T}}<0.25;Mass (GeV) #frac{#sigma_{pT}}{p_{T}^{2}}<0.001;tracks/bin",80,0,4000,80,0,4000);
   TH2F* h_mass_noCut_pTerrCut = new TH2F("mass_noCut_pTerrCut",";Mass (GeV);Mass (GeV) #frac{#sigma_{pT}}{p_{T}}<0.25;tracks/bin",80,0,4000,80,0,4000);

   TH2F* h_pTerrOverpT_pTerrOverpT2 = new TH2F("pTerrOverpT_pTerrOverpT2",";#frac{#sigma_{pT}}{p_{T}};#frac{#sigma_{pT}}{p_{T}^{2}};tracks/bin",100,0,1,100,0,0.01);

   TH1F* h_pT = new TH1F("pT",";p_{T} (GeV);tracks/bin",1000,0,6000);
   TH1F* h_pTerrOverpT = new TH1F("pTerrOverpT",";#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,1);
   TH1F* h_pTerrOverpT2 = new TH1F("pTerrOverpT2",";#frac{#sigma_{pT}}{p_{T}^{2}};tracks/bin",1000,0,0.005);
   TH2F* h_pT_pTerrOverpT = new TH2F("pT_pTerrOverpT",";p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   TH2F* h_pT_pTerrOverpT2 = new TH2F("pT_pTerrOverpT2",";p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}^{2}};tracks/bin",1000,0,6000,1000,0,0.005);

   TH1F* h_pTerrOverpT_0_pT_100 = new TH1F("pTerrOverpT_0_pT_100","0<p_{T}<100;#frac{#sigma_{pT}}{p_{T}}",1000,0,1);
   TH1F* h_pTerrOverpT_100_pT_200 = new TH1F("pTerrOverpT_100_pT_200","100<p_{T}<200;#frac{#sigma_{pT}}{p_{T}}",1000,0,1);
   TH1F* h_pTerrOverpT_200_pT_300 = new TH1F("pTerrOverpT_200_pT_300","200<p_{T}<300;#frac{#sigma_{pT}}{p_{T}}",1000,0,1);
   TH1F* h_pTerrOverpT_300_pT_400 = new TH1F("pTerrOverpT_300_pT_400","300<p_{T}<400;#frac{#sigma_{pT}}{p_{T}}",1000,0,1);
   TH1F* h_pTerrOverpT_400_pT_800 = new TH1F("pTerrOverpT_400_pT_800","400<p_{T}<800;#frac{#sigma_{pT}}{p_{T}}",1000,0,1);
   TH1F* h_pTerrOverpT_800_pT = new TH1F("pTerrOverpT_800_pT","800<p_{T};#frac{#sigma_{pT}}{p_{T}}",1000,0,1);

   TH1F* h_pTerrOverpT2_0_pT_100 = new TH1F("pTerrOverpT2_0_pT_100","0<p_{T}<100;#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);
   TH1F* h_pTerrOverpT2_100_pT_200 = new TH1F("pTerrOverpT2_100_pT_200","100<p_{T}<200;#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);
   TH1F* h_pTerrOverpT2_200_pT_300 = new TH1F("pTerrOverpT2_200_pT_300","200<p_{T}<300;#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);
   TH1F* h_pTerrOverpT2_300_pT_400 = new TH1F("pTerrOverpT2_300_pT_400","300<p_{T}<400;#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);
   TH1F* h_pTerrOverpT2_400_pT_800 = new TH1F("pTerrOverpT2_400_pT_800","400<p_{T}<800;#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);
   TH1F* h_pTerrOverpT2_800_pT = new TH1F("pTerrOverpT2_800_pT","800<p_{T};#frac{#sigma_{pT}}{p_{T}^{2}}",1000,0,0.005);


   TH2F* h_pT_pTerrOverpT_eta_0p8 = new TH2F("pT_pTerrOverpT_eta_0p8","|#eta|<0.8;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   TH2F* h_pT_pTerrOverpT_0p8_eta_1p3 = new TH2F("pT_pTerrOverpT_0p8_eta_1p3","0.8<|#eta|<1.3;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   TH2F* h_pT_pTerrOverpT_1p3_eta_1p7 = new TH2F("pT_pTerrOverpT_1p3_eta_1p7","1.3<|#eta|<1.7;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   TH2F* h_pT_pTerrOverpT_1p7_eta_2p1 = new TH2F("pT_pTerrOverpT_1p7_eta_2p1","1.7<|#eta|<2.1;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   
   
   h_massObs->Sumw2();
   h_mT->Sumw2();
   h_mT_regD->Sumw2();
   h_probQ->Sumw2();
   h_isoPFTk->Sumw2();
   h_ias_outer->Sumw2();
   h_ias->Sumw2();
   h_ECAL->Sumw2();
   h_HCAL->Sumw2();

   h_pT_pTerrOverpT_eta_0p8->Sumw2();
   h_pT_pTerrOverpT_0p8_eta_1p3->Sumw2();
   h_pT_pTerrOverpT_1p3_eta_1p7->Sumw2();
   h_pT_pTerrOverpT_1p7_eta_2p1->Sumw2();

   h_massObs->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);

   //float lumi=7.04; //fb-1
   float lumi=9.69; //fb-1

   float weightGluino2600 = (5.0e-2*lumi)/100281.0; //Gluino at 2.6 TeV
   float weightGluino2000 = (9.7e-1*lumi)/97151.0; //Gluino at 2 TeV
   float weightGluino1400 = (2.5e1*lumi)/100061.0;
   float weightGluino1000 = (3.2e2*lumi)/99564.0;

   float weightStau1599 = (1.4e-4*lumi)/7200.0;
   
   float weightgmStau871 = (6.9e-2*lumi)/25000.0;
   float weightgmStau1029 = (2.2e-2*lumi)/25000.0;

   float weightppStau871 = (9.9e-3*lumi)/25000.0;
   float weightppStau1029 = (3.5e-3*lumi)/25000.0;
   

   

   float weightWJets = (52940e3*lumi)/(6.988236e7);
   float weightTTTo2L2Nu = (88.29e3*lumi)/(7.5653e7);
   float weightTTToSemiLeptonic = (365.35e3*lumi)/(1.29985e8);
   float weightTTToHadronic = (377.96e3*lumi)/(1.0117e8);
   float weightQCD = (239000.0e3*lumi)/8994317.0;
   float weightDYJetsToLL = (15343.0e3*lumi)/(4.89014e7);

   float weightMass = weightppStau1029;
   
   weightMass = 1.;

   
   //cutindex=19 pt=60 ias=0.025
   
   double pt_cut=ptcut_;
   double ias_cut=iascut_;
  

   int n1=0, n2=0, n3=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //if(jentry%2==0) continue;
      //if(jentry%2==1) continue;


      if(jentry%1000000==0) std::cout << "entry: " << jentry << std::endl;

      if(Hscp<1) continue;
      if(Pt->size()<1) continue;

      if(invMET_ && RecoPFMET>100) continue;
      
      if(!HLT_Mu50) continue;

      h_njets->Fill(njets);

      int count_hscp_passPre = 0;

      int i=0;
      for(int i=0; i<Mass->size(); i++)
      {
          //if(invIso_ && !passPreselection_noIsolation_noIh->at(i)) continue;
          //if(!invIso_ && !passPreselection->at(i)) continue;

          //if(!passPreselection->at(i)) continue;

          float pt = Pt->at(i);

          //float ih = Ih_StripOnly->at(i);
          float ih = Ih_noL1->at(i);
          float ias = Ias->at(i);
          //float ias = Ias_noTIBnoTIDno3TEC->at(i);
          //float ias = Ias_PixelOnly->at(i);
          float Eta = eta->at(i); 
          float iso = iso_TK->at(i);
          float iso_r = iso/pt;
          float nhits = NOM->at(i);
          float p = pt*cosh(Eta);
//          float massT = Mass->at(i);
          float massT = GetMass(p,ih,K,C);
          float isotk = iso_TK->at(i);
          float isocalo = iso_ECAL->at(i)+iso_HCAL->at(i);
          float tof = TOF->at(i);
          float dz = dZ->at(i);
          float dxy = dXY->at(i);
          float probChi2 = TMath::Prob(Chi2->at(i),Ndof->at(i));

          float IsoRel = PFIsolationMuon(i);
          float IsoRelTk = PFIsolationTrack(i);
          h_isoPFTk->Fill(IsoRelTk,weightMass);

          h_ias_outer->Fill(Ias_noTIBnoTIDno3TEC->at(i),weightMass);
          h_ias->Fill(ias,weightMass);
         
          if(pt<50) continue;
          if(abs(Eta)>2.1) continue;
          if(NOPH->at(i)<2) continue;
          if(NOH->at(i)<8) continue;
          if(FOVH->at(i)<0.8) continue;
          if(NOM->at(i)<6) continue;
          if(Chi2->at(i)/Ndof->at(i)>5) continue;
          
          if(isotk>50) continue;
          if(isocalo/p>0.3) continue;
          
          if(abs(dz)>0.5) continue;
          if(abs(dxy)>0.02) continue;

          if(ih<ihcut_) continue;
          //if(ih>3.17) continue;
          if(p>pcut_ && pcut_>0) continue;

          bool pTerrCut1 = false;
          bool pTerrCut2 = false;

          if(PtErr->at(i)/pt<0.25) pTerrCut1=true;
          if(PtErr->at(i)/pow(pt,2)<0.0005) pTerrCut2=true;

          //if(mT->at(i)<60 || mT->at(i)>100) continue;
          //if(!isMuon->at(i)) continue;

          //if(PtErr->at(i)/pt>0.25) continue;

          float pT_qualMask=0, pT_highpur=0;
          float mass_qualMask=0, mass_highpur=0;

          /*if(QualityMask->at(i)>2){
              h_pT_QualityMask->Fill(pt,weightMass);
              h_mass_QualityMask->Fill(massT,weightMass);
              pT_qualMask=pt;
              mass_qualMask=massT;
          }

          if(isHighPurity->at(i)){
              h_pT_isHighPurity->Fill(pt,weightMass);
              h_mass_isHighPurity->Fill(massT,weightMass);
              pT_highpur=pt;
              mass_highpur=massT;
          }
          h_pT_isHighPurity_QualityMask->Fill(pT_highpur,pT_qualMask,weightMass);
          h_mass_isHighPurity_QualityMask->Fill(mass_highpur,mass_qualMask,weightMass);

          //if(QualityMask->at(i)<3) continue;*/

          if(!isHighPurity->at(i)) continue;
          //if(!isMuon->at(i)) continue;
          //if(probChi2<0.1) continue;


         
          float pT_cut1=0, pT_cut2=0, mass_cut1=0, mass_cut2=0;

          if(pTerrCut1){
              h_pT_pTerrCut->Fill(pt,weightMass);
              h_mass_pTerrCut->Fill(massT,weightMass);
              pT_cut1=pt;
              mass_cut1=massT;
          }
          if(pTerrCut2){
              h_pT_pTerrCut2->Fill(pt,weightMass);
              h_mass_pTerrCut2->Fill(massT,weightMass);
              pT_cut2=pt;
              mass_cut2=massT;
          }
          h_pT_noCut_pTerrCut->Fill(pt,pT_cut1,weightMass);
          h_pT_pTerrCut_pTerrCut2->Fill(pT_cut1,pT_cut2,weightMass);
          h_mass_noCut_pTerrCut->Fill(massT,mass_cut1,weightMass);
          h_mass_pTerrCut_pTerrCut2->Fill(mass_cut1,mass_cut2,weightMass);


          
          h_pTerrOverpT_pTerrOverpT2->Fill(PtErr->at(i)/pt,PtErr->at(i)/pow(pt,2),weightMass);
          
          h_pTerrOverpT2->Fill(PtErr->at(i)/pow(pt,2),weightMass);



          h_pT->Fill(pt,weightMass);
          h_pTerrOverpT->Fill(PtErr->at(i)/pt,weightMass);
          h_pTerrOverpT2->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          h_pT_pTerrOverpT->Fill(pt,PtErr->at(i)/pt,weightMass);
          h_pT_pTerrOverpT2->Fill(pt,PtErr->at(i)/pow(pt,2),weightMass);

          if(abs(Eta)<0.8) h_pT_pTerrOverpT_eta_0p8->Fill(pt,PtErr->at(i)/pt,weightMass);
          if(abs(Eta)<1.3 && abs(Eta)>=0.8) h_pT_pTerrOverpT_0p8_eta_1p3->Fill(pt,PtErr->at(i)/pt,weightMass);
          if(abs(Eta)<1.7 && abs(Eta)>=1.3) h_pT_pTerrOverpT_1p3_eta_1p7->Fill(pt,PtErr->at(i)/pt,weightMass);
          if(abs(Eta)<2.1 && abs(Eta)>=1.7) h_pT_pTerrOverpT_1p7_eta_2p1->Fill(pt,PtErr->at(i)/pt,weightMass);


          if(pt>=0 && pt<100) h_pTerrOverpT_0_pT_100->Fill(PtErr->at(i)/pt,weightMass);
          if(pt>=100 && pt<200) h_pTerrOverpT_100_pT_200->Fill(PtErr->at(i)/pt,weightMass);
          if(pt>=200 && pt<300) h_pTerrOverpT_200_pT_300->Fill(PtErr->at(i)/pt,weightMass);
          if(pt>=300 && pt<400) h_pTerrOverpT_300_pT_400->Fill(PtErr->at(i)/pt,weightMass);
          if(pt>=400 && pt<800) h_pTerrOverpT_400_pT_800->Fill(PtErr->at(i)/pt,weightMass);
          if(pt>=800) h_pTerrOverpT_800_pT->Fill(PtErr->at(i)/pt,weightMass);


          if(pt>=0 && pt<100) h_pTerrOverpT2_0_pT_100->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          if(pt>=100 && pt<200) h_pTerrOverpT2_100_pT_200->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          if(pt>=200 && pt<300) h_pTerrOverpT2_200_pT_300->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          if(pt>=300 && pt<400) h_pTerrOverpT2_300_pT_400->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          if(pt>=400 && pt<800) h_pTerrOverpT2_400_pT_800->Fill(PtErr->at(i)/pow(pt,2),weightMass);
          if(pt>=800) h_pTerrOverpT2_800_pT->Fill(PtErr->at(i)/pow(pt,2),weightMass);


          n1++;
          if(pTerrCut1) n3++;
          
          if(PtErr->at(i)/pt>pTerrOverpT_forPt(pt)) continue;
 
          n2++;

          //if(!pTerrCut2) continue;
          

   /*
          if(!(etacutinf_<Eta && Eta<etacutsup_)) continue;
          if(isocalo/p>0.3)continue;
          if((invIso_ && (isotk<50 || isotk>100))) continue;
*/
          //if(ProbQ_noL1->at(i)>0.1) continue;

          //if(mT->at(i)<90) continue;

          h_mT->Fill(mT->at(i),weightMass);
          h_probQ->Fill(ProbQ_noL1->at(i),weightMass);
          h_isoPFMuon->Fill(PFIsolationMuon(i),weightMass);
          h_ECAL->Fill(ECAL_energy->at(i),weightMass);
          h_HCAL->Fill(HCAL_energy->at(i),weightMass);
         
          count_hscp_passPre++;

          //std::cout << "mT: " << mT->at(i) << std::endl;

          rAll.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
          ntot++;




          if(pt<=pt_cut)
          {
              if(ias<=ias_cut) //regionA
              {
                  rA.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nA++;

              }
              else //regionB
              {
                  rB.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nB++;

                  if(ias<0.1)
                  {
                      //rB_boundedIas.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
                      nB_boundedIas++;
                  }
              }
              if(ias<quan40) rA_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
              if(ias<quan50) rA_med.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
              if(ias>=quan40 && ias<quan50) rB_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan50 && ias<quan60) rB_50.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan60 && ias<quan70) rB_60.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan70 && ias<quan80) rB_70.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan80 && ias<quan90) rB_80.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan90)              rB_90.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(ias>=quan50 && ias<quan90)rB_50_90.fill(Eta,nhits,p,pt,ih,ias,massT,tof);

          }
          else
          {
              if(ias<=ias_cut) //regionC
              {
                  rC.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nC++;

                  if(pt<70)
                  {
                      //rC_boundedPt.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                      nC_boundedPt++;
                  }

              }
              else //regionD
              {
                  rD.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nD++;
                  h_mT_regD->Fill(mT->at(i),weightMass);
                  h_massObs->Fill(massT,weightMass);
                  
                  if(ias<0.1)
                  {
                      //rD_boundedIas.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
                      nD_boundedIas++;
                  } 
                  
                  if(pt<70)
                  {
                      //rD_boundedPt.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                      nD_boundedPt++;
                  }

              }
              if(ias<quan50) rC_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
              if(ias<quan50) rC_med.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
              if(ias>=quan40 && ias<quan50) {rD_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q40->Fill(massT,weightMass);}
              if(ias>=quan50 && ias<quan60) {rD_50.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q50->Fill(massT,weightMass);}
              if(ias>=quan60 && ias<quan70) {rD_60.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q60->Fill(massT,weightMass);}
              if(ias>=quan70 && ias<quan80) {rD_70.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q70->Fill(massT,weightMass);}
              if(ias>=quan80 && ias<quan90) {rD_80.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q80->Fill(massT,weightMass);}
              if(ias>=quan90)              {rD_90.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q90->Fill(massT,weightMass);}
              if(ias>=quan50 && ias<quan90){rD_50_90.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q50_90->Fill(massT,weightMass);}

          }

          if(ias<=ias_cut){
              if(pt<quan40_pt) rA_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt<quan50_pt) rA_med_pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan40_pt && pt<quan50_pt) rC_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan50_pt && pt<quan60_pt) rC_50pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan60_pt && pt<quan70_pt) rC_60pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan70_pt && pt<quan80_pt) rC_70pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan80_pt && pt<quan90_pt) rC_80pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan90_pt) rC_90pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
          }
          else{
              if(pt<quan40_pt) rB_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt<quan50_pt) rB_med_pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
              if(pt>=quan40_pt && pt<quan50_pt) {rD_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q40_pt->Fill(massT,weightMass);}
              if(pt>=quan50_pt && pt<quan60_pt) {rD_50pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q50_pt->Fill(massT,weightMass);}
              if(pt>=quan60_pt && pt<quan70_pt) {rD_60pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q60_pt->Fill(massT,weightMass);}
              if(pt>=quan70_pt && pt<quan80_pt) {rD_70pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q70_pt->Fill(massT,weightMass);}
              if(pt>=quan80_pt && pt<quan90_pt) {rD_80pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q80_pt->Fill(massT,weightMass);}
              if(pt>=quan90_pt) {rD_90pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof);h_massObs_q90_pt->Fill(massT,weightMass);}

          }

      }

   h_nhscp->Fill(count_hscp_passPre);
   
   }

   ntot*=weightMass;
   nA*=weightMass;
   nB*=weightMass;
   nC*=weightMass;
   nD*=weightMass;

   /*crossHistos(ih_p_from2D,(TH1F*)rD.ih_p->ProjectionX(),(TH1F*)rD.ih_p->ProjectionY());
   crossHistos(ih_p_fromBC,(TH1F*)rC.ih_p->ProjectionX(),(TH1F*)rB.ih_p->ProjectionY());
   mapOfDifferences(diff_D_BC,ih_p_fromD,ih_p_fromBC);
   mapOfDifferences(diff_D_2D,ih_p_fromD,ih_p_from2D);
   mapOfDifferences(diff_BC_2D,ih_p_fromBC,ih_p_from2D);*/

   /*rAll.fillStdDev();
   rA.fillStdDev();
   rB.fillStdDev();
   rC.fillStdDev();
   rD.fillStdDev();
   
   rAll.fillQuantile();
   rA.fillQuantile();
   rB.fillQuantile();
   rC.fillQuantile();
   rD.fillQuantile();

   rAll.fillMassFrom1DTemplatesEtaBinning();
   rA.fillMassFrom1DTemplatesEtaBinning();
   rB.fillMassFrom1DTemplatesEtaBinning();
   rC.fillMassFrom1DTemplatesEtaBinning();
   rD.fillMassFrom1DTemplatesEtaBinning();
*/

      std::cout << "n1: " << n1 << " n2: " << n2 << " n2/n1: " << (float)n2/(float)n1 << " n3: " << n3 << " n3/n1: " << (float)n3/(float)n1 << std::endl;
    
      ofstream ofile((outfilename_+"_normalisations.txt").c_str());
   
    TFile* outfile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving..." << std::endl;
      rAll.write();
    std::cout << " region all saved " << std::endl;
      rA.write();
    std::cout << " region A saved " << std::endl;
      rB.write();
    std::cout << " region B saved " << std::endl;
      rC.write();
    std::cout << " region C saved " << std::endl;
      rD.write();
    std::cout << " region D saved " << std::endl;
/*      rB_boundedIas.write();
    std::cout << " region B_boundedIas saved " << std::endl;
      rD_boundedIas.write();
    std::cout << " region D_boundedIas saved " << std::endl;
      rC_boundedPt.write();
    std::cout << " region C_boundedPt saved " << std::endl;
      rD_boundedPt.write();
    std::cout << " region D_boundedPt saved " << std::endl;
*/

    rB_40.write();
    rB_50.write();
    rB_60.write();
    rB_70.write();
    rB_80.write();
    rB_90.write();
    rB_50_90.write();

    rD_40.write();
    rD_50.write();
    rD_60.write();
    rD_70.write();
    rD_80.write();
    rD_90.write();
    rD_50_90.write();

    rC_40pt.write();
    rC_50pt.write();
    rC_60pt.write();
    rC_70pt.write();
    rC_80pt.write();
    rC_90pt.write();

    rD_40pt.write();
    rD_50pt.write();
    rD_60pt.write();
    rD_70pt.write();
    rD_80pt.write();
    rD_90pt.write();

    rA_med.write();
    rA_med_pt.write();
    rB_med_pt.write();
    rC_med.write();

    rA_40.write();
    rA_40pt.write();
    rB_40pt.write();
    rC_40.write();

    h_mT->Write();
    h_mT_regD->Write();
    h_probQ->Write();
    h_isoPFMuon->Write();
    h_njets->Write();
    h_nhscp->Write();
    h_massObs->Write();
    h_isoPFTk->Write();
    h_ias_outer->Write();
    h_ias->Write();

    h_ECAL->Write();
    h_HCAL->Write();

    h_massObs_q40->Write();
    h_massObs_q50->Write();
    h_massObs_q60->Write();
    h_massObs_q70->Write();
    h_massObs_q80->Write();
    h_massObs_q90->Write();
    h_massObs_q50_90->Write();

    h_pT_isHighPurity->Write();
    h_pT_QualityMask->Write();
    h_pT_isHighPurity_QualityMask->Write();

    h_mass_isHighPurity->Write();
    h_mass_QualityMask->Write();
    h_mass_isHighPurity_QualityMask->Write();

    h_pTerrOverpT_pTerrOverpT2->Write();
    h_pT_pTerrCut->Write();
    h_pT_pTerrCut2->Write();
    h_pT_pTerrCut_pTerrCut2->Write();
    h_pT_noCut_pTerrCut->Write();
    h_mass_pTerrCut->Write();
    h_mass_pTerrCut2->Write();
    h_mass_pTerrCut_pTerrCut2->Write();
    h_mass_noCut_pTerrCut->Write();

    h_pTerrOverpT2->Write();

    h_pT->Write();
    h_pTerrOverpT->Write();
    h_pTerrOverpT2->Write();
    h_pT_pTerrOverpT->Write();
    h_pT_pTerrOverpT2->Write();

    h_pTerrOverpT_0_pT_100->Write();
    h_pTerrOverpT_100_pT_200->Write();
    h_pTerrOverpT_200_pT_300->Write();
    h_pTerrOverpT_300_pT_400->Write();
    h_pTerrOverpT_400_pT_800->Write();
    h_pTerrOverpT_800_pT->Write();


    h_pTerrOverpT2_0_pT_100->Write();
    h_pTerrOverpT2_100_pT_200->Write();
    h_pTerrOverpT2_200_pT_300->Write();
    h_pTerrOverpT2_300_pT_400->Write();
    h_pTerrOverpT2_400_pT_800->Write();
    h_pTerrOverpT2_800_pT->Write();

    h_pT_pTerrOverpT_eta_0p8->Write();
    h_pT_pTerrOverpT_0p8_eta_1p3->Write();
    h_pT_pTerrOverpT_1p3_eta_1p7->Write();
    h_pT_pTerrOverpT_1p7_eta_2p1->Write();



      outfile->Write();
      outfile->Close();


      ofile << "ntot: " << ntot << std::endl;
      ofile << "nA: " << nA << " " << 100*(float)nA/(float)ntot << " %" << std::endl;
      ofile << "nB: " << nB << " " << 100*(float)nB/(float)ntot << " %" << std::endl;
      ofile << "nC: " << nC << " " << 100*(float)nC/(float)ntot << " %" << std::endl;
      ofile << "nD: " << nD << " " << 100*(float)nD/(float)ntot << " %" << std::endl;

      /*ofile << "nB_boundedIas: " << nB_boundedIas << " " << 100*(float)nB_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nB_boundedIas/(float)nB << " % (/regionB)" << std::endl;
      ofile << "nD_boundedIas: " << nD_boundedIas << " " << 100*(float)nD_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedIas/(float)nD << " % (/regionD)" << std::endl;
      ofile << "nC_boundedPt: " << nC_boundedPt << " " << 100*(float)nC_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nC_boundedPt/(float)nC << " % (/regionC)" << std::endl;
      ofile << "nD_boundedPt: " << nD_boundedPt << " " << 100*(float)nD_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedPt/(float)nD << " % (/regionD)" << std::endl;
        */
      ofile.close();

}
