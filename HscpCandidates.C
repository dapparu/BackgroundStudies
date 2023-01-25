#define HscpCandidates_cxx
#include "HscpCandidates.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
//#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"


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

bool isLepton(float pID){
    int id = abs((int)pID);
    if( id == 11 ||
        id == 13 ||
        id == 15 ) return true;
    else return false;
}

bool isRHadron(float pID){
    int id = abs((int)pID);
    if( id == 1000993 ||
        id == 1009113 ||
        id == 1009213 ||
        id == 1009223 ||
        id == 1009313 ||
        id == 1009323 ||
        id == 1009333 ||
        id == 1091114 ||
        id == 1092114 ||
        id == 1092214 ||
        id == 1092224 ||
        id == 1093114 ||
        id == 1093214 ||
        id == 1093224 ||
        id == 1093314 ||
        id == 1093324 ||
        id == 1093334 ||
        id == 1000612 ||
        id == 1000622 ||
        id == 1000632 ||
        id == 1000642 ||
        id == 1000652 ||
        id == 1006113 ||
        id == 1006211 ||
        id == 1006213 ||
        id == 1006223 ||
        id == 1006311 ||
        id == 1006313 ||
        id == 1006321 ||
        id == 1006323 ||
        id == 1006333 ) return true;
    else return false;
}

double deltaR(double eta1,double phi1,double eta2,double phi2){
    double deta=eta1-eta2;
    double dphi=phi1-phi2;
    while(dphi>M_PI)
        dphi-=2*M_PI;
    while(dphi<=-M_PI)
        dphi+=2*M_PI;
    return sqrt(deta*deta+dphi*dphi);
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
   
   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   std::cout << "entries: " << nentries << std::endl;
   nof_event = nentries;

   //nentries /= 100;
   //std::cout << "entries: " << nentries << std::endl;

   int ntot=0, nA=0, nB=0, nC=0, nD=0, nB_boundedIas=0, nD_boundedIas=0, nC_boundedPt=0, nD_boundedPt=0;


   TH1::SetDefaultSumw2(kTRUE);
   TH2::SetDefaultSumw2(kTRUE);

   region rAll("_all",etabins_,ihbins_,pbins_,massbins_);
   //region rA("_regionA",etabins_,ihbins_,pbins_,massbins_);
   //region rB("_regionB",etabins_,ihbins_,pbins_,massbins_);
   //region rC("_regionC",etabins_,ihbins_,pbins_,massbins_);
   //region rD("_regionD",etabins_,ihbins_,pbins_,massbins_);

   region rA_med("_regionA_50",etabins_,ihbins_,pbins_,massbins_);
//   region rA_med_pt("_regionA_med_pt",etabins_,ihbins_,pbins_,massbins_);
//   region rB_med_pt("_regionB_med_pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_med("_regionC_50",etabins_,ihbins_,pbins_,massbins_);

   region rA_80("_regionA_80",etabins_,ihbins_,pbins_,massbins_);
   region rA_90("_regionA_90",etabins_,ihbins_,pbins_,massbins_);
   
   region rC_80("_regionC_80",etabins_,ihbins_,pbins_,massbins_);   
   region rC_90("_regionC_90",etabins_,ihbins_,pbins_,massbins_);   


/*   region rA_40("_regionA_40",etabins_,ihbins_,pbins_,massbins_);
   region rA_40pt("_regionA_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rB_40pt("_regionB_40pt",etabins_,ihbins_,pbins_,massbins_);
   region rC_40("_regionC_40",etabins_,ihbins_,pbins_,massbins_);
*/
   
//   region rB_40("_regionB_40",etabins_,ihbins_,pbins_,massbins_);
   region rB_50("_regionB_50",etabins_,ihbins_,pbins_,massbins_);
   region rB_60("_regionB_60",etabins_,ihbins_,pbins_,massbins_);
   region rB_70("_regionB_70",etabins_,ihbins_,pbins_,massbins_);
   region rB_80("_regionB_80",etabins_,ihbins_,pbins_,massbins_);
   region rB_90("_regionB_90",etabins_,ihbins_,pbins_,massbins_);
   region rB_99("_regionB_99",etabins_,ihbins_,pbins_,massbins_);
   region rB_999("_regionB_999",etabins_,ihbins_,pbins_,massbins_);
   region rB_50_90("_regionB_50_90",etabins_,ihbins_,pbins_,massbins_);
   region rB_50_99("_regionB_50_99",etabins_,ihbins_,pbins_,massbins_);
   region rB_50_999("_regionB_50_999",etabins_,ihbins_,pbins_,massbins_);
   region rB_50_100("_regionB_50_100",etabins_,ihbins_,pbins_,massbins_);

//   region rD_40("_regionD_40",etabins_,ihbins_,pbins_,massbins_);
   region rD_50("_regionD_50",etabins_,ihbins_,pbins_,massbins_);
   region rD_60("_regionD_60",etabins_,ihbins_,pbins_,massbins_);
   region rD_70("_regionD_70",etabins_,ihbins_,pbins_,massbins_);
   region rD_80("_regionD_80",etabins_,ihbins_,pbins_,massbins_);
   region rD_90("_regionD_90",etabins_,ihbins_,pbins_,massbins_);
   region rD_99("_regionD_99",etabins_,ihbins_,pbins_,massbins_);
   region rD_999("_regionD_999",etabins_,ihbins_,pbins_,massbins_);
   region rD_50_90("_regionD_50_90",etabins_,ihbins_,pbins_,massbins_);
   region rD_50_99("_regionD_50_99",etabins_,ihbins_,pbins_,massbins_);
   region rD_50_999("_regionD_50_999",etabins_,ihbins_,pbins_,massbins_);
   region rD_50_100("_regionD_50_100",etabins_,ihbins_,pbins_,massbins_);

/*   region rC_40pt("_regionC_40pt",etabins_,ihbins_,pbins_,massbins_);
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
   region rD_90pt("_regionD_90pt",etabins_,ihbins_,pbins_,massbins_);*/

   //region r_mass800("_reg_mass800",etabins_,ihbins_,pbins_,massbins_);

 /*  region rA_005ias01("_regA_005ias01",etabins_,ihbins_,pbins_,massbins_);
   region rB_005ias01("_regB_005ias01",etabins_,ihbins_,pbins_,massbins_);
   region rC_005ias01("_regC_005ias01",etabins_,ihbins_,pbins_,massbins_);
   region rD_005ias01("_regD_005ias01",etabins_,ihbins_,pbins_,massbins_);
   
   region rA_005ias015("_regA_005ias015",etabins_,ihbins_,pbins_,massbins_);
   region rB_005ias015("_regB_005ias015",etabins_,ihbins_,pbins_,massbins_);
   region rC_005ias015("_regC_005ias015",etabins_,ihbins_,pbins_,massbins_);
   region rD_005ias015("_regD_005ias015",etabins_,ihbins_,pbins_,massbins_);

   region rA_ias50_eta08("_regA_ias50_eta08",etabins_,ihbins_,pbins_,massbins_);
   region rC_ias50_eta08("_regC_ias50_eta08",etabins_,ihbins_,pbins_,massbins_);
   region rB_50ias90_eta08("_regB_50ias90_eta08",etabins_,ihbins_,pbins_,massbins_);
   region rD_50ias90_eta08("_regD_50ias90_eta08",etabins_,ihbins_,pbins_,massbins_);
   TH1F* h_massObs_50ias90_eta08 = new TH1F("massObs_50ias90_eta08","|#eta|<0.8;Mass (GeV)",80,0,4000);
   h_massObs_50ias90_eta08->Sumw2();
   h_massObs_50ias90_eta08->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);

   region rA_ias50_08eta17("_regA_ias50_08eta17",etabins_,ihbins_,pbins_,massbins_);
   region rC_ias50_08eta17("_regC_ias50_08eta17",etabins_,ihbins_,pbins_,massbins_);
   region rB_50ias90_08eta17("_regB_50ias90_08eta17",etabins_,ihbins_,pbins_,massbins_);
   region rD_50ias90_08eta17("_regD_50ias90_08eta17",etabins_,ihbins_,pbins_,massbins_);
   TH1F* h_massObs_50ias90_08eta17 = new TH1F("massObs_50ias90_08eta17","0.8<|#eta|<1.7;Mass (GeV)",80,0,4000);
   h_massObs_50ias90_08eta17->Sumw2();
   h_massObs_50ias90_08eta17->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);


   region rA_ias50_17eta21("_regA_ias50_17eta21",etabins_,ihbins_,pbins_,massbins_);
   region rC_ias50_17eta21("_regC_ias50_17eta21",etabins_,ihbins_,pbins_,massbins_);
   region rB_50ias90_17eta21("_regB_50ias90_17eta21",etabins_,ihbins_,pbins_,massbins_);
   region rD_50ias90_17eta21("_regD_50ias90_17eta21",etabins_,ihbins_,pbins_,massbins_);
   TH1F* h_massObs_50ias90_17eta21 = new TH1F("massObs_50ias90_17eta21","1.7<|#eta|<2.1;Mass (GeV)",80,0,4000);
   h_massObs_50ias90_17eta21->Sumw2();
   h_massObs_50ias90_17eta21->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);*/

   /*region rA_sc1("_regA_sc1",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc1("_regB_sc1",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc1("_regC_sc1",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc1("_regD_sc1",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc2("_regA_sc2",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc2("_regB_sc2",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc2("_regC_sc2",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc2("_regD_sc2",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc3("_regA_sc3",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc3("_regB_sc3",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc3("_regC_sc3",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc3("_regD_sc3",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc4("_regA_sc4",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc4("_regB_sc4",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc4("_regC_sc4",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc4("_regD_sc4",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc5("_regA_sc5",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc5("_regB_sc5",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc5("_regC_sc5",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc5("_regD_sc5",etabins_,ihbins_,pbins_,massbins_);
   
   region rA_sc6("_regA_sc6",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc6("_regB_sc6",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc6("_regC_sc6",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc6("_regD_sc6",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc7("_regA_sc7",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc7("_regB_sc7",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc7("_regC_sc7",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc7("_regD_sc7",etabins_,ihbins_,pbins_,massbins_);

   region rA_sc8("_regA_sc8",etabins_,ihbins_,pbins_,massbins_);
   region rB_sc8("_regB_sc8",etabins_,ihbins_,pbins_,massbins_);
   region rC_sc8("_regC_sc8",etabins_,ihbins_,pbins_,massbins_);
   region rD_sc8("_regD_sc8",etabins_,ihbins_,pbins_,massbins_);
*/

   //region rA_1o4("_regA_1o4",etabins_,ihbins_,pbins_,massbins_);
   //region rA_2o4("_regA_2o4",etabins_,ihbins_,pbins_,massbins_);
   //region rA_3o4("_regA_3o4",etabins_,ihbins_,pbins_,massbins_);
   //region rA_4o4("_regA_4o4",etabins_,ihbins_,pbins_,massbins_);
   //region rB_1o2("_regB_1o2",etabins_,ihbins_,pbins_,massbins_);
   //region rB_2o2("_regB_2o2",etabins_,ihbins_,pbins_,massbins_);
   //region rC_1o2("_regC_1o2",etabins_,ihbins_,pbins_,massbins_);
   //region rC_2o2("_regC_2o2",etabins_,ihbins_,pbins_,massbins_);
   //region rD_1o1("_regD_1o1",etabins_,ihbins_,pbins_,massbins_);

   //region rB_true("_regB_true",etabins_,ihbins_,pbins_,massbins_);
   //region rC_true("_regC_true",etabins_,ihbins_,pbins_,massbins_);
   
   //region ("",etabins_,ihbins_,pbins_,massbins_);



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
   //double q[6]={ 0.033863945, 0.039237098, 0.045422508, 0.053164483, 0.063867635, 0.082214175 }; //data or signal
   //double q[6]={ 0.032398762, 0.037014346, 0.041986481, 0.047712094, 0.054792688, 0.065430914 }; //data or signal -- IAS STRIP ONLY
   
   //double q[6]={ 0.032558433, 0.037216634, 0.042322464, 0.048092840, 0.055574009, 0.066387113 }; //data or signal -- IAS STRIP ONLY NO FSTRIP CUT
   //double q[6] = { 0.030369465, 0.035503496, 0.041515836, 0.048954613, 0.059305365, 0.077425112 }; //MC TTTo2L2Nu
   double q[8]={ 0.014565036, 0.017987774, 0.022399569, 0.028518069, 0.038047370, 0.056746799, 0.13331622, 0.22018057 }; //data or signal -- IAS STRIP ONLY NO FSTRIP CUT
   //double q[5] = { 0.038786767, 0.045297877, 0.054489588, 0.069096781, 0.097055567 }; //MC QCD
   double q_pt[6]={ 57.013944, 59.648194, 63.123116, 67.930817, 75.501106, 90.672101 }; //data or signal
   
   float quan40 = q[0];
   float quan50 = q[1];
   float quan60 = q[2];
   float quan70 = q[3];
   float quan80 = q[4];
   float quan90 = q[5];
   float quan99 = q[6];
   float quan999 = q[7];


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
   TH1F* h_massObs_q99 = new TH1F("massObs_q99",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q999 = new TH1F("massObs_q999",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_90 = new TH1F("massObs_q50_90",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_99 = new TH1F("massObs_q50_99",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_999 = new TH1F("massObs_q50_999",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_100 = new TH1F("massObs_q50_100",";Mass (GeV)",80,0,4000);


   TH1F* h_massObs_q40_pt = new TH1F("massObs_q40_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q50_pt = new TH1F("massObs_q50_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q60_pt = new TH1F("massObs_q60_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q70_pt = new TH1F("massObs_q70_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q80_pt = new TH1F("massObs_q80_pt",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_q90_pt = new TH1F("massObs_q90_pt",";Mass (GeV)",80,0,4000);

   TH1F* h_massObs_005ias01 = new TH1F("massObs_005ias01",";Mass (GeV)",80,0,4000);
   TH1F* h_massObs_005ias015 = new TH1F("massObs_005ias015",";Mass (GeV)",80,0,4000);

   if(dataset_ == "data") {
   h_massObs_q40->Sumw2();
   h_massObs_q50->Sumw2();
   h_massObs_q60->Sumw2();
   h_massObs_q70->Sumw2();
   h_massObs_q80->Sumw2();
   h_massObs_q90->Sumw2();
   h_massObs_q50_90->Sumw2();
   h_massObs_q50_100->Sumw2();
   
   h_massObs_q40_pt->Sumw2();
   h_massObs_q50_pt->Sumw2();
   h_massObs_q60_pt->Sumw2();
   h_massObs_q70_pt->Sumw2();
   h_massObs_q80_pt->Sumw2();
   h_massObs_q90_pt->Sumw2();
   }

   h_massObs_005ias01->Sumw2();
   h_massObs_005ias015->Sumw2();

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
   TH2F* h_pT_pTerrOverpT_0p8_eta_1p7 = new TH2F("pT_pTerrOverpT_0p8_eta_1p7","0.8<|#eta|<1.7;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
   TH2F* h_pT_pTerrOverpT_1p7_eta_2p1 = new TH2F("pT_pTerrOverpT_1p7_eta_2p1","1.7<|#eta|<2.1;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);
  
   TH2F* h_mass_pTerrOverpT = new TH2F("mass_pTerrOverpT",";Mass (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",80,0,2000,1000,0,1);
   TH2F* h_ias_pTerrOverpT = new TH2F("ias_pTerrOverpT",";I_{as};#frac{#sigma_{pT}}{p_{T}};tracks/bin",100,0,1,1000,0,1);
   TH2F* h_pT_pTerrOverpT_rhad = new TH2F("pT_pTerrOverpT_rhad","isRHadron;p_{T} (GeV);#frac{#sigma_{pT}}{p_{T}};tracks/bin",1000,0,6000,1000,0,1);

   TH2F* h_mass_fracLepton = new TH2F("mass_fracLepton",";Mass (GeV);is leptons;tracks/bin",80,0,2000,2,0,2);
   
   TH2F* h_pT_PULLpT = new TH2F("pT_PULLpT","DeltaR(RECO,GEN)<0.1;p_{T} (GeV);#frac{RECO pT-GEN pT}{#sigma_{pT}}",100,0,1000,100,-5,5);
   TH1F* h_PULLpT = new TH1F("PULLpT","DeltaR(RECO,GEN)<0.1;#frac{RECO pT-GEN pT}{#sigma_{pT}}",1000,-100,100);
  
   TH1F* h_N_1_miniIso = new TH1F("h_N-1_miniIso","N-1;Relative mini iso;tracks/bin",100,0,10);
   TH1F* h_N_1_miniIso_wMuon = new TH1F("h_N-1_miniIso_wMuon","N-1;Relative mini iso with muons contribution;tracks/bin",100,0,10);
   TH1F* h_N_1_TkBased_absolute = new TH1F("h_N-1_TkBased_absolute","N-1;Absolute tracker based isolation (GeV);tracks/bin",200,0,100);
   TH1F* h_N_1_TkBased_relative = new TH1F("h_N-1_TkBased_relative","N-1;Relative tracker based isolation;tracks/bin",100,0,10);
   TH1F* h_Ias_miniIso = new TH1F("h_Ias_miniIso","mini iso cut;I_{as};tracks/bin",100,0,1);
   TH1F* h_Ias_miniIso_wMuon = new TH1F("h_Ias_miniIso_wMuon","mini iso + muon cut;I_{as};tracks/bin",100,0,1);
   TH1F* h_Ias_TkBased_absolute = new TH1F("h_Ias_TkBased_absolute","Absolute Tk-based iso cut;I_{as};tracks/bin",100,0,1);
   TH1F* h_Ias_TkBased_relative = new TH1F("h_Ias_TkBased_relative","Relative Tk-based iso cut;I_{as};tracks/bin",100,0,1);
   TH1F* h_Ias_preSIso = new TH1F("h_Ias_preSIso","Absolute Tk-based iso cut + relative calo iso cut;I_{as};tracks/bin",100,0,1);
   TH1F* h_Ih_miniIso = new TH1F("h_Ih_miniIso","mini iso cut;I_{h} (MeV/cm);tracks/bin",100,0,10);
   TH1F* h_Ih_miniIso_wMuon = new TH1F("h_Ih_miniIso_wMuon","mini iso + muon cut;I_{h} (MeV/cm);tracks/bin",100,0,10);
   TH1F* h_Ih_TkBased_absolute = new TH1F("h_Ih_TkBased_absolute","Absolute Tk-based iso cut;I_{h} (MeV/cm);tracks/bin",100,0,10);
   TH1F* h_Ih_TkBased_relative = new TH1F("h_Ih_TkBased_relative","Relative Tk-based iso cut;I_{h} (MeV/cm);tracks/bin",100,0,10);
   TH1F* h_Ih_preSIso = new TH1F("h_Ih_preSIso","Absolute Tk-based iso cut + relative calo iso cut;I_{h} (MeV/cm);tracks/bin",100,0,10);
   TH1F* h_pT_miniIso = new TH1F("h_pT_miniIso","mini iso cut;p_{T} (GeV);tracks/bin",100,0,500);
   TH1F* h_pT_miniIso_wMuon = new TH1F("h_pT_miniIso_wMuon","mini iso + muon cut;p_{T} (GeV);tracks/bin",100,0,500);
   TH1F* h_pT_TkBased_absolute = new TH1F("h_pT_TkBased_absolute","Absolute Tk-based iso cut;p_{T} (GeV);tracks/bin",100,0,500);
   TH1F* h_pT_TkBased_relative = new TH1F("h_pT_TkBased_relative","Relative Tk-based iso cut;p_{T} (GeV);tracks/bin",100,0,500);
   TH1F* h_pT_preSIso = new TH1F("h_pT_preSIso","Absolute Tk-based iso cut + relative calo iso cut;p_{T} (GeV);tracks/bin",100,0,500);

   TH1F* h_FstripBef = new TH1F("h_FStripBef","Fraction of strips passing ClusterCleaning (before preselection);#frac{Nof strip measurements}{Nof valid strip hits};tracks/bin",20,0,1);
   TH1F* h_Fstrip = new TH1F("h_FStrip","Fraction of strips passing ClusterCleaning;#frac{Nof strip measurements}{Nof valid strip hits};tracks/bin",20,0,1);

   TH1F* h_vetoJetBef = new TH1F("h_vetoJetBef","Before preselection;veto jet (p_{T}>30 GeV) dR<0.3;tracks/bin",2,0,2);
   TH1F* h_vetoJet = new TH1F("h_vetoJet",";veto jet (p_{T}>30 GeV) dR<0.3;tracks/bin",2,0,2);

   TH1F* h_dRjetMinBef = new TH1F("h_dRjetMinBef","Minimum deltaR(HSCP,jet) - Before preselection;min dR(hscp,jet[p_{T}>30 GeV]);tracks/bin",100,0,2);
   TH1F* h_dRjetMin = new TH1F("h_dRjetMin","Minimum deltaR(HSCP,jet);min dR(hscp,jet[p_{T}>30 GeV]);tracks/bin",100,0,2);
   
   TH1F* h_dRjetHighestPtBef = new TH1F("h_dRjetHighestPtBef","Jet with highest p_{T} - Before preselection;dR(hscp,jet[p_{T} max]);tracks/bin",100,0,2);
   TH1F* h_dRjetHighestPt = new TH1F("h_dRjetHighestPt","Jet with highest p_{T};dR(hscp,jet[p_{T} max]);tracks/bin",100,0,2);


   TH2F* h_ias_FstripBef = new TH2F("h_ias_FstripBef","Before preselection;#frac{Nof strip measurements}{Nof valid strip hits};I_{as}",20,0,1,100,0,1);
   TH2F* h_ias_Fstrip = new TH2F("h_ias_Fstrip",";#frac{Nof strip measurements}{Nof valid strip hits};I_{as}",20,0,1,100,0,1);

   TH2F* h_ias_dRjetMinBef = new TH2F("h_ias_dRjetMinBef","Before preselection;min dR(hscp,jet[p_{T}>30 GeV]);I_{as}",100,0,2,100,0,1);
   TH2F* h_ias_dRjetMin = new TH2F("h_ias_dRjetMin",";min dR(hscp,jet[p_{T}>30 GeV]);I_{as}",100,0,2,100,0,1);
   
   TH2F* h_ias_dRjetHighestPtBef = new TH2F("h_ias_dRjetHighestPtBef","Before preselection;dR(hscp,jet[p_{T} max]);I_{as}",100,0,2,100,0,1);
   TH2F* h_ias_dRjetHighestPt = new TH2F("h_ias_dRjetHighestPt",";dR(hscp,jet[p_{T} max]);I_{as}",100,0,2,100,0,1);

   TH2F* h_PtErrOverPt_Fstrip = new TH2F("h_PtErrOverPt_Fstrip",";F^{strip}=#frac{Nof strip measurements}{Nof valid strip hits};#frac{#sigma_{p_{T}}}{p_{T}}",50,0,1,100,0,0.5);
   TH2F* h_DiffGenPt_Fstrip = new TH2F("h_DiffGenPt_Fstrip",";F^{strip}=#frac{Nof strip measurements}{Nof valid strip hits};#frac{GEN p_{T} - RECO p_{T}}{GEN p_{T}}",50,0,1,100,0,2);

   TH1F* h_Fstrip_regionB = new TH1F("h_Fstrip_regionB",";F^{strip}=#frac{Nof strip measurements}{Nof valid strip hits};",50,0,1);

   TH1F* h_Is = new TH1F("h_Is",";I_{S};",100,0,1);
   TH2F* h_Is_Ias = new TH2F("h_Is_Ias",";I_{as};I_{S}",100,0,1,100,0,1);

    TH1F* nCandidates = new TH1F("nCandidates",";Number of candidates;",10,0,10);
    TH1F* nCandidates_BefPreS = new TH1F("nCandidates_BefPreS","Before PreS;Number of candidates;",10,0,10);
    TH1F* nCandidates_AfterPreS = new TH1F("nCandidates_AfterPreS","After PreS;Number of candidates;",10,0,10);
    TH1F* nCandidates_AfterPreS_nMinus1 = new TH1F("nCandidates_AfterPreS_nMinus1","After PreS n-1;Number of candidates;",10,0,10);

    TH1F* h_EP = new TH1F("h_EP","",100,0,2);

    TH1F* h_DeltaEtaTwoCandidates = new TH1F("DeltaEtaTwoCandidates",";#Delta#eta(candidate1, candidate2);",500,-2.5,2.5);
    TH1F* h_DeltaPhiTwoCandidates = new TH1F("DeltaPhiTwoCandidates",";#Delta#Phi(candidate1, candidate2);",500,-5,5);
    TH1F* h_AbsDeltaPhiTwoCandidates = new TH1F("AbsDeltaPhiTwoCandidates",";|#Delta#Phi(candidate1, candidate2)|;",500,0,5);
    TH1F* h_DeltaRTwoCandidates = new TH1F("DeltaRTwoCandidates",";#DeltaR(candidate1, candidate2);",500,0,5);
    TH1F* h_InvMassTwoCandidates = new TH1F("InvMassTwoCandidates",";Invariant mass (candidate1, candidate2);",500,0,500);
    TH1F* h_PtTwoCandidates = new TH1F("PtTwoCandidates",";p_{T} (candidate1, candidate2);",500,0,500);
    TH1F* h_isMuonCandidate1 = new TH1F("isMuonCandidate1",";is muon;",2,0,2);
    TH1F* h_isMuonCandidate2 = new TH1F("isMuonCandidate2",";is muon;",2,0,2);
    TH1F* h_MET_TwoCandidates = new TH1F("MET_TwoCandidates",";MET;",300,0,300);
    TH1F* h_MET_TwoCandidates_cutInvM110 = new TH1F("MET_TwoCandidates_cutInvM110",";MET;",300,0,300);

   h_FstripBef->Sumw2();
   h_Fstrip->Sumw2();

   h_vetoJetBef->Sumw2();
   h_vetoJet->Sumw2();

   h_dRjetMinBef->Sumw2();
   h_dRjetMin->Sumw2();

   h_dRjetHighestPtBef->Sumw2();
   h_dRjetHighestPt->Sumw2();

   h_ias_FstripBef->Sumw2();

   h_N_1_miniIso->Sumw2();
   h_N_1_miniIso_wMuon->Sumw2();
   h_N_1_TkBased_absolute->Sumw2();
   h_N_1_TkBased_relative->Sumw2();
   h_Ias_miniIso->Sumw2();
   h_Ias_miniIso_wMuon->Sumw2();
   h_Ias_TkBased_absolute->Sumw2();
   h_Ias_TkBased_relative->Sumw2();
   h_Ias_preSIso->Sumw2();
   h_Ih_miniIso->Sumw2();
   h_Ih_miniIso_wMuon->Sumw2();
   h_Ih_TkBased_absolute->Sumw2();
   h_Ih_TkBased_relative->Sumw2();
   h_Ih_preSIso->Sumw2();
   h_pT_miniIso->Sumw2();
   h_pT_miniIso_wMuon->Sumw2();
   h_pT_TkBased_absolute->Sumw2();
   h_pT_TkBased_relative->Sumw2();
   h_pT_preSIso->Sumw2();

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
   h_pT_pTerrOverpT_0p8_eta_1p7->Sumw2();
   h_pT_pTerrOverpT_1p7_eta_2p1->Sumw2();
   h_pT_pTerrOverpT_rhad->Sumw2();

   h_mass_fracLepton->Sumw2();
   h_mass_pTerrOverpT->Sumw2();
   h_ias_pTerrOverpT->Sumw2();

   h_massObs->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   
   h_massObs_q40->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q60->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q70->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q80->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q90->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q99->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q999->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50_90->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50_99->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50_999->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50_100->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   
   h_massObs_q40_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q50_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q60_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q70_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q80_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_q90_pt->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   
   h_massObs_005ias01->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h_massObs_005ias015->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);

   //mini_region MR_Bef("_BefPreS");
   //mini_region MR_PostPreS_noFstrip("_PostPreS_noFstrip");
   //mini_region MR_PostPreS_Eta2p1("_PostPreS_noFstrip_eta2p1");
   //mini_region MR_all("_PostPreS");
   //mini_region MR_ih_l_5("_Ih_l_5");
   //mini_region MR_ih_ge_5("_Ih_ge_5");
   //mini_region MR_HCAL_pic6("_HCAL_pic6");
   //mini_region MR_pterr0p25("_pterr0p25");
   //mini_region MR_cutQ99("_cutQ99");
   //mini_region MR_isotk50("_isotk50");
   //mini_region MR_isotk15("_isotk15");
   //mini_region MR_mini_iso_03("_mini_iso_03");
   //mini_region MR_ias_l_0p6("_Ias_l_0p6");
   //mini_region MR_ias_ge_0p6("_Ias_ge_0p6");MR_ias_ge_0p6.setOutputFileName(dataset_);

   //preS_analysis preS_BefPreS("_BefPreS");
   //preS_BefPreS.loadCuts(.25,15.,0.3,3.47,1.,0.5,0.3);
  //
   //preS_analysis preS_default_noMiniIso("_default_noMiniIso");
   //preS_default_noMiniIso.loadCuts(.25,15.,999,3.47,1.,0.5,0.3);
//
   //preS_analysis preS_default_noTKIso("_default_noTKIso");
   //preS_default_noTKIso.loadCuts(.25,999,0.3,3.47,1.,0.5,0.3);
//
   //preS_analysis preS_default_cutQ99("_default_cutQ99");
   //preS_default_cutQ99.loadCuts(.25,15.,999,3.47,1.,0.5,0.3);
//
   //preS_analysis preS_default_cutQ999("_default_cutQ999");
   //preS_default_cutQ999.loadCuts(.25,15.,999,3.47,1.,0.5,0.3);

   //float lumi=7.04; //fb-1
   float lumi=51.2; //fb-1  //2018


   std::string year_xsec = "th";

   float weightGluino2600 = (5.0e-2*lumi)/nof_event; 
   float weightGluino2400 = (1.8e-2*lumi)/nof_event; 
   float weightGluino2200 = (4.9e-1*lumi)/nof_event; 
   float weightGluino2000 = (9.7e-1*lumi)/nof_event; 
   float weightGluino1800 = (7.5e-3*1e3*lumi)/nof_event; //Gluino xsec from 2015 paper --> not in 2017 AN
   float weightGluino1600 = (7.5e-3*1e3*lumi)/nof_event; //Gluino xsec from 2015 paper --> not in 2017 AN
   float weightGluino1400 = (2.5e1*lumi)/nof_event;
   float weightGluino1000 = (3.2e2*lumi)/nof_event;
   float weightGluino1600on = (8.8e-3*1e3*lumi)/nof_event;//Gluino xsec from 2015 paper --> not in 2017 AN

   //Analysis Note 2017 
   if(year_xsec=="2017"){
   weightGluino2600 = (5.0e-2*lumi)/nof_event; 
   weightGluino2000 = (9.7e-1*lumi)/nof_event; 
   weightGluino1600 = (7.5e-3*1e3*lumi)/nof_event; 
   weightGluino1400 = (2.5e1*lumi)/nof_event;
   weightGluino1000 = (3.2e2*lumi)/nof_event;
   }

   //2015 paper obs xsec
   if(year_xsec=="2015"){
   weightGluino2000 = (1.1e-2*1e3*lumi)/nof_event; 
   weightGluino1600 = (7.5e-3*1e3*lumi)/nof_event; 
   weightGluino1400 = (6.55e-3*1e3*lumi)/nof_event;
   weightGluino1000 = (5.55e-3*1e3*lumi)/nof_event;
   }

   //2016 PAS obs xsec
   if(year_xsec=="2016"){
   weightGluino2000 = (2.3e-3*1e3*lumi)/nof_event; 
   weightGluino1600 = (1.6e-3*1e3*lumi)/nof_event; 
   weightGluino1400 = (1.4e-3*1e3*lumi)/nof_event;
   weightGluino1000 = (1.15e-3*1e3*lumi)/nof_event;
   }


   float weightStau1599 = (1.4e-4*lumi)/nof_event;
   
   float weightgmStau200 = (2.8e-1*lumi)/nof_event;
   float weightgmStau432 = (3.9e-1*lumi)/nof_event;
   float weightgmStau651 = (6.9e-2*lumi)/nof_event;
   float weightgmStau745 = (6.9e-2*lumi)/nof_event;
   float weightgmStau871 = (6.9e-2*lumi)/nof_event;
   float weightgmStau1029 = (2.2e-2*lumi)/nof_event;


   float weightppStau871 = (9.9e-3*lumi)/nof_event;
   float weightppStau1029 = (3.5e-3*lumi)/nof_event;
   
   //Analysis Note 2017
   if(year_xsec=="2017"){
   weightppStau871 = (9.9e-3*lumi)/nof_event;
   weightppStau1029 = (3.5e-3*lumi)/nof_event;
   }

   //2015 paper obs xsec
   if(year_xsec=="2015"){
   weightppStau1029 = (2.0e-3*1e3*lumi)/nof_event;
   }

   //2016 PAS obs xsec
   if(year_xsec=="2016"){
   weightppStau1029 = (4.9e-4*1e3*lumi)/nof_event;
   }
   


   float weightWJets = (52940*1e3*lumi)/nof_event;
   float weightTTToSemiLeptonic = (365.35*1e3*lumi)/nof_event;
   float weightTTToHadronic = (377.96*1e3*lumi)/nof_event;
   float weightTTTo2L2Nu = (88.29*1e3*lumi)/(nof_event);
   float weightQCD = (239000.0*1e3*lumi)/nof_event;
   float weightDYJetsToLL = (15343.0*1e3*lumi)/nof_event;

   float weightQCD_pt15_7000 = (1375000000.0*1e3*lumi)/nof_event;

   //float weightQCD_Pt50to80 = (377800.0*1e3*lumi)/40377957.;
   float weightQCD_Pt30to50 = (1361000.0*1e3*lumi)/nof_event;
   float weightQCD_Pt50to80 = (377800.0*1e3*lumi)/nof_event;
   float weightQCD_Pt80to120 = (88620.0*1e3*lumi)/nof_event;
   float weightQCD_Pt120to170 = (21070.0*1e3*lumi)/nof_event;
   float weightQCD_Pt170to300 = (7019.0*1e3*lumi)/nof_event;
   float weightQCD_Pt300to470 = (622.4*1e3*lumi)/nof_event;
   float weightQCD_Pt470to600 = (58.86*1e3*lumi)/nof_event;
   float weightQCD_Pt600to800 = (18.22*1e3*lumi)/nof_event;
   float weightQCD_Pt800to1000 = (3.25*1e3*lumi)/nof_event;
   //float weightQCD_Pt1000 = (1.068*1e3*lumi)/nof_event; //??? xsec != 1 
   float weightQCD_Pt1000 = (1.613*1e3*lumi)/nof_event; //??? xsec != 1 

   /*
   //if(year_xsec=="th"){
   weightGluino1600 = (8.87e-3*1e3*lumi)/nof_event; 
   weightGluino1800 = (2.93e-3*1e3*lumi)/nof_event; 
   weightGluino2000 = (9.7e-4*1e3*lumi)/nof_event; 
   weightGluino2200 = (3.6e-4*1e3*lumi)/nof_event; 
   weightGluino2400 = (1.3e-4*1e3*lumi)/nof_event; 
   weightGluino2600 = (1.3e-4*1e3*lumi)/nof_event; 
   weightgmStau651 = (4.1e-4*1e3*lumi)/nof_event; 
   weightgmStau745 = (1.9e-4*1e3*lumi)/nof_event; 
   weightgmStau871 = (6.9e-5*1e3*lumi)/nof_event; 
   weightgmStau1029 = (2.2e-5*1e3*lumi)/nof_event; 
   //weight = (*1e3*lumi)/nof_event; 
   //}
   */

   weightGluino1600 = (8.87e-3*1e3*lumi)/nof_event; 
   weightGluino1800 = (2.93e-3*1e3*lumi)/nof_event; 
   weightGluino2000 = (1.01e-3*1e3*lumi)/nof_event; 
   weightGluino2200 = (3.56e-4*1e3*lumi)/nof_event; 
   weightGluino2400 = (1.28e-4*1e3*lumi)/nof_event; 
   weightGluino2600 = (4.62e-5*1e3*lumi)/nof_event; 
   
   weightgmStau200 = (2.8e-1*1e3*lumi)/nof_event;
   weightgmStau432 = (3.9e-3*1e3*lumi)/nof_event;
   weightgmStau651 = (4.1e-4*1e3*lumi)/nof_event; 
   weightgmStau745 = (1.9e-4*1e3*lumi)/nof_event; 
   weightgmStau871 = (6.9e-5*1e3*lumi)/nof_event; 
   weightgmStau1029 = (2.2e-5*1e3*lumi)/nof_event; 

   
   float weightMass = 1;
  
   if(dataset_ == "data") weightMass = 1;
   if(dataset_ == "gluino1000") weightMass = weightGluino1000;
   if(dataset_ == "gluino1400") weightMass = weightGluino1400;
   if(dataset_ == "gluino1600") weightMass = weightGluino1600;
   if(dataset_ == "gluino1800") weightMass = weightGluino1800;
   if(dataset_ == "gluino2000") weightMass = weightGluino2000;
   if(dataset_ == "gluino2200") weightMass = weightGluino2200;
   if(dataset_ == "gluino2400") weightMass = weightGluino2400;
   if(dataset_ == "gluino2600") weightMass = weightGluino2600;
   if(dataset_ == "gluino1600on") weightMass = weightGluino1600on;
   if(dataset_ == "gmStau200") weightMass = weightgmStau200;
   if(dataset_ == "gmStau432") weightMass = weightgmStau432;
   if(dataset_ == "gmStau651") weightMass = weightgmStau651;
   if(dataset_ == "gmStau745") weightMass = weightgmStau745;
   if(dataset_ == "gmStau871") weightMass = weightgmStau871;
   if(dataset_ == "gmStau1029") weightMass = weightgmStau1029;
   //if(dataset_ == "") weightMass = weight;

   if(dataset_ == "WJets") weightMass = weightWJets;
   if(dataset_ == "TTTo2L2Nu") weightMass = weightTTTo2L2Nu;
   if(dataset_ == "TTToSemiLept") weightMass = weightTTToSemiLeptonic;
   if(dataset_ == "TTToHadr") weightMass = weightTTToHadronic;
   if(dataset_ == "TTToFullLept") weightMass = weightTTTo2L2Nu+weightTTToSemiLeptonic;
   if(dataset_ == "QCD") weightMass = weightQCD;
   if(dataset_ == "QCD_pt_15_7000") weightMass = weightQCD_pt15_7000;
   if(dataset_ == "DY") weightMass = weightDYJetsToLL;

   if(dataset_ == "QCD_Pt30to50") weightMass = weightQCD_Pt30to50;
   if(dataset_ == "QCD_Pt50to80") weightMass = weightQCD_Pt50to80;
   if(dataset_ == "QCD_Pt80to120") weightMass = weightQCD_Pt80to120;
   if(dataset_ == "QCD_Pt120to170") weightMass = weightQCD_Pt120to170;
   if(dataset_ == "QCD_Pt170to300") weightMass = weightQCD_Pt170to300;
   if(dataset_ == "QCD_Pt300to470") weightMass = weightQCD_Pt300to470;
   if(dataset_ == "QCD_Pt470to600") weightMass = weightQCD_Pt470to600;
   if(dataset_ == "QCD_Pt600to800") weightMass = weightQCD_Pt600to800;
   if(dataset_ == "QCD_Pt800to1000") weightMass = weightQCD_Pt800to1000;
   if(dataset_ == "QCD_Pt1000") weightMass = weightQCD_Pt1000;

   std::cout << "weight: " << weightMass << std::endl;
   std::cout << "K: " << K << " C: " << C << std::endl;
   
   //cutindex=19 pt=60 ias=0.025
  
   const TrackerTopology* tTopo;

   double pt_cut=ptcut_;
   double ias_cut=iascut_;
  

   int n1=0, n2=0, n3=0, n4=0;
   int n50=0, n95=0, n99=0, n999=0;

   int n_print = 0;

   int n_passPre=0;
   int n_bef=0, n_aft=0;

   //fChain->SetBranchStatus("*",1);
   //fChain->SetBranchStatus("*HLT*",0);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       Long64_t ientryGen = LoadTreeGen(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       if(fChainGen)fChainGen->GetEntry(jentry);

       //if(jentry%2==0) continue;
       //if(jentry%2==1) continue;

       float percentEntry = jentry*100/(float)nentries;

       if(percentEntry>=n_print){
               std::cout << "percent: " << percentEntry << "%" << std::endl;
               //std::cout << "entry: " << jentry << std::endl;
               n_print+=5;
       }
       bool passpre_ev = false;


       //if(jentry%1000000==0) std::cout << "entry: " << jentry << " --> " << percentEntry << std::endl;

       //if(Pt->size()<1) continue;

       //if(invMET_ && RecoPFMET>100) continue;


       if(!HLT_Mu50) continue;
       //if(!HLT_Mu50 && !HLT_PFMET120_PFMHT120_IDTight && !HLT_PFHT500_PFMET100_PFMHT100_IDTight && !HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && !HLT_MET105_IsoTrk50) continue;
       //if(!HLT_PFMET120_PFMHT120_IDTight && !HLT_PFHT500_PFMET100_PFMHT100_IDTight && !HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && !HLT_MET105_IsoTrk50) continue;

       //if(dataset_ == "QCD_pt_15_7000") weightMass *= GeneratorWeight;

       //h_njets->Fill(njets);
       //fChain->SetBranchStatus("*Hscp*",1);
       //fChain->SetBranchStatus("*njets*",1);
       int ntracks = Hscp;
       int nofjets = njets; 

       int count_hscp_passPre = 0;

       int run = Run;
       int ev = Event;
       int lumi = Lumi;
       int pu = PileUp;
       float NPV = nofVtx; 

       int indiceToKeep=-1;
       float maxIh=-1;
        
        //fChain->SetBranchStatus("*Mass*",1);

       /*if(Mass->size()<2) continue;

       nCandidates->Fill(Mass->size());

       for(int cand=0; cand<Mass->size(); cand++){
        float ih = Ih_noL1->at(cand);
        if(ih>maxIh){
            maxIh=ih;
            indiceToKeep=cand;
        }
       }*/
       nCandidates->Fill(Mass->size());
       int nHSCPpassingPreS=0;
       std::vector<candidat> vCand;
        for(int i=0; i<Mass->size(); i++){ 
           float pt = Pt->at(i);
           float sigmapt = PtErr->at(i);
           float ih = Ih_noL1->at(i);
           float ismirnov = -1;
           float ias = Ias_StripOnly->at(i);
           float ias_all = Ias->at(i);
           float ias_stripOnly = Ias_StripOnly->at(i);
           float ias_outer = Ias_noTIBnoTIDno3TEC->at(i);
           float ias_pixelOnly = Ias_PixelOnly->at(i);
           float ias_pixelOnly_noL1 = Ias_PixelOnly_noL1->at(i);
           float probQ = ProbQ_noL1->at(i);
           float Eta = eta->at(i); 
           float Phi = phi->at(i);
           float iso = iso_TK->at(i);
           float iso_r = iso/pt;
           float nhits = NOM->at(i);
           float p = pt*cosh(Eta);
           float massT = GetMass(p,ih,K,C);
           p = 10000./p;
           float isotk = iso_TK->at(i);
           float isocalo = iso_ECAL->at(i)+iso_HCAL->at(i);
           float tof = TOF->at(i);
           float dz = dZ->at(i);
           float dxy = dXY->at(i);
           float probChi2 = TMath::Prob(Chi2->at(i),Ndof->at(i));
           float miniIso = iso_MiniIso->at(i);
           float miniIsoWithMuon = iso_MiniIso_wMuon->at(i);
           miniIso = miniIsoWithMuon;
           float IsoRel = PFIsolationMuon(i);
           float IsoRelTk = PFIsolationTrack(i);
           float EP = EoverP->at(i);
           int nom = NOM->at(i);
           int noh = NOH->at(i);
           int nstrp = 0; 
           float sat254 = 0;
           float sat255 = 0;
           float clusterCleaning = 0;
           int NOSH = NOH->at(i)-NOPH->at(i); 
           int nPix = 0;
           int nPixNoL1_mask = 0;
           int nPixNoL1 = 0;
           int nStrip = 0;
           nom+=nPixNoL1_mask;
           float FStrip = 0;
           if(NOSH>0) FStrip=(nom-nPixNoL1_mask)/(float)NOSH;
           float sigmaptsurpt = (float)sigmapt/(float)pt; 
   
           if(pt<55) continue;
           if(NOPH->at(i)<2) continue;
           if(FOVH->at(i)<0.8) continue;
           if(nom<=9) continue;
           if(Chi2->at(i)/Ndof->at(i)>5) continue;
           if(!isHighPurity->at(i)) continue;
           if(abs(dz)>0.1) continue;
           if(abs(dxy)>0.02) continue;
           if(probQ>0.7) continue;
           if(EP>0.3) continue;
           if(abs(Eta)>1.0) continue;
           //if(ih<C) continue;

           bool bool_isotk = isotk<15 ? true : false;
           bool bool_miniiso = miniIso<0.02 ? true : false;
           bool bool_pterr1 = sigmapt/pt < pTerr_over_pT_etaBin(pt,Eta,99) ? true : false;
           bool bool_pterr2 = sigmapt/pt < pTerr_over_pT_etaBin(pt,Eta,999) ? true : false;
           bool bool_pterr3 = sigmapt/(pt*pt) < 0.0008 ? true : false; 
           
            if(!bool_isotk) continue;
            if(!bool_miniiso) continue;
            if(!bool_pterr3) continue;

            if(ih>maxIh){
                maxIh=ih;
                indiceToKeep=i;
            }

            candidat cand(pt,Eta,Phi,isMuon->at(i));
            vCand.push_back(cand);
            
            count_hscp_passPre++;
            nHSCPpassingPreS++;

        }
        //if(nHSCPpassingPreS>1)cout<<vCand.size()<<" "<<nHSCPpassingPreS<<endl;
        nCandidates_AfterPreS->Fill(nHSCPpassingPreS);
        bool afterZ=false;
        if(vCand.size()==2){
            float phi1=std::fmod(abs(vCand[0].getPhi()),TMath::TwoPi()), phi2=std::fmod(abs(vCand[1].getPhi()),TMath::TwoPi());
            float deltaEta=vCand[0].getEta()-vCand[1].getEta();
            double deltaPhi=phi1-phi2;
            float deltaR1=deltaR(vCand[0].getEta(),vCand[0].getPhi(),vCand[1].getEta(),vCand[1].getPhi());
            TLorentzVector v1;
            TLorentzVector v2;
            TLorentzVector v3;
            float mass_hyp=0.105;
            v1.SetPtEtaPhiM(vCand[0].getPt(),vCand[0].getEta(),vCand[0].getPhi(),mass_hyp);
            v2.SetPtEtaPhiM(vCand[1].getPt(),vCand[1].getEta(),vCand[1].getPhi(),mass_hyp);
            v3=v1+v2;
            float invMass=v3.M();
            h_DeltaEtaTwoCandidates->Fill(deltaEta);
            h_DeltaPhiTwoCandidates->Fill(v1.DeltaPhi(v2));
            h_AbsDeltaPhiTwoCandidates->Fill(abs(v1.DeltaPhi(v2)));
            h_DeltaRTwoCandidates->Fill(v1.DeltaR(v2));
            h_isMuonCandidate1->Fill(vCand[0].getIsMuon());
            h_isMuonCandidate2->Fill(vCand[1].getIsMuon());
            h_MET_TwoCandidates->Fill(RecoPFMET);
            h_InvMassTwoCandidates->Fill(invMass);
            h_PtTwoCandidates->Fill(v3.Pt());
            if(invMass>110) afterZ=true;
            if(afterZ) h_MET_TwoCandidates_cutInvM110->Fill(RecoPFMET);
            //cout<<deltaR1<<" "<<v1.DeltaR(v2)<<" "<<invMass<<endl;
        } 
        
        //if(vCand.size()!=1) continue;
        //if(vCand.size()<2) continue;
        //if(!afterZ) continue;


       for(int i=0; i<Mass->size(); i++){    
            if(i!=indiceToKeep) continue;
           n_bef++;

           float pt = Pt->at(i);
           float sigmapt = PtErr->at(i);
           float ih = Ih_noL1->at(i);
           float ismirnov = -1;
           //float ismirnov = Is->at(i);
           float ias = Ias_StripOnly->at(i);
           float ias_all = Ias->at(i);
           float ias_stripOnly = Ias_StripOnly->at(i);
           float ias_outer = Ias_noTIBnoTIDno3TEC->at(i);
           float ias_pixelOnly = Ias_PixelOnly->at(i);
           float ias_pixelOnly_noL1 = Ias_PixelOnly_noL1->at(i);
           float probQ = ProbQ_noL1->at(i);
           float Eta = eta->at(i); 
           float Phi = phi->at(i);
           float iso = iso_TK->at(i);
           float iso_r = iso/pt;
           float nhits = NOM->at(i);
           float p = pt*cosh(Eta);
           float massT = GetMass(p,ih,K,C);
           p = 10000./p;
           float isotk = iso_TK->at(i);
           float isocalo = iso_ECAL->at(i)+iso_HCAL->at(i);
           float tof = TOF->at(i);
           float dz = dZ->at(i);
           float dxy = dXY->at(i);
           float probChi2 = TMath::Prob(Chi2->at(i),Ndof->at(i));

           float miniIso = iso_MiniIso->at(i);
           float miniIsoWithMuon = iso_MiniIso_wMuon->at(i);

           miniIso = miniIsoWithMuon;

           float IsoRel = PFIsolationMuon(i);
           float IsoRelTk = PFIsolationTrack(i);


           float EP = EoverP->at(i);

           int nom = NOM->at(i);
           int noh = NOH->at(i);
           int nstrp = 0; 
           float sat254 = 0;
           float sat255 = 0;
           float clusterCleaning = 0;

           int NOSH = NOH->at(i)-NOPH->at(i); 

           int nPix = 0;
           int nPixNoL1_mask = 0;
           int nPixNoL1 = 0;
           int nStrip = 0;
           /*for(int j=0; j<clust_charge->at(i).size(); j++){
               float ch = clust_charge->at(i)[j];
               bool cleaned = true;
               //bool cleaned = clust_ClusterCleaning->at(i)[j];
               DetId detid(clust_detid->at(i)[j]);
               bool strip = clust_isStrip->at(i)[j];
               bool pixel = clust_isPixel->at(i)[j];

               if(detid.subdetId() == PixelSubdetector::PixelEndcap || detid.subdetId() == PixelSubdetector::PixelBarrel) nPix++;
               if((detid.subdetId() == PixelSubdetector::PixelEndcap) || (detid.subdetId() == 1 && ((detid >> 20) & 0xF) != 1)) nPixNoL1_mask++;
               //if((detid.subdetId() == PixelSubdetector::PixelEndcap) || (detid.subdetId() == PixelSubdetector::PixelBarrel && tTopo->pxbLayer(detid)!=1)) nPixNoL1++;
               if(detid.subdetId()>=3) nStrip++;


               //if(pixel)std::cout << "ch: " << ch << " cleaned: " << cleaned << " detid: " << detid << " strip: " << strip << " pixel: " << pixel << std::endl;

           }*/
           //std::cout << "nPix: " << nPix << " NOPH: " << NOPH->at(i) << " pixel hits no L1: " << nPixNoL1 << " pixel hits no L1 via mask: " << nPixNoL1_mask << std::endl;

           //std::cout << "nom: " << nom << " noh: " << noh << " nstrip: " << nStrip << " NOSH: " << NOSH << " npixel: " << nPix << " npixelnoL1: " << nPixNoL1_mask << std::endl;

           nom+=nPixNoL1_mask;
          
           //int NOSH = NOH->at(i); 
           float FStrip = 0;
           if(NOSH>0) FStrip=(nom-nPixNoL1_mask)/(float)NOSH;

           //if(FStrip>1)std::cout << "FStrip: " << FStrip << " NOH: " << NOH->at(i) << " NOSH: " << NOSH << " NOM: " << NOM->at(i) << std::endl;

           h_isoPFTk->Fill(IsoRelTk,weightMass);

           h_ias_outer->Fill(Ias_noTIBnoTIDno3TEC->at(i),weightMass);
           h_ias->Fill(ias,weightMass);

           //if(abs(Eta)<=1.7 || abs(Eta)>2.1) continue;

           //if(pt>500)continue;
           
           //std::cout << "miniIso: " << miniIso*pt << std::endl;


           
           bool vetoJet = false;
           float drMin = 9999.;
           float ptjetMin = 0.;
           float drHighestPt = 9999.;
           for(int j=0;j<njets;j++){
               float jetPt = jet_pt->at(j);
               float jetEta = jet_eta->at(j);
               float jetPhi = jet_phi->at(j);
               float jetM = jet_mass->at(j);
               float jetEnergy = jet_energy->at(j);
               float jetPdgId = jet_pdgId->at(j);
               float jetEt = jet_et->at(j);
               float jetChEmEnFraction = jet_chargedEmEnergyFraction->at(j);
               float jetNtEnEnFraction = jet_neutralEmEnergyFraction->at(j);
               
               if(jetPt<30) continue;
               if(jetPt>ptjetMin) {
                   ptjetMin = jetPt;
                   drHighestPt = deltaR(Eta,Phi,jetEta,jetPhi);
               }
               float dR_jet_hscp = deltaR(Eta,Phi,jetEta,jetPhi);
               if(dR_jet_hscp<0.3) vetoJet = true;
               if(dR_jet_hscp<drMin) drMin = dR_jet_hscp;
           }

           h_FstripBef->Fill(FStrip,weightMass);
           h_vetoJetBef->Fill(vetoJet,weightMass);
           h_dRjetMinBef->Fill(drMin,weightMass);
           h_dRjetHighestPtBef->Fill(drHighestPt,weightMass);
           
           //MR_Bef.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all, ias_pixelOnly_noL1,ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh, NOSH,nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //preselection base
           //if(!passPreselection) continue;
           /*if(pt<50) continue;
           if(abs(Eta)>2.1) continue;
           if(NOPH->at(i)<2) continue;
           if(NOH->at(i)<8) continue;
           if(FOVH->at(i)<0.8) continue;
           if(NOM->at(i)<6) continue;
           if(Chi2->at(i)/Ndof->at(i)>5) continue;
           //if(!isMuon->at(i)) continue;
           if(!isHighPurity->at(i)) continue;
           if(isotk>50) continue;*/

           //preselection tighter --> minimal selection 

           float sigmaptsurpt = (float)sigmapt/(float)pt; 
   
           //preS_BefPreS.fill(weightMass,pt,sigmaptsurpt,iso,miniIso,ih,abs(Eta),FStrip,EP);

           //if(!isPFMuon) continue;
           if(pt<55) continue;
           //if(nPixNoL1_mask<=1) continue;
           if(NOPH->at(i)<2) continue;
           //if(NOH->at(i)<10) continue;
           if(FOVH->at(i)<0.8) continue;
           if(nom<=9) continue;
           if(Chi2->at(i)/Ndof->at(i)>5) continue;
           if(!isHighPurity->at(i)) continue;
           if(abs(dz)>0.1) continue;
           if(abs(dxy)>0.02) continue;
           //if(probQ>0.15) continue;
           if(probQ>0.7) continue;
           if(EP>0.3) continue;
           

           //preS_default_noMiniIso.fill(weightMass,pt,sigmaptsurpt,iso,miniIso,ih,abs(Eta),FStrip,EP);
           //preS_default_noTKIso.fill(weightMass,pt,sigmaptsurpt,iso,miniIso,ih,abs(Eta),FStrip,EP);
           //preS_default_cutQ99.changePtErrCut(pt,Eta,99);
           //preS_default_cutQ99.fill(weightMass,pt,sigmaptsurpt,iso,miniIso,ih,abs(Eta),FStrip,EP);
           //preS_default_cutQ999.changePtErrCut(pt,Eta,999);
           //preS_default_cutQ999.fill(weightMass,pt,sigmaptsurpt,iso,miniIso,ih,abs(Eta),FStrip,EP);

           if(abs(Eta)>1.0) continue;
          // if(abs(Eta)<=1.0 || abs(Eta)>2.1) continue;
           
           //if(!isMuon->at(i)) continue;

           //if(isocalo/p>0.3) continue;

           //if(ih<C) continue;
           //if(ih<ihcut_) continue;
           //if(ih>3.17) continue;
           //if(p>pcut_ && pcut_>0) continue;
           //if(vetoJet) continue;

           //if(isotk>15) continue;
           //if(sigmapt/pt > pTerr_over_pT_etaBin(pt,Eta,99)) continue;

           bool bool_isotk = isotk<15 ? true : false;
           bool bool_miniiso = miniIso<0.02 ? true : false;
           bool bool_pterr1 = sigmapt/pt < pTerr_over_pT_etaBin(pt,Eta,99) ? true : false;
           bool bool_pterr2 = sigmapt/pt < pTerr_over_pT_etaBin(pt,Eta,999) ? true : false;
           bool bool_pterr3 = sigmapt/(pt*pt) < 0.0008 ? true : false; 
           
            if(!bool_isotk) continue;
            if(!bool_miniiso) continue;
            if(!bool_pterr3) continue;


            if(i!=indiceToKeep) continue;

            //nCandidates_AfterPreS->Fill(Mass->size());

            h_EP->Fill(EP);


           h_Is->Fill(ismirnov,weightMass);
           h_Is_Ias->Fill(ias,ismirnov,weightMass);

           /*if(pt<65){
               if(pt<60){
                   if(ias<0.064){
                       if(ias<0.053){
                           rA_1o4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                       }
                       else{
                           rA_3o4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                       }
                   }
                   else{
                        rB_1o2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                   }
               }
               else{ 
                   if(ias<0.064){
                       if(ias<0.053){
                           rA_2o4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                       }
                       else{
                           rA_4o4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                       }
                   }
                   else{
                        rB_2o2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                   }
               }
           }
           else{
               if(ias<0.064){
                   if(ias<0.053){
                       rC_1o2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                   }
                   else{
                       rC_2o2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                   }
               }
               else{
                   rD_1o1.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
               }
           }*/


           /*if(pt<=pt_cut){
                if(ias<quan50) {
                    if(bool_isotk && bool_pterr1) rA_sc1.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr2) rA_sc2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr3) rA_sc3.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr2) rA_sc4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr3) rA_sc5.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr1) rA_sc6.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk) rA_sc7.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso) rA_sc8.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                }
                if(ias>=quan50 && ias<quan80) {
                //if(ias>=quan90) {
                    if(bool_isotk && bool_pterr1) rB_sc1.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr2) rB_sc2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr3) rB_sc3.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr2) rB_sc4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr3) rB_sc5.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);    
                    if(bool_miniiso && bool_pterr1) rB_sc6.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk) rB_sc7.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso) rB_sc8.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                }   
           }
           else{

               if(ias<quan50) {
                    if(bool_isotk && bool_pterr1) rC_sc1.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr2) rC_sc2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr3) rC_sc3.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr2) rC_sc4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr3) rC_sc5.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr1) rC_sc6.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk) rC_sc7.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso) rC_sc8.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                }
                if(ias>=quan50 && ias<quan80) {
                //if(ias>=quan90) {
                    if(bool_isotk && bool_pterr1) rD_sc1.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr2) rD_sc2.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk && bool_pterr3) rD_sc3.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr2) rD_sc4.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr3) rD_sc5.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso && bool_pterr1) rD_sc6.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_isotk) rD_sc7.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                    if(bool_miniiso) rD_sc8.fill(Eta,nhits,p,pt,sigmapt,ih,ias,massT,tof,weightMass);
                }
            }*/

           
           //MR_PostPreS_Eta2p1.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh, NOSH,nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           
           //if(abs(Eta)>1) continue;
           
           //MR_PostPreS_noFstrip.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh, NOSH,nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           
           //if(FStrip<0.7) continue;
           
           h_Fstrip->Fill(FStrip,weightMass);
           h_vetoJet->Fill(vetoJet,weightMass);
           h_dRjetMin->Fill(drMin,weightMass);
           h_dRjetHighestPt->Fill(drHighestPt,weightMass);

           h_PtErrOverPt_Fstrip->Fill(FStrip,sigmapt/pt,weightMass);

           /*if(GenMass!= NULL){
               float drminGen = 0.1;
               float drgen=9999.;
               int indexgenmindr=-1;
               for(unsigned int g=0;g<GenMass->size();g++){
                   float pID = GenId->at(g);
                if(GenCharge->at(g)==0) continue;
               h_mass_fracLepton->Fill(massT,isLepton(pID),weightMass);
                //if(abs((int)pID) == 13) continue; // is muon
               float dR=deltaR(Eta,phi->at(i),GenEta->at(g),GenPhi->at(g));
               if(dR>0.1) continue;
               if(dR<drgen) {drgen=drgen; indexgenmindr=g;}
               //std::cout << "deltaR: " << dR << " pT: " << pt << " genPt: " << GenPt->at(g) << std::endl;
               h_pT_PULLpT->Fill(pt,(pt-GenPt->at(g))/PtErr->at(i));
               h_PULLpT->Fill((pt-GenPt->at(g))/PtErr->at(i));
               if(isRHadron(pID)) h_pT_pTerrOverpT_rhad->Fill(pt,PtErr->at(i)/pt,weightMass);
               }
               if(indexgenmindr!=-1){
                    h_DiffGenPt_Fstrip->Fill(FStrip,(GenPt->at(indexgenmindr)-pt)/GenPt->at(indexgenmindr),weightMass);
               }
           }*/

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


/*
           h_pT->Fill(pt,weightMass);
           h_pTerrOverpT->Fill(PtErr->at(i)/pt,weightMass);
           h_pTerrOverpT2->Fill(PtErr->at(i)/pow(pt,2),weightMass);
           h_pT_pTerrOverpT->Fill(pt,PtErr->at(i)/pt,weightMass);
           h_pT_pTerrOverpT2->Fill(pt,PtErr->at(i)/pow(pt,2),weightMass);

           if(abs(Eta)<0.8) h_pT_pTerrOverpT_eta_0p8->Fill(pt,PtErr->at(i)/pt,weightMass);
           if(abs(Eta)<1.3 && abs(Eta)>=0.8) h_pT_pTerrOverpT_0p8_eta_1p3->Fill(pt,PtErr->at(i)/pt,weightMass);
           if(abs(Eta)<1.7 && abs(Eta)>=1.3) h_pT_pTerrOverpT_1p3_eta_1p7->Fill(pt,PtErr->at(i)/pt,weightMass);
           if(abs(Eta)<2.1 && abs(Eta)>=1.7) h_pT_pTerrOverpT_1p7_eta_2p1->Fill(pt,PtErr->at(i)/pt,weightMass);
           if(abs(Eta)<1.7 && abs(Eta)>=0.8) h_pT_pTerrOverpT_0p8_eta_1p7->Fill(pt,PtErr->at(i)/pt,weightMass);


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

           h_mass_pTerrOverpT->Fill(massT,PtErr->at(i)/pt,weightMass);
           h_ias_pTerrOverpT->Fill(ias,PtErr->at(i)/pt,weightMass);

           n1++;
           if(pTerrCut1) n3++;

           if(PtErr->at(i)/pt < pTerr_over_pT_etaBin(pt,Eta,50)) n50++;
           if(PtErr->at(i)/pt < pTerr_over_pT_etaBin(pt,Eta,95)) n95++;
           if(PtErr->at(i)/pt < pTerr_over_pT_etaBin(pt,Eta,99)) n99++;
           if(PtErr->at(i)/pt < pTerr_over_pT_etaBin(pt,Eta,999)) n999++;
           n2++;
*/
//---------------- pT cut --------------------------------------------------------

           //if(PtErr->at(i)/pt > pTerr_over_pT_etaBin(pt,Eta,99)) continue;
           //if(sigmapt/pt > pTerr_over_pT_etaBin(pt,Eta,99)) continue;
           //if(PtErr->at(i)/pt > 0.25) continue;

//--------------------------------------------------------------------------------

           h_mT->Fill(mT->at(i),weightMass);
           h_probQ->Fill(ProbQ_noL1->at(i),weightMass);
           h_isoPFMuon->Fill(PFIsolationMuon(i),weightMass);
           h_ECAL->Fill(ECAL_energy->at(i),weightMass);
           h_HCAL->Fill(HCAL_energy->at(i),weightMass);

           h_N_1_miniIso->Fill(miniIso,weightMass);
           h_N_1_miniIso_wMuon->Fill(miniIsoWithMuon,weightMass);
           h_N_1_TkBased_absolute->Fill(isotk,weightMass);
           h_N_1_TkBased_relative->Fill(isotk/pt,weightMass);


           //if(clust_nstrip) nstrp=clust_nstrip->at(i).size();
           //if(clust_sat254) {for(int j=0; j<clust_sat254->at(i).size(); j++){sat254+=clust_sat254->at(i)[j];} if(clust_sat254->at(i).size()>0) sat254/=clust_sat254->at(i).size();}
           //if(clust_sat255) {for(int j=0; j<clust_sat255->at(i).size(); j++){sat255+=clust_sat255->at(i)[j];} if(clust_sat255->at(i).size()>0) sat255/=clust_sat255->at(i).size();}
           //if(clust_ClusterCleaning) {for(int j=0; j<clust_ClusterCleaning->at(i).size(); j++){clusterCleaning+=clust_ClusterCleaning->at(i)[j];} if(clust_ClusterCleaning->at(i).size()>0) clusterCleaning/=clust_ClusterCleaning->at(i).size();}

           //std::cout << sat255 << std::endl;

           //if(ih<5) MR_ih_l_5.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer, ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(ih>=5) MR_ih_ge_5.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(HCAL_energy->at(i)>=5.5 && HCAL_energy->at(i)<=7) MR_HCAL_pic6.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(sigmapt/pt < 0.25) MR_pterr0p25.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt), FStrip, drMin, drHighestPt;
           //if(sigmapt/pt < pTerr_over_pT_etaBin(pt,Eta,99)) MR_cutQ99.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1,ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(isotk < 50) MR_isotk50.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(isotk < 15) MR_isotk15.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(miniIso<0.3) MR_mini_iso_03.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(ias < 0.6) MR_ias_l_0p6.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if(ias >= 0.6) {MR_ias_ge_0p6.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt); MR_ias_ge_0p6.saveEvent(run,ev,lumi,pu);}
           
           //MR_all.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias_all,ias_pixelOnly_noL1, ias_stripOnly,ias_outer,ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh,NOSH, nstrp, sat254, sat255, clusterCleaning, FStrip, drMin, drHighestPt);
           //if() {MR_.fill(weightMass, Eta, phi->at(i), p, ih, pt, sigmapt, ias, ProbQ_noL1->at(i), isotk, miniIso, mT->at(i), massT, ntracks, nofjets, ECAL_energy->at(i), HCAL_energy->at(i), iso_ECAL->at(i), iso_HCAL->at(i), nom, noh, nstrp, sat254, sat255, clusterCleaning);}


           if(miniIso < 0.3){
               h_Ias_miniIso->Fill(ias,weightMass);
               h_Ih_miniIso->Fill(ih,weightMass);
               h_pT_miniIso->Fill(pt,weightMass);
           }
           if(miniIsoWithMuon < 0.3){
               h_Ias_miniIso_wMuon->Fill(ias,weightMass);
               h_Ih_miniIso_wMuon->Fill(ih,weightMass);
               h_pT_miniIso_wMuon->Fill(pt,weightMass);
           }
           if(isotk < 15){
               h_Ias_TkBased_absolute->Fill(ias,weightMass);
               h_Ih_TkBased_absolute->Fill(ih,weightMass);
               h_pT_TkBased_absolute->Fill(pt,weightMass);
           }
           if(isotk/pt < 0.3){
               h_Ias_TkBased_relative->Fill(ias,weightMass);
               h_Ih_TkBased_relative->Fill(ih,weightMass);
               h_pT_TkBased_relative->Fill(pt,weightMass);
           }
           if(isotk<50 && isocalo/p<0.3){
               h_Ias_preSIso->Fill(ias,weightMass);
               h_Ih_preSIso->Fill(ih,weightMass);
               h_pT_preSIso->Fill(pt,weightMass);
           }




           count_hscp_passPre++;

           //if(massT>800) r_mass800.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);

           //std::cout << "mT: " << mT->at(i) << std::endl;

           rAll.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
           ntot++;

           
           if(pt<=pt_cut)
           {
                /*if(ias<=ias_cut) //regionA
                {
                    rA.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                    nA++;
                }
                else //regionB
                {
                    rB.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                    nB++;
                    h_Fstrip_regionB->Fill(FStrip,weightMass);

                    if(ias<0.1)
                    {
                        //rB_boundedIas.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
                        nB_boundedIas++;
                    }
                }*/
                //if(ias<quan40) rA_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass); 
                if(ias<quan50) rA_med.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                if(ias<quan80) rA_80.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                if(ias<quan90) rA_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                //if(ias>=quan40 && ias<quan50) rB_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(ias>=quan50 && ias<quan60) rB_50.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan60 && ias<quan70) rB_60.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan70 && ias<quan80) rB_70.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan80 && ias<quan90) rB_80.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan90)              rB_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan99)              rB_99.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan999)              rB_999.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan50 && ias<quan90)rB_50_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan50 && ias<quan99)rB_50_99.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan50 && ias<quan999)rB_50_999.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>=quan50)rB_50_100.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                
                /*if(ias<=0.05) rA_005ias01.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>0.05 && ias<=0.1) rB_005ias01.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                
                if(ias<=0.05) rA_005ias015.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>0.05 && ias<=0.15) rB_005ias015.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                
                if(ias<quan50){
                    if(abs(Eta)<0.8) rA_ias50_eta08.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=0.8 && abs(Eta)<1.7) rA_ias50_08eta17.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=1.7 && abs(Eta)<2.1) rA_ias50_17eta21.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                }
                if(ias>=quan50 && ias<quan90) {
                    if(abs(Eta)<0.8) rB_50ias90_eta08.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=0.8 && abs(Eta)<1.7) rB_50ias90_08eta17.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=1.7 && abs(Eta)<2.1) rB_50ias90_17eta21.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                }*/

           }
           else
           {
                /*if(ias<=ias_cut) //regionC
                {
                    rC.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                    nC++;

                    if(pt<70)
                    {
                        //rC_boundedPt.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                        nC_boundedPt++;
                    }
                }
                else //regionD
                {
                    rD.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
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

                }*/
                //if(ias<quan40) rC_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass); 
                if(ias<quan50) rC_med.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                if(ias<quan80) rC_80.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                if(ias<quan90) rC_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); 
                //if(ias>=quan40 && ias<quan50) {rD_40.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q40->Fill(massT,weightMass);}
                if(ias>=quan50 && ias<quan60) {rD_50.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q50->Fill(massT,weightMass);}
                if(ias>=quan60 && ias<quan70) {rD_60.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q60->Fill(massT,weightMass);}
                if(ias>=quan70 && ias<quan80) {rD_70.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q70->Fill(massT,weightMass);}
                if(ias>=quan80 && ias<quan90) {rD_80.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q80->Fill(massT,weightMass);}
                if(ias>=quan90)              {rD_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q90->Fill(massT,weightMass);}
                if(ias>=quan99)              {rD_99.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q99->Fill(massT,weightMass);}
                if(ias>=quan999)              {rD_999.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q999->Fill(massT,weightMass);}
                if(ias>=quan50 && ias<quan90){rD_50_90.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q50_90->Fill(massT,weightMass);}
                if(ias>=quan50 && ias<quan99){rD_50_99.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q50_99->Fill(massT,weightMass);}
                if(ias>=quan50 && ias<quan999){rD_50_999.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q50_999->Fill(massT,weightMass);}
                if(ias>=quan50){rD_50_100.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);h_massObs_q50_100->Fill(massT,weightMass);}
                
                /*if(ias<=0.05) rC_005ias01.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>0.05 && ias<=0.1) {rD_005ias01.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); h_massObs_005ias01->Fill(massT,weightMass);}
                
                if(ias<=0.05) rC_005ias015.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                if(ias>0.05 && ias<=0.15) {rD_005ias015.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); h_massObs_005ias015->Fill(massT,weightMass);}

                if(ias<quan50){
                    if(abs(Eta)<0.8) rC_ias50_eta08.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=0.8 && abs(Eta)<1.7) rC_ias50_08eta17.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                    if(abs(Eta)>=1.7 && abs(Eta)<2.1) rC_ias50_17eta21.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass);
                }
                if(ias>=quan50 && ias<quan90) {
                    if(abs(Eta)<0.8) {rD_50ias90_eta08.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); h_massObs_50ias90_eta08->Fill(massT,weightMass);}
                    if(abs(Eta)>=0.8 && abs(Eta)<1.7) {rD_50ias90_08eta17.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); h_massObs_50ias90_08eta17->Fill(massT,weightMass);}
                    if(abs(Eta)>=1.7 && abs(Eta)<2.1) {rD_50ias90_17eta21.fill(Eta,nhits,p,pt,sigmapt,ih,ias,ismirnov,massT,tof,NPV,weightMass); h_massObs_50ias90_17eta21->Fill(massT,weightMass);}
                }*/
            }
           n_aft++;
           passpre_ev = true;

            /*if(ias<=ias_cut){
                if(pt<quan40_pt) rA_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt<quan50_pt) rA_med_pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan40_pt && pt<quan50_pt) rC_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan50_pt && pt<quan60_pt) rC_50pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan60_pt && pt<quan70_pt) rC_60pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan70_pt && pt<quan80_pt) rC_70pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan80_pt && pt<quan90_pt) rC_80pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan90_pt) rC_90pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
            }
            else{
                if(pt<quan40_pt) rB_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt<quan50_pt) rB_med_pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);
                if(pt>=quan40_pt && pt<quan50_pt) {rD_40pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q40_pt->Fill(massT,weightMass);}
                if(pt>=quan50_pt && pt<quan60_pt) {rD_50pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q50_pt->Fill(massT,weightMass);}
                if(pt>=quan60_pt && pt<quan70_pt) {rD_60pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q60_pt->Fill(massT,weightMass);}
                if(pt>=quan70_pt && pt<quan80_pt) {rD_70pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q70_pt->Fill(massT,weightMass);}
                if(pt>=quan80_pt && pt<quan90_pt) {rD_80pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q80_pt->Fill(massT,weightMass);}
                if(pt>=quan90_pt) {rD_90pt.fill(Eta,nhits,p,pt,ih,ias,massT,tof,weightMass);h_massObs_q90_pt->Fill(massT,weightMass);}

            }*/

  //     }
  // }

        }// end loop mass --> loop over the candidates of the event 'ientry'
           
       h_nhscp->Fill(count_hscp_passPre);
       if(passpre_ev) n_passPre++;


   }// end loop HscpCandidateTree entry

   /*preS_BefPreS.write();
   preS_BefPreS.computeEfficiencies();
   preS_default_noMiniIso.write();
   preS_default_noMiniIso.computeEfficiencies();
   preS_default_noTKIso.write();
   preS_default_noTKIso.computeEfficiencies();
   preS_default_cutQ99.write();
   preS_default_cutQ99.computeEfficiencies();
   preS_default_cutQ999.write();
   preS_default_cutQ999.computeEfficiencies();  */

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

      std::cout << "n1: " << n1 << " n50/n1: " << (float)n50/(float)n1 << " n95/n1: " << (float)n95/(float)n1 << " n99/n1: " << (float)n99/(float)n1 << " n999/n1: " << (float)n999/(float)n1 << std::endl;    
     
      ofstream ofile((outfilename_+"_normalisations.txt").c_str());
   
    TFile* outfile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving..." << std::endl;


    h_EP->Write();

    nCandidates->Write();
    nCandidates_AfterPreS->Write();
    nCandidates_AfterPreS_nMinus1->Write();

    h_DeltaEtaTwoCandidates->Write();
    h_AbsDeltaPhiTwoCandidates->Write();
    h_DeltaPhiTwoCandidates->Write();
    h_DeltaRTwoCandidates->Write();
    h_InvMassTwoCandidates->Write();
    h_PtTwoCandidates->Write();
    h_isMuonCandidate1->Write();
    h_isMuonCandidate2->Write();
    h_MET_TwoCandidates->Write();
    h_MET_TwoCandidates_cutInvM110->Write();

    h_Is->Write();
    h_Is_Ias->Write();

    //rA_1o4.write();
    //rA_2o4.write();
    //rA_3o4.write();
    //rA_4o4.write();
    //rB_1o2.write();
    //rB_2o2.write();
    //rC_1o2.write();
    //rC_2o2.write();
    //rD_1o1.write();
    //rB_true.write();
    //rC_true.write();

    h_vetoJetBef->Write();
    h_vetoJet->Write();
    h_FstripBef->Write();
    h_Fstrip->Write();
    h_dRjetMinBef->Write();
    h_dRjetMin->Write();
    h_dRjetHighestPtBef->Write();
    h_dRjetHighestPt->Write();
    h_PtErrOverPt_Fstrip->Write();
    h_DiffGenPt_Fstrip->Write();
    h_Fstrip_regionB->Write();
   
/*    rA_sc1.write();
    rA_sc2.write();
    rA_sc3.write();
    rA_sc4.write();
    rA_sc5.write();
    rA_sc6.write();
    rA_sc7.write();
    rA_sc8.write();

    rB_sc1.write();
    rB_sc2.write();
    rB_sc3.write();
    rB_sc4.write();
    rB_sc5.write();
    rB_sc6.write();
    rB_sc7.write();
    rB_sc8.write();

    
    rC_sc1.write();
    rC_sc2.write();
    rC_sc3.write();
    rC_sc4.write();
    rC_sc5.write();
    rC_sc6.write();
    rC_sc7.write();
    rC_sc8.write();

    
    rD_sc1.write();
    rD_sc2.write();
    rD_sc3.write();
    rD_sc4.write();
    rD_sc5.write();
    rD_sc6.write();
    rD_sc7.write();
    rD_sc8.write();
*/

    //r_mass800.write();

      rAll.write();
//    std::cout << " region all saved " << std::endl;
      //rA.write();
//    std::cout << " region A saved " << std::endl;
      //rB.write();
//    std::cout << " region B saved " << std::endl;
      //rC.write();
//    std::cout << " region C saved " << std::endl;
      //rD.write();
//    std::cout << " region D saved " << std::endl;
/*      rB_boundedIas.write();
    std::cout << " region B_boundedIas saved " << std::endl;
      rD_boundedIas.write();
    std::cout << " region D_boundedIas saved " << std::endl;
      rC_boundedPt.write();
    std::cout << " region C_boundedPt saved " << std::endl;
      rD_boundedPt.write();
    std::cout << " region D_boundedPt saved " << std::endl;
*/

/*      rA_ias50_eta08.write();
      rC_ias50_eta08.write();
      rB_50ias90_eta08.write();
      rD_50ias90_eta08.write();
      h_massObs_50ias90_eta08->Write();

      rA_ias50_08eta17.write();
      rC_ias50_08eta17.write();
      rB_50ias90_08eta17.write();
      rD_50ias90_08eta17.write();
      h_massObs_50ias90_08eta17->Write();

      rA_ias50_17eta21.write();
      rC_ias50_17eta21.write();
      rB_50ias90_17eta21.write();
      rD_50ias90_17eta21.write();
      h_massObs_50ias90_17eta21->Write();

      rA_005ias01.write();
      rB_005ias01.write();
      rC_005ias01.write();
      rD_005ias01.write();
      h_massObs_005ias01->Write();

      rA_005ias015.write();
      rB_005ias015.write();
      rC_005ias015.write();
      rD_005ias015.write();
      h_massObs_005ias015->Write();*/

//    rB_40.write();
    rB_50.write();
    rB_60.write();
    rB_70.write();
    rB_80.write();
    rB_90.write();
    rB_99.write();
    rB_999.write();
    rB_50_90.write();
    rB_50_99.write();
    rB_50_999.write();
    rB_50_100.write();

//    rD_40.write();
    rD_50.write();
    rD_60.write();
    rD_70.write();
    rD_80.write();
    rD_90.write();
    rD_99.write();
    rD_999.write();
    rD_50_90.write();
    rD_50_99.write();
    rD_50_999.write();
    rD_50_100.write();

/*    rC_40pt.write();
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
    rD_90pt.write();*/

    rA_med.write();
//    rA_med_pt.write();
//    rB_med_pt.write();
    rC_med.write();

    rA_80.write();
    rA_90.write();
    rC_80.write();
    rC_90.write();

/*    rA_40.write();
    rA_40pt.write();
    rB_40pt.write();
    rC_40.write();*/


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
    h_massObs_q99->Write();
    h_massObs_q999->Write();
    h_massObs_q50_90->Write();
    h_massObs_q50_99->Write();
    h_massObs_q50_999->Write();
    h_massObs_q50_100->Write();
    
    h_massObs_q40_pt->Write();
    h_massObs_q50_pt->Write();
    h_massObs_q60_pt->Write();
    h_massObs_q70_pt->Write();
    h_massObs_q80_pt->Write();
    h_massObs_q90_pt->Write();

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
    h_pT_pTerrOverpT_0p8_eta_1p7->Write();
    h_pT_pTerrOverpT_1p7_eta_2p1->Write();

    h_pT_pTerrOverpT_rhad->Write();
    h_mass_pTerrOverpT->Write();
    h_ias_pTerrOverpT->Write();
    h_mass_fracLepton->Write();

    h_pT_PULLpT->Write();
    h_PULLpT->Write();
    
   h_N_1_miniIso->Write();
   h_N_1_miniIso_wMuon->Write();
   h_N_1_TkBased_absolute->Write();
   h_N_1_TkBased_relative->Write();
   h_Ias_miniIso->Write();
   h_Ias_miniIso_wMuon->Write();
   h_Ias_TkBased_absolute->Write();
   h_Ias_TkBased_relative->Write();
   h_Ias_preSIso->Write();
   h_Ih_miniIso->Write();
   h_Ih_miniIso_wMuon->Write();
   h_Ih_TkBased_absolute->Write();
   h_Ih_TkBased_relative->Write();
   h_Ih_preSIso->Write();
   h_pT_miniIso->Write();
   h_pT_miniIso_wMuon->Write();
   h_pT_TkBased_absolute->Write();
   h_pT_TkBased_relative->Write();
   h_pT_preSIso->Write();

   //MR_Bef.write();
   //MR_PostPreS_noFstrip.write();
   //MR_PostPreS_Eta2p1.write();
   //MR_all.write();
   //MR_ih_l_5.write();
   //MR_ih_ge_5.write();
   //MR_ias_l_0p6.write();
   //MR_ias_ge_0p6.write();
   //MR_HCAL_pic6.write();
   //MR_pterr0p25.write();
   //MR_cutQ99.write();
   //MR_isotk50.write();
   //MR_isotk15.write();
   //MR_mini_iso_03.write();

      outfile->Write();
      outfile->Close();

      /*std::cout << "nevent: " << nentries << " n_passPre/nevent: " << (float)n_passPre/(float)nentries << std::endl;
      std::cout << "n_aft: " << n_aft << " n_bef: " << n_bef << " n_aft/n_bef: " << (float)n_aft/(float)n_bef << std::endl;

      ofile << "nevent: " << nentries << " n_passPre/nevent: " << (float)n_passPre/(float)nentries << std::endl;
      ofile << "n_aft: " << n_aft << " n_bef: " << n_bef << " n_aft/n_bef: " << (float)n_aft/(float)n_bef << std::endl;
      
      ofile << "ntot: " << ntot << std::endl;
      ofile << "nA: " << nA << " " << 100*(float)nA/(float)ntot << " %" << std::endl;
      ofile << "nB: " << nB << " " << 100*(float)nB/(float)ntot << " %" << std::endl;
      ofile << "nC: " << nC << " " << 100*(float)nC/(float)ntot << " %" << std::endl;
      ofile << "nD: " << nD << " " << 100*(float)nD/(float)ntot << " %" << std::endl;*/

      /*ofile << "nB_boundedIas: " << nB_boundedIas << " " << 100*(float)nB_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nB_boundedIas/(float)nB << " % (/regionB)" << std::endl;
      ofile << "nD_boundedIas: " << nD_boundedIas << " " << 100*(float)nD_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedIas/(float)nD << " % (/regionD)" << std::endl;
      ofile << "nC_boundedPt: " << nC_boundedPt << " " << 100*(float)nC_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nC_boundedPt/(float)nC << " % (/regionC)" << std::endl;
      ofile << "nD_boundedPt: " << nD_boundedPt << " " << 100*(float)nD_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedPt/(float)nD << " % (/regionD)" << std::endl;
        */
      ofile.close();

}

