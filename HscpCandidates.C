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



void crossHistos(TH2F* res, TH1F* h1, TH1F* h2)
{
    scale(h1); scale(h2);
    for(int i=0;i<h1->GetNbinsX();i++)
    {
        for(int j=0;j<h2->GetNbinsX();j++)
        {   
            double prob = h1->GetBinContent(i)*h2->GetBinContent(j);
            res->SetBinContent(i,j,prob);
        }
    }
}

void crossHistosEtaBinning(TH2F* res, TH2F* eta_p, TH2F* ih_eta)
{
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    for(int i=1;i<eta->GetNbinsX();i++)
    {
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i);
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i);
        scale(p);
        for(int j=1;j<p->GetNbinsX()+1;j++)
        {
            for(int k=1;k<ih->GetNbinsX()+1;k++)
            {
                float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                if(prob<=0 || isnan((float)prob)) continue;
                //std::cout << p->GetBinContent(j) << std::endl;
                res->SetBinContent(i,j,prob);
            }
        }
        delete p;
        delete ih;
    }
    res->Sumw2();
}

void mapOfDifferences(TH2F* res, TH2F* h1, TH2F* h2)
{
    for(int i=0;i<h1->GetNbinsX();i++)
    {
        for(int j=0;j<h1->GetNbinsY();j++)
        {
            double diff = h1->GetBinContent(i,j)>0 ? (h2->GetBinContent(i,j))/h1->GetBinContent(i,j) : 0;
            res->SetBinContent(i,j,diff);
        }
    }
}

void predMass(TH1F* res, TH1F* h_p, TH1F* h_ih, float norm=0)
{
    //double norm = h_p->Integral();
//    std::cout << "norm: " << norm << std::endl;
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
    //std::cout << "entries: " << res->GetEntries() << " intgral: " << res->Integral() << std::endl;
}



float chi2testWith1(TH1F* h)
{
    float res=0.;
    for(int i=1;i<h->GetNbinsX();i++)
    {
        res += h->GetBinError(i)>0 ? pow((1-h->GetBinContent(i)),2)/h->GetBinError(i) : 0;
    }
    return res;
}

void ratioHist(TH1F* res, TH1F* h1, TH1F* h2)
{
    for(int i=0;i<h1->GetNbinsX();i++)
    {
        double ratio = h2->GetBinContent(i)>0 ? h1->GetBinContent(i)/h2->GetBinContent(i) : 0;
        res->SetBinContent(i,ratio);
    }
    res = (TH1F*) h1->Clone();
    res->Divide(h2);
    res->Sumw2();
    //chi2test(h1,h2);
    //chi2testWith1(res);
}



void ratioIntegral(TH1F* res, TH1F* h1, TH1F* h2)
{
    for(int i=0;i<h1->GetNbinsX()+1;i++)
    {   
        double ratio = h1->Integral(i,h1->GetNbinsX()+1)/h2->Integral(i,h2->GetNbinsX()+1);
        res->SetBinContent(i,ratio);
    }
}





void pseudoHisto(TH1F* h)
{
    TRandom3* RNG = new TRandom3();
    for(int i=0;i<h->GetNbinsX();i++)
    {
        h->SetBinContent(i,RNG->Poisson(h->GetBinContent(i)));
    }
}

TH2F* BetheBlochForMass(float mass)
{
    std::string strmass = to_string(mass);
    TH2F* tmp = new TH2F(strmass.c_str(),";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
    for(int i=1;i<tmp->GetNbinsX();i++)
    {
        float mom = tmp->GetXaxis()->GetBinCenter(i);
        float dedx = K*pow(mass/mom,2)+C;
        for(int j=1;j<tmp->GetNbinsY();j++)
        {
            if(tmp->GetYaxis()->GetBinLowEdge(j)<=dedx && tmp->GetYaxis()->GetBinLowEdge(j+1)>dedx)
            {
                tmp->SetBinContent(i,j,0.001);
            }   
        }
    }
    return tmp;
}

TH1F* rebinning(TH1F* h1, float min, std::vector<float>& vect)
{
    std::vector<float> v_val;
    int i=0;
    //v_val.push_back(0);
    while(i<h1->GetNbinsX()+1)
    {
        float cont=0;
        while(cont<min && i<h1->GetNbinsX()+1)
        {
            cont += h1->GetBinContent(i);
            i++;
        }
        v_val.push_back(h1->GetBinLowEdge(i));
    }



    //v_val.pop_back();


    TH1F* h2 = new TH1F("Rebinned","",v_val.size()-1,v_val.data());
    for(int j=0;j<h1->GetNbinsX()+1;j++)
    {
        for(int c=0;c<h1->GetBinContent(j);c++)
        {
            h2->Fill(h1->GetBinLowEdge(j));
        }
    }
    h2->Sumw2();
    vect = v_val;
    return h2;
}



TH1F* rebinningGraph(TH1F* h1, float min, std::vector<float>& vect1, std::vector<double>& vect2)
{
    std::vector<float> v_val;
    std::vector<float> v_cont;
    std::vector<float> v_err;
    std::vector<double> v_bin;
    int i=0;
    while(i<=h1->GetNbinsX()+1)
    {
        float cont=0;
        float div=0;
        float err=0;
        while(div<min && i<=h1->GetNbinsX()+1)
        {
            cont += h1->GetBinCenter(i)*h1->GetBinContent(i);
            div += h1->GetBinContent(i);
            err += pow(h1->GetBinError(i),2);
            i++;
        }
        v_val.push_back(div>0?cont/div:0);
        v_cont.push_back(div);
        v_err.push_back(sqrt(err));
        v_bin.push_back(h1->GetBinLowEdge(i));
    }

    TH1F* h2 = new TH1F("Rebinned","",v_bin.size()-1,v_bin.data());

    for(int j=0;j<h1->GetNbinsX()+1;j++)
    {
        for(int c=0;c<h1->GetBinContent(j);c++)
        {
            h2->Fill(h1->GetBinLowEdge(j));
        }
    }

    vect1 = v_val;
    vect2 = v_bin;

    TGraphErrors* gr = new TGraphErrors(v_val.size()-1,v_val.data(),v_cont.data(),0,v_err.data());
    return h2;
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

   TH1F* distribPt = new TH1F("pt",";pt [GeV]",500,0,2000);
   TH1F* distribPtPreselection = new TH1F("ptPreselection","Preselection;pt [GeV]",500,0,2000);
   TH1F* distribIas = new TH1F("ias",";I_{as}",100,0,1);
   TH1F* distribIasPreselection = new TH1F("iasPreselection","Preselection;I_{as}",100,0,1);

   TH2F* ih_p_fromD = new TH2F("ih_p_fromD",";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
   TH2F* ih_p_from2D = new TH2F("ih_p_from2D",";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
   TH2F* ih_p_fromBC = new TH2F("ih_p_fromBC",";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
   TH2F* ih_p_fromBCrew = new TH2F("ih_p_fromBCrew",";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
   TH2F* ih_p_fromBCrewbinning = new TH2F("ih_p_fromBCrewbinning",";p [GeV];I_{h} [MeV/cm]",200,0,2000,100,0,10);
   TH2F* diff_D_BC = new TH2F("diff_D_BC","",200,0,2000,100,0,10);
   TH2F* diff_D_2D = new TH2F("diff_D_2D","",200,0,2000,100,0,10);
   TH2F* diff_BC_2D = new TH2F("diff_BC_2D","",200,0,2000,100,0,10);

   
   TH1F* etaA = new TH1F("etaA","etaA",100,-3,3);
   TH1F* etaB = new TH1F("etaB","etaB",100,-3,3);
   TH1F* etaC = new TH1F("etaC","etaC",100,-3,3);
   TH1F* etaD = new TH1F("etaD","etaD",100,-3,3);

   TH2F* ih_ptA = new TH2F("ih_ptA","ih_ptA",100,0,1000,100,0,15);
   TH2F* ih_ptB = new TH2F("ih_ptB","ih_ptB",100,0,1000,100,0,15);
   TH2F* ih_ptC = new TH2F("ih_ptC","ih_ptC",100,0,1000,100,0,15);
   TH2F* ih_ptD = new TH2F("ih_ptD","ih_ptD",100,0,1000,100,0,15);

   TProfile* prof_ih_ptA = new TProfile("prof_ih_ptA","prof_ih_ptA",1000,0,1000);
   TProfile* prof_ih_ptB = new TProfile("prof_ih_ptB","prof_ih_ptB",1000,0,1000);
   TProfile* prof_ih_ptC = new TProfile("prof_ih_ptC","prof_ih_ptC",1000,0,1000);
   TProfile* prof_ih_ptD = new TProfile("prof_ih_ptD","prof_ih_ptD",1000,0,1000);

   TH1F* iso_abs = new TH1F("iso_abs","iso_abs",1000,0,100);
   TH1F* iso_rel = new TH1F("iso_rel","iso_rel",1000,0,100);

   TH1F* h_eta = new TH1F("eta","eta",100,-3,3);
   TH1F* eta_moderateIso = new TH1F("eta_moderateIso","eta_moderateIso",100,-3,3);
   TH1F* eta_moderateIsoD = new TH1F("eta_moderateIsoD","eta_moderateIsoD",100,-3,3);

   int CutIndex = 19;

   TH1F* regD_ih_proj = (TH1F*) regD_ih->ProjectionY("ih_proj",CutIndex+1,CutIndex+1);
   TH1F* regD_p_proj = (TH1F*) regD_p->ProjectionY("p_proj",CutIndex+1,CutIndex+1);
   TH1F* regD_mass_proj = (TH1F*) regD_mass->ProjectionY("regD_mass_proj",CutIndex+1,CutIndex+1);

   TH2F* regD_ih_p = new TH2F("regD_ih_p_fromTH1",";p [GeV];I_{h} [MeV/cm]",regD_p_proj->GetNbinsX(),0,regD_p_proj->GetBinLowEdge(regD_p_proj->GetNbinsX()+1),regD_ih_proj->GetNbinsX(),0,regD_ih_proj->GetBinLowEdge(regD_ih_proj->GetNbinsX()+1));


   std::vector<float> v_bins;
   std::vector<double> p_bins;
   std::vector<float> vector_bin_gr;

    regD_mass_proj = rebinning(regD_mass_proj,100,v_bins);
    //regD_p_proj = rebinning(regD_p_proj,50,p_bins);

    regD_p_proj = rebinningGraph(regD_p_proj,100,vector_bin_gr,p_bins);
   
   TH1F* pred_mass = (TH1F*)regD_mass_proj->Clone("pred_mass"); pred_mass->Reset();
   TH1F* pred_massBC = (TH1F*)regD_mass_proj->Clone("pred_massBC"); pred_massBC->Reset();
   TH1F* pred_massBCrew = (TH1F*)regD_mass_proj->Clone("pred_massBCrew"); pred_massBCrew->Reset();
   TH1F* pred_massBCrewbinning = (TH1F*)regD_mass_proj->Clone("pred_massBCrewbinning"); pred_massBCrewbinning->Reset();
   TH1F* pred_massDCrewbinning = (TH1F*)regD_mass_proj->Clone("pred_massDCrewbinning"); pred_massDCrewbinning->Reset();
   TH1F* pred_massBDbinning = (TH1F*)regD_mass_proj->Clone("pred_massBDbinning"); pred_massBDbinning->Reset();
   TH1F* pred_massBCrewbinningCutI = (TH1F*)regD_mass_proj->Clone("pred_massBCrewbinning"); pred_massBCrewbinning->Reset();
   TH1F* pred_massTemplatesFromD = (TH1F*)regD_mass_proj->Clone("pred_massTemplatesFromD"); pred_massTemplatesFromD->Reset();

   TH1F* ratio_mass = (TH1F*)regD_mass_proj->Clone("ratio_mass"); ratio_mass->Reset();

   TH1F* massTree = (TH1F*)regD_mass_proj->Clone("massTree"); massTree->Reset();
   TH1F* massTreeAll = (TH1F*)regD_mass_proj->Clone("massTreeAll"); massTreeAll->Reset();

   /*regD_ih_proj->Scale(1./regD_ih_proj->Integral(0,regD_ih_proj->GetNbinsX()+1));
   regD_p_proj->Scale(1./regD_p_proj->Integral(0,regD_p_proj->GetNbinsX()+1));

   double norma = regD_mass_proj->Integral();

   for(int i=0;i<regD_p_proj->GetNbinsX();i++)
   {
       double p = regD_p_proj->GetBinCenter(i);

       for(int j=0;j<regD_ih_proj->GetNbinsX();j++)
       {
           double ih = regD_ih_proj->GetBinCenter(j);
           double proba = regD_p_proj->GetBinContent(i) * regD_ih_proj->GetBinContent(j);
           double mass = GetMass(p,ih,K,C);
           pred_mass->Fill(mass,proba*norma);
           regD_ih_p->SetBinContent(i,j,proba);
           //std::cout << "mass: " << mass << " proba: " << proba << " norma: " << norma << std::endl;
   
       }
   }
   //pred_mass->Scale(norma);
*/

   //predMass(pred_mass,regD_p_proj,regD_ih_proj);
   //ratioIntegral(ratio_mass,regD_mass_proj,pred_mass);

   Long64_t nentries = fChain->GetEntriesFast();
   
   std::cout << "nentries: " << nentries << std::endl;

   //nentries=1000000;

   int ntot=0, nA=0, nB=0, nC=0, nD=0;

   int nbins = v_bins.size()-1;
   float* xbins = v_bins.data();

   /*region rAll("_all",nbins,xbins,p_bins,vector_bin_gr);
   region rA("_regionA",nbins,xbins,p_bins,vector_bin_gr);
   region rB("_regionB",nbins,xbins,p_bins,vector_bin_gr);
   region rC("_regionC",nbins,xbins,p_bins,vector_bin_gr);
   region rD("_regionD",nbins,xbins,p_bins,vector_bin_gr);
   region rCutIh("_cutIh",nbins,xbins,p_bins,vector_bin_gr);
   region rEta1("_eta1",nbins,xbins,p_bins,vector_bin_gr);
   region rEta2("_eta2",nbins,xbins,p_bins,vector_bin_gr);
   region rEta10("_eta10",nbins,xbins,p_bins,vector_bin_gr);
   region rDcutIh("_regD_CutIh4p5",nbins,xbins,p_bins,vector_bin_gr);*/

   region rAll("_all",etabins_,ihbins_,pbins_,massbins_);
   region rA("_regionA",etabins_,ihbins_,pbins_,massbins_);
   region rB("_regionB",etabins_,ihbins_,pbins_,massbins_);
   region rC("_regionC",etabins_,ihbins_,pbins_,massbins_);
   region rD("_regionD",etabins_,ihbins_,pbins_,massbins_);


   
   //cutindex=19 pt=60 ias=0.025

   //double pt_cut=60.;
   //double ias_cut=0.025;
   
   
   double pt_cut=ptcut_;
   double ias_cut=iascut_;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
        
      if(jentry%1000000==0) std::cout << "entry: " << jentry << std::endl;

      if(Hscp<1) continue;
      if(Pt->size()<1) continue;

      if(invMET_ && RecoPFMET>100) continue;

      int i=0;
      for(int i=0; i<Mass->size(); i++)
      {
          if(invIso_ && !passPreselection_noIsolation_noIh->at(i)) continue;
          if(!invIso_ && !passPreselection->at(i)) continue;

          float pt = Pt->at(i);
          if(pt<50) continue;
          float ih = Ih->at(i);
          float ias = Ias->at(i);
          float Eta = eta->at(i); 
          float iso = iso_TK->at(i);
          float iso_r = iso/pt;
          float nhits = NOM->at(i);
          float p = pt*cosh(Eta);
          float massT = Mass->at(i);

          float isotk = iso_TK->at(i);
          float isocalo = iso_ECAL->at(i)+iso_HCAL->at(i);


          if(invIso_ && (isotk<50 || isocalo<0.3)) continue;

          //if(abs(Eta)>0.4) continue;
          
          if(ih<ihcut_) continue;

          if(p>pcut_ && pcut_>0) continue;

          distribPt->Fill(pt);
          distribIas->Fill(ias);


          /*massTreeAll->Fill(massT);

          distribPtPreselection->Fill(pt);
          distribIasPreselection->Fill(ias);

          //if(pt<51) std::cout << pt << std::endl;
          
          iso_abs->Fill(iso);
          iso_rel->Fill(iso_r);

          h_eta->Fill(Eta);
          if(iso>20) eta_moderateIso->Fill(Eta);


          if(abs(Eta)>=0.0 && abs(Eta)<0.1) rEta1.fill(Eta,nhits,p,pt,ih,ias,massT);
          if(abs(Eta)>=0.1 && abs(Eta)<0.2) rEta2.fill(Eta,nhits,p,pt,ih,ias,massT);
          if(abs(Eta)>=1.0 && abs(Eta)<1.1) rEta10.fill(Eta,nhits,p,pt,ih,ias,massT);

          if(ih<3.4) rCutIh.fill(Eta,nhits,p,pt,ih,ias,massT);*/

          rAll.fill(Eta,nhits,p,pt,ih,ias,massT);
          ntot++;

          if(pt<=pt_cut)
          {
              if(ias<=ias_cut) //regionA
              {
                  /*etaA->Fill(Eta);
                  ih_ptA->Fill(pt,ih);
                  prof_ih_ptA->Fill(pt,ih);*/
                  rA.fill(Eta,nhits,p,pt,ih,ias,massT); 
                  nA++;

              }
              else //regionB
              {
                  /*etaB->Fill(Eta);
                  ih_ptB->Fill(pt,ih);
                  prof_ih_ptB->Fill(pt,ih);*/
                  rB.fill(Eta,nhits,p,pt,ih,ias,massT); 
                  nB++;

              }
          }
          else
          {
              if(ias<=ias_cut) //regionC
              {
                  /*etaC->Fill(Eta);
                  ih_ptC->Fill(pt,ih);
                  prof_ih_ptC->Fill(pt,ih);*/
                  rC.fill(Eta,nhits,p,pt,ih,ias,massT); 
                  nC++;

              }
              else //regionD
              {
                  /*massTree->Fill(massT);
                  ih_p_fromD->Fill(p,ih);
                  etaD->Fill(Eta);
                  ih_ptD->Fill(pt,ih);
                  prof_ih_ptD->Fill(pt,ih);
                  if(iso>20) eta_moderateIsoD->Fill(Eta);
                  //if(p<600 && ih<4.5 && ih>3) */
                  
                  rD.fill(Eta,nhits,p,pt,ih,ias,massT); 
                      //rDcutIh.fill(Eta,nhits,p,pt,ih,ias,massT);
                  
                  nD++;
              }
          }

      }

   
   }

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
   rCutIh.fillStdDev();
   rEta1.fillStdDev();
   rEta2.fillStdDev();
   rEta10.fillStdDev();
   rDcutIh.fillStdDev();
   
   rAll.fillQuantile();
   rA.fillQuantile();
   rB.fillQuantile();
   rC.fillQuantile();
   rD.fillQuantile();
   rCutIh.fillQuantile();
   rEta1.fillQuantile();
   rEta2.fillQuantile();
   rEta10.fillQuantile();
   rDcutIh.fillQuantile();

   rAll.fillMassFrom1DTemplates();
   rA.fillMassFrom1DTemplates();
   rB.fillMassFrom1DTemplates();
   rC.fillMassFrom1DTemplates();
   rD.fillMassFrom1DTemplates();
   rCutIh.fillMassFrom1DTemplates();
   rDcutIh.fillMassFrom1DTemplates();

   rAll.fillMassFrom1DTemplatesEtaBinning();
   rA.fillMassFrom1DTemplatesEtaBinning();
   rB.fillMassFrom1DTemplatesEtaBinning();
   rC.fillMassFrom1DTemplatesEtaBinning();
   rD.fillMassFrom1DTemplatesEtaBinning();
   rDcutIh.fillMassFrom1DTemplatesEtaBinning();*/

    double normalisation = nB*nC/nA;

    std::cout << "normalisation: " << normalisation << " nD: " << nD << std::endl;

   /*predMass(pred_mass,(TH1F*)rD.ih_p->ProjectionX(),(TH1F*)rD.ih_p->ProjectionY());
   predMass(pred_massBC,(TH1F*)rC.ih_p->ProjectionX(),(TH1F*)rB.ih_p->ProjectionY(),normalisation);
   //predMass(pred_massTemplatesFromD,(TH1F*)regD_p_proj,(TH1F*)rD.ih_p->ProjectionY(),normalisation);

   etaReweighingP(rC.eta_p,rB.eta_p);
   TH1F* prew = (TH1F*)rC.eta_p->ProjectionX(); prew->Sumw2();
   
   //pseudoHisto(prew);

   predMass(pred_massBCrew,prew,(TH1F*)rB.ih_p->ProjectionY(),normalisation);
   predMassEtaBinning(pred_massBCrewbinning,rC.eta_p,rB.ih_eta,normalisation); 
   predMassEtaBinning(pred_massDCrewbinning,rC.eta_p,rD.ih_eta,normalisation);
   predMassEtaBinning(pred_massBDbinning,rD.eta_p,rB.ih_eta,normalisation);

   crossHistos(ih_p_fromBCrew,(TH1F*)rC.eta_p->ProjectionX(),(TH1F*)rB.ih_eta->ProjectionY());
   crossHistosEtaBinning(ih_p_fromBCrewbinning,rC.eta_p,rB.ih_eta);*/
   
   //TFile* outfile = new TFile("outfile_nocut.root","RECREATE");
   ofstream ofile((outfilename_+"_normalisations.txt").c_str());

   TFile* outfile = new TFile((outfilename_+".root").c_str(),"RECREATE");

   
   /*
//      plotting(regD_mass_proj,rD.mass,true,"obs")->Write();
      plotting(rD.mass,pred_mass,false,"predMasstree")->Write();
      plotting(regD_mass_proj,pred_mass,false,"predMassD")->Write();
      //plotting(regD_mass_proj,pred_massBC,false,"predMassBC")->Write();
      plotting(rD.mass,pred_massBC,false,"predMassBCtree")->Write();
      //plotting(regD_mass_proj,pred_massBCrew,true,"predMassBCrew")->Write();
      plotting(rD.mass,pred_massBCrew,true,"predMassBCrewtree")->Write();
      plotting(pred_mass,pred_massBCrew,false,"predMass_predMassBCrew")->Write();
      plotting(pred_massBC,pred_massBCrew,false,"predMassBC_predMassBCrew")->Write();
      plotting(regD_mass_proj,pred_mass,true,"RatioSimple")->Write();
      plotting((TH1F*)rB.ih_p->ProjectionY(),(TH1F*)rD.ih_p->ProjectionY(),true,"Ih_BD")->Write();
      plotting((TH1F*)rC.ih_p->ProjectionX(),(TH1F*)rD.ih_p->ProjectionX(),true,"P_CD")->Write();
      plotting(prew,(TH1F*)rD.ih_p->ProjectionX(),true,"P_CD_rew")->Write();
*/

   //rebinning(massTree,100)->Write();

   /*scale(massTree); scale(massTreeAll);
   scale(pred_massBC);
   scale(pred_massBCrew);
   scale(pred_massBCrewbinning);
      massTree->Write();

      regD_mass_proj->Write();
      regD_mass->Write();

      pred_mass->Write();
      pred_massBC->Write();
      pred_massBCrew->Write();
      pred_massBCrewbinning->Write();
      pred_massDCrewbinning->Write();
      pred_massBDbinning->Write();

      ratioIntegral(pred_massBCrew,massTree)->Write();
      ratioIntegral(pred_massBCrewbinning,massTree)->Write();

      plotting(massTree,pred_massBCrew,false,"BCrew")->Write();
      plotting(massTree,pred_massBCrew,true,"BCrew_simple")->Write();
      plotting(massTree,pred_massBCrewbinning,false,"BCrewbinning")->Write();
      plotting(massTree,pred_massBCrewbinning,true,"BCrewbinning_simple")->Write();
      plotting(massTree,pred_massDCrewbinning,false,"DCrewbinning")->Write();
      plotting(massTree,pred_massDCrewbinning,true,"DCrewbinning_simple")->Write();
      plotting(massTree,pred_massBDbinning,false,"BDbinning")->Write();
      plotting(massTree,pred_massBDbinning,true,"BDbinning_simple")->Write();
      plotting(massTree,pred_massTemplatesFromD,true,"templatesFromD")->Write();

      plotting((TH1F*)rB.eta_p->ProjectionY(),(TH1F*)rD.eta_p->ProjectionY(),true,"eta_b_d")->Write();
      plotting((TH1F*)rC.eta_p->ProjectionX(),(TH1F*)rD.eta_p->ProjectionX(),true,"p_c_d")->Write();
      plotting((TH1F*)rB.ih_eta->ProjectionY(),(TH1F*)rD.ih_eta->ProjectionY(),true,"ih_b_d")->Write();



      ih_p_fromBCrew->Write();
      ih_p_fromBCrewbinning->Write();
      
      BetheBlochForMass(600)->Write();
      BetheBlochForMass(800)->Write();
      BetheBlochForMass(1000)->Write();
      BetheBlochForMass(1200)->Write();
      BetheBlochForMass(1400)->Write();

      rD.ih_p->Write();
      ih_p_fromD->Write();
      ih_p_from2D->Write();
      ih_p_fromBC->Write();
      diff_D_BC->Write();
      diff_D_2D->Write();
      diff_BC_2D->Write();


      ratio_mass->Write();

      regD_p_proj->Write();
      regD_ih_proj->Write();

      regD_ih_p->Write();

      distribPt->Write();
      distribPtPreselection->Write();

      distribIas->Write();
      distribIasPreselection->Write();
      
      iso_abs->Write();
      iso_rel->Write();
      etaA->Write();
      etaB->Write();
      etaC->Write();
      etaD->Write();
      ih_ptA->Write();
      ih_ptB->Write();
      ih_ptC->Write();
      ih_ptD->Write();
      prof_ih_ptA->Write();
      prof_ih_ptB->Write();
      prof_ih_ptC->Write();
      prof_ih_ptD->Write();

      h_eta->Scale(1./h_eta->Integral());
      eta_moderateIso->Scale(1./eta_moderateIso->Integral());
      eta_moderateIsoD->Scale(1./eta_moderateIsoD->Integral());
      TCanvas* c1 = new TCanvas("c1","c1");
      etaA->Scale(1./etaA->Integral());
      etaB->Scale(1./etaB->Integral());
      etaC->Scale(1./etaC->Integral());
      etaD->Scale(1./etaD->Integral());
      etaA->SetLineColor(1);
      etaB->SetLineColor(2);
      etaC->SetLineColor(3);
      etaD->SetLineColor(4);
      etaA->Draw();
      etaB->Draw("same");
      etaC->Draw("same");
      etaD->Draw("same");

      h_eta->Write();
      eta_moderateIso->Write();
      eta_moderateIsoD->Write();
      c1->Write();
*/

      rAll.write();
      rA.write();
      rB.write();
      rC.write();
      rD.write();
      //rDcutIh.write();
      //rCutIh.write();
      //rEta1.write();
      //rEta2.write();
      //rEta10.write();



      outfile->Write();
      outfile->Close();

      std::cout << "ntot: " << ntot << std::endl;
      std::cout << "nA: " << nA << " " << 100*(float)nA/(float)ntot << " %" << std::endl;
      std::cout << "nB: " << nB << " " << 100*(float)nB/(float)ntot << " %" << std::endl;
      std::cout << "nC: " << nC << " " << 100*(float)nC/(float)ntot << " %" << std::endl;
      std::cout << "nD: " << nD << " " << 100*(float)nD/(float)ntot << " %" << std::endl;


      ofile << "ntot: " << ntot << std::endl;
      ofile << "nA: " << nA << " " << 100*(float)nA/(float)ntot << " %" << std::endl;
      ofile << "nB: " << nB << " " << 100*(float)nB/(float)ntot << " %" << std::endl;
      ofile << "nC: " << nC << " " << 100*(float)nC/(float)ntot << " %" << std::endl;
      ofile << "nD: " << nD << " " << 100*(float)nD/(float)ntot << " %" << std::endl;


      ofile.close();

}
