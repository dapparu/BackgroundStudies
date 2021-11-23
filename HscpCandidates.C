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
   region rB_boundedIas("_regionB_boundedIas",etabins_,ihbins_,pbins_,massbins_);
   region rD_boundedIas("_regionD_boundedIas",etabins_,ihbins_,pbins_,massbins_);
   region rC_boundedPt("_regionC_boundedPt",etabins_,ihbins_,pbins_,massbins_);
   region rD_boundedPt("_regionD_boundedPt",etabins_,ihbins_,pbins_,massbins_);

   
   //cutindex=19 pt=60 ias=0.025
   
   double pt_cut=ptcut_;
   double ias_cut=iascut_;
   
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
          float tof = TOF->at(i);
          if(!(etacutinf_<Eta && Eta<etacutsup_)) continue;
          if(isocalo/p>0.3)continue;
          if((invIso_ && (isotk<50 || isotk>100))) continue;
          if(ih<ihcut_) continue;
          if(p>pcut_ && pcut_>0) continue;


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
                      rB_boundedIas.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
                      nB_boundedIas++;
                  }
              }
          }
          else
          {
              if(ias<=ias_cut) //regionC
              {
                  rC.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nC++;

                  if(pt<70)
                  {
                      rC_boundedPt.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                      nC_boundedPt++;
                  }

              }
              else //regionD
              {
                  rD.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                  nD++;
                  
                  if(ias<0.1)
                  {
                      rD_boundedIas.fill(Eta,nhits,p,pt,ih,ias,massT,tof);
                      nD_boundedIas++;
                  } 
                  
                  if(pt<70)
                  {
                      rD_boundedPt.fill(Eta,nhits,p,pt,ih,ias,massT,tof); 
                      nD_boundedPt++;
                  }
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
      rB_boundedIas.write();
    std::cout << " region B_boundedIas saved " << std::endl;
      rD_boundedIas.write();
    std::cout << " region D_boundedIas saved " << std::endl;
      rC_boundedPt.write();
    std::cout << " region C_boundedPt saved " << std::endl;
      rD_boundedPt.write();
    std::cout << " region D_boundedPt saved " << std::endl;



      outfile->Write();
      outfile->Close();

      ofile << "ntot: " << ntot << std::endl;
      ofile << "nA: " << nA << " " << 100*(float)nA/(float)ntot << " %" << std::endl;
      ofile << "nB: " << nB << " " << 100*(float)nB/(float)ntot << " %" << std::endl;
      ofile << "nC: " << nC << " " << 100*(float)nC/(float)ntot << " %" << std::endl;
      ofile << "nD: " << nD << " " << 100*(float)nD/(float)ntot << " %" << std::endl;
      ofile << "nB_boundedIas: " << nB_boundedIas << " " << 100*(float)nB_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nB_boundedIas/(float)nB << " % (/regionB)" << std::endl;
      ofile << "nD_boundedIas: " << nD_boundedIas << " " << 100*(float)nD_boundedIas/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedIas/(float)nD << " % (/regionD)" << std::endl;
      ofile << "nC_boundedPt: " << nC_boundedPt << " " << 100*(float)nC_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nC_boundedPt/(float)nC << " % (/regionC)" << std::endl;
      ofile << "nD_boundedPt: " << nD_boundedPt << " " << 100*(float)nD_boundedPt/(float)ntot << " % (/total)" << " " << 100*(float)nD_boundedPt/(float)nD << " % (/regionD)" << std::endl;

      ofile.close();

}
