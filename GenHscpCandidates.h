//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 21 11:14:29 2022 by ROOT version 6.14/09
// from TTree GenHscpCandidates/GenHscpCandidates
// found on file: histoGluino_M-1600_TuneCP5_20UL17_v16.root
//////////////////////////////////////////////////////////

#ifndef GenHscpCandidates_h
#define GenHscpCandidates_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class GenHscpCandidates {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Run;
   ULong64_t       Event;
   UInt_t          Lumi;
   Float_t         Weight;
   Float_t         GeneratorWeight;
   vector<float>   *GenId;
   vector<float>   *GenCharge;
   vector<float>   *GenMass;
   vector<float>   *GenPt;
   vector<float>   *GenEta;
   vector<float>   *GenPhi;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_GeneratorWeight;   //!
   TBranch        *b_GenId;   //!
   TBranch        *b_GenCharge;   //!
   TBranch        *b_GenMass;   //!
   TBranch        *b_GenPt;   //!
   TBranch        *b_GenEta;   //!
   TBranch        *b_GenPhi;   //!

   GenHscpCandidates(TTree *tree=0);
   virtual ~GenHscpCandidates();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GenHscpCandidates_cxx
GenHscpCandidates::GenHscpCandidates(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   std::cout << "here " << std::endl;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histoGluino_M-1600_TuneCP5_20UL17_v16.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("histoGluino_M-1600_TuneCP5_20UL17_v16.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("histoGluino_M-1600_TuneCP5_20UL17_v16.root:/analyzer/BaseName");
      dir->GetObject("GenHscpCandidates",tree);

   }
   std::cout << "here " << std::endl;
   Init(tree);
}

GenHscpCandidates::~GenHscpCandidates()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GenHscpCandidates::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   std::cout << "here " << std::endl;
   return fChain->GetEntry(entry);
}
Long64_t GenHscpCandidates::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GenHscpCandidates::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

    std::cout << "init genTree" << std::endl;
   // Set object pointer
   GenId = 0;
   GenCharge = 0;
   GenMass = 0;
   GenPt = 0;
   GenEta = 0;
   GenPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
   fChain->SetBranchAddress("GenId", &GenId, &b_GenId);
   fChain->SetBranchAddress("GenCharge", &GenCharge, &b_GenCharge);
   fChain->SetBranchAddress("GenMass", &GenMass, &b_GenMass);
   fChain->SetBranchAddress("GenPt", &GenPt, &b_GenPt);
   fChain->SetBranchAddress("GenEta", &GenEta, &b_GenEta);
   fChain->SetBranchAddress("GenPhi", &GenPhi, &b_GenPhi);
   Notify();
}

Bool_t GenHscpCandidates::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GenHscpCandidates::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GenHscpCandidates::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GenHscpCandidates_cxx
