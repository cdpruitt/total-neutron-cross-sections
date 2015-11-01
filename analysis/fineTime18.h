//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct 31 14:55:32 2015 by ROOT version 5.34/32
// from TTree tree/
// found on file: run18-0000.root
//////////////////////////////////////////////////////////

#ifndef fineTime18_h
#define fineTime18_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class fineTime18 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          macroNo;
   UInt_t          evtNo;
   UInt_t          chNo;
   UInt_t          evtType;
   UInt_t          extTime;
   UInt_t          timetag;
   UInt_t          fineTime;
   UInt_t          sgQ;
   UInt_t          lgQ;

   // List of branches
   TBranch        *b_macroNo;   //!
   TBranch        *b_evtNo;   //!
   TBranch        *b_chNo;   //!
   TBranch        *b_evtType;   //!
   TBranch        *b_extTime;   //!
   TBranch        *b_timetag;   //!
   TBranch        *b_fineTime;   //!
   TBranch        *b_sgQ;   //!
   TBranch        *b_lgQ;   //!

   fineTime18(TTree *tree=0);
   virtual ~fineTime18();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fineTime18_cxx
fineTime18::fineTime18(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run23-0000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("run23-0000.root","UPDATE");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

fineTime18::~fineTime18()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fineTime18::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fineTime18::LoadTree(Long64_t entry)
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

void fineTime18::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("macroNo", &macroNo, &b_macroNo);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evtNo);
   fChain->SetBranchAddress("chNo", &chNo, &b_chNo);
   fChain->SetBranchAddress("evtType", &evtType, &b_evtType);
   fChain->SetBranchAddress("extTime", &extTime, &b_extTime);
   fChain->SetBranchAddress("timetag", &timetag, &b_timetag);
   fChain->SetBranchAddress("fineTime", &fineTime, &b_fineTime);
   fChain->SetBranchAddress("sgQ", &sgQ, &b_sgQ);
   fChain->SetBranchAddress("lgQ", &lgQ, &b_lgQ);
   Notify();
}

Bool_t fineTime18::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fineTime18::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fineTime18::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fineTime18_cxx
