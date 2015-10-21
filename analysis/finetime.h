//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 21 11:31:20 2015 by ROOT version 5.34/32
// from TTree tempTree/
// found on file: tempTree.root
//////////////////////////////////////////////////////////

#ifndef finetime_h
#define finetime_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class finetime {
public :
   TFile          *f;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UShort_t        runNo;
   UShort_t        macroNo;
   UShort_t        evtNo;
   UShort_t        chNo;
   UInt_t          evtType;
   UInt_t          timetag;
   UShort_t        fineTime;
   UShort_t        sgQ;
   UShort_t        lgQ;
   UShort_t        waveform;

   // List of branches
   TBranch        *b_runNo;   //!
   TBranch        *b_macroNo;   //!
   TBranch        *b_evtNo;   //!
   TBranch        *b_chNo;   //!
   TBranch        *b_evtType;   //!
   TBranch        *b_timetag;   //!
   TBranch        *b_fineTime;   //!
   TBranch        *b_sgQ;   //!
   TBranch        *b_lgQ;   //!
   TBranch        *b_waveform;   //!

   finetime(TTree *tree=0);
   virtual ~finetime();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef finetime_cxx
finetime::finetime(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      f = (TFile*)gROOT->GetListOfFiles()->FindObject("tempTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tempTree.root","UPDATE");
      }
      f->GetObject("tempTree",tree);

   }
   Init(tree);
}

finetime::~finetime()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t finetime::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t finetime::LoadTree(Long64_t entry)
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

void finetime::Init(TTree *tree)
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

   fChain->SetBranchAddress("runNo", &runNo, &b_runNo);
   fChain->SetBranchAddress("macroNo", &macroNo, &b_macroNo);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evtNo);
   fChain->SetBranchAddress("chNo", &chNo, &b_chNo);
   fChain->SetBranchAddress("evtType", &evtType, &b_evtType);
   fChain->SetBranchAddress("timetag", &timetag, &b_timetag);
   fChain->SetBranchAddress("fineTime", &fineTime, &b_fineTime);
   fChain->SetBranchAddress("sgQ", &sgQ, &b_sgQ);
   fChain->SetBranchAddress("lgQ", &lgQ, &b_lgQ);
   fChain->SetBranchAddress("waveform", &waveform, &b_waveform);
   Notify();
}

Bool_t finetime::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void finetime::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t finetime::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef finetime_cxx
