#define finetime_cxx
#include "finetime.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void finetime::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L finetime.C
//      Root > finetime t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

    //TCanvas *c1 = f.Get("c1");
    //if (f->GetObject("FTLR")) {
    //    f->Delete("FTLR")
    //}

    TFile* f = new TFile("tempTree.root","UPDATE");

    f->ls();

    if (f->Get("FTLR")) {
        TH2S *FTLR = f->Get("FTLR");
        FTLR->Reset();
    }

    else {
        TH2S *FTLR = new TH2S("FTLR","FTLR", 1023, 0, 1023, 1023, 0, 1023);
    }

    TH1S *FTLR1D = new TH1S("FTLR1D","FTLR1D", 1023, 0, 1023);
    FTLR->SetMarkerStyle(20);

    UInt_t leftTime;
    UInt_t leftFine;
    UInt_t leftMacro;
    UInt_t leftRun;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        b_chNo->GetEntry(jentry);

        if (chNo == 6) { // left detector
            b_runNo->GetEntry(jentry);
            b_macroNo->GetEntry(jentry);

            b_timetag->GetEntry(jentry);
            b_fineTime->GetEntry(jentry);

            leftRun = runNo;
            leftMacro = macroNo;
            leftTime = timetag;
            leftFine = fineTime;

            for (Long64_t kentry=0; kentry<nentries; kentry++) {
                b_chNo->GetEntry(kentry);

                if (chNo == 7) { // right detector
                    b_runNo->GetEntry(jentry);
                    b_macroNo->GetEntry(jentry);
                    b_timetag->GetEntry(kentry);

                    if (leftTime == timetag && (leftMacro == macroNo && leftRun == runNo)) {
                        // run, macro, and time match; same event
                        b_fineTime->GetEntry(kentry);
                        FTLR->Fill(fineTime,leftFine); // left det on abscissa
                        FTLR1D->Fill(fineTime-leftFine); // right-left diff
                        //cout << "finetime is " << fineTime << " leftTime is " << leftTime << endl;
                    }
                }
            }
        }

        //Long64_t ientry = LoadTree(jentry);
        //if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }

    FTLR->Write();
    //FTLR1D->Write();

}
