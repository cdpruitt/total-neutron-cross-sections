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

    //TFile *f = new TFile("tempTree.root","UPDATE");

    f->cd();

    if(f->Get("FTLR"))
    {
        f->Get("FTLR")->Delete();
    }

    if(f->Get("FTLR1D"))
    {
        f->Get("FTLR1D")->Delete();
    }

    TH2S *FTLR = new TH2S("FTLR","FTLR", 500, -500, 2500, 2500, -2500, 2500);
    TH1S *FTLR1D = new TH1S("FTLR1D","FTLR1D", 50, -2500, 2500);

    FTLR->GetXaxis()->SetTitle("Right fine time (ps)");
    FTLR->GetXaxis()->CenterTitle();
    FTLR->GetYaxis()->SetTitle("Left fine time (ps)");
    FTLR->GetYaxis()->CenterTitle();
    FTLR->SetMarkerStyle(7);

    FTLR1D->GetXaxis()->SetTitle("Right-Left difference (ps)");
    FTLR1D->GetXaxis()->CenterTitle();

    UInt_t leftTime;
    Int_t leftFine;
    UInt_t leftMacro;
    UInt_t leftRun;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        Long64_t kentry, lentry; // loop variables below

        // clear the left detector's info
        leftRun = -1;
        leftMacro = -1;
        leftTime = -1;
        leftFine = -1;

        b_chNo->GetEntry(jentry);

        if (chNo == 6) { // get left detector info
            b_runNo->GetEntry(jentry);
            b_macroNo->GetEntry(jentry);

            b_timetag->GetEntry(jentry);
            b_fineTime->GetEntry(jentry);

            leftRun = runNo;
            leftMacro = macroNo;
            leftTime = timetag;
            leftFine = fineTime;

            if (jentry<500)
            {
                kentry = 0;
            }

            else
            {
                kentry = jentry-500;
            }

            if (jentry+500>nentries)
            {
                lentry = nentries;
            }

            else
            {
                lentry = jentry+500;
            }

            for (kentry; kentry<lentry; kentry++) {
                b_chNo->GetEntry(kentry);

                if (chNo == 7) { // get right detector info
                    b_runNo->GetEntry(kentry);
                    b_macroNo->GetEntry(kentry);
                    b_timetag->GetEntry(kentry);

                    if ((leftTime == timetag || leftTime == timetag+2 || leftTime == timetag-2) && (leftMacro == macroNo && leftRun == runNo)) {
                        // run, macro, and time match; left/right detected same event
                        b_fineTime->GetEntry(kentry);
                        int diff = (timetag-leftTime)*1000; // adjust fine time by coarse time difference
                        //cout << leftFine*2000/1024-diff << endl;
                        FTLR->Fill(fineTime*2000/1024,(leftFine*2000/1024-diff)); // left det on abscissa
                        FTLR1D->Fill((fineTime-leftFine)*2000/1024); // right-left diff

                        //cout << "finetime is " << fineTime << " leftFine is " << leftFine << endl;
                        //cout << "entry " << jentry << " entry2 " << kentry << endl;
                    }
                    /*if ((leftTime == timetag+2 || leftTime == timetag-2) && (leftMacro == macroNo && leftRun == runNo)) {
                        // run, macro, and time match; left/right detected same event
                        b_fineTime->GetEntry(kentry);
                        FTLR->Fill(fineTime*2000/1024,leftFine*2000/1024); // left det on abscissa
                        FTLR1D->Fill((fineTime-leftFine)*2000/1024); // right-left diff

                        //cout << "finetime is " << fineTime << " leftFine is " << leftFine << endl;
                        //cout << "entry " << jentry << " entry2 " << kentry << endl;
                    }
                    */
                }
            }
        }

        //Long64_t ientry = LoadTree(jentry);
        //if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if(jentry%1000 == 0)
        {
            cout << jentry << endl;
        }
    }

    //f->Delete("tempTree");
    //f->Write();
    //FTLR->Write();
    //FTLR1D->Write();
}

int main()
{
    Loop();
    return 0;
}
