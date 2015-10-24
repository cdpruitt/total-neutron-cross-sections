#define analyze_cxx
#include "analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <vector>
#include <string>

#pragma link C++ class vector<short>+;
#pragma link C++ class vector<TH1S*>+;
//#pragma link C++ class vector<vector<TH1S*>>+;
#pragma link C++ class vector<TDirectory*>+;

void analyze::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L analyze.C
//      Root > analyze t
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

    TFile *histograms = new TFile("histograms.root","RECREATE");
    histograms->cd();
    
    UInt_t temptime;
    UInt_t tempfine;
    UInt_t temprunNo = -1;
    vector<TDirectoryFile*> direct;

    TDirectory *targetChanger;
    TDirectory *monitor;
    TDirectory *detectorT;
    TDirectory *detectorL;
    TDirectory *detectorR;

    vector<TH1S*> listWaveforms;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<1000;jentry++) {

        if(jentry%1000 == 0)
        {
            cout << jentry << endl;
        }

        GetEntry(jentry);
        if(temprunNo != runNo) // new run; new directory; new histos
        {
            stringstream temp;
            temp << "run " << runNo;

            direct.push_back(new TDirectoryFile(temp.str().c_str(),temp.str().c_str()));
            direct.back()->cd();

            string subDirs[8] = {"targetChanger","","monitor","","detT","","detL","detR"}

            vector<TH1S*> histos; // holds all histograms for the run
            // split into sub-vectors on a per-channel basis

            for(int i=0; i<8; i++)
            {
                //histos.push_back(new vector<TH1S>); // create sub-vector for this channel
                if (subDirs[i].compare("") != 0)
                {
                    direct.back()->mkdir(subDirs[i].c_str(),subDirs[i].c_str());
                    direct.back()->GetDirectory(subDirs[i].c_str())->cd();

                    // instantiate histograms
                    histos.push_back(new TH1S("outMacro","outMacro",100000,0,10000000));
                    histos.push_back(new TH1S("outEvt","outEvt",500,0,1000));
                    histos.push_back(new TH1S("outTime","outTime",100000,0,100000000));
                    histos.push_back(new TH1S("outSGQ","outSGQ",1024,0,70000));
                    histos.push_back(new TH1S("outLGQ","outLGQ",1024,0,70000));
                    histos.push_back(new TH1S("outFT","outFT",1023,0,1023));

                    direct.back()->GetDirectory(subDirs[i].c_str())->mkdir("waveforms","waveforms");
                    listWaveforms.clear();
                }
            }
        }

        temprunNo = runNo; // update the runNo counter to prep for new runNo

        direct.back()->GetDirectory(subDirs[chNo].c_str())->cd();

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outMacro"));
        outMacro->Fill(macroNo);

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outEvt"));
        outEvt->Fill(evtNo);

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outTime"));
        outTime->Fill(timetag);

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outSGQ"));
        outSGQ->Fill(sgQ);

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outLGQ"));
        outLGQ->Fill(lgQ);

        TH1S* outMacro = (TH1S*)(gDirectory->FindObject("outFT"));
        outFT->Fill(fineTime);

        stringstream temp;
        temp << "event " << evtNo;

        waveforms->cd();
        listWaveforms.push_back(new TH1S(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2));

        for(int i=0;i<waveform->size();i++)
        {
            listWaveforms.back()->SetBinContent(i,waveform->at(i));
        }

        //Long64_t ientry = LoadTree(jentry);
        //if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }

    outMacro->Write();
}
