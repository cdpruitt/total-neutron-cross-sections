#define analyze_cxx
#include "analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <vector>
#include <string>

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
    std::vector<TDirectoryFile*> direct;

    TDirectory *targetChangerDir;
    TDirectory *monitorDir;
    TDirectory *detectorTDir;
    TDirectory *detectorLDir;
    TDirectory *detectorRDir;

    TH1S* outMacro;
    TH1S* outEvt;
    TH1S* outTime;
    TH1S* outFT;
    TH1S* outSGQ;
    TH1S* outLGQ;
    vector<TH1S*> listWaveforms;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        fChain->GetEntry(jentry);
        if(temprunNo != runNo) // new run; new directory; new histos
        {
            //std::string s = std::to_string(runNo);
            direct.push_back(new TDirectoryFile("test","test"));
            direct.back()->cd();

            targetChangerDir = new TDirectory("targetChanger","Target Changer");
            monitorDir = new TDirectory("monitor","Monitor");
            detectorTDir = new TDirectory("detT","Detector sum");
            detectorLDir = new TDirectory("detL","Detector left)");
            detectorRDir = new TDirectory("detR","Detector right");

            // instantiate histograms
            outMacro = new TH1S("outMacro","outMacro",100000,0,10000000);
            outEvt = new TH1S("outEvt","outEvt",500,0,1000);
            outTime = new TH1S("outTime","outTime",100000,0,100000000);
            outSGQ = new TH1S("outSGQ","outSGQ",1024,0,70000);
            outLGQ = new TH1S("outLGQ","outLGQ",1024,0,70000);
            outFT = new TH1S("outFT","outFT",1023,0,1023);
            listWaveforms.clear();
        }

        temprunNo = runNo;
            
        outMacro->Fill(macroNo);
        outEvt->Fill(evtNo);
        outTime->Fill(timetag);
        outSGQ->Fill(sgQ);
        outLGQ->Fill(lgQ);
        outFT->Fill(fineTime);

        listWaveforms.push_back(new TH1S("waveform","waveform",waveform.size(),0,waveform.size()*2));
        for(int i=0;i<waveform.size();i++)
        {
            listWaveforms.back()->SetBinContent(i,waveform[i])
        }

        //Long64_t ientry = LoadTree(jentry);
        //if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }
}
