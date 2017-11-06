#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

#include "../include/dataStructures.h"
#include "../include/branches.h"

using namespace std;

const double VETO_WINDOW = 8; // in ns

int vetoEvents(string detectorFileName, string outputFileName, ofstream& logFile, string detTreeName, string vetoTreeName)
{
    // check to see if output file already exists; if so, exit
    ifstream f(outputFileName);

    if(f.good())
    {
        cout << outputFileName << " already exists; skipping vetoing of events." << endl;
        logFile << outputFileName << " already exists; skipping vetoing of events." << endl;
        return 0;
    }

    f.close();

    TFile* detectorFile = new TFile(detectorFileName.c_str(),"READ");

    TTree* vetoTree = (TTree*)detectorFile->Get(vetoTreeName.c_str());

    DetectorEvent veto;
    vetoTree->SetBranchAddress("cycleNumber",&veto.cycleNumber);
    vetoTree->SetBranchAddress("completeTime",&veto.completeTime);

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    DetectorEvent event;
    vector<int>* waveformPointer = 0;

    TTree* detTree = (TTree*)detectorFile->Get(detTreeName.c_str());
    detTree->SetBranchAddress("cycleNumber",&event.cycleNumber);
    detTree->SetBranchAddress("macroNo",&event.macroNo);
    detTree->SetBranchAddress("macroTime",&event.macroTime);
    detTree->SetBranchAddress("fineTime",&event.fineTime);
    detTree->SetBranchAddress("eventNo",&event.eventNo);
    detTree->SetBranchAddress("completeTime",&event.completeTime);
    detTree->SetBranchAddress("targetPos",&event.targetPos);
    detTree->SetBranchAddress("sgQ",&event.sgQ);
    detTree->SetBranchAddress("lgQ",&event.lgQ);
    detTree->SetBranchAddress("waveform",&waveformPointer);

    // cleanTree holds events that have survived the veto
    string clean = detTreeName;
    TTree* cleanTree = new TTree(clean.c_str(),clean.c_str());
    cleanTree->Branch("cycleNumber",&event.cycleNumber, "cycleNumber/i");
    cleanTree->Branch("macroNo",&event.macroNo, "macroNo/i");
    cleanTree->Branch("macroTime",&event.macroTime, "macroTime/d");
    cleanTree->Branch("fineTime",&event.fineTime, "fineTime/d");
    cleanTree->Branch("eventNo",&event.eventNo, "eventNo/i");
    cleanTree->Branch("completeTime",&event.completeTime, "completeTime/d");
    cleanTree->Branch("targetPos",&event.targetPos, "targetPos/I");
    cleanTree->Branch("sgQ",&event.sgQ, "sgQ/i");
    cleanTree->Branch("lgQ",&event.lgQ, "lgQ/i");
    cleanTree->Branch("waveform",&waveformPointer);

    // dirtyTree holds events that have been vetoed
    string dirty = detTreeName + "Dirty";
    TTree* dirtyTree = new TTree(dirty.c_str(),dirty.c_str());
    dirtyTree->Branch("cycleNumber",&event.cycleNumber, "cycleNumber/i");
    dirtyTree->Branch("macroNo",&event.macroNo, "macroNo/i");
    dirtyTree->Branch("macroTime",&event.macroTime, "macroTime/d");
    dirtyTree->Branch("fineTime",&event.fineTime, "fineTime/d");
    dirtyTree->Branch("eventNo",&event.eventNo, "eventNo/i");
    dirtyTree->Branch("completeTime",&event.completeTime, "completeTime/d");
    dirtyTree->Branch("targetPos",&event.targetPos, "targetPos/I");
    dirtyTree->Branch("sgQ",&event.sgQ, "sgQ/i");
    dirtyTree->Branch("lgQ",&event.lgQ, "lgQ/i");
    dirtyTree->Branch("waveform",&waveformPointer);

    long detTreeEntries = detTree->GetEntries();
    long vetoTreeEntries = vetoTree->GetEntries();

    // load the first veto event; j is the veto event counter
    int j=0;
    vetoTree->GetEntry(j);

    for(int i=0; i<detTreeEntries; i++)
    {
        detTree->GetEntry(i);

        // shift vetoEvent up to macropulse of current event
        while(veto.cycleNumber < event.cycleNumber)
        {
            j++;

            if(j<vetoTreeEntries)
            {
                vetoTree->GetEntry(j);
            }

            else
            {
                cout << "Reached end of veto tree - allowing all remaining events." << endl;
                unsigned int k=i;
                while(k<detTreeEntries)
                {
                    detTree->GetEntry(k);
                    cleanTree->Fill();
                    k++;

                    if(k%10000==0)
                    {
                        cout << "Processed " << k << " events on " << detTreeName << " through veto \r";
                        fflush(stdout);
                    }
                }

                cleanTree->Write();
                dirtyTree->Write();

                double numberCleanEvents = cleanTree->GetEntries();
                logFile << "Fraction of events surviving veto: " << numberCleanEvents/detTreeEntries << endl;

                detectorFile->Close();
                outputFile->Close();

                return 0;
            }
        }

        while(veto.cycleNumber == event.cycleNumber &&
                veto.completeTime < event.completeTime-VETO_WINDOW)
        {
            j++;

            if(j<vetoTreeEntries)
            {
                vetoTree->GetEntry(j);
            }

            else
            {
                cout << "Reached end of veto tree - allowing all remaining events." << endl;
                unsigned int k=i;
                while(k<detTreeEntries)
                {
                    detTree->GetEntry(k);
                    cleanTree->Fill();
                    k++;

                    if(k%10000==0)
                    {
                        cout << "Processed " << k << " events on " << detTreeName << " through veto \r";
                        fflush(stdout);
                    }
                }

                cleanTree->Write();
                dirtyTree->Write();

                double numberCleanEvents = cleanTree->GetEntries();
                logFile << "Fraction of events surviving veto: " << numberCleanEvents/detTreeEntries << endl;

                detectorFile->Close();
                outputFile->Close();

                return 0;
            }
        }

        // test for coincidence, within VETO_WINDOW
        if(abs(event.completeTime-veto.completeTime)<VETO_WINDOW)
        {
            // coincidence found - throw out event
            dirtyTree->Fill();
            continue;
        }

        // event survived veto, so add to cleanTree
        cleanTree->Fill();

        if(i%10000==0)
        {
            cout << "Processed " << i << " events on " << detTreeName << " through veto \r";
            fflush(stdout);
        }
    }

    cleanTree->Write();
    dirtyTree->Write();

    double numberCleanEvents = cleanTree->GetEntries();
    logFile << "Fraction of events surviving veto: " << numberCleanEvents/detTreeEntries << endl;

    detectorFile->Close();
    outputFile->Close();

    return 0;
}
