#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "../include/config.h"

#include "../include/dataStructures.h"
#include "../include/branches.h"

using namespace std;

extern Config config;

const double VETO_WINDOW = 3; // in ns

int vetoEvents(string detectorFileName, string outputFileName, ofstream& logFile, string vetoTreeName)
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

    // create output file
    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");

    TFile* detectorFile = new TFile(detectorFileName.c_str(),"READ");
    TTree* vetoTree = (TTree*)detectorFile->Get(vetoTreeName.c_str());

    DetectorEvent veto;
    vetoTree->SetBranchAddress("cycleNumber",&veto.cycleNumber);
    vetoTree->SetBranchAddress("completeTime",&veto.completeTime);

    for(string detTreeName : config.cs.DETECTOR_NAMES)
    {
        // read input trees to prepare for vetoing
        TTree* detTree = (TTree*)detectorFile->Get(detTreeName.c_str());

        DetectorEvent event;
        vector<int>* waveformPointer = 0;

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

        outputFile->cd();

        // output tree holds events that have survived the veto
        TTree* tree = new TTree(detTreeName.c_str(),detTreeName.c_str());
        tree->Branch("cycleNumber",&event.cycleNumber, "cycleNumber/i");
        tree->Branch("macroNo",&event.macroNo, "macroNo/i");
        tree->Branch("macroTime",&event.macroTime, "macroTime/d");
        tree->Branch("fineTime",&event.fineTime, "fineTime/d");
        tree->Branch("eventNo",&event.eventNo, "eventNo/i");
        tree->Branch("completeTime",&event.completeTime, "completeTime/d");
        tree->Branch("targetPos",&event.targetPos, "targetPos/I");
        tree->Branch("sgQ",&event.sgQ, "sgQ/i");
        tree->Branch("lgQ",&event.lgQ, "lgQ/i");
        tree->Branch("vetoed",&event.vetoed,"vetoed/O");

        tree->Branch("waveform",&waveformPointer);

        TH1D* vetoedEventHisto = new TH1D("vetoed event time diff",
            "vetoed event time diff", 100*VETO_WINDOW, -10*VETO_WINDOW, 10*VETO_WINDOW);

        long detTreeEntries = detTree->GetEntries();
        long vetoTreeEntries = vetoTree->GetEntries();
        double numberVetoedEvents = 0;

        bool endVeto = false;

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
                        tree->Fill();
                        k++;

                        if(k%10000==0)
                        {
                            cout << "Processed " << k << " events on " << detTreeName << " through veto \r";
                            fflush(stdout);
                        }
                    }

                    vetoedEventHisto->Write();

                    tree->Write();

                    logFile << "Fraction of events surviving veto: " << (detTreeEntries-numberVetoedEvents)/detTreeEntries << endl;

                    endVeto = true;
                    break;
                }
            }

            if(endVeto)
            {
                break;
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
                        tree->Fill();
                        k++;

                        if(k%10000==0)
                        {
                            cout << "Processed " << k << " events on " << detTreeName << " through veto \r";
                            fflush(stdout);
                        }
                    }

                    vetoedEventHisto->Write();

                    tree->Write();

                    logFile << "Fraction of events surviving veto: " << (detTreeEntries-numberVetoedEvents)/detTreeEntries << endl;

                    endVeto = true;
                    break;
                }
            }

            if(endVeto)
            {
                break;
            }

            // test for coincidence, within VETO_WINDOW
            if(abs(event.completeTime-veto.completeTime)<10*VETO_WINDOW)
            {
                vetoedEventHisto->Fill(event.completeTime-veto.completeTime);
            }

            if(abs(event.completeTime-veto.completeTime)<VETO_WINDOW)
            {
                // coincidence found - mark event as "vetoed"
                event.vetoed = true;
                numberVetoedEvents++;
            }

            else
            {
                event.vetoed = false;
            }

            tree->Fill();

            if(i%10000==0)
            {
                cout << "Processed " << i << " events on " << detTreeName << " through veto \r";
                fflush(stdout);
            }
        }

        vetoedEventHisto->Write();

        tree->Write();

        logFile << "Fraction of events surviving veto: " << (detTreeEntries-numberVetoedEvents)/detTreeEntries << endl;
    }

    detectorFile->Close();
    outputFile->Close();

    return 0;
}
