#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "../include/dataStructures.h"
#include "../include/branches.h"

using namespace std;

const double VETO_WINDOW = 30; // in ns

void vetoEvents(string sortedFileName, string vetoedFileName, vector<string> eventTreeNames, string vetoTreeName)
{
    TFile* sortedFile = new TFile(sortedFileName.c_str(),"READ");

    TTree* vetoTree = (TTree*)sortedFile->Get(vetoTreeName.c_str());
    ProcessedEvent vetoEvent;
    setBranchesVeto(vetoTree, vetoEvent);

    TFile* vetoedFile = new TFile(vetoedFileName.c_str(),"RECREATE");

    for(string s : eventTreeNames)
    {
        TTree* eventTree = (TTree*)sortedFile->Get(s.c_str());
        ProcessedEvent procEvent;
        setBranchesProcessed(eventTree, procEvent);

        // cleanTree holds events that have survived the veto
        string clean = s + "Clean";
        TTree* cleanTree = new TTree(clean.c_str(),clean.c_str());
        branchProc(cleanTree, procEvent);

        // dirtyTree holds events that have been vetoed
        string dirty = s + "Dirty";
        TTree* dirtyTree = new TTree(dirty.c_str(),dirty.c_str());
        branchProc(dirtyTree, procEvent);

        long eventTreeEntries = eventTree->GetEntries();
        long vetoTreeEntries = vetoTree->GetEntries();

        // load the first veto event; j is the veto event counter
        int j=0;
        vetoTree->GetEntry(j);

        // if we run off the end of the vetoTree, accept all remaining events
        bool endOfVetoTree = false;

        for(int i=0; i<eventTreeEntries; i++)
        {
            eventTree->GetEntry(i);

            if(!endOfVetoTree)
            {
                // shift vetoEvent up to macropulse of current event
                while(vetoEvent.macroNo<procEvent.macroNo)
                {
                    j++;
                    if(j<vetoTreeEntries)
                    {
                        vetoTree->GetEntry(j);
                    }

                    else
                    {
                        cout << "Reached end of veto tree - allowing all remaining events." << endl;
                        endOfVetoTree = true;
                        break;
                    }
                }

                if(endOfVetoTree)
                {
                    continue;
                }

                // shift vetoEvent up to timestamp of current event
                while(vetoEvent.macroNo==procEvent.macroNo && vetoEvent.completeTime<procEvent.completeTime-VETO_WINDOW)
                {
                    j++;
                    if(j<vetoTreeEntries)
                    {
                        vetoTree->GetEntry(j);
                    }

                    else
                    {
                        cout << "Reached end of veto tree - allowing all remaining events." << endl;
                        endOfVetoTree = true;
                        break;
                    }
                }

                // test for coincidence, within VETO_WINDOW
                if(procEvent.completeTime-VETO_WINDOW<vetoEvent.completeTime &&
                        procEvent.completeTime+VETO_WINDOW>vetoEvent.completeTime)
                {
                    // coincidence found - throw out event
                    dirtyTree->Fill();
                    continue;
                }
            }

            // event survived veto, so add to cleanTree
            cleanTree->Fill();

            if(i%10000==0)
            {
                cout << "Processed " << i << " events on " << s << " through veto \r";
                fflush(stdout);
            }
        }

        cleanTree->Write();
        dirtyTree->Write();
    }

    vetoedFile->Close();
}
