#include <iostream>

#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/branches.h"

extern ProcessedEvent procEvent;
extern ProcessedEvent vetoEvent;

using namespace std;

double VETO_WINDOW = 10; // in ns

void vetoEvents(TTree*& eventTree, TTree*& vetoTree, string name)
{
    setBranchesProcessed(eventTree);
    setBranchesVeto(vetoTree);

    // cleanTree holds events that have survived the veto
    string clean = name + "Clean";
    TTree* cleanTree = new TTree(clean.c_str(),clean.c_str());
    branchProc(cleanTree);

    // dirtyTree holds events that have been vetoed
    string dirty = name + "Dirty";
    TTree* dirtyTree = new TTree(dirty.c_str(),dirty.c_str());
    branchProc(dirtyTree);

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
            cout << "Processed " << i << " events on " << name << " through veto \r";
            fflush(stdout);
        }
    }

    cleanTree->Write();
    dirtyTree->Write();
}
