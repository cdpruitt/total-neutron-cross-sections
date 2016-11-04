#include <iostream>

#include "TTree.h"

#include "../include/dataStructures.h"
#include "../include/branches.h"

extern ProcessedEvent procEvent;
extern ProcessedEvent vetoEvent;

using namespace std;

double VETO_WINDOW = 100; // in ns

void vetoEvents(TTree*& eventTree, TTree*& vetoTree)
{
    setBranchesProcessed(eventTree);
    setBranchesVeto(vetoTree);

    // cleanTree holds events that have survived the veto
    TTree* cleanTree = new TTree("cleanTree","cleanTree");
    branchProc(cleanTree);

    // vetoTree holds events that have been vetoed
    TTree* dirtyTree = new TTree("dirtyTree","dirtyTree");
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
    }
}
