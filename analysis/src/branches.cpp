/******************************************************************************/
// Define tree structures and functions to fill/read from trees
/******************************************************************************/

#include "TTree.h"
#include "../include/dataStructures.h"

RawEvent rawEvent;
SortedEvent sortedEvent;
ProcessedEvent procEvent;
TargetChangerEvent tcEvent;

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchRaw(TTree*& tree)
{
    tree->Branch("evtNo",&sortedEvent.evtNo,"evtNo/i");
    tree->Branch("chNo",&sortedEvent.chNo,"chNo/i");
    tree->Branch("timetag",&sortedEvent.timetag,"timetag/d");
    tree->Branch("extTime",&sortedEvent.extTime,"extTime/i");
    tree->Branch("fineTime",&sortedEvent.fineTime,"fineTime/i");
    tree->Branch("sgQ",&sortedEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&sortedEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&sortedEvent.waveform);
}

// Used to connect a channel-specific tree to waveform event variables so we can
// start populating it with waveform events
void branchRawW(TTree*& tree)
{
    tree->Branch("timetag",&sortedEvent.timetag,"timetag/d");
    tree->Branch("extTime",&sortedEvent.extTime,"extTime/i");
    tree->Branch("chNo",&sortedEvent.chNo,"chNo/i");
    tree->Branch("evtNo",&sortedEvent.evtNo,"evtNo/i");
    tree->Branch("waveform",&sortedEvent.waveform);
}

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchProc(TTree*& tree)
{
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("macroTime",&procEvent.macroTime,"macroTime/d");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/d");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("sgQ",&procEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&procEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&procEvent.waveform);
}

// Used to connect a channel-specific tree to waveform procEvent variables so we can
// start populating it with waveform procEvents
void branchProcW(TTree*& tree)
{
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/d");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("waveform",&procEvent.waveform);
}

void branchTargetChanger(TTree*& tree)
{
    tree->Branch("macroNo",&tcEvent.macroNo,"macroNo/i");
    tree->Branch("macroTime",&tcEvent.macroTime,"macroTime/d");
    tree->Branch("modeChange",&tcEvent.modeChange,"modeChange/i");
    tree->Branch("targetPos",&tcEvent.targetPos,"targetPos/i");
}


// Re-link to an already-existing tree's data so we can read the tree
void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("evtNo",&procEvent.evtNo);
    tree->SetBranchAddress("macroTime",&procEvent.macroTime);
    tree->SetBranchAddress("completeTime",&procEvent.completeTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("sgQ",&procEvent.sgQ);
    tree->SetBranchAddress("lgQ",&procEvent.lgQ);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}

// Re-link to an already-existing tree's data so we can read the tree
void setBranchesW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("evtNo",&procEvent.evtNo);
    tree->SetBranchAddress("completeTime",&procEvent.completeTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}

void setTCBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("macroTime",&procEvent.macroTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}
