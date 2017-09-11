/******************************************************************************/
// Define tree structures and functions to fill/read from trees
/******************************************************************************/

#include "TTree.h"
#include "../include/dataStructures.h"

RawEvent rawEvent;
SeparatedEvent separatedEvent;
ProcessedEvent procEvent;
ProcessedEvent vetoEvent;
TargetChangerEvent tcEvent;

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchRaw(TTree*& tree)
{
    tree->Branch("fineTime",&rawEvent.fineTime,"fineTime/d");
    tree->Branch("evtType",&rawEvent.evtType,"evtType/i");
    tree->Branch("chNo",&rawEvent.chNo,"chNo/i");
    tree->Branch("extTime",&rawEvent.extTime,"extTime/i");
    tree->Branch("timetag",&rawEvent.timetag,"timetag/i");
    tree->Branch("sgQ",&rawEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&rawEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&rawEvent.waveform);
}

// Used to connect a channel-specific tree to waveform event variables so we can
// start populating it with waveform events
void branchSplitW(TTree*& tree)
{
    tree->Branch("evtNo",&separatedEvent.evtNo,"evtNo/i");
    tree->Branch("timetag",&separatedEvent.timetag,"timetag/i");
    tree->Branch("extTime",&separatedEvent.extTime,"extTime/i");
    tree->Branch("chNo",&separatedEvent.chNo,"chNo/i");
    tree->Branch("waveform",&separatedEvent.waveform);
}

// Used to connect a channel-specific tree to DPP event variables so we can
// start populating it with DPP events
void branchProc(TTree*& tree)
{
    tree->Branch("macroTime",&procEvent.macroTime,"macroTime/D");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/D");
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("sgQ",&procEvent.sgQ,"sgQ/i");
    tree->Branch("lgQ",&procEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&procEvent.waveform);
}

// Used to connect a channel-specific tree to waveform procEvent variables so we can
// start populating it with waveform procEvents
void branchProcW(TTree*& tree)
{
    tree->Branch("macroTime",&procEvent.macroTime,"macroTime/D");
    tree->Branch("completeTime",&procEvent.completeTime,"completeTime/D");
    tree->Branch("macroNo",&procEvent.macroNo,"macroNo/i");
    tree->Branch("evtNo",&procEvent.evtNo,"evtNo/i");
    tree->Branch("targetPos",&procEvent.targetPos,"targetPos/i");
    tree->Branch("waveform",&procEvent.waveform);
}

void branchTargetChanger(TTree*& tree)
{
    tree->Branch("macroTime",&tcEvent.macroTime,"macroTime/D");
    tree->Branch("fineTime",&tcEvent.fineTime,"fineTime/d");
    tree->Branch("macroNo",&tcEvent.macroNo,"macroNo/i");
    tree->Branch("modeChange",&tcEvent.modeChange,"modeChange/i");
    tree->Branch("targetPos",&tcEvent.targetPos,"targetPos/i");
    tree->Branch("lgQ",&tcEvent.lgQ,"lgQ/i");
    tree->Branch("waveform",&tcEvent.waveform);
}

void setBranchesSeparated(TTree* tree)
{
   tree->SetBranchAddress("evtType",&separatedEvent.evtType);
   tree->SetBranchAddress("chNo",&separatedEvent.chNo);
   tree->SetBranchAddress("extTime",&separatedEvent.extTime);
   tree->SetBranchAddress("timetag",&separatedEvent.timetag);
   tree->SetBranchAddress("fineTime",&separatedEvent.fineTime);
   tree->SetBranchAddress("sgQ",&separatedEvent.sgQ);
   tree->SetBranchAddress("lgQ",&separatedEvent.lgQ);
   tree->SetBranchAddress("waveform",&separatedEvent.waveform);
}

void setBranchesSeparatedW(TTree* tree)
{
    tree->SetBranchAddress("timetag",&separatedEvent.timetag);
    tree->SetBranchAddress("extTime",&separatedEvent.extTime);
    tree->SetBranchAddress("evtNo",&separatedEvent.evtNo);
    tree->SetBranchAddress("waveform",&separatedEvent.waveform);
}

void setBranchesProcessed(TTree* tree)
{
   tree->SetBranchAddress("macroNo",&procEvent.macroNo);
   tree->SetBranchAddress("macroTime",&procEvent.macroTime);
   tree->SetBranchAddress("evtNo",&procEvent.evtNo);
   tree->SetBranchAddress("completeTime",&procEvent.completeTime);
   tree->SetBranchAddress("fineTime",&procEvent.fineTime);
   tree->SetBranchAddress("targetPos",&procEvent.targetPos);
   tree->SetBranchAddress("sgQ",&procEvent.sgQ);
   tree->SetBranchAddress("lgQ",&procEvent.lgQ);
   tree->SetBranchAddress("waveform",&procEvent.waveform);
}

void setBranchesProcessedW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("evtNo",&procEvent.evtNo);
    tree->SetBranchAddress("completeTime",&procEvent.completeTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}

void setBranchesHistos(TTree* tree)
{
   tree->SetBranchAddress("macroNo",&procEvent.macroNo);
   tree->SetBranchAddress("macroTime",&procEvent.macroTime);
   tree->SetBranchAddress("evtNo",&procEvent.evtNo);
   tree->SetBranchAddress("completeTime",&procEvent.completeTime);
   tree->SetBranchAddress("targetPos",&procEvent.targetPos);
   tree->SetBranchAddress("sgQ",&procEvent.sgQ);
   tree->SetBranchAddress("lgQ",&procEvent.lgQ);
   tree->SetBranchAddress("waveform",&procEvent.waveform);
}

void setBranchesHistosW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("evtNo",&procEvent.evtNo);
    tree->SetBranchAddress("completeTime",&procEvent.completeTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}

void setBranchesVeto(TTree* tree)
{
   tree->SetBranchAddress("macroNo",&vetoEvent.macroNo);
   tree->SetBranchAddress("macroTime",&vetoEvent.macroTime);
   tree->SetBranchAddress("evtNo",&vetoEvent.evtNo);
   tree->SetBranchAddress("completeTime",&vetoEvent.completeTime);
   tree->SetBranchAddress("targetPos",&vetoEvent.targetPos);
   tree->SetBranchAddress("sgQ",&vetoEvent.sgQ);
   tree->SetBranchAddress("lgQ",&vetoEvent.lgQ);
   tree->SetBranchAddress("waveform",&vetoEvent.waveform);
}

void setBranchesProcessedTC(TTree* tree)
{
    tree->SetBranchAddress("macroTime",&tcEvent.macroTime);
    tree->SetBranchAddress("fineTime",&tcEvent.fineTime);
    tree->SetBranchAddress("macroNo",&tcEvent.macroNo);
    tree->SetBranchAddress("modeChange",&tcEvent.modeChange);
    tree->SetBranchAddress("targetPos",&tcEvent.targetPos);
    tree->SetBranchAddress("lgQ",&tcEvent.lgQ);
    tree->SetBranchAddress("waveform",&tcEvent.waveform);
}

void setBranchesTC(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&procEvent.macroNo);
    tree->SetBranchAddress("macroTime",&procEvent.macroTime);
    tree->SetBranchAddress("fineTime",&procEvent.fineTime);
    tree->SetBranchAddress("targetPos",&procEvent.targetPos);
    tree->SetBranchAddress("waveform",&procEvent.waveform);
}
