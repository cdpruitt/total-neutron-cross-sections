#ifndef BRANCHES_H
#define BRANCHES_H

void branchRaw(TTree*& tree, RawEvent& rawEvent);
void branchProc(TTree*& tree, ProcessedEvent& procEvent);
void branchProcW(TTree*& tree, ProcessedEvent& procEvent);
void branchTargetChanger(TTree*& tree, TargetChangerEvent& tcEvent);

void setBranchesSeparated(TTree* tree, SeparatedEvent& separatedEvent);
void setBranchesSeparatedW(TTree* tree, SeparatedEvent& separatedEvent);

void setBranchesTC(TTree* tree, ProcessedEvent& procEvent);

void setBranchesProcessed(TTree* tree, ProcessedEvent& procEvent);
void setBranchesProcessedW(TTree* tree, ProcessedEvent& procEvent);
void setBranchesProcessedTC(TTree* tree, TargetChangerEvent& tcEvent);

void setBranchesVeto(TTree* tree, ProcessedEvent& vetoEvent);

void setBranchesHistos(TTree* tree, ProcessedEvent& procEvent);
void setBranchesHistosW(TTree* tree, ProcessedEvent& procEvent);

#endif
