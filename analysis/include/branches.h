#ifndef BRANCHES_H
#define BRANCHES_H

void branchRaw(TTree*& tree);
void branchSplit(TTree*& tree);
void branchSplitW(TTree*& tree);
void branchProc(TTree*& tree);
void branchProcW(TTree*& tree);
void branchTargetChanger(TTree*& tree);
void setBranchesSeparated(TTree* tree);
void setBranchesSeparatedW(TTree* tree);
void setBranchesProcessed(TTree* tree);
void setBranchesProcessedW(TTree* tree);
void setBranchesVeto(TTree* tree);
void setBranchesProcessedTC(TTree* tree);
void setBranchesHistos(TTree* tree);
void setBranchesHistosW(TTree* tree);
void setBranchesTC(TTree* tree);

#endif
