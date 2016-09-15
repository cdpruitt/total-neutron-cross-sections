#ifndef BRANCHES_H
#define BRANCHES_H

void branchRaw(TTree*& tree);
void branchRawW(TTree*& tree);
void branchProc(TTree*& tree);
void branchProcW(TTree*& tree);
void branchTargetChanger(TTree*& tree);
void setBranches(TTree* tree);
void setBranchesW(TTree* tree);

#endif
