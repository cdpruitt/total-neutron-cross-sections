#ifndef SEPARATE_H
#define SEPARATE_H

#include "TTree.h"

void separateByChannel(std::string rawFileName, std::string tempFileName, std::vector<TTree*>& orchardRaw, std::vector<TTree*>& orchardRawW);

void processDPPEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRaw, std::vector<TTree*>& orchardProcessed);

void processWaveformEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRawW, std::vector<TTree*>& orchardProcessedW);

void processTargetChanger(std::string rawFileName, TFile*& sortedFile);

#endif
