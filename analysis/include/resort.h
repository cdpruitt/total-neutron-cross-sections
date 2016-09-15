#ifndef RESORT_H
#define RESORT_H

#include "TTree.h"

void separateByChannel(std::string rawFileName, std::string tempFileName, std::vector<TTree*>& orchardRaw, std::vector<TTree*>& orchardRawW);

void processDPPEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRaw, std::vector<TTree*>& orchardProcessed, std::ofstream& error);

void processWaveformEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRawW, std::vector<TTree*>& orchardProcessedW, std::ofstream& error);

void processTargetChanger(std::string rawFileName, TFile*& sortedFile, std::ofstream& error);

#endif
