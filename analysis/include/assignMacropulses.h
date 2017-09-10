#ifndef SEPARATE_H
#define SEPARATE_H

#include "TTree.h"

void assignMacropulses(std::string rawFileName, std::string sortedFileName, std::vector<std::string> channelMap);

void processDPPEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRaw, std::vector<TTree*>& orchardProcessed);

void processWaveformEvents(TFile*& sortedFile, std::vector<TTree*>& orchardRawW, std::vector<TTree*>& orchardProcessedW);

void processTargetChanger(std::string rawFileName, TFile*& sortedFile);

double calculateFineTime(std::vector<int>* waveform, unsigned int threshold, bool isPositiveSignal);

#endif
