#ifndef RAW_H
#define RAW_H

#include <string>
#include <fstream>

void readRawData(std::string inFileName, std::string outFileName, std::vector<std::string> channelMap);
bool readEvent(std::ifstream& evtfile);
bool readEventHeader(std::ifstream& evtfile);
bool readDPPEventBody(std::ifstream& evtfile);
bool readWaveformEventBody(std::ifstream& evtfile);
double calculateCFDTime(std::vector<int>* waveform,
        double baseline,
        double fraction,
        unsigned int delay,
        bool isPositiveSignal);

#endif
