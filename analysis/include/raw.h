#ifndef RAW_H
#define RAW_H

#include <string>
#include <fstream>

#include "../include/dataStructures.h"

int readRawData(std::string inFileName, std::string outFileName, std::vector<std::string> channelMap);
bool readEvent(std::ifstream& evtfile, RawEvent& rawEvent);
bool readEventHeader(std::ifstream& evtfile, RawEvent& rawEvent);
bool readDPPEventBody(std::ifstream& evtfile, RawEvent& rawEvent);
bool readWaveformEventBody(std::ifstream& evtfile, RawEvent& rawEvent);

#endif /* RAW_H */
