#ifndef IDENTIFY_MACROPULSES_H
#define IDENTIFY_MACROPULSES_H

#include <string>
#include "../include/dataStructures.h"

int identifyMacropulses(
        std::string inputFileName,
        std::string outputFileName,
        std::ofstream& logFile,
        std::vector<MacropulseEvent>& macropulseList);

#endif /* IDENTIFY_MACROPULSES_H */
