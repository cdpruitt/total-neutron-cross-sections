#ifndef IDENTIFY_MACROPULSES_H
#define IDENTIFY_MACROPULSES_H

#include <string>

int identifyMacropulses(
        std::string inputFileName,
        std::string macropulseTimeTreeName,
        std::string targetPositionTreeName,
        std::string outputFileName,
        std::string outputTreeName);

#endif /* IDENTIFY_MACROPULSES_H */
