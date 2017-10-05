#ifndef ASSIGN_EVENTS_TO_MACROPULSES_H
#define ASSIGN_EVENTS_TO_MACROPULSES_H

#include <string>

int assignEventsToMacropulses(
        std::string inputFileName,
        std::string inputTreeName,
        std::string outputFileName,
        std::string macropulseTreeName,
        unsigned int channelNo,
        std::string outputTreeName);

#endif /* ASSIGN_EVENTS_TO_MACROPULSES_H */

