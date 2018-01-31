#ifndef ASSIGN_EVENTS_TO_MACROPULSES_H
#define ASSIGN_EVENTS_TO_MACROPULSES_H

#include <string>
#include <utility>
#include "dataStructures.h"

int assignEventsToMacropulses(
        std::string inputFileName,
        std::string outputFileName,
        std::ofstream& log,
        std::vector<MacropulseEvent>& macropulseList);

#endif /* ASSIGN_EVENTS_TO_MACROPULSES_H */

