#ifndef ASSIGN_EVENTS_TO_MACROPULSES_H
#define ASSIGN_EVENTS_TO_MACROPULSES_H

#include <string>
#include <utility>

int assignEventsToMacropulses(
        std::string inputFileName,
        std::string outputFileName,
        std::ofstream& log,
        std::pair<unsigned int, std::string> channel);

#endif /* ASSIGN_EVENTS_TO_MACROPULSES_H */

