#ifndef IDENTIFY_GOOD_MACROS_H
#define IDENTIFY_GOOD_MACROS_H

#include <fstream>
#include <vector>
#include <string>

#include "dataStructures.h"

int identifyGoodMacros(std::string macropulseFileName, std::vector<MacropulseEvent>& macropulseList, std::ofstream& logFile);

#endif /* IDENTIFY_GOOD_MACROS_H */
