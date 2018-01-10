#ifndef IDENTIFY_GOOD_MACROS_H
#define IDENTIFY_GOOD_MACROS_H

#include <fstream>
#include <vector>
#include <string>

#include "TTree.h"

#include "dataStructures.h"

int identifyGoodMacros(std::string inputFileName, std::vector<MacropulseEvent>& macropulseList, std::ofstream& logFile);

#endif /* IDENTIFY_GOOD_MACROS_H */
