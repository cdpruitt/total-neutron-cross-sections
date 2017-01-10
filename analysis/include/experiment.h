#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include <string>
#include <utility>

std::vector<std::string> getChannelMap(std::string expName, int runNumber);
std::vector<std::string> getTargetOrder(std::string expName, int runNumber);
std::vector<std::pair<std::string,std::string>> getRelativePlotNames(std::string expName,std::string fileName);
std::vector<std::string> readExperimentConfig(std::string expName, std::string fileName);

#endif