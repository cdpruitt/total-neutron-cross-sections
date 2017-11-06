#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "config.h"

#include <vector>
#include <string>
#include <utility>

std::vector<std::pair<unsigned int, std::string>> getChannelMap(std::string expName, unsigned int runNumber);
std::vector<std::string> getTargetOrder(std::string expName, int runNumber);
std::vector<std::pair<std::string,std::string>> getRelativePlotNames(std::string expName,std::string fileName);
std::vector<std::string> readExperimentConfig(std::string expName, std::string fileName);

AnalysisConfig readAnalysisConfig(std::string expName);
FacilityConfig readFacilityConfig(std::string expName, int runNumber);
SoftwareCFDConfig readSoftwareCFDConfig(std::string expName, int runNumber);
TimeOffsetsConfig readTimeOffsetsConfig(std::string expName, int runNumber);
DigitizerConfig readDigitizerConfig(std::string expName, int runNumber);
CSConfig readCSConfig(std::string expName, int runNumber);
PlotConfig readPlotConfig(std::string expName, int runNumber);
TargetConfig readTargetConfig(std::string expName, int runNumber);


#endif
