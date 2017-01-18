#ifndef HISTOS_H
#define HISTOS_H

#include "TH1I.h"
#include "plots.h"

int histos(std::string sortedFileName, std::string vetoedFileName, std::string histoFileName, std::vector<std::string> channelMap);
TH1I* convertTOFtoEn(TH1I* tof, std::string name);

#endif
